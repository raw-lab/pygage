#!/usr/bin/env python3
"""
Data Processing and Expression Analysis Utilities (corrected).

The critical fix here is :meth:`DataTransformer.row_normalize`, which now
z-scores each *gene* (row) across samples, matching its name/docstring and the
behaviour GAGE and heatmap conventions expect.  The original implementation
z-scored each *column* (sample) across genes -- a transposed operation that
silently produced misleading heatmaps and per-gene statistics.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from .visualization_utils import HeatmapPlotter


class DataTransformer:
    """Data transformation utilities."""

    @staticmethod
    def row_normalize(data: pl.DataFrame, gene_col: str = "gene_id") -> pl.DataFrame:
        """Z-score each gene (row) across its sample columns.

        Args:
            data: DataFrame with an optional ``gene_col`` plus numeric samples.
            gene_col: Name of the ID column to leave untouched (if present).

        Returns:
            DataFrame with the same shape; numeric columns replaced by
            per-row z-scores.  Rows with zero variance are centred only.
        """
        num_cols = [
            c
            for c in data.columns
            if c != gene_col and data[c].dtype.is_numeric()
        ]
        if not num_cols:
            return data

        mat = data.select(num_cols).to_numpy().astype(float)
        row_mean = np.nanmean(mat, axis=1, keepdims=True)
        row_std = np.nanstd(mat, axis=1, ddof=0, keepdims=True)
        safe_std = np.where(row_std == 0, 1.0, row_std)
        z = (mat - row_mean) / safe_std

        result = data.clone()
        for j, c in enumerate(num_cols):
            result = result.with_columns(pl.Series(c, z[:, j]))
        return result

    @staticmethod
    def column_normalize(data: pl.DataFrame, gene_col: str = "gene_id") -> pl.DataFrame:
        """Z-score each sample (column) across genes.

        Provided explicitly for the cases where per-sample standardisation is
        actually wanted (e.g. quantile-style preprocessing), so the intent is
        never ambiguous at the call site.
        """
        num_cols = [
            c
            for c in data.columns
            if c != gene_col and data[c].dtype.is_numeric()
        ]
        result = data.clone()
        for c in num_cols:
            v = data[c].to_numpy().astype(float)
            mu, sd = np.nanmean(v), np.nanstd(v, ddof=0)
            result = result.with_columns(
                pl.Series(c, (v - mu) / sd if sd > 0 else v - mu)
            )
        return result

    @staticmethod
    def prepare_paired_data(
        data: pl.DataFrame,
        ref_indices: List[int],
        samp_indices: List[int],
        gene_col: str = "gene_id",
        comparison: str = "paired",
        use_fold: bool = True,
        input_logged: bool = True,
        log_base: float = 2.0,
        pseudocount: float = 1.0,
    ) -> pl.DataFrame:
        """Prepare paired/unpaired fold-change data, preserving ``gene_col``.

        Fold change is a difference on log data (GAGE convention).  Set
        ``input_logged=False`` to log-transform raw values first.  Indices are
        interpreted against the non-ID columns.
        """
        value_cols = [c for c in data.columns if c != gene_col]
        ref_cols = [value_cols[i] for i in ref_indices]
        samp_cols = [value_cols[i] for i in samp_indices]

        def _prep(col: str) -> pl.Expr:
            e = pl.col(col).cast(pl.Float64)
            if use_fold and not input_logged:
                e = (e + pseudocount).log(base=log_base)
            return e

        exprs = []
        if comparison == "paired":
            if len(ref_cols) != len(samp_cols):
                raise ValueError("Paired comparison requires equal ref/samp counts")
            for r, s in zip(ref_cols, samp_cols):
                exprs.append((_prep(s) - _prep(r)).alias(f"{s}_vs_{r}"))
        elif comparison == "unpaired":
            for s in samp_cols:
                for r in ref_cols:
                    exprs.append((_prep(s) - _prep(r)).alias(f"{s}_vs_{r}"))
        else:
            raise ValueError(f"Unknown comparison type: {comparison}")

        keep = [pl.col(gene_col)] if gene_col in data.columns else []
        return data.select(keep + exprs)


class GeneExtractor:
    """Gene extraction utilities."""

    @staticmethod
    def extract_essential_genes(
        gene_set: List[str],
        expression_data: pl.DataFrame,
        gene_col: str = "gene_id",
        threshold: float = 1.0,
        rank_by_abs: bool = False,
    ) -> pl.DataFrame:
        """Return genes from ``gene_set`` whose mean fold change is an outlier.

        The z-score is computed relative to the per-gene mean distribution of
        all measured genes; genes with ``|z| > threshold`` are retained.
        """
        filtered = expression_data.filter(pl.col(gene_col).is_in(gene_set))
        if filtered.height == 0:
            return pl.DataFrame()

        num_cols = [
            c for c in filtered.columns if c != gene_col and filtered[c].dtype.is_numeric()
        ]
        if not num_cols:
            return filtered

        filtered = filtered.with_columns(
            pl.mean_horizontal(num_cols).alias("mean_expression")
        )
        all_means = expression_data.select(
            pl.mean_horizontal(num_cols)
        ).to_numpy().ravel()
        global_mean = float(np.nanmean(all_means))
        global_std = float(np.nanstd(all_means, ddof=1))
        if global_std == 0:
            global_std = 1.0

        filtered = filtered.with_columns(
            ((pl.col("mean_expression") - global_mean) / global_std).alias("z_score")
        )
        if rank_by_abs:
            filtered = filtered.sort(pl.col("z_score").abs(), descending=True)
        else:
            filtered = filtered.sort("mean_expression", descending=True)

        essential = filtered.filter(pl.col("z_score").abs() > threshold)
        return essential.drop(["mean_expression", "z_score"], strict=False)


class GeneDataExporter:
    """Gene data export and visualization."""

    @staticmethod
    def export_gene_data(
        genes: List[str],
        expression_data: pl.DataFrame,
        gene_col: str = "gene_id",
        output_file: Optional[Path] = None,
        create_heatmap: bool = False,
        heatmap_output: Optional[Path] = None,
        normalize: bool = True,
    ) -> None:
        """Write selected genes to CSV/TSV and optionally draw a heatmap."""
        gene_data = expression_data.filter(pl.col(gene_col).is_in(genes))
        if gene_data.height == 0:
            print("Warning: No genes found in expression data")
            return

        if output_file is not None:
            sep = "," if output_file.suffix == ".csv" else "\t"
            gene_data.write_csv(output_file, separator=sep)
            print(f"Gene data exported to {output_file}")

        if create_heatmap:
            num_cols = [
                c for c in gene_data.columns if c != gene_col and gene_data[c].dtype.is_numeric()
            ]
            if not num_cols:
                return
            heatmap_data = gene_data.select([gene_col] + num_cols)
            if normalize:
                heatmap_data = DataTransformer.row_normalize(heatmap_data, gene_col=gene_col)
            HeatmapPlotter.plot_heatmap(
                heatmap_data.select(num_cols),
                row_labels=heatmap_data[gene_col].to_list(),
                col_labels=num_cols,
                output_file=heatmap_output,
                title="Gene Expression Heatmap",
                center=0 if normalize else None,
            )

    @staticmethod
    def create_scatterplot(
        expression_data: pl.DataFrame,
        ref_col: str,
        samp_col: str,
        gene_col: str = "gene_id",
        genes: Optional[List[str]] = None,
        output_file: Optional[Path] = None,
        title: Optional[str] = None,
    ) -> None:
        """Scatter reference vs. sample with an identity line (NaN-safe)."""
        plot_data = (
            expression_data.filter(pl.col(gene_col).is_in(genes))
            if genes
            else expression_data
        )
        x = plot_data[ref_col].to_numpy().astype(float)
        y = plot_data[samp_col].to_numpy().astype(float)

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(x, y, alpha=0.5, s=20)

        finite = np.isfinite(x) & np.isfinite(y)
        if finite.any():
            lo = float(min(x[finite].min(), y[finite].min()))
            hi = float(max(x[finite].max(), y[finite].max()))
            ax.plot([lo, hi], [lo, hi], "r--", alpha=0.5, linewidth=2)

        ax.set_xlabel(f"{ref_col} (Control)", fontsize=12)
        ax.set_ylabel(f"{samp_col} (Experiment)", fontsize=12)
        ax.set_title(title or "Expression Comparison", fontsize=14, weight="bold")
        ax.set_aspect("equal")
        plt.tight_layout()

        if output_file is not None:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            print(f"Scatterplot saved to {output_file}")
        else:
            plt.show()
        plt.close()
