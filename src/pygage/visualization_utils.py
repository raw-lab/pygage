#!/usr/bin/env python3
"""
Visualization Utilities (corrected).

Fixes vs. the original draft
----------------------------
* ``plot_clustered_heatmap`` no longer calls a hard-coded
  ``set_index("gene_id")`` that raised a KeyError whenever the matrix had no
  such column; row labels are applied directly.
* ``plot_heatmap`` no longer defaults row labels to the *column* names.
"""

from __future__ import annotations

import logging
from itertools import product
from pathlib import Path
from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle

logger = logging.getLogger("pygage.viz")


def _as_polars_frame(data) -> pl.DataFrame:
    """Coerce a result table to a Polars frame at the plotting boundary.

    Accepts a Polars or pandas ``DataFrame``, a PyArrow ``Table``, or a column
    dict, and returns a Polars frame. This is the single conversion point for
    the visualization layer, so plotters never assume a particular dataframe
    backend regardless of what an upstream pipeline supplies.
    """
    if isinstance(data, pl.DataFrame):
        return data
    try:
        import pandas as pd

        if isinstance(data, pd.DataFrame):
            return pl.from_pandas(data)
    except ImportError:  # pragma: no cover
        pass
    try:
        import pyarrow as pa

        if isinstance(data, pa.Table):
            return pl.from_arrow(data)
    except ImportError:  # pragma: no cover
        pass
    if isinstance(data, dict):
        return pl.DataFrame(data)
    raise TypeError(
        f"expected a Polars/pandas/Arrow frame or column dict, got {type(data).__name__}"
    )


class ColorUtils:
    @staticmethod
    def create_colormap(
        low: str, mid: Optional[str], high: str, n: int = 256
    ) -> LinearSegmentedColormap:
        colors = [low, high] if mid is None else [low, mid, high]
        return LinearSegmentedColormap.from_list("custom", colors, N=n)

    @staticmethod
    def greenred(n: int = 256) -> LinearSegmentedColormap:
        return ColorUtils.create_colormap("green", "black", "red", n)


class VennDiagram:
    @staticmethod
    def venn_counts(data: pl.DataFrame, include: str = "both") -> pl.DataFrame:
        if include == "up":
            data = data.select([(pl.col(c) > 0).cast(pl.Int32) for c in data.columns])
        elif include == "down":
            data = data.select([(pl.col(c) < 0).cast(pl.Int32) for c in data.columns])
        else:
            data = data.select([(pl.col(c).abs() > 0).cast(pl.Int32) for c in data.columns])

        n_sets = data.width
        if n_sets > 3:
            raise ValueError("Can't create Venn diagram for more than 3 sets")

        counts = []
        for combo in product([0, 1], repeat=n_sets):
            mask = pl.lit(True)
            for i, val in enumerate(combo):
                mask = mask & (pl.col(data.columns[i]) == val)
            counts.append(list(combo) + [data.filter(mask).height])
        return pl.DataFrame(counts, schema=data.columns + ["Counts"], orient="row")

    @staticmethod
    def plot_venn2(
        counts: pl.DataFrame,
        names: List[str],
        output_file: Optional[Path] = None,
        figsize: Tuple[int, int] = (8, 8),
    ) -> None:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_aspect("equal"); ax.axis("off")
        ax.add_patch(Circle((-1, 0), 1.5, fill=False, edgecolor="blue", linewidth=2))
        ax.add_patch(Circle((1, 0), 1.5, fill=False, edgecolor="red", linewidth=2))
        ax.text(-1.5, 2, names[0], fontsize=14, ha="center", weight="bold")
        ax.text(1.5, 2, names[1], fontsize=14, ha="center", weight="bold")
        c = counts["Counts"].to_list()
        ax.text(-1.5, 0, str(c[2]), fontsize=12, ha="center")   # A only  (1,0)
        ax.text(1.5, 0, str(c[1]), fontsize=12, ha="center")    # B only  (0,1)
        ax.text(0, 0, str(c[3]), fontsize=12, ha="center")      # both    (1,1)
        ax.text(0, -2.5, str(c[0]), fontsize=10, ha="center")   # neither (0,0)
        plt.tight_layout()
        _save_or_show(output_file, "Venn diagram")

    @staticmethod
    def plot_venn3(
        counts: pl.DataFrame,
        names: List[str],
        output_file: Optional[Path] = None,
        figsize: Tuple[int, int] = (10, 10),
    ) -> None:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_aspect("equal"); ax.axis("off")
        r, s3 = 1.5, np.sqrt(3)
        ax.add_patch(Circle((-1, 1 / s3), r, fill=False, edgecolor="blue", linewidth=2))
        ax.add_patch(Circle((1, 1 / s3), r, fill=False, edgecolor="red", linewidth=2))
        ax.add_patch(Circle((0, -2 / s3), r, fill=False, edgecolor="green", linewidth=2))
        ax.text(-1.5, 2.5, names[0], fontsize=14, ha="center", weight="bold")
        ax.text(1.5, 2.5, names[1], fontsize=14, ha="center", weight="bold")
        ax.text(0, -3.2, names[2], fontsize=14, ha="center", weight="bold")
        c = counts["Counts"].to_list()
        for x, y, v in [
            (0, -3, c[0]), (-1.5, 1, c[4]), (1.5, 1, c[2]), (0, -1.8, c[1]),
            (-0.7, -0.3, c[5]), (0.7, -0.3, c[3]), (0, 1, c[6]), (0, 0, c[7]),
        ]:
            ax.text(x, y, str(v), fontsize=11, ha="center")
        plt.tight_layout()
        _save_or_show(output_file, "Venn diagram")


class HeatmapPlotter:
    @staticmethod
    def plot_heatmap(
        data: Union[pl.DataFrame, np.ndarray],
        row_labels: Optional[List[str]] = None,
        col_labels: Optional[List[str]] = None,
        cmap: str = "RdYlGn_r",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        center: Optional[float] = 0,
        figsize: Tuple[int, int] = (10, 8),
        output_file: Optional[Path] = None,
        title: Optional[str] = None,
        **kwargs,
    ) -> None:
        if isinstance(data, pl.DataFrame):
            if col_labels is None:
                col_labels = data.columns
            arr = data.to_numpy()
        else:
            arr = data

        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            arr, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, center=center,
            xticklabels=col_labels if col_labels else False,
            yticklabels=row_labels if row_labels else False,
            cbar_kws={"label": "Value"}, **kwargs,
        )
        if title:
            ax.set_title(title, fontsize=14, weight="bold")
        plt.tight_layout()
        _save_or_show(output_file, "Heatmap")

    @staticmethod
    def plot_clustered_heatmap(
        data: Union[pl.DataFrame, np.ndarray],
        row_labels: Optional[List[str]] = None,
        col_labels: Optional[List[str]] = None,
        gene_col: str = "gene_id",
        cmap: str = "RdYlGn_r",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        figsize: Tuple[int, int] = (12, 10),
        output_file: Optional[Path] = None,
        title: Optional[str] = None,
        **kwargs,
    ) -> None:
        import pandas as pd

        if isinstance(data, pl.DataFrame):
            if gene_col in data.columns and row_labels is None:
                row_labels = data[gene_col].to_list()
                data = data.drop(gene_col)
            df = data.to_pandas()
        else:
            df = pd.DataFrame(data)

        if row_labels is not None:
            df.index = row_labels
        if col_labels is not None:
            df.columns = col_labels

        g = sns.clustermap(
            df, cmap=cmap, vmin=vmin, vmax=vmax, figsize=figsize,
            cbar_kws={"label": "Value"}, **kwargs,
        )
        if title:
            g.fig.suptitle(title, fontsize=14, weight="bold", y=0.98)
        if output_file is not None:
            g.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info("Clustered heatmap saved to %s", output_file)
        else:
            plt.show()
        plt.close()


def _save_or_show(output_file: Optional[Path], label: str) -> None:
    if output_file is not None:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info("%s saved to %s", label, output_file)
    else:
        plt.show()
    plt.close()


# --------------------------------------------------------------------------- #
# Enrichment-oriented plots (added)
# --------------------------------------------------------------------------- #
class EnrichmentPlots:
    """Publication-style views of GAGE results."""

    @staticmethod
    def bubble_plot(
        results: pl.DataFrame,
        top_n: int = 20,
        stat_col: str = "stat_mean",
        q_col: str = "q_val",
        size_col: str = "set_size",
        name_col: str = "gene_set",
        title: str = "Enriched gene sets",
        output_file: Optional[Path] = None,
    ) -> None:
        """Dot plot: x = stat.mean, colour = -log10(q), size = set_size.

        Accepts a Polars, pandas, or Arrow results table (coerced internally).
        """
        df = _as_polars_frame(results).sort(q_col).head(top_n)
        stat = df[stat_col].to_numpy().astype(float)
        q = df[q_col].to_numpy().astype(float)
        size = df[size_col].to_numpy().astype(float)
        labels = [s[:48] for s in df[name_col].to_list()]
        nlq = -np.log10(np.clip(q, 1e-300, 1.0))
        y = np.arange(len(labels))[::-1]

        fig, ax = plt.subplots(figsize=(9, max(3, 0.42 * len(labels))))
        sizes = 40 + 260 * (size - size.min()) / (np.ptp(size) + 1e-9)
        sc = ax.scatter(stat, y, s=sizes, c=nlq, cmap="viridis", edgecolor="k", linewidth=0.4)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color="grey", lw=0.8, ls="--")
        ax.set_xlabel("stat.mean (mean per-sample statistic)")
        ax.set_title(title, weight="bold")
        cb = fig.colorbar(sc, ax=ax, pad=0.01)
        cb.set_label(r"$-\log_{10}(q)$")
        # size legend
        for sv in (size.min(), np.median(size), size.max()):
            ax.scatter([], [], s=40 + 260 * (sv - size.min()) / (np.ptp(size) + 1e-9),
                       c="grey", edgecolor="k", label=f"{int(sv)}")
        ax.legend(title="set size", loc="lower right", fontsize=7, framealpha=0.8)
        plt.tight_layout()
        _save_or_show(output_file, "Bubble plot")

    @staticmethod
    def enrichment_heatmap(
        results_by_condition: dict,
        value: str = "stat_mean",
        top_n: int = 25,
        q_col: str = "q_val",
        name_col: str = "gene_set",
        output_file: Optional[Path] = None,
        title: str = "Gene-set enrichment across conditions",
    ) -> None:
        """Heatmap of a statistic per gene set (rows) x condition (cols).

        ``results_by_condition`` maps condition name -> GAGE greater/less frame.
        Rows are the union of the top-N sets per condition.
        """
        import pandas as pd
        frames = {cond: _as_polars_frame(df) for cond, df in results_by_condition.items()}
        top_sets = []
        for df in frames.values():
            top_sets += df.sort(q_col).head(top_n)[name_col].to_list()
        rows = list(dict.fromkeys(top_sets))
        data = {}
        for cond, df in frames.items():
            lut = dict(zip(df[name_col].to_list(), df[value].to_numpy()))
            data[cond] = [lut.get(r, np.nan) for r in rows]
        pdf = pd.DataFrame(data, index=[r[:48] for r in rows])
        fig, ax = plt.subplots(figsize=(1.6 + 1.1 * len(data), max(4, 0.4 * len(rows))))
        sns.heatmap(pdf, cmap="RdBu_r", center=0, linewidths=0.4,
                    cbar_kws={"label": value}, ax=ax)
        ax.set_title(title, weight="bold")
        plt.tight_layout()
        _save_or_show(output_file, "Enrichment heatmap")

    @staticmethod
    def running_enrichment(
        ranked: pl.DataFrame,
        gene_set: List[str],
        gene_col: str = "gene_id",
        score_col: str = "score",
        weight: float = 1.0,
        title: Optional[str] = None,
        output_file: Optional[Path] = None,
    ) -> dict:
        """GSEA-style weighted running-enrichment plot for one gene set.

        Ranks genes by ``score_col`` (descending), walks the list accumulating a
        weighted hit / miss running sum, and returns the enrichment score (ES)
        and leading-edge genes.  Complements pygage's KS mode (gs.KSTest).
        """
        df = ranked.select([gene_col, score_col]).drop_nulls().sort(score_col, descending=True)
        genes = df[gene_col].to_list()
        scores = np.abs(df[score_col].to_numpy().astype(float)) ** weight
        in_set = np.array([g in set(gene_set) for g in genes])
        n = len(genes)
        n_hit = in_set.sum()
        if n_hit == 0:
            raise ValueError("no gene-set members found in the ranked list")
        n_miss = n - n_hit
        hit_inc = np.where(in_set, scores / scores[in_set].sum(), 0.0)
        miss_dec = np.where(in_set, 0.0, 1.0 / n_miss)
        running = np.cumsum(hit_inc - miss_dec)
        es_idx = int(np.argmax(np.abs(running)))
        es = float(running[es_idx])
        leading = [genes[i] for i in range(es_idx + 1) if in_set[i]] if es >= 0 \
            else [genes[i] for i in range(es_idx, n) if in_set[i]]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 5), height_ratios=[3, 1], sharex=True)
        ax1.plot(running, color="green", lw=1.5)
        ax1.axhline(0, color="grey", lw=0.8)
        ax1.axvline(es_idx, color="red", ls="--", lw=1, label=f"ES = {es:.3f}")
        ax1.set_ylabel("running enrichment score")
        ax1.legend(loc="upper right", fontsize=8)
        ax1.set_title(title or "Running enrichment", weight="bold")
        ax2.vlines(np.where(in_set)[0], 0, 1, color="black", lw=0.5)
        ax2.set_yticks([])
        ax2.set_xlabel("gene rank")
        plt.tight_layout()
        _save_or_show(output_file, "Running enrichment")
        return {"ES": es, "leading_edge": leading, "n_hits": int(n_hit)}

    @staticmethod
    def pathway_gene_colors(
        pathway_genes: List[str],
        fold_changes: dict,
        title: str = "Pathway member fold changes",
        vmax: Optional[float] = None,
        output_file: Optional[Path] = None,
    ) -> dict:
        """Pathview-style colouring of a pathway's member genes by fold change.

        Produces a colour per gene (RdBu_r on log2FC) and a compact grid figure.
        Returns the gene->hex-colour map so the same colours can drive a KGML
        overlay in your existing Pathview/SBGNview pipeline.
        """
        import matplotlib
        genes = [g for g in pathway_genes]
        vals = np.array([fold_changes.get(g, np.nan) for g in genes], dtype=float)
        finite = vals[np.isfinite(vals)]
        vm = vmax or (np.max(np.abs(finite)) if finite.size else 1.0)
        norm = matplotlib.colors.Normalize(-vm, vm)
        cmap = matplotlib.cm.get_cmap("RdBu_r")
        colors = {g: matplotlib.colors.to_hex(cmap(norm(v)) if np.isfinite(v) else (0.85, 0.85, 0.85))
                  for g, v in zip(genes, vals)}

        ncol = int(np.ceil(np.sqrt(len(genes)))) or 1
        nrow = int(np.ceil(len(genes) / ncol))
        fig, ax = plt.subplots(figsize=(0.9 * ncol + 1, 0.6 * nrow + 1))
        for i, g in enumerate(genes):
            r, c = divmod(i, ncol)
            ax.add_patch(plt.Rectangle((c, nrow - r - 1), 0.95, 0.95,
                         facecolor=colors[g], edgecolor="k", lw=0.4))
            ax.text(c + 0.475, nrow - r - 0.525, str(g)[:8], ha="center", va="center", fontsize=6)
        ax.set_xlim(0, ncol); ax.set_ylim(0, nrow); ax.axis("off")
        ax.set_title(title, weight="bold")
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
        fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="log2 fold change")
        plt.tight_layout()
        _save_or_show(output_file, "Pathway gene colours")
        return colors
