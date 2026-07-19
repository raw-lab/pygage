#!/usr/bin/env python3
"""
GAGE Results Analysis and Comparison (corrected).

Fixes vs. the original draft
----------------------------
* ``compare_results`` used ``how='outer'`` (deprecated -> removed in polars);
  now uses ``how='full'`` consistently.
* ``SignificanceFilter.filter_significant`` (dual_sig == 1) tried to vstack the
  greater/less frames whose p/q columns have different names, raising a
  ShapeError.  Columns are now aligned to a common ``score`` before concat.
* ``group_gene_sets`` restricts the hypergeometric universe and set sizes to
  genes actually measured, so the overlap test can't receive n > N.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import polars as pl

from .visualization_utils import VennDiagram

logger = logging.getLogger("pygage.results")


def esset_grp(
    results: pl.DataFrame,
    expression_data: pl.DataFrame,
    gene_sets: Dict[str, List[str]],
    gene_col: str = "gene_id",
    test4up: bool = True,
    same_dir: bool = True,
    cutoff: float = 0.01,
    use_q: bool = False,
    pc: float = 1e-10,
) -> Dict[str, object]:
    """Faithful port of gage's ``esset.grp``: collapse redundant gene sets.

    Redundancy is defined on **core genes** (set genes whose mean fold-change
    exceeds one SD of the gene-mean distribution in the direction of interest),
    tested for overlap against the pool of **essential genes** (mean fc beyond
    2 SD) via the hypergeometric ``phyper`` upper tail.  Sets whose overlap
    p-value < ``pc`` are edges; connected components are the redundant groups.

    Returns ``groups`` (representative -> members), ``core_genes`` per set, and
    ``essential_genes``.  Mirrors gage's use of the p-value (not q-value) by
    default for stable non-redundant selection.
    """
    from scipy.stats import hypergeom

    pcol = "p_val" if not use_q else "q_val"
    if pcol not in results.columns:
        pcol = next((c for c in ("p_val", "p_greater", "q_val", "q_greater") if c in results.columns), None)
        if pcol is None:
            raise ValueError("results needs a p_val/q_val column")
    ssp = results.filter(pl.col(pcol).is_not_null() & (pl.col(pcol) < cutoff))
    if ssp.height < 2:
        return {"groups": {}, "core_genes": {}, "essential_genes": [],
                "note": "fewer than 2 significant gene sets"}

    genes = expression_data[gene_col].to_list()
    gpos = {g: i for i, g in enumerate(genes)}
    mat = expression_data.select([c for c in expression_data.columns if c != gene_col]).to_numpy().astype(float)
    gene_means = np.nanmean(mat, axis=1)              # per-gene mean fold change
    s = float(np.nanstd(gene_means, ddof=0))
    m = float(np.nanmean(gene_means))
    sign = test4up

    names = ssp["gene_set"].to_list()
    core: Dict[str, List[str]] = {}
    for nm in names:
        idx = [gpos[g] for g in gene_sets.get(nm, []) if g in gpos]
        b = gene_means[idx] - m
        order = np.argsort(-b if (sign or not same_dir) else b)
        keep = [genes[idx[k]] for k in order
                if (b[k] > s if (sign or not same_dir) else b[k] < -s)]
        core[nm] = keep
    len_core = {nm: len(v) for nm, v in core.items()}

    ball = gene_means - m
    ess_mask = (ball > 2 * s) if (sign or not same_dir) else (ball < -2 * s)
    ess_genes = [genes[i] for i in np.where(ess_mask)[0]]
    n_ess = len(ess_genes)

    k = len(names)
    adj = np.zeros((k, k), dtype=int)
    for i in range(k):
        ci = set(core[names[i]])
        for j in range(i + 1, k):
            overlap = len(ci & set(core[names[j]]))
            if overlap < 1:
                continue
            # gage: phyper(overlap, len.core[i], n_ess-len.core[i], len.core[j],
            # lower.tail=FALSE) == P(X > overlap) == scipy hypergeom.sf(overlap, ...)
            p_ov = hypergeom.sf(overlap, n_ess, len_core[names[i]], len_core[names[j]])
            if p_ov < pc:
                adj[i, j] = adj[j, i] = 1

    visited, groups = set(), []
    for start in range(k):
        if start in visited:
            continue
        stack, comp = [start], []
        while stack:
            nd = stack.pop()
            if nd in visited:
                continue
            visited.add(nd)
            comp.append(nd)
            stack.extend(nb for nb in range(k) if adj[nd, nb] and nb not in visited)
        groups.append(comp)

    # representative = most significant (first in ssp order); members labelled
    group_dict = {}
    for g in groups:
        member_names = [names[x] for x in g]
        rep = member_names[0]
        group_dict[rep] = member_names
    return {"groups": group_dict, "core_genes": core,
            "essential_genes": ess_genes, "n_essential": n_ess}


class ResultsComparator:
    """Compare enrichment results across multiple datasets."""

    @staticmethod
    def _read(path) -> pl.DataFrame:
        path = Path(path)
        sep = "," if path.suffix == ".csv" else "\t"
        return pl.read_csv(path, separator=sep)

    @staticmethod
    def compare_results(
        result_files: List[Path],
        sample_names: List[str],
        q_cutoff: float = 0.1,
        output_file: Optional[Path] = None,
    ) -> pl.DataFrame:
        if len(result_files) != len(sample_names):
            raise ValueError("Number of files must match number of sample names")

        all_results = []
        for path, name in zip(result_files, sample_names):
            df = ResultsComparator._read(path)
            rename: Dict[str, str] = {}
            if "stat_mean" in df.columns:
                rename["stat_mean"] = f"{name}_stat"
            if "q_val" in df.columns:
                rename["q_val"] = f"{name}_q"
            elif "q_greater" in df.columns:
                rename["q_greater"] = f"{name}_q"
            elif "p_greater" in df.columns:
                rename["p_greater"] = f"{name}_p"
            all_results.append(df.rename(rename))

        combined = all_results[0]
        gene_set_col = "gene_set" if "gene_set" in combined.columns else combined.columns[0]
        for df in all_results[1:]:
            combined = combined.join(
                df, on=gene_set_col, how="full", coalesce=True
            ).fill_null(0)
        combined = combined.sort(gene_set_col)

        sig_cols = [c for c in combined.columns if c.endswith("_q") or c.endswith("_p")]
        if sig_cols:
            hits = pl.lit(0)
            for c in sig_cols:
                hits = hits + (pl.col(c) < q_cutoff).cast(pl.Int32)
            combined = combined.with_columns(hits.alias("hits")).sort(
                "hits", descending=True
            )

        if output_file is not None:
            output_file = Path(output_file)
            sep = "," if output_file.suffix == ".csv" else "\t"
            combined.write_csv(output_file, separator=sep)
            logger.info("Combined results written to %s", output_file)
        return combined

    @staticmethod
    def create_venn_comparison(
        result_files: List[Path],
        sample_names: List[str],
        q_cutoff: float = 0.1,
        output_file: Optional[Path] = None,
    ) -> None:
        if len(result_files) > 3:
            raise ValueError("Venn diagrams support 2-3 comparisons only")

        sig_frames: Dict[str, pl.DataFrame] = {}
        for path, name in zip(result_files, sample_names):
            df = ResultsComparator._read(path)
            q_col = next(
                (c for c in ("q_val", "q_greater", "p_greater") if c in df.columns),
                None,
            )
            if q_col is None:
                raise ValueError(f"No q/p-value column found in {path}")
            gcol = "gene_set" if "gene_set" in df.columns else df.columns[0]
            sig_frames[name] = df.select(
                [pl.col(gcol).alias("gene_set"),
                 (pl.col(q_col) < q_cutoff).cast(pl.Int32).alias(name)]
            )

        combined = sig_frames[sample_names[0]]
        for name in sample_names[1:]:
            combined = combined.join(
                sig_frames[name], on="gene_set", how="full", coalesce=True
            ).fill_null(0)

        venn = VennDiagram()
        counts = venn.venn_counts(combined.select(sample_names), include="both")
        if len(sample_names) == 2:
            venn.plot_venn2(counts, sample_names, output_file)
        else:
            venn.plot_venn3(counts, sample_names, output_file)


class GeneSetGrouper:
    """Group overlapping significant gene sets by shared membership."""

    @staticmethod
    def group_gene_sets(
        results: pl.DataFrame,
        gene_sets: Dict[str, List[str]],
        expression_data: pl.DataFrame,
        gene_col: str = "gene_id",
        p_cutoff: float = 0.01,
        overlap_cutoff: float = 1e-10,
        output_file: Optional[Path] = None,
    ) -> Dict[str, List[str]]:
        from scipy.stats import hypergeom

        p_col = next(
            (c for c in ("p_val", "p_greater", "q_greater") if c in results.columns),
            None,
        )
        if p_col is None:
            raise ValueError("results needs a p_val/p_greater/q_greater column")
        sig = results.filter(pl.col(p_col) < p_cutoff)
        if sig.height < 2:
            logger.info("Fewer than 2 significant gene sets; nothing to group")
            return {}

        names = sig["gene_set"].to_list()
        universe = set(expression_data[gene_col].to_list())
        N = len(universe)
        members = {n: set(gene_sets.get(n, [])) & universe for n in names}

        k = len(names)
        adjacency = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(i + 1, k):
                a, b = members[names[i]], members[names[j]]
                ov = len(a & b)
                if ov > 0 and a and b:
                    p_ov = hypergeom.sf(ov - 1, N, len(a), len(b))
                    if p_ov < overlap_cutoff:
                        adjacency[i, j] = adjacency[j, i] = 1

        visited = set()
        groups: List[List[int]] = []
        for start in range(k):
            if start in visited:
                continue
            stack, comp = [start], []
            while stack:                       # iterative DFS (no recursion limit)
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.append(node)
                stack.extend(
                    nb for nb in range(k) if adjacency[node, nb] and nb not in visited
                )
            groups.append(comp)

        group_dict = {f"Group_{i+1}": [names[x] for x in g] for i, g in enumerate(groups)}
        if output_file is not None:
            Path(output_file).write_text(json.dumps(group_dict, indent=2))
            logger.info("Gene set groups written to %s", output_file)
        return group_dict


class SignificanceFilter:
    """Filter significant gene sets from GAGE results."""

    @staticmethod
    def filter_significant(
        results: Dict[str, pl.DataFrame],
        cutoff: float = 0.1,
        use_q: bool = True,
        dual_sig: int = 2,
    ) -> Dict[str, pl.DataFrame]:
        """Filter results; ``dual_sig`` controls two-direction handling.

        0 = exclusive (drop sets significant in both directions)
        1 = keep the better-scoring direction per gene set
        2 = keep both directions (default)
        """
        def pick(df: pl.DataFrame, direction: str) -> str:
            for cand in (
                ("q_val" if use_q else "p_val"),
                (f"q_{direction}" if use_q else f"p_{direction}"),
            ):
                if cand in df.columns:
                    return cand
            raise ValueError(f"No score column in {direction} results")

        filtered: Dict[str, pl.DataFrame] = {}
        has_less = "less" in results and results["less"] is not None

        gcol = pick(results["greater"], "greater")
        greater_sig = results["greater"].filter(pl.col(gcol) < cutoff)

        if not has_less:
            filtered["greater"] = greater_sig
            if "stats" in results:
                filtered["stats"] = results["stats"]
            return filtered

        lcol = pick(results["less"], "less")
        less_sig = results["less"].filter(pl.col(lcol) < cutoff)

        if dual_sig == 0:
            filtered["greater"] = greater_sig.filter(
                ~pl.col("gene_set").is_in(less_sig["gene_set"])
            )
            filtered["less"] = less_sig.filter(
                ~pl.col("gene_set").is_in(greater_sig["gene_set"])
            )
        elif dual_sig == 1:
            # Align to a common 'score' column, then keep the best per gene set.
            g = greater_sig.select(
                pl.col("gene_set"),
                pl.col(gcol).alias("score"),
                pl.lit("greater").alias("direction"),
            )
            l = less_sig.select(
                pl.col("gene_set"),
                pl.col(lcol).alias("score"),
                pl.lit("less").alias("direction"),
            )
            best = pl.concat([g, l]).sort("score").unique(
                subset=["gene_set"], keep="first", maintain_order=True
            )
            keep_g = set(best.filter(pl.col("direction") == "greater")["gene_set"].to_list())
            keep_l = set(best.filter(pl.col("direction") == "less")["gene_set"].to_list())
            filtered["greater"] = greater_sig.filter(pl.col("gene_set").is_in(keep_g))
            filtered["less"] = less_sig.filter(pl.col("gene_set").is_in(keep_l))
        else:
            filtered["greater"] = greater_sig
            filtered["less"] = less_sig

        if "stats" in results:
            filtered["stats"] = results["stats"]
        return filtered
