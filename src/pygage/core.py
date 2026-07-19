"""GAGE core: a faithful Python port of the gage R stat engine.

Ported directly from datapplab/gage (gagePrep.R, gs.tTest.R, gs.zTest.R,
gs.KSTest.R, gageSum.R).  The per-gene-set statistic and the meta
summarisation reproduce gage's numbers (see tests/regression)::

    gs.tTest per column j, set S (size n, present genes):
        a  = var(S_j)/n ;  b = var(all_j)/n        # note: BOTH over n (gage)
        df = (a+b)^2 / (a^2/(n-1) + b^2/(n-1))
        stat = (mean(S_j) - mean(all_j)) * (a+b)^(-1/2)
        p_up = P(T_df > stat) ;  p_down = P(T_df < stat)

    meta (default Stouffer): p.val = Phi( sum_j qnorm(p_j) / sqrt(nc) )
         (Fisher/gamma alt) : p.val = pgamma(sum_j -log(p_j), shape=nc, lower=F)
    p.geomean = exp( -sum_j -log(p_j) / nc )
    stat.mean = mean_j stat_j ;  q.val = BH(p.val) within each direction
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, Optional, Sequence, Tuple

if TYPE_CHECKING:
    from ._types import GeneSetsLike

import numpy as np
import polars as pl
from scipy import stats

logger = logging.getLogger("pygage.core")

VALID_TESTS = ("t-test", "z-test", "ks-test")
VALID_META = ("stouffer", "fisher")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """BH FDR, matching R's p.adjust(method='BH'). NaNs are preserved."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    m = int(mask.sum())
    if m == 0:
        return out
    pv = p[mask]
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    adj = np.empty(m)
    adj[order] = q
    out[mask] = adj
    return out


@dataclass
class GAGEResult:
    """Typed GAGE result (tidy access without leaving polars behind)."""

    greater: pl.DataFrame
    less: Optional[pl.DataFrame]
    stats: pl.DataFrame
    meta: Dict[str, object] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, pl.DataFrame]:
        d = {"greater": self.greater, "stats": self.stats}
        if self.less is not None:
            d["less"] = self.less
        return d

    def significant(self, cutoff: float = 0.1, use_q: bool = True) -> Dict[str, pl.DataFrame]:
        col = ("q_val" if use_q else "p_val")
        out = {"greater": self.greater.filter(pl.col(col) < cutoff)}
        if self.less is not None:
            out["less"] = self.less.filter(pl.col(col) < cutoff)
        return out


# --------------------------------------------------------------------------- #
# preparation (gagePrep)
# --------------------------------------------------------------------------- #
class GAGEPreparation:
    """Per-gene fold-change / statistic preparation (gagePrep port)."""

    @staticmethod
    def prepare_expression(
        expression_data: pl.DataFrame,
        ref_indices: Optional[Sequence[int]] = None,
        samp_indices: Optional[Sequence[int]] = None,
        gene_col: str = "gene_id",
        comparison: str = "paired",
        same_dir: bool = True,
        use_fold: bool = True,
        input_logged: bool = True,
        rank_test: bool = False,
    ) -> pl.DataFrame:
        """Form the per-gene fold-change / statistic matrix GAGE tests on.

        ``ref_indices``/``samp_indices`` are 0-based positions among the value
        columns (i.e. excluding ``gene_col``), matching gage's ``ref``/``samp``.

        Parameters
        ----------
        comparison:
            How reference and sample columns are contrasted:

            - ``"paired"`` (default): element-wise ``sample - ref`` per pair;
              ``len(samp)`` must be a multiple of ``len(ref)``. Use when samples
              and references are matched (e.g. tumour/normal from one patient).
            - ``"unpaired"``: every sample-vs-every-reference difference (an
              ``n_samp x n_ref`` fan-out); use for unmatched group designs.
            - ``"as.group"``: a single column, ``mean(samples) - mean(refs)``.
            - ``"1ongroup"``: each sample minus the mean of the references
              (one-vs-group).
        same_dir:
            If ``False``, take ``|fold change|`` so up- and down-regulation both
            count as "changed" (pairs with ``run_gage(same_dir=False)``).
        use_fold:
            Must be ``True`` here; to run on precomputed per-gene statistics
            (e.g. moderated t), pass them directly to ``run_gage`` instead.
        input_logged:
            Set ``False`` to ``log2(x + 1)`` raw (non-log) inputs first.
        rank_test:
            Rank-transform each column (for the KS ``saaTest``).

        Returns
        -------
        polars.DataFrame
            ``gene_col`` plus one fold-change column per resulting comparison.
        """
        if gene_col not in expression_data.columns:
            raise ValueError(
                f"gene_col {gene_col!r} not in columns {expression_data.columns}. "
                "Provide the gene identifier column name via gene_col=."
            )
        value_cols = [c for c in expression_data.columns if c != gene_col]
        genes = expression_data[gene_col]
        mat = expression_data.select(value_cols).to_numpy().astype(float)
        if not input_logged:
            mat = np.log2(mat + 1.0)

        if ref_indices is None:
            fc = mat  # already fold-changes / statistics
            names = value_cols
        else:
            ref = list(ref_indices)
            samp = list(samp_indices) if samp_indices is not None else [
                i for i in range(mat.shape[1]) if i not in ref
            ]
            if not use_fold:
                raise ValueError("use_fold=False (t-stat prep) unsupported; pass prepared stats instead")
            if comparison == "as.group":
                fc = (mat[:, samp].mean(axis=1) - mat[:, ref].mean(axis=1)).reshape(-1, 1)
                names = ["mean.fc"]
            elif comparison == "paired":
                if len(samp) % len(ref) != 0:
                    raise ValueError("paired: len(samp) must be a multiple of len(ref)")
                fc = mat[:, samp] - mat[:, ref]
                names = [value_cols[j] for j in samp]
            elif comparison == "unpaired":
                cols = []
                names = []
                for sj in samp:
                    for rj in ref:
                        cols.append(mat[:, sj] - mat[:, rj])
                        names.append(f"{value_cols[sj]}_vs_{value_cols[rj]}")
                fc = np.column_stack(cols)
            elif comparison in ("1ongroup", "as.ref"):
                ref_mean = mat[:, ref].mean(axis=1)
                fc = mat[:, samp] - ref_mean[:, None]
                names = [value_cols[j] for j in samp]
            else:
                raise ValueError(f"comparison must be paired/unpaired/as.group/1ongroup, got {comparison!r}")

        if not same_dir:
            fc = np.abs(fc)
        if rank_test:
            fc = np.apply_along_axis(lambda c: stats.rankdata(c), 0, fc)

        out = {gene_col: genes}
        for k, name in enumerate(names):
            out[name] = fc[:, k]
        return pl.DataFrame(out)


# --------------------------------------------------------------------------- #
# analysis (gs.tTest / gs.zTest / gs.KSTest + gageSum)
# --------------------------------------------------------------------------- #
class GAGEAnalysis:
    def __init__(self) -> None:
        self.results: Optional[Dict[str, pl.DataFrame]] = None
        self.result_obj: Optional[GAGEResult] = None

    # ---- per-column gene-set statistics (returns stat, p_up, p_down) ---- #
    @staticmethod
    def _stats_ttest(set_mat, mu, s, n):
        a = np.nanvar(set_mat, axis=0, ddof=1) / n
        b = s / n
        with np.errstate(divide="ignore", invalid="ignore"):
            df = (a + b) ** 2 / (a ** 2 / (n - 1) + b ** 2 / (n - 1))
            mod = (a + b) ** (-0.5)
        stat = (np.nanmean(set_mat, axis=0) - mu) * mod
        p_up = stats.t.sf(stat, df)
        p_down = stats.t.cdf(stat, df)
        return stat, p_up, p_down

    @staticmethod
    def _stats_ztest(set_mat, mu, s, n):
        with np.errstate(divide="ignore", invalid="ignore"):
            mod = np.sqrt(n / s)
        stat = (np.nanmean(set_mat, axis=0) - mu) * mod
        return stat, stats.norm.sf(stat), stats.norm.cdf(stat)

    @staticmethod
    def _stats_kstest(set_rank_cols, all_rank_cols):
        nc = len(all_rank_cols)
        stat = np.empty(nc)
        p_up = np.empty(nc)
        p_down = np.empty(nc)
        for j in range(nc):
            x = set_rank_cols[j]
            comp = np.setdiff1d(all_rank_cols[j], x, assume_unique=False)
            r_less = stats.ks_2samp(x, comp, alternative="less")
            r_greater = stats.ks_2samp(x, comp, alternative="greater")
            p_up[j] = r_less.pvalue
            p_down[j] = r_greater.pvalue
            stat[j] = max(r_less.statistic, r_greater.statistic)
        return stat, p_up, p_down

    # ---- meta summarisation (gageSum) ---- #
    @staticmethod
    def _meta_pval(p_cols: np.ndarray, method: str) -> float:
        p = np.clip(p_cols[~np.isnan(p_cols)], 1e-300, 1.0)
        if p.size == 0:
            return np.nan
        nc = p.size
        if method == "stouffer":
            return float(stats.norm.cdf(np.sum(stats.norm.ppf(p)) / np.sqrt(nc)))
        # fisher / gamma
        sg = np.sum(-np.log(p))
        return float(stats.gamma.sf(sg, a=nc, scale=1.0))

    @staticmethod
    def _p_geomean(p_cols: np.ndarray) -> float:
        p = np.clip(p_cols[~np.isnan(p_cols)], 1e-300, 1.0)
        if p.size == 0:
            return np.nan
        return float(np.exp(-np.sum(-np.log(p)) / p.size))

    def run_gage(
        self,
        expression_data: pl.DataFrame,
        gene_sets: "GeneSetsLike",
        gene_col: str = "gene_id",
        set_size_range: Tuple[int, int] = (10, 500),
        same_dir: bool = True,
        test_method: str = "t-test",
        meta_method: str = "stouffer",
        fdr_method: str = "BH",
        control_genes: Optional[Sequence[str]] = None,
        global_bh: bool = False,
        compute_effect: bool = True,
        leading_edge: bool = False,
        permutations: int = 0,
        n_jobs: int = 1,
        random_state: int = 0,
    ) -> Dict[str, pl.DataFrame]:
        """Run GAGE and return per-direction result tables.

        The default configuration (``test_method="t-test"``, ``meta_method="stouffer"``,
        ``fdr_method="BH"``) reproduces the GAGE R package to machine precision.

        Parameters
        ----------
        expression_data:
            A prepared per-gene fold-change / statistic frame (``gene_col`` plus
            one column per sample). Produce it with
            :meth:`GAGEPreparation.prepare_expression`, or pass a single-column
            DE / pre-ranked frame for one-sample analysis.
        gene_sets:
            Gene sets as a mapping ``{name: [genes]}`` or a
            :class:`pygage.gene_sets.GeneSetCollection`; either is normalised to
            the canonical schema internally.
        gene_col:
            Name of the gene-identifier column (default ``"gene_id"``).
        set_size_range:
            ``(min, max)`` number of set genes that must be *present* in the data
            for a set to be tested; others are dropped (GAGE default ``(10, 500)``).
        same_dir:
            If ``True`` (default) return separate ``greater`` (up) and ``less``
            (down) tables. If ``False``, test the magnitude of change only
            (``|fold change|``) and return a single directionless ranking.
        test_method:
            Per-sample statistic: ``"t-test"`` (gage ``gs.tTest``, default),
            ``"z-test"`` (PAGE-style ``gs.zTest``), or ``"ks-test"`` (rank-based
            ``gs.KSTest``; parity is algorithmic — see the docs).
        meta_method:
            Cross-sample p-value combination: ``"stouffer"`` (gage default) or
            ``"fisher"`` (gamma/Fisher).
        fdr_method:
            ``"BH"`` for Benjamini–Hochberg q-values (default), or any other
            value to report the raw p-value as ``q_val``.
        control_genes:
            Optional gene set used as the background reference; if omitted the
            background is all genes (standard GAGE).
        global_bh:
            If ``True`` (and ``same_dir``), apply BH across the union of the
            greater and less p-values instead of per direction.
        compute_effect:
            Add an ``effect`` column (mean fold change across the set).
        leading_edge:
            Add a ``leading_edge`` column listing the member genes driving the
            signal.
        permutations:
            If ``> 0``, add a ``p_perm`` sample-label permutation p-value
            (non-parametric alternative to the analytic meta-test).
        n_jobs:
            Worker count for the per-set loop (``-1`` = all cores). See
            :func:`pygage.config.default_n_jobs` and
            :func:`pygage.config.set_thread_limits` for HPC scaling.
        random_state:
            Seed for the permutation null.

        Returns
        -------
        dict of str to polars.DataFrame
            ``{"greater": ..., "less": ..., "stats": ...}`` (``less`` only when
            ``same_dir``). Each result row has ``gene_set``, ``set_size``,
            ``stat_mean``, ``p_geomean``, ``p_val``, ``q_val`` (plus any optional
            columns). See :class:`GAGEResult` for a typed view via ``result_obj``.
        """
        from .gene_sets import normalize_gene_sets

        if test_method not in VALID_TESTS:
            raise ValueError(f"test_method must be one of {VALID_TESTS}, got {test_method!r}")
        if meta_method not in VALID_META:
            raise ValueError(f"meta_method must be one of {VALID_META}, got {meta_method!r}")
        gene_sets = normalize_gene_sets(gene_sets)
        if gene_col not in expression_data.columns:
            raise ValueError(f"gene_col {gene_col!r} not found. Did prepare_expression drop it?")
        value_cols = [c for c in expression_data.columns if c != gene_col]
        if not value_cols:
            raise ValueError("No fold-change/statistic columns found besides gene_col")

        genes = expression_data[gene_col].to_list()
        gene_pos = {g: i for i, g in enumerate(genes)}
        mat = expression_data.select(value_cols).to_numpy().astype(float)
        nc = mat.shape[1]

        # background moments (over all genes, or a control set)
        if control_genes is not None:
            cidx = [gene_pos[g] for g in set(control_genes) if g in gene_pos]
            if len(cidx) < 2:
                raise ValueError("control_genes: <2 present in data")
            bg = mat[cidx, :]
        else:
            bg = mat
        mu = np.nanmean(bg, axis=0)
        s = np.nanvar(bg, axis=0, ddof=1)

        rank_cols = None
        if test_method == "ks-test":
            rank_cols = [stats.rankdata(mat[:, j]) for j in range(nc)]

        def eval_set(item):
            name, glist = item
            idx = [gene_pos[g] for g in dict.fromkeys(glist) if g in gene_pos]  # ordered-unique
            present = len(idx)
            if present < set_size_range[0] or present > set_size_range[1]:
                return None
            if test_method == "t-test":
                stat, p_up, p_down = self._stats_ttest(mat[idx, :], mu, s, present)
            elif test_method == "z-test":
                stat, p_up, p_down = self._stats_ztest(mat[idx, :], mu, s, present)
            else:
                set_rc = [rank_cols[j][idx] for j in range(nc)]
                stat, p_up, p_down = self._stats_kstest(set_rc, rank_cols)
            row = {
                "gene_set": name,
                "set_size": present,
                "stat_mean": float(np.nanmean(stat)),
                "p_geomean": self._p_geomean(p_up),
                "p_geomean_down": self._p_geomean(p_down) if same_dir else np.nan,
                "p_greater": self._meta_pval(p_up, meta_method),
                "p_less": self._meta_pval(p_down, meta_method) if same_dir else np.nan,
            }
            if compute_effect:
                row["effect"] = float(np.nanmean(mat[idx, :]))  # mean fold-change
            if leading_edge:
                gene_mean = np.nanmean(mat[idx, :], axis=1)
                thr = np.nanmean(mat[idx, :])
                lead = [genes[idx[k]] for k in np.argsort(-gene_mean)
                        if gene_mean[k] > thr][:25]
                row["leading_edge"] = ";".join(lead)
            if permutations > 0:
                row["_stat_mean_obs"] = row["stat_mean"]
                row["_idx"] = idx
            return row

        items = list(gene_sets.items())
        if n_jobs and n_jobs != 1:
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=(None if n_jobs < 0 else n_jobs)) as ex:
                rows = [r for r in ex.map(eval_set, items) if r is not None]
        else:
            rows = [r for r in map(eval_set, items) if r is not None]

        # optional sample-label permutation null (on stat.mean)
        if permutations > 0 and rows:
            rng = np.random.default_rng(random_state)
            perm_stats = {r["gene_set"]: np.empty(permutations) for r in rows}
            for p in range(permutations):
                perm = mat[:, rng.permutation(nc)]
                pmu = np.nanmean(perm if control_genes is None else perm[cidx, :], axis=0)
                ps = np.nanvar(perm if control_genes is None else perm[cidx, :], axis=0, ddof=1)
                for r in rows:
                    idx = r["_idx"]
                    if test_method == "t-test":
                        st, _, _ = self._stats_ttest(perm[idx, :], pmu, ps, len(idx))
                    elif test_method == "z-test":
                        st, _, _ = self._stats_ztest(perm[idx, :], pmu, ps, len(idx))
                    else:
                        prc = [stats.rankdata(perm[:, j]) for j in range(nc)]
                        st, _, _ = self._stats_kstest([prc[j][idx] for j in range(nc)], prc)
                    perm_stats[r["gene_set"]][p] = np.nanmean(st)
            for r in rows:
                obs = r.pop("_stat_mean_obs")
                r.pop("_idx")
                null = perm_stats[r["gene_set"]]
                r["p_perm"] = float((np.sum(null >= obs) + 1) / (permutations + 1))

        if not rows:
            final = {"gene_set": pl.Utf8, "set_size": pl.Int64, "stat_mean": pl.Float64,
                     "p_geomean": pl.Float64, "p_val": pl.Float64, "q_val": pl.Float64}
            empty = pl.DataFrame(schema=final)
            out = {"greater": empty, "stats": pl.DataFrame(schema={"gene_set": pl.Utf8, "stat_mean": pl.Float64, "set_size": pl.Int64})}
            if same_dir:
                out["less"] = empty
            self.results = out
            self.result_obj = GAGEResult(greater=out["greater"], less=out.get("less"),
                                         stats=out["stats"], meta={"test_method": test_method,
                                         "meta_method": meta_method, "n_samples": nc, "same_dir": same_dir})
            return out
        df = pl.DataFrame(rows, schema=None)

        # FDR: BH within each direction (gage), or global across the union
        if df.height and fdr_method.upper() == "BH":
            if global_bh and same_dir:
                allp = np.concatenate([df["p_greater"].to_numpy(), df["p_less"].to_numpy()])
                allq = benjamini_hochberg(allp)
                df = df.with_columns(
                    pl.Series("q_greater", allq[: df.height]),
                    pl.Series("q_less", allq[df.height:]),
                )
            else:
                df = df.with_columns(pl.Series("q_greater", benjamini_hochberg(df["p_greater"].to_numpy())))
                if same_dir:
                    df = df.with_columns(pl.Series("q_less", benjamini_hochberg(df["p_less"].to_numpy())))
        else:
            df = df.with_columns(pl.col("p_greater").alias("q_greater"))
            if same_dir:
                df = df.with_columns(pl.col("p_less").alias("q_less"))

        extra = [c for c in ("effect", "leading_edge", "p_perm") if c in df.columns]
        gcols = ["gene_set", "set_size", "stat_mean", "p_geomean", "p_greater", "q_greater"] + extra
        greater = df.select([c for c in gcols if c in df.columns]).rename(
            {"p_greater": "p_val", "q_greater": "q_val"}
        ).sort("p_val")
        out: Dict[str, pl.DataFrame] = {
            "greater": greater,
            "stats": df.select(["gene_set", "stat_mean", "set_size"]),
        }
        if same_dir:
            lcols = ["gene_set", "set_size", "stat_mean", "p_geomean_down", "p_less", "q_less"] + extra
            out["less"] = df.select([c for c in lcols if c in df.columns]).rename(
                {"p_geomean_down": "p_geomean", "p_less": "p_val", "q_less": "q_val"}
            ).sort("p_val")

        self.results = out
        self.result_obj = GAGEResult(
            greater=out["greater"], less=out.get("less"), stats=out["stats"],
            meta={"test_method": test_method, "meta_method": meta_method,
                  "n_samples": nc, "same_dir": same_dir},
        )
        return out

    def filter_significant(self, cutoff: float = 0.1, use_q: bool = True) -> Dict[str, pl.DataFrame]:
        if self.results is None:
            raise ValueError("No results available. Run run_gage() first.")
        col = "q_val" if use_q else "p_val"
        out = {}
        for key in ("greater", "less"):
            if key in self.results and col in self.results[key].columns:
                out[key] = self.results[key].filter(pl.col(col) < cutoff)
        if "stats" in self.results:
            out["stats"] = self.results["stats"]
        return out
