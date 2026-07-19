#!/usr/bin/env python3
"""Per-test gene-set statistics, delegating to the gage-faithful core engine.

``GeneSetTests`` is a thin, lower-level view over ``core.GAGEAnalysis`` so the
three saaTests (t / z / KS) return **exactly** the numbers gage R produces
(verified to ~1e-14 in tests/test_regression_gage.py).  Historically this module
carried its own t-test whose background variance term used the set size for one
term and the background size for the other; that has been replaced by
delegation, because gage's definition uses the *set* size for the background
term (``b = s/length(ix)`` in gs.tTest.R) -- and matching gage is the point.

Output schema per set: gene_set, set_size, statistic (= mean per-sample stat),
p_greater, p_less.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import polars as pl

from .core import GAGEAnalysis

_SCHEMA = {
    "gene_set": pl.Utf8,
    "set_size": pl.Int64,
    "statistic": pl.Float64,
    "p_greater": pl.Float64,
    "p_less": pl.Float64,
}


def _empty() -> pl.DataFrame:
    return pl.DataFrame(schema=_SCHEMA)


def _reshape(res: Dict[str, pl.DataFrame], same_dir: bool) -> pl.DataFrame:
    """Map core's greater/less/stats output to the GeneSetTests schema."""
    g = res["greater"].select(["gene_set", "set_size", "stat_mean", "p_val"]).rename(
        {"stat_mean": "statistic", "p_val": "p_greater"}
    )
    if same_dir and "less" in res:
        l = res["less"].select(["gene_set", "p_val"]).rename({"p_val": "p_less"})
        g = g.join(l, on="gene_set", how="left")
    else:
        g = g.with_columns(pl.lit(None, dtype=pl.Float64).alias("p_less"))
    return g.select(list(_SCHEMA.keys())).cast(_SCHEMA)


class GeneSetTests:
    """Set-vs-array-background tests (t / z / KS), numerically identical to gage."""

    @staticmethod
    def t_test(
        expression_data: pl.DataFrame,
        gene_sets: Dict[str, List[str]],
        gene_col: str = "gene_id",
        set_size_range: Tuple[int, int] = (10, 500),
        same_dir: bool = True,
    ) -> Dict[str, object]:
        res = GAGEAnalysis().run_gage(
            expression_data, gene_sets, gene_col=gene_col,
            set_size_range=set_size_range, same_dir=same_dir,
            test_method="t-test", compute_effect=False,
        )
        out = _reshape(res, same_dir) if res["greater"].height else _empty()
        return {"results": out, "method": "t-test"}

    @staticmethod
    def z_test(
        expression_data: pl.DataFrame,
        gene_sets: Dict[str, List[str]],
        gene_col: str = "gene_id",
        set_size_range: Tuple[int, int] = (10, 500),
        same_dir: bool = True,
    ) -> Dict[str, object]:
        res = GAGEAnalysis().run_gage(
            expression_data, gene_sets, gene_col=gene_col,
            set_size_range=set_size_range, same_dir=same_dir,
            test_method="z-test", compute_effect=False,
        )
        out = _reshape(res, same_dir) if res["greater"].height else _empty()
        return {"results": out, "method": "z-test"}

    @staticmethod
    def kolmogorov_smirnov_test(
        expression_data: pl.DataFrame,
        gene_sets: Dict[str, List[str]],
        gene_col: str = "gene_id",
        set_size_range: Tuple[int, int] = (10, 500),
        same_dir: bool = True,
    ) -> Dict[str, object]:
        res = GAGEAnalysis().run_gage(
            expression_data, gene_sets, gene_col=gene_col,
            set_size_range=set_size_range, same_dir=same_dir,
            test_method="ks-test", compute_effect=False,
        )
        out = _reshape(res, same_dir) if res["greater"].height else _empty()
        return {"results": out, "method": "ks-test"}
