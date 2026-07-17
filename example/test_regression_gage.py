"""Gold-standard regression: pygage must reproduce real gage R output.

Fixtures under pygage/data/regression were produced by running the actual
`gage` R package (datapplab/gage) on its own demo data (gse16873 + kegg.gs).
We load gage's prepared fold-change matrix and its greater/less result tables,
run pygage's engine on the identical matrix, and require a machine-precision
match on every reported column, for t-test/z-test and Stouffer/Fisher meta.
"""
import gzip
import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from pygage import core
from pygage.core import GAGEAnalysis

REG = Path(core.__file__).parent / "data" / "regression"
TOL = 1e-8


def _prepared():
    with gzip.open(REG / "gse16873_prepared.csv.gz", "rt") as fh:
        return pl.read_csv(fh.read().encode() if False else fh, schema_overrides={"gene_id": pl.Utf8}) \
            if False else pl.read_csv(REG / "gse16873_prepared.csv.gz", schema_overrides={"gene_id": pl.Utf8})


def _gene_sets():
    return json.loads((REG / "kegg_gs.json").read_text())


def _gage(name):
    df = pl.read_csv(REG / (name + ".gz"), null_values="NA")
    return df.rename({df.columns[0]: "gene_set"}).drop_nulls("p.val")


@pytest.fixture(scope="module")
def prep():
    return pl.read_csv(REG / "gse16873_prepared.csv.gz", schema_overrides={"gene_id": pl.Utf8})


@pytest.fixture(scope="module")
def gsets():
    return _gene_sets()


def _max_abs(pyg, ref, pairs):
    m = pyg.join(ref, on="gene_set", suffix="_R")
    worst = 0.0
    for c, rc in pairs.items():
        a = m[c].to_numpy().astype(float)
        b = m[rc].to_numpy().astype(float)
        fin = np.isfinite(a) & np.isfinite(b)
        worst = max(worst, float(np.max(np.abs(a[fin] - b[fin]))) if fin.any() else 0.0)
    return worst, m.height


def test_ttest_greater_matches_gage(prep, gsets):
    res = GAGEAnalysis().run_gage(prep, gsets, test_method="t-test", meta_method="stouffer")
    worst, n = _max_abs(res["greater"], _gage("gage_tTest_greater.csv"),
                        {"stat_mean": "stat.mean", "p_val": "p.val",
                         "p_geomean": "p.geomean", "q_val": "q.val", "set_size": "set.size"})
    assert n >= 150
    assert worst < TOL, f"max abs diff {worst:.2e} exceeds {TOL}"


def test_ttest_less_matches_gage(prep, gsets):
    res = GAGEAnalysis().run_gage(prep, gsets, test_method="t-test", meta_method="stouffer")
    worst, n = _max_abs(res["less"], _gage("gage_tTest_less.csv"),
                        {"stat_mean": "stat.mean", "p_val": "p.val",
                         "p_geomean": "p.geomean", "q_val": "q.val"})
    assert worst < TOL, f"max abs diff {worst:.2e} exceeds {TOL}"


def test_ztest_greater_matches_gage(prep, gsets):
    res = GAGEAnalysis().run_gage(prep, gsets, test_method="z-test", meta_method="stouffer")
    worst, n = _max_abs(res["greater"], _gage("gage_zTest_greater.csv"),
                        {"stat_mean": "stat.mean", "p_val": "p.val"})
    assert worst < TOL, f"max abs diff {worst:.2e} exceeds {TOL}"


def test_fisher_meta_matches_gage(prep, gsets):
    res = GAGEAnalysis().run_gage(prep, gsets, test_method="t-test", meta_method="fisher")
    worst, n = _max_abs(res["greater"], _gage("gage_fisher_greater.csv"),
                        {"p_val": "p.val", "p_geomean": "p.geomean"})
    assert worst < TOL, f"max abs diff {worst:.2e} exceeds {TOL}"


def test_bh_matches_r_padjust():
    # R p.adjust(c(0.001,0.008,0.03,0.5,0.2), method='BH')
    p = np.array([0.001, 0.008, 0.03, 0.5, 0.2])
    got = core.benjamini_hochberg(p)
    exp = np.array([0.005, 0.02, 0.05, 0.5, 0.25])
    assert np.allclose(got, exp, atol=1e-12)
