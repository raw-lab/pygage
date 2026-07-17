#!/usr/bin/env python3
"""Compare pygage to gage R on the IDENTICAL prepared fold-change matrix.

Reads gage's outputs from ./gage_out (produced by 02_run_gage_reference.R),
runs pygage's engine on the same prepared matrix + gene sets, prints a
head-to-head table and the max deviation per column, writes a parity plot,
and exits non-zero if any column exceeds the tolerance.

Usage:  python3 03_compare.py [--gage-out gage_out] [--tol 1e-8]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pygage.core import GAGEAnalysis   # pip install -e . from the pygage repo


def load_gage(path: Path) -> pl.DataFrame:
    # gage writes NA for size-filtered sets; gene_ids can be mixed Entrez/AFFX strings
    df = pl.read_csv(path, null_values="NA")
    return df.rename({df.columns[0]: "gene_set"}).drop_nulls("p.val")


def max_abs(m: pl.DataFrame, a: str, b: str) -> float:
    x = m[a].to_numpy().astype(float); y = m[b].to_numpy().astype(float)
    fin = np.isfinite(x) & np.isfinite(y)
    return float(np.nanmax(np.abs(x[fin] - y[fin]))) if fin.any() else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gage-out", default="gage_out")
    ap.add_argument("--tol", type=float, default=1e-8)
    args = ap.parse_args()
    out = Path(args.gage_out)

    prep = pl.read_csv(out / "prepared_matrix.csv", schema_overrides={"gene_id": pl.Utf8})
    gs = json.loads((out / "kegg_gs.json").read_text())

    # ---- run pygage on the identical input ----
    ga = GAGEAnalysis()
    py_t = ga.run_gage(prep, gs, test_method="t-test", meta_method="stouffer")
    py_z = GAGEAnalysis().run_gage(prep, gs, test_method="z-test", meta_method="stouffer")
    py_f = GAGEAnalysis().run_gage(prep, gs, test_method="t-test", meta_method="fisher")

    checks = {
        "t-test greater": (py_t["greater"], load_gage(out / "gage_tTest_greater.csv"),
                           {"stat_mean": "stat.mean", "p_val": "p.val",
                            "p_geomean": "p.geomean", "q_val": "q.val", "set_size": "set.size"}),
        "t-test less":    (py_t["less"], load_gage(out / "gage_tTest_less.csv"),
                           {"stat_mean": "stat.mean", "p_val": "p.val",
                            "p_geomean": "p.geomean", "q_val": "q.val"}),
        "z-test greater": (py_z["greater"], load_gage(out / "gage_zTest_greater.csv"),
                           {"stat_mean": "stat.mean", "p_val": "p.val"}),
        "fisher greater": (py_f["greater"], load_gage(out / "gage_fisher_greater.csv"),
                           {"p_val": "p.val", "p_geomean": "p.geomean"}),
    }

    ok = True
    for label, (pyg, gage, cols) in checks.items():
        m = pyg.join(gage, on="gene_set", suffix="_R")
        worst = max(max_abs(m, a, b) for a, b in cols.items())
        status = "PASS" if worst < args.tol else "FAIL"
        ok = ok and worst < args.tol
        print(f"[{status}] {label:16s} matched {m.height:3d} sets   max|Δ| = {worst:.2e}")
        for a, b in cols.items():
            print(f"           {a:11s} vs {b:11s}  max|Δ| = {max_abs(m, a, b):.2e}")

    # ---- head-to-head table (top 12 by pygage p) for the default t-test ----
    m = py_t["greater"].join(load_gage(out / "gage_tTest_greater.csv"), on="gene_set", suffix="_R").sort("p_val")
    print(f"\nTop 12 KEGG sets (default t-test + Stouffer):\n")
    print(f"{'pathway':44s} {'pygage p':>11s} {'gage p':>11s} {'pygage stat':>11s} {'gage stat':>10s}")
    print("-" * 92)
    for r in m.head(12).iter_rows(named=True):
        print(f"{r['gene_set'][:44]:44s} {r['p_val']:11.3e} {r['p.val']:11.3e} {r['stat_mean']:11.4f} {r['stat.mean']:10.4f}")

    # ---- parity plot ----
    pg = -np.log10(np.clip(m["p_val"].to_numpy(), 1e-300, 1))
    gg = -np.log10(np.clip(m["p.val"].to_numpy(), 1e-300, 1))
    sa, sb = m["stat_mean"].to_numpy(), m["stat.mean"].to_numpy()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2))
    ax1.scatter(gg, pg, s=28, alpha=0.7, edgecolor="k", linewidth=0.3, color="#2c7fb8")
    lim = [0, max(pg.max(), gg.max()) * 1.05]; ax1.plot(lim, lim, "r--", label="y = x")
    ax1.set_xlabel(r"gage R  $-\log_{10}(p)$"); ax1.set_ylabel(r"pygage  $-\log_{10}(p)$")
    ax1.set_title("p-value parity"); ax1.legend(); ax1.set_xlim(lim); ax1.set_ylim(lim)
    ax2.scatter(sb, sa, s=28, alpha=0.7, edgecolor="k", linewidth=0.3, color="#31a354")
    sl = [min(sa.min(), sb.min()) - 0.2, max(sa.max(), sb.max()) + 0.2]; ax2.plot(sl, sl, "r--", label="y = x")
    ax2.set_xlabel("gage R  stat.mean"); ax2.set_ylabel("pygage  stat.mean")
    ax2.set_title("statistic parity"); ax2.legend(); ax2.set_xlim(sl); ax2.set_ylim(sl)
    fig.suptitle("pygage vs. gage R — identical input", weight="bold")
    plt.tight_layout(); fig.savefig("pygage_vs_gageR.png", dpi=150, bbox_inches="tight")
    print(f"\nPearson r: -log10(p) {np.corrcoef(pg, gg)[0,1]:.8f} | stat.mean {np.corrcoef(sa, sb)[0,1]:.8f}")
    print("wrote pygage_vs_gageR.png")

    print("\nRESULT:", "ALL PASS ✅" if ok else "FAILURES ❌")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
