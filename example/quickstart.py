#!/usr/bin/env python3
"""PyGAGE quickstart: raw matrix, DE table, gene-set sourcing, and a plot.

Run from the repo root:  python examples/quickstart.py
Uses the shipped gse16873 regression fixture (real gage demo data).
"""
import gzip
import json
from pathlib import Path

import polars as pl

from pygage import core, gage, GAGEAnalysis, GAGEPreparation, read_de_table, load_gmt
from pygage.visualization_utils import EnrichmentPlots

REG = Path(core.__file__).parent / "data" / "regression"


def main() -> None:
    prepared = pl.read_csv(REG / "gse16873_prepared.csv.gz", schema_overrides={"gene_id": pl.Utf8})
    gene_sets = json.loads((REG / "kegg_gs.json").read_text())

    # 1) one-call GAGE on a prepared fold-change matrix -> tidy frame
    tidy = gage(prepared, gene_sets, prepared=True)
    print("Top enriched (greater):")
    print(tidy.filter(pl.col("direction") == "greater").sort("p_val").head(5))

    # 2) explicit two-step API with extras (effect size + leading-edge genes)
    res = GAGEAnalysis().run_gage(
        prepared, gene_sets, test_method="t-test", meta_method="stouffer",
        compute_effect=True, leading_edge=True,
    )
    print("\nWith effect size + leading edge:")
    print(res["greater"].select(["gene_set", "stat_mean", "p_val", "q_val", "effect"]).head(3))

    # 3) DE-table path (DESeq2/edgeR/limma) — auto-detect columns
    #    de = read_de_table("deseq2_results.csv", value="log2FC")
    #    res = gage(de, gene_sets)

    # 4) MSigDB / GMT sourcing (offline, versioned)
    #    coll = load_gmt("h.all.v2023.2.Hs.symbols.gmt", source="MSigDB", release="2023.2")
    #    res = gage(matrix, coll.gene_sets, ref_indices=[...], samp_indices=[...])

    # 5) a bubble plot of the enriched sets
    EnrichmentPlots.bubble_plot(res["greater"], top_n=15, output_file="enriched_bubble.png")
    print("\nWrote enriched_bubble.png")


if __name__ == "__main__":
    main()
