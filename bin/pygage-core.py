#!/usr/bin/env python3
"""Command-line interface for GAGE analysis (corrected)."""

import argparse
import json
from pathlib import Path

import polars as pl

from pygage.core import GAGEPreparation, GAGEAnalysis


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GAGE (Generally Applicable Gene set Enrichment) analysis"
    )
    parser.add_argument("--expression", type=Path, required=True,
                        help="Expression data file (CSV/TSV)")
    parser.add_argument("--gene-sets", type=Path, required=True,
                        help="Gene sets JSON file (dict or {'gene_sets': {...}})")
    parser.add_argument("--gene-col", default="gene_id", help="Gene ID column name")
    parser.add_argument("--ref-indices", nargs="+", type=int,
                        help="Reference column indices (0-based, excluding gene col)")
    parser.add_argument("--samp-indices", nargs="+", type=int,
                        help="Sample column indices (0-based, excluding gene col)")
    parser.add_argument("--comparison",
                        choices=["paired", "unpaired", "1ongroup", "as.group"],
                        default="paired")
    parser.add_argument("--test-method",
                        choices=["t-test", "z-test"], default="t-test")
    parser.add_argument("--set-size-min", type=int, default=10)
    parser.add_argument("--set-size-max", type=int, default=500)
    parser.add_argument("--same-dir", action="store_true", default=True)
    parser.add_argument("--not-log", action="store_true",
                        help="Input is raw (not log-scaled); log2-transform first")
    parser.add_argument("--cutoff", type=float, default=0.1,
                        help="Q-value cutoff for significance")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output directory for results")
    args = parser.parse_args()

    sep = "," if args.expression.suffix == ".csv" else "\t"
    expr = pl.read_csv(args.expression, separator=sep)
    print(f"Loaded expression data: {expr.shape}")

    payload = json.loads(args.gene_sets.read_text())
    gene_sets = payload.get("gene_sets", payload)
    print(f"Loaded {len(gene_sets)} gene sets")

    prepared = GAGEPreparation.prepare_expression(
        expr,
        ref_indices=args.ref_indices,
        samp_indices=args.samp_indices,
        gene_col=args.gene_col,
        comparison=args.comparison,
        same_dir=args.same_dir,
        input_logged=not args.not_log,
    )

    gage = GAGEAnalysis()
    gage.run_gage(
        prepared,
        gene_sets,
        gene_col=args.gene_col,
        set_size_range=(args.set_size_min, args.set_size_max),
        same_dir=args.same_dir,
        test_method=args.test_method,
    )
    significant = gage.filter_significant(cutoff=args.cutoff)

    args.output.mkdir(parents=True, exist_ok=True)
    for key, df in significant.items():
        out = args.output / f"{key}.tsv"
        df.write_csv(out, separator="\t")
        print(f"Wrote {df.height} rows to {out}")

    print("\nSummary:")
    print(f"  Up-regulated gene sets:   {significant['greater'].height}")
    if "less" in significant:
        print(f"  Down-regulated gene sets: {significant['less'].height}")


if __name__ == "__main__":
    main()
