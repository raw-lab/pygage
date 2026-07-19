#!/usr/bin/env python3
"""Unified pygage command-line interface.

Subcommands
-----------
  pygage run      Run GAGE on an expression matrix or DE table + a gene-set file
  pygage kegg     Download KEGG pathway / KO / module gene sets to JSON
  pygage go       Build GO gene sets from a GAF (+ optional OBO)
  pygage compare  Combine several GAGE result tables across conditions

Replaces the previous collection of single-purpose scripts.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path


def _setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )


def _load_gene_sets(path: str) -> dict:
    p = Path(path)
    if p.suffix == ".json" or p.suffix == ".gz":
        from .gene_sets import _open
        with _open(p) as fh:
            obj = json.load(fh)
        return obj["gene_sets"] if isinstance(obj, dict) and "gene_sets" in obj else obj
    from .gene_sets import load_gmt
    return load_gmt(p).gene_sets


# --------------------------------------------------------------------------- #
def cmd_run(args: argparse.Namespace) -> int:
    from .io_loaders import gage, read_de_table, read_matrix
    gene_sets = _load_gene_sets(args.gene_sets)

    if args.de_table:
        frame = read_de_table(args.input, value=args.value)
        res = gage(frame, gene_sets, same_dir=not args.one_direction,
                   test_method=args.test, meta_method=args.meta,
                   set_size_range=(args.min_size, args.max_size), tidy=True)
    else:
        frame = read_matrix(args.input)
        ref = [int(x) for x in args.ref.split(",")] if args.ref else None
        samp = [int(x) for x in args.samp.split(",")] if args.samp else None
        res = gage(frame, gene_sets, ref_indices=ref, samp_indices=samp,
                   comparison=args.compare, same_dir=not args.one_direction,
                   test_method=args.test, meta_method=args.meta,
                   set_size_range=(args.min_size, args.max_size),
                   prepared=args.prepared, tidy=True)

    out = Path(args.output)
    (res.write_csv(out) if out.suffix == ".csv" else res.write_csv(out, separator="\t"))
    sys.stderr.write(f"Wrote {res.height} rows -> {out}\n")
    print(res.head(args.top))
    return 0


def cmd_kegg(args: argparse.Namespace) -> int:
    from .pathway_database_utils import KEGGPathwayRetriever
    r = KEGGPathwayRetriever(cache_dir=args.cache)
    out = Path(args.output)
    if args.kind == "ko":
        r.download_ko_gene_sets(out, reference=args.reference)
    elif args.kind == "pathway":
        res = r.get_pathway_genes(args.species, id_type=args.id_type)
        out.write_text(json.dumps({"gene_sets": res["gene_sets"],
                                   "set_names": res["pathway_names"],
                                   "categories": res["categories"]}))
        sys.stderr.write(f"Wrote {len(res['gene_sets'])} pathway gene sets -> {out}\n")
    else:  # module
        res = r.get_module_gene_sets(args.species)
        out.write_text(json.dumps({"gene_sets": res["gene_sets"],
                                   "set_names": res["module_names"]}))
        sys.stderr.write(f"Wrote {len(res['gene_sets'])} module gene sets -> {out}\n")
    return 0


def cmd_go(args: argparse.Namespace) -> int:
    from .gene_sets import load_go
    coll = load_go(args.gaf, obo_path=args.obo, aspect=args.aspect,
                   id_field=args.id_field, include_iea=not args.no_iea,
                   propagate=args.propagate)
    out = Path(args.output)
    out.write_text(json.dumps({"gene_sets": coll.gene_sets, "metadata": coll.metadata()}))
    sys.stderr.write(f"Wrote {coll.n_sets} GO gene sets ({coll.extra.get('aspect')}) -> {out}\n")
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    from .results_analysis import ResultsComparator
    files = args.results
    names = args.names.split(",") if args.names else [Path(f).stem for f in files]
    combined = ResultsComparator.compare_results(files, names, q_cutoff=args.q_cutoff,
                                                 output_file=args.output)
    sys.stderr.write(f"Compared {len(files)} result files -> {args.output}\n")
    print(combined.head(args.top))
    return 0


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="pygage", description="GAGE gene-set enrichment (Python port of gage R)")
    p.add_argument("-v", "--verbose", action="store_true", help="log progress to stderr")
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run GAGE on a matrix or DE table")
    r.add_argument("input", help="expression matrix (genes x samples) or DE table")
    r.add_argument("-g", "--gene-sets", required=True, help="GMT or JSON gene-set file")
    r.add_argument("-o", "--output", required=True, help="output CSV/TSV")
    r.add_argument("--de-table", action="store_true", help="input is a DESeq2/edgeR/limma table")
    r.add_argument("--value", default="log2FC", choices=["log2FC", "stat"], help="DE column to rank on")
    r.add_argument("--ref", help="0-based reference column indices, comma-separated")
    r.add_argument("--samp", help="0-based sample column indices, comma-separated")
    r.add_argument("--compare", default="paired", choices=["paired", "unpaired", "as.group", "1ongroup"])
    r.add_argument("--prepared", action="store_true", help="input is already a fold-change matrix")
    r.add_argument("--test", default="t-test", choices=["t-test", "z-test", "ks-test"])
    r.add_argument("--meta", default="stouffer", choices=["stouffer", "fisher"])
    r.add_argument("--one-direction", action="store_true", help="same.dir=FALSE (abs changes)")
    r.add_argument("--min-size", type=int, default=10)
    r.add_argument("--max-size", type=int, default=500)
    r.add_argument("--top", type=int, default=15)
    r.set_defaults(func=cmd_run)

    k = sub.add_parser("kegg", help="download KEGG gene sets")
    k.add_argument("kind", choices=["pathway", "ko", "module"])
    k.add_argument("-o", "--output", required=True)
    k.add_argument("-s", "--species", default="hsa", help="KEGG organism code (pathway/module)")
    k.add_argument("--reference", default="pathway", choices=["pathway", "module"], help="KO reference type")
    k.add_argument("--id-type", default="entrez", choices=["kegg", "entrez"])
    k.add_argument("--cache", help="cache directory for KEGG responses")
    k.set_defaults(func=cmd_kegg)

    g = sub.add_parser("go", help="build GO gene sets from a GAF")
    g.add_argument("gaf", help="GAF 2.x annotation file")
    g.add_argument("-o", "--output", required=True)
    g.add_argument("--obo", help="go-basic.obo (names + propagation)")
    g.add_argument("--aspect", choices=["BP", "MF", "CC"], help="restrict to one domain")
    g.add_argument("--id-field", default="symbol", choices=["symbol", "object_id"])
    g.add_argument("--no-iea", action="store_true", help="drop electronic (IEA) annotations")
    g.add_argument("--propagate", action="store_true", help="propagate up the GO DAG (needs --obo)")
    g.set_defaults(func=cmd_go)

    c = sub.add_parser("compare", help="combine GAGE result tables")
    c.add_argument("results", nargs="+", help="GAGE result files (CSV/TSV)")
    c.add_argument("-o", "--output", required=True)
    c.add_argument("--names", help="comma-separated condition names")
    c.add_argument("--q-cutoff", type=float, default=0.1)
    c.add_argument("--top", type=int, default=15)
    c.set_defaults(func=cmd_compare)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    _setup_logging(getattr(args, "verbose", False))
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
