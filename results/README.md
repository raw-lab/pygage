# pygage ↔ gage R parity — reproducibility bundle

Runs the real GAGE R package and pygage on the same data and shows they agree
to machine precision. See METHODS.md for the full writeup, formulas, and the
gage-R↔pygage code mapping.

## Quickstart

```bash
# 1. R + gage (source or Bioconductor) and Python deps
bash 01_setup.sh

# 2. run the real gage R package on its demo data -> ./gage_out/
Rscript 02_run_gage_reference.R

# 3. run pygage on the identical input and diff every column
#    point PYTHONPATH at your pygage checkout (or `pip install -e .` it first)
PYTHONPATH=/path/to/pygage-1.2.0 python3 03_compare.py
```

Expect `RESULT: ALL PASS`, max |Δ| ~1e-15 per column, and `pygage_vs_gageR.png`.

## What it proves

pygage is fed gage's *exact* prepared fold-change matrix, so the only thing
being compared is the gene-set test + meta step. Under the default t-test +
Stouffer, and also under z-test and Fisher meta, both directions, pygage
reproduces gage's `stat.mean`, `p.val`, `p.geomean`, and `q.val` to
floating-point noise. It is the same statistic, not an approximation.
