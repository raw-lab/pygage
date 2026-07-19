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

## Reproducing the pygage ↔ gage R parity comparison

This bundle reproduces the claim that **pygage's statistic engine reproduces the
GAGE R package to machine precision** on real data. It runs the actual `gage`
R package and pygage on the *same* input and diffs every reported column.

## TL;DR

```bash
bash 01_setup.sh                       # R + gage, Python deps
Rscript 02_run_gage_reference.R        # runs real gage R -> ./gage_out/*.csv
PYTHONPATH=/path/to/pygage python3 03_compare.py   # runs pygage, diffs, plots
```

Expected: every column matches within ~1e-15 (`RESULT: ALL PASS`), a parity
plot `pygage_vs_gageR.png`, and Pearson r = 1.00000000.

## Environment used

- R version 4.3.3 (2024-02-29) (installed via `apt-get install r-base-core` from the Ubuntu archive)
- Python 3.12.3 — polars 1.42.1, numpy 2.4.4, scipy 1.17.1, matplotlib
- gage source: `github.com/datapplab/gage` (the maintainer's repo; ships the demo
  data as `.rda`). The **core** gage path used here — `gagePrep` → `gs.tTest` /
  `gs.zTest` → `gageSum` — depends only on base R + `stats`, so no Bioconductor
  install is required to run it. `02_run_gage_reference.R` prefers the installed
  `gage` Bioconductor package if present, and otherwise `source()`s the repo's
  `.R` files (that is what was done in the sandbox, where Bioconductor was
  unreachable).

## Data

- `gse16873` — the GAGE demo expression set: 11,979 genes × 12 samples
  (6 head-and-neck "HN" and 6 ductal carcinoma in situ "DCIS"), already
  log2-transformed. Shipped inside the gage package.
- `kegg.gs` — 177 human KEGG signaling/metabolism gene sets (Entrez IDs).
  Shipped inside the gage package.
- Comparison: **6 HN vs 6 DCIS, paired** (the standard vignette call
  `gage(gse16873, gsets=kegg.gs, ref=hn, samp=dcis)`).

## Comparison design (why it is a fair, tight test)

The single biggest source of spurious "differences" between two GSEA
implementations is the **preprocessing** (how per-gene fold changes are formed).
To isolate the *statistic engine*, we do not let pygage do its own prep: we
export gage's **exact prepared fold-change matrix** (`gagePrep` output,
11,979 × 6) and feed that identical matrix to pygage. Any remaining difference
is therefore purely the gene-set test + meta-summarisation, not prep. We also
export gage's gene sets to JSON so set memberships are byte-identical.

We compare four configurations:

| configuration | gage R call | pygage call |
|---|---|---|
| t-test + Stouffer (default) | `gage(...)` | `run_gage(test_method="t-test", meta_method="stouffer")` |
| t-test, `less` direction | `res$less` | `run_gage(...)["less"]` |
| z-test (PAGE-style) | `saaTest=gs.zTest` | `test_method="z-test"` |
| Fisher/gamma meta | `use.stouffer=FALSE` | `meta_method="fisher"` |

For each we join on gene-set name and take `max |Δ|` over
`stat.mean`, `p.val`, `p.geomean`, `q.val`, and `set.size`.

## The algorithm, and how pygage maps to it

Ported 1:1 from the gage R sources (file names below are in `gage_R/R/`).
pygage code is in `pygage/core.py`.

### 1. Preparation — `gagePrep.R` → `GAGEPreparation.prepare_expression`
Paired fold changes: `exprs[,samp] - exprs[,ref]` per pair (also supports
`unpaired`, `as.group`, `1ongroup`; `abs()` when `same.dir=FALSE`; optional
rank transform for KS).

### 2. Per-sample gene-set statistic — `gs.tTest.R` → `GAGEAnalysis._stats_ttest`
For column *j* and set *S* of *n* present genes, with background mean `mu_j` and
variance `s_j` over **all** genes:

```
a  = var(S_j) / n            # set variance / set size
b  = s_j       / n           # BACKGROUND variance / SET size   <-- gage's definition
df = (a + b)^2 / ( a^2/(n-1) + b^2/(n-1) )
stat_j = ( mean(S_j) - mu_j ) * (a + b)^(-1/2)
p_up_j   = P(T_df > stat_j)        # pt(stat, df, lower.tail=FALSE)
p_down_j = P(T_df < stat_j)
```

**Key subtlety (worth knowing if you audit this):** gage divides *both* variance
terms by the **set** size `n` (`b = s/length(ix)` in `gs.tTest.R`), and uses
Welch df on `n-1`. That is GAGE's statistic by definition — it is *not* a
textbook two-sample Welch test (which would divide the background term by the
background gene count). pygage matches gage exactly.

z-test (`gs.zTest.R` → `_stats_ztest`): `stat_j = (mean(S_j) - mu_j) * sqrt(n / s_j)`, p via the normal.

KS-test (`gs.KSTest.R` → `_stats_kstest`): rank genes within each column; per column, KS of set-ranks vs complement-ranks (`alternative="less"` and `"greater"`), statistic = max of the two D's.

### 3. Cross-sample meta-summarisation — `gageSum.R` → `_meta_pval`, `_p_geomean`
Default is **Stouffer's Z** over the `nc` per-sample one-sided p-values:

```
p.val = Phi( sum_j qnorm(p_j) / sqrt(nc) )          # use.stouffer=TRUE (default)
```

Fisher/gamma alternative (`use.stouffer=FALSE`):

```
p.val = pgamma( sum_j -log(p_j), shape=nc, rate=1, lower.tail=FALSE )
```

Other reported columns:

```
p.geomean = exp( -sum_j -log(p_j) / nc )    # geometric mean of per-sample p's; direction-specific
stat.mean = mean_j stat_j
q.val     = p.adjust(p.val, method="BH")    # per direction
```

`p.geomean` is **direction-specific**: gage's `less` table uses the down-tail
p-values, so pygage reports `p_geomean_down` in its `less` table to match.

### 4. BH FDR — `benjamini_hochberg`
Matches R's `p.adjust(method="BH")` (verified separately to 1e-12), applied
within each direction (or optionally across the greater∪less union via
`global_bh=True`, which is a pygage extension, not gage default).

## Results (this bundle, rerun live)

Matched 160 of 177 sets (the other 17 fall outside the 10–500 size window in
both tools). `max |Δ|` across all 160 sets:

| column | max abs Δ |
|---|---|
| stat.mean | 4.9e-15 |
| p.val (Stouffer) | 1.4e-15 |
| p.geomean | 8.9e-16 |
| q.val (BH) | 2.9e-15 |
| set.size | 0 |

z-test p.val 1.9e-15; Fisher p.val 1.8e-15; `less` direction identical. Pearson
r = 1.00000000 on both −log10(p) and stat.mean. These are floating-point-noise
magnitudes, i.e. the same computation, not a re-implementation that "agrees
closely".

## Caveats / honest scope

- **KS mode** parity is *algorithmic*, not machine-precision: R's `ks.test` and
  SciPy's `ks_2samp` differ slightly on tie handling and exact-vs-asymptotic
  p-values. The t-test and z-test are the tightly-validated paths; those are the
  gage defaults.
- **Background** is "all genes" by default in both. pygage adds an optional
  control-gene-set reference (`control_genes=`), which gage does not have.
- **Live KEGG** (`rest.kegg.jp`) was not reachable inside the sandbox, so the
  KEGG *retrieval* code was validated against wire-format mock payloads, not the
  parity comparison above (which uses gage's bundled `kegg.gs`). In your
  environment `pygage kegg` hits the live API.
- Numbers depend only on R/SciPy CDF implementations; different SciPy/R
  point-releases could shift the last digit or two but not the conclusion.

## Files

- `01_setup.sh` — install R + gage, Python deps
- `02_run_gage_reference.R` — run real gage R, export prepared matrix + results + gene sets to `gage_out/`
- `03_compare.py` — run pygage on the identical input, diff every column, plot, exit non-zero on any failure
- `gage_out/` — gage R outputs (created by step 2)
- `pygage_vs_gageR.png` — parity plot (created by step 3)

> Note: the intermediate `gage_out/` tables are **generated** by
> `02_run_gage_reference.R` and are intentionally not committed (see the repo
> `.gitignore`). Run the three steps above to reproduce them and the figure.
