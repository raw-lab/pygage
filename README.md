# PyGAGE — Generally Applicable Gene-set Enrichment in Python

[![Docs](https://img.shields.io/badge/docs-readthedocs-blue)](https://pygage.readthedocs.io)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/pygage?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/pygage)
[![Downloads PePy monthly](https://static.pepy.tech/badge/pygage/month)](https://pepy.tech/project/pygage)
[![Downloads PePy weekly](https://static.pepy.tech/badge/pygage/week)](https://pepy.tech/project/pygage)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org)
[![Validated vs gage R](https://img.shields.io/badge/gage%20R%20parity-~1e--15-brightgreen)](#validation--parity-with-gage-r)

**PyGAGE** is a fast, dependency-light Python implementation of **GAGE**
(*Generally Applicable Gene-set Enrichment*, Luo *et al.* 2009) for pathway
analysis. 

It reproduces the GAGE R package **to machine precision** on real data
(~1e-15 across every reported column, both directions, t-test / z-test / Fisher
meta), and adds first-class support for the inputs people actually have
(DESeq2/edgeR/limma tables, pre-ranked vectors, AnnData), broad gene-set sourcing
(KEGG, KEGG Orthology, GO, Reactome, MSigDB), extra statistical rigor, a unified
command-line interface, and publication-ready plots.

---

## Table of contents

- [Highlights](#highlights)
- [Installation](#installation)
- [60-second quickstart](#60-second-quickstart)
- [How GAGE works](#how-gage-works)
- [Inputs](#inputs)
- [Gene-set sourcing](#gene-set-sourcing)
- [Running GAGE](#running-gage)
- [The `gage()` convenience function](#the-gage-convenience-function)
- [Understanding the results](#understanding-the-results)
- [Collapsing redundant sets (`esset_grp`)](#collapsing-redundant-sets-esset_grp)
- [Visualization](#visualization)
- [Command-line interface](#command-line-interface)
- [Validation / parity with gage R](#validation--parity-with-gage-r)
- [Data format requirements](#data-format-requirements)
- [Performance](#performance)
- [API overview](#api-overview)
- [Migrating from 1.0.0](#migrating-from-100)
- [Citation](#citation)
- [License · Support · Contributing](#license)

---

## Highlights

- **It *is* GAGE.** The two-level design is intact: a per-sample test of each gene
  set against the array-wide background, combined across samples by a meta-test
  (Stouffer's Z by default, Fisher/gamma optional). Rankings and p-values match
  GAGE — validated numerically against the R package (see [below](#validation--parity-with-gage-r)).
- **Inputs people have.** `read_de_table` auto-detects DESeq2 / edgeR / limma
  columns; pre-ranked (fgsea-style) vectors and **pandas / AnnData** are accepted
  alongside polars.
- **Gene sets from everywhere.** KEGG pathways (any organism), **KEGG Orthology**
  (species-agnostic, for metagenomes / viromes / phages), GO (GAF + OBO with DAG
  propagation), Reactome, and MSigDB/GMT — with a versioned offline cache and
  provenance stamped into every collection.
- **More rigor.** Optional control-gene-set background, sample-label permutation
  null, effect size + leading-edge genes, global BH across greater∪less,
  multi-core over gene sets, NaN-robust throughout.
- **One CLI.** A single `pygage` command with `run` / `kegg` / `go` / `compare`
  subcommands replaces the seven per-module scripts.
- **Plots.** Bubble/dot plot, cross-condition enrichment heatmap, GSEA-style
  running-enrichment plot, and pathview-style KEGG member colouring.

---

## Installation

### From PyPI

```bash
pip install pygage
```

### From Bioconda

```bash
conda install -c bioconda pygage
```

### From source

```bash
git clone https://github.com/raw-lab/pygage
cd pygage
pip install .
# with optional extras:
pip install ".[anndata]"     # AnnData input
pip install ".[test]"        # run the test suite (incl. the gage-R regression)
pip install ".[docs]"        # build the documentation
```

**Requirements:** Python ≥ 3.8; `polars`, `numpy`, `scipy`, `matplotlib`,
`seaborn`, `pandas`, `pyarrow`, `requests`. `anndata` is optional.

---

## 60-second quickstart

PyGAGE ships the real GAGE demo data (`gse16873`, the 6 HN vs 6 DCIS breast-cancer
set) so you can run a full analysis with zero downloads:

```python
import json
from pathlib import Path
import polars as pl
from pygage import core, gage

reg = Path(core.__file__).parent / "data" / "regression"
prepared  = pl.read_csv(reg / "gse16873_prepared.csv.gz", schema_overrides={"gene_id": pl.Utf8})
gene_sets = json.loads((reg / "kegg_gs.json").read_text())

# one call -> a tidy result frame (direction-labelled)
result = gage(prepared, gene_sets, prepared=True)
print(result.filter(pl.col("direction") == "greater").sort("p_val").head(5))
```

```
┌───────────────────────────────────────────────┬──────────┬───────────┬────────────┬───────────┬───────────┐
│ gene_set                                       ┆ set_size ┆ stat_mean ┆ p_val      ┆ q_val     ┆ direction │
╞════════════════════════════════════════════════╪══════════╪═══════════╪════════════╪═══════════╪═══════════╡
│ hsa04141 Protein processing in endoplasmic ret ┆ 164      ┆ 3.5171    ┆ 9.2371e-18 ┆ 1.48e-15  ┆ greater   │
│ hsa00190 Oxidative phosphorylation             ┆ 118      ┆ 2.8488    ┆ 3.2788e-12 ┆ 2.62e-10  ┆ greater   │
│ hsa03050 Proteasome                            ┆ 44       ┆ 2.6311    ┆ 2.1089e-10 ┆ 1.12e-08  ┆ greater   │
│ …                                              ┆ …        ┆ …         ┆ …          ┆ …         ┆ …         │
└────────────────────────────────────────────────┴──────────┴───────────┴────────────┴───────────┴───────────┘
```

Those are the published GAGE vignette hits (secretory / degradation machinery up
in DCIS vs HN) — and they match the R package to ~1e-15.

Or from the shell:

```bash
pygage run expression.csv -g gene_sets.json -o results.csv --ref 0,1,2 --samp 3,4,5
```

---

## How GAGE works

GAGE's defining feature is a **two-level** test:

1. **Preparation.** Per-gene fold changes are formed from reference vs sample
   columns (paired / unpaired / as-group / one-on-group).
2. **Per-sample test.** For each sample (column) *j* and gene set *S* of *n*
   present genes, a two-sample-style statistic compares the set to the whole-array
   background, using background mean `mu_j` and variance `s_j`:

   ```
   a  = var(S_j) / n          # set variance / set size
   b  = s_j       / n          # background variance / SET size   (GAGE's definition)
   df = (a + b)^2 / (a^2/(n-1) + b^2/(n-1))
   stat_j   = (mean(S_j) - mu_j) * (a + b)^(-1/2)
   p_up_j   = P(T_df > stat_j) ;  p_down_j = P(T_df < stat_j)
   ```

3. **Cross-sample meta.** The per-sample p-values are combined into one p-value
   per set — **Stouffer's Z** by default:

   ```
   p.val     = Phi( sum_j qnorm(p_j) / sqrt(nc) )       # meta_method="stouffer"
   p.geomean = exp( -sum_j -log(p_j) / nc )
   stat.mean = mean_j stat_j
   q.val     = BH(p.val)                                 # per direction
   ```

   A Fisher/gamma alternative is available (`meta_method="fisher"`).

Separate **greater** and **less** tables report up- and down-regulated sets.
PyGAGE also offers the PAGE-style **z-test** (`test_method="z-test"`) and a
rank-based **KS test** (`test_method="ks-test"`).

> The three `saaTest`s and both meta-combinations are ported 1:1 from the gage R
> sources; the formulas above match `gs.tTest.R` / `gs.zTest.R` / `gageSum.R`
> exactly. See the [validation section](#validation--parity-with-gage-r).

---

## Inputs

### Raw expression matrix (genes × samples)

```python
from pygage import read_matrix, GAGEPreparation, GAGEAnalysis

expr = read_matrix("expression.csv")            # gene_id + sample columns
prepared = GAGEPreparation.prepare_expression(
    expr, ref_indices=[0, 1, 2], samp_indices=[3, 4, 5],
    comparison="paired",       # paired | unpaired | as.group | 1ongroup
    input_logged=True,         # set False to log2(x+1) first
)
result = GAGEAnalysis().run_gage(prepared, gene_sets)
```

### Converting between gene symbols and Entrez IDs
sym2eg: Converts gene symbols to Entrez IDs

eg2sym: Converts Entrez IDs to gene symbols
```python
from pygage.gene_id_utils import GeneIDConverter
from pygage import read_matrix

matrix = read_matrix("expression.csv", gene_col='ID')

symbols = matrix['ID']
eg = GeneIDConverter().sym2eg(symbols)
matrix = matrix.with_columns(pl.Series("ID", eg))
```

### DE tables (DESeq2 / edgeR / limma)

Most users arrive with a differential-expression table, not a raw matrix.
`read_de_table` auto-detects the gene, log2-fold-change, and statistic columns
by common aliases (`log2FoldChange`, `logFC`, `coef`; `stat`, `t`, `LR`; …):

```python
from pygage import read_de_table, gage

de = read_de_table("deseq2_results.csv", value="log2FC")   # or value="stat"
result = gage(de, gene_sets)                                # single-sample ranked GAGE
```

### Pre-ranked vector (fgsea-style)

```python
from pygage import read_preranked, gage

ranked = read_preranked({"TP53": 3.1, "BRCA1": -2.4, "EGFR": 1.8, ...})
result = gage(ranked, gene_sets)
```

### pandas / AnnData

```python
import anndata as ad
from pygage import gage

adata = ad.read_h5ad("counts.h5ad")     # obs × var (samples × genes)
result = gage(adata, gene_sets, ref_indices=[...], samp_indices=[...])
```

AnnData is transposed to genes × samples automatically; pandas DataFrames are
accepted anywhere a polars frame is.

---

## Gene-set sourcing

All loaders live in `pygage.gene_sets` (except KEGG/GO retrieval, which is in
`pygage.pathway_database_utils`) and return either a plain `{name: [genes]}` dict
or a `GeneSetCollection` carrying provenance (source, release, retrieval date,
checksum).

### KEGG pathways (any organism)

```python
from pygage.pathway_database_utils import KEGGPathwayRetriever

kegg = KEGGPathwayRetriever()
res = kegg.get_pathway_genes("hsa", id_type="entrez")   # mmu, eco, ath, dme, sce, ...
gene_sets, names, categories = res["gene_sets"], res["pathway_names"], res["categories"]
```

### KEGG Orthology (species-agnostic — metagenomes, viromes, phages)

For non-model organisms and metagenomic data, use KO-keyed gene sets. Annotate
your genes to KO IDs (KofamScan / MetaCerberus / eggNOG) and run GAGE on the
KO-level matrix:

```python
kegg = KEGGPathwayRetriever()
all_kos = kegg.list_all_kos()                       # the full KO namespace (~26k) {K00001: "..."}
ko_sets = kegg.get_ko_gene_sets(reference="pathway")# {map00010: [K00844, K12407, ...]}

# or download once, reproducibly, to disk (sets + names + categories + full KO catalog + provenance)
kegg.download_ko_gene_sets("ko_gene_sets.json", reference="pathway")
```

### Gene Ontology (GAF + optional OBO, with DAG propagation)

```python
from pygage.gene_sets import load_go

go = load_go(
    "goa_human.gaf",
    obo_path="go-basic.obo",   # enables propagation and term names
    aspect="BP",               # BP | MF | CC | None
    include_iea=True,          # keep electronic annotations
    propagate=True,            # roll annotations up the is_a / part_of DAG
)
print(go.n_sets, "GO sets;", go.metadata())
```

### Reactome

```python
from pygage.gene_sets import load_reactome

# a ReactomePathways.gmt export...
rx = load_reactome("ReactomePathways.gmt", id_type="gmt")
# ...or the NCBI2Reactome mapping, filtered to a species
rx = load_reactome("NCBI2Reactome_All_Levels.txt", id_type="ncbi2reactome",
                   species="Homo sapiens")
```

### MSigDB / any GMT

```python
from pygage.gene_sets import load_gmt, load_msigdb

hallmark = load_gmt("h.all.v2023.2.Hs.symbols.gmt", source="MSigDB", release="2023.2")
c2       = load_msigdb("c2.cp.v2023.2.Hs.symbols.gmt", collection="C2")
```

### Versioned offline cache

```python
from pygage.gene_sets import GeneSetCache

cache = GeneSetCache()                       # ~/.cache/pygage/gene_sets
cache.save("hallmark_2023.2", hallmark)      # gzip on disk
same = cache.load("hallmark_2023.2")         # offline-reproducible; checksum-verified
```

---

## Running GAGE

`GAGEAnalysis.run_gage` is the full engine. Defaults reproduce gage R exactly;
the remaining arguments are opt-in extensions.

```python
from pygage import GAGEAnalysis

ga  = GAGEAnalysis()
res = ga.run_gage(
    prepared, gene_sets,
    gene_col="gene_id",
    set_size_range=(10, 500),   # min/max genes present per set

    # --- statistic & combination (defaults = gage R) ---
    test_method="t-test",       # t-test | z-test | ks-test
    meta_method="stouffer",     # stouffer | fisher
    same_dir=True,              # separate greater/less tables (False = |changes|)
    fdr_method="BH",            # BH q-values (per direction)

    # --- extensions (opt-in) ---
    control_genes=None,         # background from a control gene set instead of all genes
    global_bh=False,            # BH across the greater∪less union
    compute_effect=True,        # add per-set mean fold change ("effect")
    leading_edge=False,         # add member genes driving the signal
    permutations=0,             # sample-label permutation null p-value (e.g. 1000)
    n_jobs=1,                   # parallelise over gene sets (-1 = all cores)
)
greater, less, stats = res["greater"], res["less"], res["stats"]

# filter to significant sets
sig = ga.filter_significant(cutoff=0.05, use_q=True)   # {"greater": ..., "less": ...}
```

A typed result object is also available:

```python
r = ga.result_obj                       # GAGEResult
r.greater; r.less; r.stats; r.meta      # dataframes + run metadata
r.significant(cutoff=0.1)               # dict of filtered frames
```

---

## The `gage()` convenience function

`gage()` wraps preparation + running for any input and returns a tidy,
direction-labelled frame (or the raw dict with `tidy=False`):

```python
from pygage import gage

# raw matrix
gage(expr, gene_sets, ref_indices=[0,1,2], samp_indices=[3,4,5])
# already-prepared fold-change matrix
gage(prepared, gene_sets, prepared=True)
# DE table / pre-ranked
gage(read_de_table("deseq2.csv"), gene_sets)
# switch test / meta and request extras
gage(expr, gene_sets, ref_indices=[0,1,2], samp_indices=[3,4,5],
     test_method="z-test", meta_method="fisher", compute_effect=True)
```

---

## Understanding the results

Each result table has one row per gene set:

| column | meaning |
|---|---|
| `gene_set` | gene-set / pathway name |
| `set_size` | number of set genes present in the data |
| `stat_mean` | mean per-sample statistic (GAGE `stat.mean`) |
| `p_geomean` | geometric mean of per-sample p-values (direction-specific) |
| `p_val` | combined p-value (Stouffer or Fisher; GAGE `p.val`) |
| `q_val` | Benjamini–Hochberg FDR of `p_val` (GAGE `q.val`) |
| `effect` | *(optional)* mean fold change across the set |
| `leading_edge` | *(optional)* member genes driving the signal |
| `p_perm` | *(optional)* sample-label permutation p-value |

`greater` ranks up-regulated sets, `less` ranks down-regulated sets. Lower
`q_val` = stronger, more significant enrichment.

---

## Collapsing redundant sets (`esset_grp`)

Highly overlapping pathways (e.g. *Lysosome* and *Other glycan degradation*) can
appear together at the top. `esset_grp` is a faithful port of GAGE's `esset.grp`:
it defines each set's **core genes** (those beyond one SD of the gene-mean
distribution in the direction of interest), tests overlap against the pool of
**essential genes** (beyond 2 SD) with the hypergeometric upper tail, and merges
sets whose overlap p-value is below a threshold.

```python
from pygage.results_analysis import esset_grp

groups = esset_grp(res["greater"], prepared, gene_sets,
                   test4up=True, cutoff=0.01, pc=1e-10)
for representative, members in groups["groups"].items():
    if len(members) > 1:
        print(representative, "<=", members)
```

---

## Visualization

```python
from pygage.visualization_utils import EnrichmentPlots

# 1) bubble / dot plot: x = stat.mean, colour = -log10(q), size = set_size
EnrichmentPlots.bubble_plot(res["greater"], top_n=20, output_file="bubble.png")

# 2) cross-condition enrichment heatmap
EnrichmentPlots.enrichment_heatmap(
    {"DCIS_vs_HN": res["greater"], "reverse": res["less"]},
    output_file="enrichment_heatmap.png",
)

# 3) GSEA-style running-enrichment plot (for a ranked list)
info = EnrichmentPlots.running_enrichment(ranked, gene_sets["hsa03050 Proteasome"],
                                          output_file="running.png")
print("ES =", info["ES"], "| leading edge:", len(info["leading_edge"]))

# 4) pathview-style KEGG member colouring -> returns gene->hex for KGML reuse
colors = EnrichmentPlots.pathway_gene_colors(
    gene_sets["hsa03050 Proteasome"], fold_changes, output_file="pathway.png")
```

The `pathway_gene_colors` map plugs straight into a Pathview/SBGNview KGML
overlay, bridging PyGAGE to the RAW Lab Pathview work.

---

## Command-line interface

One command, four subcommands:

```bash
# run GAGE on a matrix or a DE table
pygage run expression.csv -g gene_sets.json -o results.csv \
    --ref 0,1,2 --samp 3,4,5 --compare paired --test t-test --meta stouffer

pygage run deseq2_results.csv -g gene_sets.json -o results.csv --de-table --value log2FC

# download gene sets from KEGG (pathway / KO / module)
pygage kegg pathway -o kegg_hsa.json -s hsa --id-type entrez
pygage kegg ko      -o ko_sets.json  --reference pathway
pygage kegg module  -o modules.json  -s hsa

# build GO gene sets from a GAF (+ optional OBO propagation)
pygage go goa_human.gaf -o go_bp.json --obo go-basic.obo --aspect BP --propagate

# combine result tables across conditions
pygage compare ctrl.csv treat.csv -o combined.tsv --names Control,Treatment

# verbose logging + help
pygage -v run ...          # progress to stderr
pygage run --help
```

**Migration from the 1.0.0 scripts** (the seven `pygage-*.py` wrappers are
replaced by subcommands):

| 1.0.0 script | 1.2.0 command |
|---|---|
| `pygage-core.py` | `pygage run` |
| `pygage-pathway_database_utils.py kegg` | `pygage kegg pathway` |
| `pygage-pathway_database_utils.py go` | `pygage go` |
| `pygage-results_analysis.py compare` | `pygage compare` |
| `pygage-tests.py` | `pygage run --test {t-test,z-test,ks-test}` |

---

## Validation / parity with gage R

PyGAGE is not merely "close" to GAGE — it reproduces it. We ran the **actual
`gage` R package** on its own demo data (`gse16873`, 11,979 genes × 12 samples;
`kegg.gs`, 177 KEGG sets) and fed PyGAGE gage's *exact* prepared fold-change
matrix, so the comparison isolates the statistic + meta step. Across all 160
size-passing sets:

| column | max &#124;Δ&#124; vs gage R |
|---|---|
| `stat.mean` | 4.9e-15 |
| `p.val` (Stouffer) | 1.4e-15 |
| `p.geomean` | 8.9e-16 |
| `q.val` (BH) | 2.9e-15 |
| `set.size` | 0 |

The `less` direction is identical, the **z-test** matches to 1.9e-15, and the
**Fisher/gamma** meta matches to 1.8e-15. Pearson *r* = 1.00000000 on both
−log10(p) and `stat.mean`. `benjamini_hochberg` matches R's
`p.adjust(method="BH")` to 1e-12.

This is locked in as a **regression test** shipped with the package
(`tests/test_regression_gage.py`, tolerance 1e-8) using gzipped gage-output
fixtures, so it runs in CI on every change:

```bash
pip install ".[test]"
pytest -v            # includes the gage-R regression
```

> **Scope note.** The t-test and z-test are the tightly-validated (default) paths.
> KS-mode parity is algorithmic (R's `ks.test` and SciPy's `ks_2samp` differ on
> tie handling / exact-vs-asymptotic p-values). A full reproducibility bundle
> (run gage R yourself and diff) is available on request.

---

## Data format requirements

**Expression matrix (CSV/TSV)** — genes × samples, gene column first:

```
gene_id,HN_1,HN_2,HN_3,DCIS_1,DCIS_2,DCIS_3
7157,5.2,5.4,5.1,8.3,8.5,8.1
672,3.1,3.3,3.2,3.4,3.5,3.3
```

**Gene sets (JSON)** — either a plain mapping or with metadata:

```json
{ "hsa03050 Proteasome": ["5683","5684","5685"], "hsa00190 Oxidative phosphorylation": ["498","513"] }
```

```json
{ "gene_sets": { "path1": ["5683","5684"] }, "pathway_names": { "path1": "Proteasome" } }
```

**Gene-ID mapping (TSV)** — for Entrez ↔ symbol conversion (the real 40,784-row
GAGE `egSymb` map ships with the package):

```
entrez_id	symbol
7157	TP53
672	BRCA1
```

---

## Performance

The polars/numpy/scipy engine handles KEGG/GO-scale collections (hundreds of
sets) comfortably, and `n_jobs` parallelises over gene sets. For very large
collections at many samples (e.g. MSigDB C2 × dozens of samples) or heavy
permutation nulls, the inner per-set loop is a clean target for a native
Rust/PyO3 kernel matching the RAW Lab pure-Rust pattern — a design note ships in
`docs/` and is deliberately gated on the regression test so correctness is never
traded for speed.

---

## API overview

| Module | Key entry points |
|---|---|
| `pygage.core` | `GAGEPreparation`, `GAGEAnalysis`, `GAGEResult`, `benjamini_hochberg` |
| `pygage.io_loaders` | `gage`, `read_matrix`, `read_de_table`, `read_preranked` |
| `pygage.gene_sets` | `load_gmt`, `load_msigdb`, `load_reactome`, `load_go`, `GeneSetCollection`, `GeneSetCache` |
| `pygage.pathway_database_utils` | `KEGGPathwayRetriever`, `GOGeneSetRetriever` |
| `pygage.results_analysis` | `esset_grp`, `ResultsComparator`, `GeneSetGrouper`, `SignificanceFilter` |
| `pygage.visualization_utils` | `EnrichmentPlots`, `HeatmapPlotter`, `VennDiagram`, `ColorUtils` |
| `pygage.gene_id_utils` | `GeneIDConverter` |
| `pygage.data_processing_utils` | `DataTransformer`, `GeneExtractor`, `GeneDataExporter` |
| `pygage.tests` | `GeneSetTests` (t / z / KS, delegating to the core engine) |

Full API docs: **https://pygage.readthedocs.io**

---

## Citation

If you publish results obtained with **PyGAGE**, please cite:

- Figueroa III JL, Brouwer CR, White III RA. 2026. *Statistically resolving
  gene-set enrichment for pathway analysis that is broadly applicable via
  PyGAGE.* bioRxiv.

If you use the original R version, please also cite:

- Luo W, Friedman MS, Shedden K, Hankenson KD, Woolf PJ. 2009. *GAGE: generally
  applicable gene set enrichment for pathway analysis.* BMC Bioinformatics 10:161.
  https://doi.org/10.1186/1471-2105-10-161

---

## License

Creative Commons Attribution-NonCommercial (**CC BY-NC 4.0**) — academic and
non-commercial use. See [LICENSE](LICENSE).

## Support

- **Issues:** https://github.com/raw-lab/pygage/issues
- **Email:** [Dr. Richard Allen White III](mailto:rwhit101@charlotte.edu)

---

## Contributing

Contributions are welcome — additional statistical tests, gene-set databases,
visualizations, and performance work (including the R and Python versions).
Please open an issue or reach out via Support.

## Development setup

```bash
git clone https://github.com/raw-lab/pygage
cd pygage
python -m pip install -e ".[dev]"     # editable install + dev tools
```

The package uses a `src/` layout: the importable code lives in `src/pygage/`.

## Running the tests

```bash
pytest                 # runs tests/ (includes the gage-R regression)
pytest --cov=pygage    # with coverage
```

All contributions must keep the test suite green, including
`tests/test_regression_gage.py`, which asserts machine-precision parity with the
GAGE R package on shipped fixtures.

## Linting / style

```bash
ruff check src tests
```

Please write NumPy/Google-style docstrings for public functions and classes;
the API documentation is generated from them via Sphinx autodoc.

## Building the documentation

```bash
pip install -e ".[docs]"
sphinx-build -b html docs/source docs/build/html
```

## Pull-request checklist

- [ ] Tests pass (`pytest`) and new behaviour is covered by a test.
- [ ] Public API changes are documented (docstrings + a note in `CHANGELOG.md`).
- [ ] No build artifacts, data dumps, or secrets are committed (see `.gitignore`).
- [ ] `ruff check` is clean.

---
