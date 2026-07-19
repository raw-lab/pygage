# Changelog

All notable changes to PyGAGE are documented here. This project follows
[Semantic Versioning](https://semver.org).

## [1.2.1] — 2026

Quality & typing hardening (backward compatible) in response to an external code
audit.

### Added
- **PEP 561 typing**: ships a ``py.typed`` marker and a ``pygage._types`` module
  of shared aliases (``FrameLike``, ``GeneSetsLike``, ``PathLike``, ``FrameT``);
  loader and engine entry points are annotated with them.
- **Canonical gene-set schema**: ``pygage.gene_sets.normalize_gene_sets`` and
  ``GeneSetCollection.as_dict`` give one internal registry format; ``run_gage``
  now accepts a mapping **or** a ``GeneSetCollection`` interchangeably.
- **Compute controls** in ``pygage.config``: ``default_n_jobs``,
  ``set_thread_limits`` (Polars/BLAS pools for HPC), and ``thread_config``.
- **Live-refreshable identifier map**: ``config.ensure_egsymb(fetch=True)`` can
  rebuild the Entrez↔symbol map from KEGG, with the packaged file as a cached
  fallback.
- Expanded parameter docstrings (``run_gage``, ``prepare_expression``) and a
  documented compute/threading section; the ``config`` module is now in the API
  reference. New tests for the schema, config, and plotting-boundary helpers.

### Changed
- The visualization layer coerces any Polars/pandas/Arrow/dict table to Polars at
  a single boundary (``_as_polars_frame``), removing implicit cross-backend
  assumptions in the plotters.

## [1.2.0] — 2026

A correctness + capability release. The engine now reproduces the GAGE R package
to machine precision, and the package gains DE-table/AnnData inputs, broad
gene-set sourcing, a unified CLI, plots, tests, and docs.

### Fixed (highlights of 28 audited issues)
- **Restored the GAGE two-level statistic.** The core no longer collapses all
  samples to one per-gene mean and runs a one-sample t-test against a scalar; it
  now performs a per-sample test of each set against the array-wide background,
  combined across samples by a meta-test (Stouffer's Z by default). Rankings and
  p-values now match GAGE.
- **Public API works without a monkey-patch.** `prepare_expression()` preserves
  the gene-ID column, so `run_gage()` and every documented example run under
  direct library use.
- **KEGG retrieval returns gene sets.** The link URL uses the correct `path:`
  entry prefix (not `pathway:`); gene sets come back for every species.
- `row_normalize` normalizes genes (rows), not samples.
- Hardened crash paths: non-BH FDR, empty result sets, dual-significance
  filtering, missing gene-ID map, NaN-containing inputs.
- GO/GAF parsing uses the Aspect column (P/F/C → BP/MF/CC), honours `NOT`, and no
  longer mislabels gene names as GO-term names.
- KS complement is computed by set membership; grouping intersects members with
  the measured universe; deprecated polars `how="outer"` replaced by `how="full"`.

### Added
- **Inputs:** `read_de_table` (DESeq2/edgeR/limma auto-detection), `read_preranked`
  (fgsea-style), pandas + AnnData ingestion, and the `gage()` one-call convenience
  function returning a tidy frame.
- **Gene-set sourcing:** GMT/MSigDB (`load_gmt`, `load_msigdb`), Reactome
  (`load_reactome`), GO with GAF+OBO propagation (`load_go`), a versioned offline
  cache (`GeneSetCache`), and provenance/checksum metadata on every collection.
- **KEGG Orthology:** `list_all_kos()` (full ~26k KO namespace) and
  `download_ko_gene_sets()` for species-agnostic (metagenome/virome/phage) analysis.
- **Method rigor:** optional control-gene-set background, sample-label permutation
  null, effect size + leading-edge genes, global BH across greater∪less, multi-core
  over gene sets, and the Fisher/gamma meta-combination.
- **`esset_grp`:** faithful port of GAGE's `esset.grp` redundancy grouping.
- **Visualization:** `EnrichmentPlots` — bubble/dot plot, cross-condition
  enrichment heatmap, GSEA-style running-enrichment plot, and pathview-style KEGG
  member colouring.
- **Ergonomics:** a single unified `pygage` CLI (`run` / `kegg` / `go` / `compare`),
  a typed `GAGEResult`, `logging` instead of `print`, actionable errors, and
  progress reporting on network fetch loops.
- **Validation & packaging:** a gold-standard gage-R regression test with gzipped
  fixtures, a pytest suite, GitHub Actions CI, a Bioconda recipe, the real 40,784-row
  `egSymb` map, ReadTheDocs documentation, and shipped example data.

### Changed
- The seven per-module `pygage-*.py` scripts are replaced by the unified `pygage`
  command (a compatibility mapping is documented in the README and CLI docs).
- `GeneSetTests` (t/z/KS) now delegates to the single gage-faithful core engine,
  eliminating a second, divergent statistical implementation.

## [1.0.0]
- Initial public release: polars-based Python port of GAGE with per-module CLI
  scripts.
