# Reference & demo data

This directory holds **optional demo and reference data** used by the examples
and documentation. It is **not** required by the installed package and is **not**
bundled into the wheel — PyGAGE ships only a minimal set of runtime assets inside
the package (`src/pygage/data/`: the Entrez↔symbol map and small regression
fixtures).

Production gene sets are obtained **live**, so database updates never require a
package version bump:

- KEGG pathways / modules / Orthology — fetched via `pygage.pathway_database_utils.KEGGPathwayRetriever`
- Gene Ontology — parsed from a user-supplied GAF/OBO via `pygage.gene_sets.load_go`
- Reactome / MSigDB — loaded from user-supplied GMT/mapping files

You can also point PyGAGE at a curated local data directory without editing code:

```bash
export PYGAGE_DATA_DIR=/path/to/my/pygage-data   # overrides packaged assets
```

## Files (derived from the GAGE R package demo data, Luo et al. 2009)

| File | Description |
|------|-------------|
| `gse16873.demo` | small text sample of the GAGE demo expression set |
| `kegg_gs.tsv` / `kegg_gs.json` | human KEGG signaling/metabolism gene sets |
| `kegg_gs_dise.tsv` / `kegg.gs.dise` | human KEGG disease gene sets |
| `go_gs.tsv` | GO gene sets |
| `carta_gs.tsv` | BioCarta gene sets |
| `c2.demo.gmt` | small MSigDB C2 GMT sample |
| `korg.tsv`, `khier.tsv`, `bods.tsv` | KEGG organism / hierarchy / species-code tables |

The full `gse16873` expression matrix is **not bundled** (it is available from the
`gage`/`gageData` R packages, and `results/02_run_gage_reference.R` regenerates the
prepared form). The package ships only the gzipped *prepared* fold-change matrix
inside `src/pygage/data/regression/` for the regression test.

Provenance: converted from the `gage`/`gageData` R packages
(https://bioconductor.org/packages/gage). See the top-level `CITATION.cff`.
