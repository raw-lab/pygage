# pygage - GAGE Analysis in Python

Python version of the GAGE (Generally Applicable Gene Set Enrichment) analysis package.

## Overview

Uses Python modules using:
- **polars** for fast dataframe operations (instead of pandas)
- **seaborn/matplotlib** for visualization
- **scipy** for statistical tests
- **argparse** for command-line interfaces
- **numpy** for numerical operations

## Installation

### Quick install

```bash
pip install pygage
```

### Custom install

```bash
# Clone repository
git clone https://github.com/raw-lab/pygage
cd pygage

# Install dependencies
pip install polars numpy scipy matplotlib seaborn requests
pip install .
```

## Modules

### 1. gene_id_utils.py
Convert between Entrez Gene IDs and gene symbols.

**Usage:**
```bash
# Convert Entrez IDs to symbols
pygage-gene_id_utils.py input_ids.txt \
    --mapping egSymb.csv \
    --direction eg2sym \
    --output output_symbols.csv

# Convert symbols to Entrez IDs
pygage-gene_id_utils.py input_symbols.txt \
    --mapping egSymb.csv \
    --direction sym2eg \
    --output output_ids.csv
```

**Python API:**
```python
from pygage.gene_id_utils import GeneIDConverter

converter = GeneIDConverter('egSymb.csv')
symbols = converter.eg2sym(['1', '2', '3'])
entrez_ids = converter.sym2eg(['TP53', 'BRCA1', 'EGFR'])
```

---

### 2. pathway_database_utils.py
Retrieve gene sets from KEGG and Gene Ontology databases.

**Usage:**
```bash
# Retrieve KEGG pathways
pygage-pathway_database_utils.py kegg \
    --species hsa \
    --id-type entrez \
    --output kegg_pathways.json

# Retrieve GO gene sets
pygage-pathway_database_utils.py go \
    --annotation-file goa_human.gaf \
    --species human \
    --output go_genesets.json
```

**Python API:**
```python
from pygage.pathway_database_utils import KEGGPathwayRetriever, GOGeneSetRetriever

# KEGG
kegg = KEGGPathwayRetriever()
pathways = kegg.get_pathway_genes(species='hsa', id_type='entrez')

# GO
go = GOGeneSetRetriever()
gene_sets = go.get_go_gene_sets(species='human', annotation_file='goa_human.gaf')
```

---

### 3. visualization_utils.py
Create heatmaps, Venn diagrams, and color palettes.

**Usage:**
```bash
# Create Venn diagram
pygage-visualization_utils.py venn \
    --input comparison_data.csv \
    --names "Sample1" "Sample2" "Sample3" \
    --include both \
    --output venn_diagram.png

# Create heatmap
pygage-visualization_utils.py heatmap \
    --input expression_data.csv \
    --output heatmap.png \
    --cluster \
    --cmap RdYlGn_r \
    --title "Gene Expression Heatmap"
```

**Python API:**
```python
from pygage.visualization_utils import HeatmapPlotter, VennDiagram, ColorUtils
import polars as pl

# Heatmap
data = pl.read_csv('expression_data.csv')
plotter = HeatmapPlotter()
plotter.plot_clustered_heatmap(data, output_file='heatmap.png')

# Venn diagram
venn = VennDiagram()
counts = venn.venn_counts(data, include='both')
venn.plot_venn2(counts, ['Set1', 'Set2'], 'venn.png')

# Color palette
cmap = ColorUtils.greenred(256)
```

---

### 4. data_processing_utils.py
Data transformation, normalization, and gene extraction.

**Usage:**
```bash
# Normalize data
pygage-data_processing_utils.py normalize \
    --input expression_data.csv \
    --output normalized_data.csv

# Extract essential genes
pygage-data_processing_utils.py extract \
    --input expression_data.csv \
    --genes gene_list.txt \
    --output essential_genes.csv \
    --threshold 1.5

# Export and visualize
pygage-data_processing_utils.py export \
    --input expression_data.csv \
    --genes gene_list.txt \
    --output gene_data.csv \
    --heatmap gene_heatmap.png \
    --normalize
```

**Python API:**
```python
from pygage.data_processing_utils import DataTransformer, GeneExtractor, GeneDataExporter
import polars as pl

# Normalize
data = pl.read_csv('expression_data.csv')
transformer = DataTransformer()
normalized = transformer.row_normalize(data)

# Extract essential genes
extractor = GeneExtractor()
essential = extractor.extract_essential_genes(
    gene_set=['TP53', 'BRCA1', 'EGFR'],
    expression_data=data,
    threshold=1.0
)

# Export with visualization
exporter = GeneDataExporter()
exporter.export_gene_data(
    genes=['TP53', 'BRCA1'],
    expression_data=data,
    output_file='output.csv',
    create_heatmap=True,
    heatmap_output='heatmap.png'
)
```

---

### 5. gage_core.py
Core GAGE analysis functions.

**Usage:**
```bash
pygage-core.py \
    --expression expression_data.csv \
    --gene-sets pathways.json \
    --gene-col gene_id \
    --ref-indices 0 1 2 \
    --samp-indices 3 4 5 \
    --comparison paired \
    --test-method t-test \
    --cutoff 0.1 \
    --output results/
```

**Python API:**
```python
from pygage.gage_core import GAGEPreparation, GAGEAnalysis
import polars as pl
import json

# Load data
expr_data = pl.read_csv('expression_data.csv')
with open('gene_sets.json') as f:
    gene_sets = json.load(f)

# Prepare data
prep = GAGEPreparation()
prepared = prep.prepare_expression(
    expr_data,
    ref_indices=[0, 1, 2],
    samp_indices=[3, 4, 5],
    comparison='paired'
)

# Run GAGE
gage = GAGEAnalysis()
results = gage.run_gage(
    prepared,
    gene_sets,
    test_method='t-test'
)

# Filter significant
significant = gage.filter_significant(cutoff=0.1)
```

---

### 6. gage_tests.py
Statistical tests for gene set analysis.

**Usage:**
```bash
pygage-tests.py \
    --expression expression_data.csv \
    --gene-sets pathways.json \
    --gene-col gene_id \
    --method t-test \
    --min-size 10 \
    --max-size 500 \
    --output test_results.tsv
```

**Python API:**
```python
from pygage.tests import GeneSetTests
import polars as pl
import json

expr_data = pl.read_csv('expression_data.csv')
with open('gene_sets.json') as f:
    gene_sets = json.load(f)

tester = GeneSetTests()

# t-test
results = tester.t_test(expr_data, gene_sets)

# z-test
results = tester.z_test(expr_data, gene_sets)

# KS test
results = tester.kolmogorov_smirnov_test(expr_data, gene_sets)
```

---

### 7. results_analysis.py
Analyze and compare GAGE results.

**Usage:**
```bash
# Compare multiple results
pygage-results_analysis.py compare \
    --inputs result1.tsv result2.tsv result3.tsv \
    --names "Control" "Treatment1" "Treatment2" \
    --cutoff 0.1 \
    --output combined_results.tsv \
    --venn comparison_venn.png

# Filter significant results
pygage-results_analysis.py filter \
    --greater greater_results.tsv \
    --less less_results.tsv \
    --cutoff 0.1 \
    --output-dir filtered_results/

# Group overlapping gene sets
pygage-results_analysis.py group \
    --results gage_results.tsv \
    --gene-sets pathways.json \
    --expression expression_data.csv \
    --output gene_set_groups.json
```

**Python API:**
```python
from pygage.results_analysis import ResultsComparator, GeneSetGrouper, SignificanceFilter
from pathlib import Path

# Compare results
comparator = ResultsComparator()
combined = comparator.compare_results(
    [Path('r1.tsv'), Path('r2.tsv')],
    ['Sample1', 'Sample2'],
    q_cutoff=0.1,
    output_file=Path('combined.tsv')
)

# Create Venn diagram
comparator.create_venn_comparison(
    [Path('r1.tsv'), Path('r2.tsv')],
    ['Sample1', 'Sample2'],
    output_file=Path('venn.png')
)

# Filter significant
filterer = SignificanceFilter()
filtered = filterer.filter_significant(
    results={'greater': greater_df, 'less': less_df},
    cutoff=0.1
)

# Group gene sets
grouper = GeneSetGrouper()
groups = grouper.group_gene_sets(
    results_df,
    gene_sets,
    expression_df,
    output_file=Path('groups.json')
)
```

---

## Complete Workflow Example

```bash
#!/bin/bash

# 1. Convert gene IDs (if needed)
pygage-gene_id_utils.py gene_list.txt \
    --mapping egSymb.csv \
    --direction eg2sym \
    --output gene_symbols.csv

# 2. Retrieve KEGG pathways
pygage-pathway_database_utils.py kegg \
    --species hsa \
    --id-type entrez \
    --output kegg_pathways.json

# 3. Prepare and normalize expression data
pygage-data_processing_utils.py normalize \
    --input raw_expression.csv \
    --output normalized_expression.csv

# 4. Run GAGE analysis
pygage-core.py \
    --expression normalized_expression.csv \
    --gene-sets kegg_pathways.json \
    --gene-col gene_id \
    --ref-indices 0 1 2 \
    --samp-indices 3 4 5 \
    --comparison paired \
    --test-method t-test \
    --cutoff 0.1 \
    --output gage_results/

# 5. Filter significant results
pygage-results_analysis.py filter \
    --greater gage_results/greater.tsv \
    --less gage_results/less.tsv \
    --cutoff 0.05 \
    --output-dir significant_results/

# 6. Create visualizations
pygage-visualization_utils.py heatmap \
    --input significant_results/greater_significant.tsv \
    --output heatmap_upregulated.png \
    --cluster \
    --title "Up-regulated Pathways"

# 7. Extract essential genes from top pathway
pygage-data_processing_utils.py extract \
    --input normalized_expression.csv \
    --genes top_pathway_genes.txt \
    --output essential_genes.csv \
    --threshold 2.0

# 8. Export with visualization
pygage-data_processing_utils.py export \
    --input normalized_expression.csv \
    --genes essential_genes.csv \
    --output essential_genes_data.csv \
    --heatmap essential_genes_heatmap.png \
    --normalize
```

---

## Data Format Requirements

### Expression Data (CSV/TSV)
```
gene_id,sample1,sample2,sample3,sample4,sample5,sample6
GENE001,5.2,5.4,5.1,8.3,8.5,8.1
GENE002,3.1,3.3,3.2,3.4,3.5,3.3
GENE003,7.8,7.9,8.0,4.2,4.1,4.3
```

### Gene Sets (JSON)
```json
{
  "pathway1": ["GENE001", "GENE002", "GENE005"],
  "pathway2": ["GENE003", "GENE004", "GENE006"],
  "pathway3": ["GENE001", "GENE003", "GENE007"]
}
```

Or with metadata:
```json
{
  "gene_sets": {
    "pathway1": ["GENE001", "GENE002"],
    "pathway2": ["GENE003", "GENE004"]
  },
  "pathway_names": {
    "pathway1": "Cell cycle regulation",
    "pathway2": "DNA repair"
  }
}
```

### Gene ID Mapping (CSV/TSV)
```
entrez_id,symbol
1,TP53
2,BRCA1
3,EGFR
```

---

## Key Differences from R Version

### 1. **Polars instead of Pandas**
- Faster performance
- More intuitive API
- Better memory efficiency
- Lazy evaluation support

```python
# Polars syntax
df.filter(pl.col('value') > 0.5)
df.select(['col1', 'col2'])
df.join(other_df, on='key', how='left')
```

### 2. **Seaborn/Matplotlib instead of base R graphics**
- More modern visualizations
- Better default aesthetics
- Easier customization

### 3. **Argparse for CLI**
- Subcommands for different operations
- Help text and validation built-in
- Type checking

### 4. **No Random Seed Issues**
- Random seed properly handled by scipy
- No hardcoded values

### 5. **JSON for Gene Sets**
- More portable than R's RData format
- Human-readable
- Easy to integrate with other tools

---

## Error Handling

All scripts include comprehensive error handling:

```python
try:
    converter = GeneIDConverter('mapping.csv')
    symbols = converter.eg2sym(['1', '2', '3'])
except ValueError as e:
    print(f"Error: {e}")
except FileNotFoundError:
    print("Mapping file not found")
```

---

## Troubleshooting

### Issue: "No numeric columns found"
**Solution:** Ensure expression data has numeric columns and correct gene ID column name

### Issue: "Mapping data not loaded"
**Solution:** Call `load_mapping()` or initialize with mapping file path

### Issue: "Can't create Venn diagram for more than 3 sets"
**Solution:** Venn diagrams limited to 2-3 sets; use heatmap for more comparisons

### Issue: "Memory error with large datasets"
**Solution:** Use polars lazy evaluation or process in chunks

---

## Contributing

Improvements welcome:
1. Add more statistical tests
2. Implement additional visualizations
3. Add support for more gene set databases
4. Improve performance with parallel processing

We welcome contributions of other experts expanding features in PyGAGE including the R and python versions. Please contact us via support. 

---

## 📄 License

Creative Commons Attribution-NonCommercial (CC BY-NC 4.0) — See LICENSE file

## 📚 Citing

If you are publishing results obtained using PyGAGE, please cite: <br />
- Pre-Print PyGAGE: Figueroa III JL, Brouwer CR, White III RA. 2026. Statistically resolving gene-set enrichment for pathway analysis that is broadly applicable via PyGAGE. bioRxiv.

If you using the R version please cite: <br />
- Luo, W., Friedman, M. S., Shedden, K., Hankenson, K. D., & Woolf, P. J. (2009). GAGE: generally applicable gene set enrichment for pathway analysis.
BMC Bioinformatics, 10, 161. [GAGE](https://doi.org/10.1186/1471-2105-10-161)

---

## 📞 Support

- **Issues:** [open an issue](https://github.com/raw-lab/pygage/issues).  
- **Email:** [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)

---
