Inputs
======

PyGAGE accepts a raw expression matrix, a differential-expression table, a
pre-ranked vector, or an AnnData object.

Raw expression matrix
---------------------

Genes × samples, gene column first:

.. code-block:: python

   from pygage import read_matrix, GAGEPreparation, GAGEAnalysis

   expr = read_matrix("expression.csv")
   prepared = GAGEPreparation.prepare_expression(
       expr, ref_indices=[0, 1, 2], samp_indices=[3, 4, 5],
       comparison="paired",   # paired | unpaired | as.group | 1ongroup
       input_logged=True,     # False -> log2(x + 1) first
   )
   result = GAGEAnalysis().run_gage(prepared, gene_sets)

DE tables (DESeq2 / edgeR / limma)
----------------------------------

:func:`pygage.io_loaders.read_de_table` auto-detects the gene, log2-fold-change,
and statistic columns by common aliases:

.. code-block:: python

   from pygage import read_de_table, gage

   de = read_de_table("deseq2_results.csv", value="log2FC")   # or value="stat"
   result = gage(de, gene_sets)

Pre-ranked vector (fgsea-style)
-------------------------------

.. code-block:: python

   from pygage import read_preranked, gage

   ranked = read_preranked({"TP53": 3.1, "BRCA1": -2.4, "EGFR": 1.8})
   result = gage(ranked, gene_sets)

pandas / AnnData
----------------

.. code-block:: python

   import anndata as ad
   from pygage import gage

   adata = ad.read_h5ad("counts.h5ad")   # obs × var (samples × genes)
   result = gage(adata, gene_sets, ref_indices=[0, 1, 2], samp_indices=[3, 4, 5])

AnnData is transposed to genes × samples automatically; pandas frames are
accepted anywhere a polars frame is.
