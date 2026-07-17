Quickstart
==========

PyGAGE ships the real GAGE demo data (``gse16873``: 6 head-and-neck vs 6 DCIS
breast-cancer samples) so you can run a full analysis with no downloads.

.. code-block:: python

   import json
   from pathlib import Path
   import polars as pl
   from pygage import core, gage

   reg = Path(core.__file__).parent / "data" / "regression"
   prepared  = pl.read_csv(reg / "gse16873_prepared.csv.gz",
                           schema_overrides={"gene_id": pl.Utf8})
   gene_sets = json.loads((reg / "kegg_gs.json").read_text())

   # one call -> a tidy, direction-labelled result frame
   result = gage(prepared, gene_sets, prepared=True)
   print(result.filter(pl.col("direction") == "greater").sort("p_val").head(5))

The top hits are the published GAGE vignette results (secretory / degradation
machinery up in DCIS vs HN) and match the R package to ~1e-15.

From the command line
---------------------

.. code-block:: bash

   pygage run expression.csv -g gene_sets.json -o results.csv \
       --ref 0,1,2 --samp 3,4,5

Next steps
----------

* :doc:`guide/inputs` — raw matrices, DE tables, pre-ranked vectors, AnnData
* :doc:`guide/genesets` — KEGG, KO, GO, Reactome, MSigDB
* :doc:`guide/running` — the full engine and its options
* :doc:`guide/visualization` — bubble plots, heatmaps, running-enrichment
