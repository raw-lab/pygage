.. raw:: html

   <div class="pg-tagline">GENE-SET ENRICHMENT · PYTHON</div>

PyGAGE
======

.. raw:: html

   <p class="pg-hero">A fast, faithful Python port of <strong>GAGE</strong>
   (Generally Applicable Gene-set Enrichment) for pathway analysis — matching the
   GAGE&nbsp;R package to <strong>machine precision</strong> (~1e-15), with modern
   inputs, broad gene-set sourcing, and publication-ready plots.</p>

.. grid:: 2
   :gutter: 2
   :margin: 0 4 0 0

   .. grid-item::

      .. button-ref:: quickstart
         :ref-type: doc
         :color: primary
         :expand:

         Get started →

   .. grid-item::

      .. button-link:: https://github.com/raw-lab/pygage
         :color: secondary
         :expand:

         View on GitHub

----

Why PyGAGE
----------

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: :octicon:`check-circle` It *is* GAGE
      :class-card: sd-border-1

      The real two-level model — a per-sample test against the array-wide
      background, combined by a cross-sample meta-test. Rankings and p-values
      match the R package.

   .. grid-item-card:: :octicon:`table` Inputs you have
      :class-card: sd-border-1

      DESeq2 / edgeR / limma tables, pre-ranked vectors, **AnnData**, pandas, or a
      raw matrix — no mandatory restructuring.

   .. grid-item-card:: :octicon:`database` Gene sets everywhere
      :class-card: sd-border-1

      KEGG, KEGG **Orthology** (metagenomes / viromes), GO (GAF + OBO), Reactome
      and MSigDB — with a versioned offline cache.

   .. grid-item-card:: :octicon:`zap` Fast by design
      :class-card: sd-border-1

      Polars + vectorised NumPy/SciPy, analytic nulls (no slow permutations by
      default), and multi-core over gene sets.

   .. grid-item-card:: :octicon:`graph` Beautiful plots
      :class-card: sd-border-1

      Bubble/dot plots, cross-condition heatmaps, GSEA running-enrichment, and
      pathview-style KEGG colouring.

   .. grid-item-card:: :octicon:`terminal` One clean CLI
      :class-card: sd-border-1

      A single ``pygage`` command — ``run`` / ``kegg`` / ``go`` / ``compare`` —
      plus a typed, importable API.

----

.. admonition:: Validated against gage R — to machine precision
   :class: tip

   Run on the GAGE demo data (``gse16873``, 11,979 genes × 12 samples), fed the
   *identical* prepared matrix, PyGAGE reproduces every reported column:

   .. list-table::
      :header-rows: 1
      :class: pg-nums
      :widths: 60 40

      * - Column
        - max \|Δ\| vs gage R
      * - ``stat.mean``
        - 4.9e-15
      * - ``p.val`` (Stouffer)
        - 1.4e-15
      * - ``q.val`` (BH)
        - 2.9e-15

   Locked in as a shipped regression test. See :doc:`validation`.

----

Install
-------

.. tab-set::

   .. tab-item:: pip

      .. code-block:: bash

         pip install pygage

   .. tab-item:: conda

      .. code-block:: bash

         conda install -c bioconda pygage

   .. tab-item:: source

      .. code-block:: bash

         git clone https://github.com/raw-lab/pygage
         cd pygage && pip install ".[anndata]"

60-second example
-----------------

.. code-block:: python

   import json
   from pathlib import Path
   import polars as pl
   from pygage import core, gage

   reg = Path(core.__file__).parent / "data" / "regression"
   prepared  = pl.read_csv(reg / "gse16873_prepared.csv.gz",
                           schema_overrides={"gene_id": pl.Utf8})
   gene_sets = json.loads((reg / "kegg_gs.json").read_text())

   result = gage(prepared, gene_sets, prepared=True)   # tidy, direction-labelled
   print(result.filter(pl.col("direction") == "greater").sort("p_val").head(5))

The top hits are the published GAGE vignette results — matching R to ~1e-15.

.. toctree::
   :hidden:
   :caption: Getting started

   installation
   quickstart
   method

.. toctree::
   :hidden:
   :caption: User guide

   guide/inputs
   guide/genesets
   guide/running
   guide/results
   guide/visualization

.. toctree::
   :hidden:
   :caption: Reference

   cli
   validation
   performance
   api
   changelog
