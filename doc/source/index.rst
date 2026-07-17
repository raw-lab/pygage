PyGAGE documentation
====================

**PyGAGE** is a fast, dependency-light Python implementation of **GAGE**
(*Generally Applicable Gene-set Enrichment*; Luo *et al.* 2009) for pathway
analysis. It reproduces the GAGE R package **to machine precision** on real data
(~1e-15 across every reported column), and adds first-class support for DE tables
(DESeq2/edgeR/limma), pre-ranked vectors and AnnData; gene-set sourcing from
KEGG, KEGG Orthology, GO, Reactome and MSigDB; extra statistical rigor; a unified
command-line interface; and publication-ready plots.

Built on polars / numpy / scipy / seaborn, consistent with the RAW Lab toolchain
(MetaCerberus, MerCat2, NFixDB, DeGenPrime).

.. note::

   PyGAGE is validated numerically against the GAGE R package. See
   :doc:`validation` for the head-to-head results and the shipped regression test.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart
   method

.. toctree::
   :maxdepth: 2
   :caption: User guide

   guide/inputs
   guide/genesets
   guide/running
   guide/results
   guide/visualization

.. toctree::
   :maxdepth: 2
   :caption: Reference

   cli
   validation
   performance
   api
   changelog

Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
