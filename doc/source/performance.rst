Performance & the Rust kernel
=============================


Scaling on shared / HPC systems
-------------------------------

The per-gene-set loop parallelises over cores. Choose the worker count per call
(``run_gage(..., n_jobs=-1)`` for all cores) or set a default and cap the
dataframe/BLAS thread pools so PyGAGE behaves well on shared nodes:

.. code-block:: python

   from pygage import config

   # default worker count for the gene-set loop (or export PYGAGE_NUM_THREADS)
   config.default_n_jobs()

   # cap Polars + BLAS thread pools (call before heavy work / before importing polars)
   config.set_thread_limits(8)

   # inspect what is active
   config.thread_config()

On a scheduler, set the environment once so every worker inherits it::

   export PYGAGE_NUM_THREADS=8
   export POLARS_MAX_THREADS=8

The polars/numpy/scipy engine handles KEGG/GO-scale collections (hundreds of
sets) comfortably, and ``n_jobs`` parallelises the per-set loop over cores. For
very large collections at many samples (e.g. MSigDB C2 × dozens of samples) or
heavy permutation nulls, the inner per-set statistic loop is a clean target for a
native Rust/PyO3 kernel, matching the RAW Lab pure-Rust pattern.

The design note below is deliberately **gated on the gage-R regression test** so
correctness is never traded for speed.

.. include:: ../RUST_KERNEL.md
   :parser: myst_parser.sphinx_
