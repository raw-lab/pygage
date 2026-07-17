Performance & the Rust kernel
=============================

The polars/numpy/scipy engine handles KEGG/GO-scale collections (hundreds of
sets) comfortably, and ``n_jobs`` parallelises the per-set loop over cores. For
very large collections at many samples (e.g. MSigDB C2 × dozens of samples) or
heavy permutation nulls, the inner per-set statistic loop is a clean target for a
native Rust/PyO3 kernel, matching the RAW Lab pure-Rust pattern.

The design note below is deliberately **gated on the gage-R regression test** so
correctness is never traded for speed.

.. include:: ../RUST_KERNEL.md
   :parser: myst_parser.sphinx_
