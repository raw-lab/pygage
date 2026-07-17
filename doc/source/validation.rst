Validation / parity with gage R
===============================

PyGAGE reproduces the GAGE R package rather than approximating it. The actual
``gage`` R package was run on its own demo data (``gse16873``, 11,979 genes × 12
samples; ``kegg.gs``, 177 sets), and PyGAGE was fed gage's *exact* prepared
fold-change matrix so the comparison isolates the statistic + meta step.

Across all 160 size-passing sets:

.. list-table::
   :header-rows: 1

   * - Column
     - max \|Δ\| vs gage R
   * - ``stat.mean``
     - 4.9e-15
   * - ``p.val`` (Stouffer)
     - 1.4e-15
   * - ``p.geomean``
     - 8.9e-16
   * - ``q.val`` (BH)
     - 2.9e-15
   * - ``set.size``
     - 0

The ``less`` direction is identical, the z-test matches to 1.9e-15, and the
Fisher/gamma meta matches to 1.8e-15. Pearson *r* = 1.00000000 on both
−log10(p) and ``stat.mean``. ``benjamini_hochberg`` matches R's
``p.adjust(method="BH")`` to 1e-12.

This is shipped as a regression test (``tests/test_regression_gage.py``, tolerance
1e-8) with gzipped gage-output fixtures, so it runs in CI on every change:

.. code-block:: bash

   pip install ".[test]"
   pytest -v

.. note::

   The t-test and z-test are the tightly-validated (default) paths. KS-mode
   parity is algorithmic — R's ``ks.test`` and SciPy's ``ks_2samp`` differ on tie
   handling and exact-vs-asymptotic p-values.
