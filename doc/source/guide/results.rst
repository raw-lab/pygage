Understanding the results
=========================

Each result table has one row per gene set:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Meaning
   * - ``gene_set``
     - gene-set / pathway name
   * - ``set_size``
     - number of set genes present in the data
   * - ``stat_mean``
     - mean per-sample statistic (GAGE ``stat.mean``)
   * - ``p_geomean``
     - geometric mean of per-sample p-values (direction-specific)
   * - ``p_val``
     - combined p-value (Stouffer or Fisher; GAGE ``p.val``)
   * - ``q_val``
     - Benjamini–Hochberg FDR of ``p_val`` (GAGE ``q.val``)
   * - ``effect``
     - *(optional)* mean fold change across the set
   * - ``leading_edge``
     - *(optional)* member genes driving the signal
   * - ``p_perm``
     - *(optional)* sample-label permutation p-value

``greater`` ranks up-regulated sets; ``less`` ranks down-regulated sets. Lower
``q_val`` means stronger, more significant enrichment.

Collapsing redundant sets
-------------------------

Highly overlapping pathways can appear together. :func:`pygage.results_analysis.esset_grp`
is a faithful port of GAGE's ``esset.grp`` that merges sets sharing core genes:

.. code-block:: python

   from pygage.results_analysis import esset_grp

   groups = esset_grp(res["greater"], prepared, gene_sets,
                      test4up=True, cutoff=0.01, pc=1e-10)
   for rep, members in groups["groups"].items():
       if len(members) > 1:
           print(rep, "<=", members)
