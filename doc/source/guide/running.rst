Running GAGE
============

:meth:`pygage.core.GAGEAnalysis.run_gage` is the full engine. Defaults reproduce
gage R exactly; the remaining arguments are opt-in extensions.

.. code-block:: python

   from pygage import GAGEAnalysis

   ga  = GAGEAnalysis()
   res = ga.run_gage(
       prepared, gene_sets,
       gene_col="gene_id",
       set_size_range=(10, 500),

       # statistic & combination (defaults = gage R)
       test_method="t-test",      # t-test | z-test | ks-test
       meta_method="stouffer",    # stouffer | fisher
       same_dir=True,             # separate greater/less (False = |changes|)
       fdr_method="BH",

       # extensions (opt-in)
       control_genes=None,        # background from a control set, not all genes
       global_bh=False,           # BH across greater ∪ less
       compute_effect=True,       # per-set mean fold change
       leading_edge=False,        # member genes driving the signal
       permutations=0,            # sample-label permutation null (e.g. 1000)
       n_jobs=1,                  # parallelise over gene sets (-1 = all cores)
   )
   greater, less, stats = res["greater"], res["less"], res["stats"]

   sig = ga.filter_significant(cutoff=0.05, use_q=True)

Typed result
------------

.. code-block:: python

   r = ga.result_obj             # pygage.core.GAGEResult
   r.greater; r.less; r.stats; r.meta
   r.significant(cutoff=0.1)

The ``gage()`` convenience function
-----------------------------------

:func:`pygage.io_loaders.gage` wraps preparation + running for any input and
returns a tidy, direction-labelled frame:

.. code-block:: python

   from pygage import gage

   gage(expr, gene_sets, ref_indices=[0,1,2], samp_indices=[3,4,5])
   gage(prepared, gene_sets, prepared=True)
   gage(read_de_table("deseq2.csv"), gene_sets)
