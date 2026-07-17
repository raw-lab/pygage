Command-line interface
======================

One command, four subcommands.

Run GAGE
--------

.. code-block:: bash

   pygage run expression.csv -g gene_sets.json -o results.csv \
       --ref 0,1,2 --samp 3,4,5 --compare paired --test t-test --meta stouffer

   pygage run deseq2_results.csv -g gene_sets.json -o results.csv \
       --de-table --value log2FC

Download KEGG gene sets
-----------------------

.. code-block:: bash

   pygage kegg pathway -o kegg_hsa.json -s hsa --id-type entrez
   pygage kegg ko      -o ko_sets.json  --reference pathway
   pygage kegg module  -o modules.json  -s hsa

Build GO gene sets
------------------

.. code-block:: bash

   pygage go goa_human.gaf -o go_bp.json --obo go-basic.obo --aspect BP --propagate

Compare result tables
---------------------

.. code-block:: bash

   pygage compare ctrl.csv treat.csv -o combined.tsv --names Control,Treatment

Use ``pygage -v ...`` for progress logging and ``pygage <command> --help`` for
full options.

Migration from the 1.0.0 scripts
--------------------------------

.. list-table::
   :header-rows: 1

   * - 1.0.0 script
     - 1.2.0 command
   * - ``pygage-core.py``
     - ``pygage run``
   * - ``pygage-pathway_database_utils.py kegg``
     - ``pygage kegg pathway``
   * - ``pygage-pathway_database_utils.py go``
     - ``pygage go``
   * - ``pygage-results_analysis.py compare``
     - ``pygage compare``
   * - ``pygage-tests.py``
     - ``pygage run --test {t-test,z-test,ks-test}``
