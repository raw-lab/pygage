Installation
============

From PyPI
---------

.. code-block:: bash

   pip install pygage

From Bioconda
-------------

.. code-block:: bash

   conda install -c bioconda pygage

From source
-----------

.. code-block:: bash

   git clone https://github.com/raw-lab/pygage
   cd pygage
   pip install .

Optional extras:

.. code-block:: bash

   pip install ".[anndata]"   # AnnData input support
   pip install ".[test]"      # test suite (incl. the gage-R regression)
   pip install ".[docs]"      # build this documentation

Requirements
------------

* Python >= 3.8
* ``polars`` >= 1.0, ``numpy``, ``scipy``, ``matplotlib``, ``seaborn``,
  ``pandas``, ``pyarrow``, ``requests``
* ``anndata`` (optional, for AnnData input)

Verify
------

.. code-block:: python

   import pygage
   print(pygage.__version__)
