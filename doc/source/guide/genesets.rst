Gene-set sourcing
=================

Loaders return either a plain ``{name: [genes]}`` dict or a
:class:`pygage.gene_sets.GeneSetCollection` carrying provenance (source, release,
retrieval date, checksum).

KEGG pathways (any organism)
----------------------------

.. code-block:: python

   from pygage.pathway_database_utils import KEGGPathwayRetriever

   kegg = KEGGPathwayRetriever()
   res = kegg.get_pathway_genes("hsa", id_type="entrez")   # mmu, eco, ath, ...

KEGG Orthology (species-agnostic)
---------------------------------

For metagenomes, viromes, phages and non-model organisms, annotate genes to KO
IDs and use KO-keyed gene sets:

.. code-block:: python

   kegg = KEGGPathwayRetriever()
   all_kos = kegg.list_all_kos()                        # full ~26k KO namespace
   ko_sets = kegg.get_ko_gene_sets(reference="pathway")

   # download once, reproducibly, with names + categories + catalog + provenance
   kegg.download_ko_gene_sets("ko_gene_sets.json", reference="pathway")

Gene Ontology (GAF + OBO, with propagation)
-------------------------------------------

.. code-block:: python

   from pygage.gene_sets import load_go

   go = load_go("goa_human.gaf", obo_path="go-basic.obo",
                aspect="BP", include_iea=True, propagate=True)

Reactome
--------

.. code-block:: python

   from pygage.gene_sets import load_reactome

   rx = load_reactome("ReactomePathways.gmt", id_type="gmt")
   rx = load_reactome("NCBI2Reactome_All_Levels.txt", id_type="ncbi2reactome",
                      species="Homo sapiens")

MSigDB / any GMT
----------------

.. code-block:: python

   from pygage.gene_sets import load_gmt, load_msigdb

   hallmark = load_gmt("h.all.v2023.2.Hs.symbols.gmt", source="MSigDB", release="2023.2")
   c2       = load_msigdb("c2.cp.v2023.2.Hs.symbols.gmt", collection="C2")

Versioned offline cache
-----------------------

.. code-block:: python

   from pygage.gene_sets import GeneSetCache

   cache = GeneSetCache()                   # ~/.cache/pygage/gene_sets
   cache.save("hallmark_2023.2", hallmark)
   same = cache.load("hallmark_2023.2")     # checksum-verified, offline
