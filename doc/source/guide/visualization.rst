Visualization
=============

:class:`pygage.visualization_utils.EnrichmentPlots` provides publication-ready
views of GAGE results.

.. code-block:: python

   from pygage.visualization_utils import EnrichmentPlots

   # bubble / dot plot: x = stat.mean, colour = -log10(q), size = set_size
   EnrichmentPlots.bubble_plot(res["greater"], top_n=20, output_file="bubble.png")

   # cross-condition enrichment heatmap
   EnrichmentPlots.enrichment_heatmap(
       {"DCIS_vs_HN": res["greater"], "reverse": res["less"]},
       output_file="enrichment_heatmap.png")

   # GSEA-style running-enrichment plot for a ranked list
   info = EnrichmentPlots.running_enrichment(
       ranked, gene_sets["hsa03050 Proteasome"], output_file="running.png")

   # pathview-style KEGG member colouring -> gene->hex map for KGML reuse
   colors = EnrichmentPlots.pathway_gene_colors(
       gene_sets["hsa03050 Proteasome"], fold_changes, output_file="pathway.png")

The ``pathway_gene_colors`` map plugs directly into a Pathview/SBGNview KGML
overlay.
