__version__ = "1.2.0"
__authors__ = "Richard Allen White III, Jose Luis Figueroa III"

from . import (
    core,
    data_processing_utils,
    gene_id_utils,
    gene_sets,
    io_loaders,
    pathway_database_utils,
    results_analysis,
    tests,
    visualization_utils,
)

# convenience re-exports
from .core import GAGEAnalysis, GAGEPreparation, GAGEResult, benjamini_hochberg
from .io_loaders import gage, read_matrix, read_de_table, read_preranked
from .gene_sets import (
    GeneSetCollection,
    GeneSetCache,
    load_gmt,
    load_msigdb,
    load_reactome,
    load_go,
)
from .results_analysis import esset_grp
from .visualization_utils import EnrichmentPlots

__all__ = [
    # modules
    "core", "data_processing_utils", "gene_id_utils", "gene_sets",
    "io_loaders", "pathway_database_utils", "results_analysis", "tests",
    "visualization_utils",
    # core API
    "GAGEAnalysis", "GAGEPreparation", "GAGEResult", "benjamini_hochberg",
    "gage", "read_matrix", "read_de_table", "read_preranked",
    "GeneSetCollection", "GeneSetCache", "load_gmt", "load_msigdb",
    "load_reactome", "load_go", "esset_grp", "EnrichmentPlots",
]
