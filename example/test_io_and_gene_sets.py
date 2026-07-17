"""Input adapters and gene-set loaders."""
import gzip, json, os, tempfile
import numpy as np
import polars as pl
import pytest

from pygage.io_loaders import read_de_table, read_preranked, gage, _to_polars
from pygage.gene_sets import load_gmt, load_reactome, load_go, GeneSetCache


def test_de_table_deseq2_aliases():
    de = pl.DataFrame({"gene": ["1", "2", "3"], "baseMean": [100, 200, 50],
                       "log2FoldChange": [2.1, -1.3, 0.2], "stat": [5.0, -3.1, 0.4]})
    t = read_de_table(de)
    assert t.columns == ["gene_id", "log2FC"] and t.height == 3
    t2 = read_de_table(de, value="stat")
    assert "stat" in t2.columns


def test_de_table_edger_aliases():
    de = pl.DataFrame({"GeneID": ["a", "b"], "logFC": [1.0, -2.0], "LR": [4.0, 9.0], "FDR": [0.1, 0.2]})
    assert read_de_table(de).columns == ["gene_id", "log2FC"]


def test_preranked_dict_and_frame():
    assert read_preranked({"g1": 2.0, "g2": -1.0}).height == 2
    fr = pl.DataFrame({"gene_id": ["a", "b"], "score": [1.0, 2.0]})
    assert read_preranked(fr).height == 2


def test_anndata_ingestion():
    ad = pytest.importorskip("anndata")
    X = np.random.default_rng(0).normal(size=(6, 40))
    a = ad.AnnData(X)
    a.var_names = [f"g{i}" for i in range(40)]
    a.obs_names = [f"s{j}" for j in range(6)]
    fr = _to_polars(a)
    assert fr.shape == (40, 7) and fr.columns[0] == "gene_id"


def test_gmt_metadata_and_checksum():
    t = tempfile.mkdtemp(); p = os.path.join(t, "h.gmt")
    open(p, "w").write("SET_A\tdesc\tG1\tG2\tG3\nSET_B\tdesc\tG4\tG5\n")
    c = load_gmt(p, source="MSigDB", release="2023.2")
    assert c.n_sets == 2 and len(c.checksum) == 16
    assert c.metadata()["source"] == "MSigDB" and c.metadata()["release"] == "2023.2"


def test_reactome_species_filter():
    t = tempfile.mkdtemp(); p = os.path.join(t, "n2r.txt")
    open(p, "w").write("\n".join([
        "1\tR-HSA-1\turl\tPathA\tTAS\tHomo sapiens",
        "2\tR-HSA-1\turl\tPathA\tTAS\tHomo sapiens",
        "3\tR-MMU-1\turl\tPathA\tTAS\tMus musculus",
    ]))
    c = load_reactome(p, id_type="ncbi2reactome", species="Homo sapiens")
    assert c.gene_sets == {"R-HSA-1": ["1", "2"]}


def test_go_obo_propagation():
    t = tempfile.mkdtemp()
    gaf = os.path.join(t, "t.gaf")
    open(gaf, "w").write("\n".join([
        "!gaf",
        "DB\tP1\tA\t\tGO:0006915\tPMID\tIDA\t\tP\tn\t\tprotein\ttaxon:9606\t2020\tX",
        "DB\tP2\tB\t\tGO:0097190\tPMID\tIDA\t\tP\tn\t\tprotein\ttaxon:9606\t2020\tX",
    ]))
    obo = os.path.join(t, "go.obo")
    open(obo, "w").write("[Term]\nid: GO:0097190\nis_a: GO:0006915 ! x\n[Term]\nid: GO:0006915\n")
    c = load_go(gaf, obo_path=obo, aspect="BP", propagate=True)
    assert set(c.gene_sets["GO:0006915"]) == {"A", "B"}


def test_cache_roundtrip():
    t = tempfile.mkdtemp(); p = os.path.join(t, "h.gmt")
    open(p, "w").write("SET_A\td\tG1\tG2\n")
    c = load_gmt(p)
    cache = GeneSetCache(cache_dir=os.path.join(t, "cache"))
    cache.save("k", c)
    back = cache.load("k")
    assert back.checksum == c.checksum and cache.list_keys() == ["k"]


def test_gage_convenience_matches_top_hit():
    from pygage import core
    from pathlib import Path
    reg = Path(core.__file__).parent / "data" / "regression"
    prep = pl.read_csv(reg / "gse16873_prepared.csv.gz", schema_overrides={"gene_id": pl.Utf8})
    gs = json.loads((reg / "kegg_gs.json").read_text())
    tidy = gage(prep, gs, prepared=True)
    top = tidy.filter(pl.col("direction") == "greater").sort("p_val")["gene_set"][0]
    assert top.startswith("hsa04141")
