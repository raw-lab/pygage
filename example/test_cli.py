"""CLI subcommand smoke tests (run / go / compare)."""
import gzip, json, os, tempfile
from pathlib import Path
import polars as pl
from pygage import core
from pygage.cli import main

REG = Path(core.__file__).parent / "data" / "regression"


def _fixture_matrix(tmp):
    prep = pl.read_csv(REG / "gse16873_prepared.csv.gz", schema_overrides={"gene_id": pl.Utf8})
    mp = os.path.join(tmp, "m.csv"); prep.write_csv(mp)
    gsp = os.path.join(tmp, "gs.json")
    json.dump({"gene_sets": json.loads((REG / "kegg_gs.json").read_text())}, open(gsp, "w"))
    return mp, gsp


def test_cli_run_prepared(tmp_path):
    tmp = str(tmp_path); mp, gsp = _fixture_matrix(tmp)
    out = os.path.join(tmp, "o.csv")
    assert main(["run", mp, "-g", gsp, "-o", out, "--prepared"]) == 0
    res = pl.read_csv(out)
    assert res.height > 0 and "p_val" in res.columns


def test_cli_go(tmp_path):
    tmp = str(tmp_path)
    gaf = os.path.join(tmp, "t.gaf")
    open(gaf, "w").write("\n".join([
        "!gaf",
        "DB\tP1\tA\t\tGO:0006915\tPMID\tIDA\t\tP\tn\t\tprotein\ttaxon:9606\t2020\tX",
        "DB\tP2\tB\t\tGO:0006915\tPMID\tIDA\t\tP\tn\t\tprotein\ttaxon:9606\t2020\tX",
    ]))
    out = os.path.join(tmp, "go.json")
    assert main(["go", gaf, "-o", out]) == 0
    assert "GO:0006915" in json.load(open(out))["gene_sets"]


def test_cli_compare(tmp_path):
    tmp = str(tmp_path); mp, gsp = _fixture_matrix(tmp)
    o1 = os.path.join(tmp, "a.csv"); o2 = os.path.join(tmp, "b.csv")
    main(["run", mp, "-g", gsp, "-o", o1, "--prepared"])
    main(["run", mp, "-g", gsp, "-o", o2, "--prepared", "--test", "z-test"])
    cmp = os.path.join(tmp, "c.tsv")
    assert main(["compare", o1, o2, "--names", "A,B", "-o", cmp]) == 0
    assert Path(cmp).exists()
