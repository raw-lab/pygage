#!/usr/bin/env python3
"""
Gene ID Conversion Utilities (corrected).

Fixes vs. the original draft
----------------------------
* ``eg2sym`` / ``sym2eg`` now return what the signature promises -- a list of
  optional strings aligned to the input -- with an opt-in ``as_frame=True`` for
  the two-column DataFrame.  (The old code always returned a DataFrame,
  breaking any caller that indexed the result as a list.)
* A missing bundled mapping file no longer raises on construction; the object
  is created empty and a clear error is raised only if a conversion is
  attempted without a mapping.  A helper can build the Entrez<->symbol map
  from KEGG for any species so the tool is not human-only.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import polars as pl


class GeneIDConverter:
    """Convert between gene IDs and symbols using a two-column mapping."""

    def __init__(self, mapping_file: Optional[Path] = None):
        self.mapping_df: Optional[pl.DataFrame] = None
        if mapping_file is None:
            # Honour an explicit path, then PYGAGE_DATA_DIR, then the packaged map.
            from . import config
            default = config.egsymb_path()
            mapping_file = default if default.exists() else None
        if mapping_file is not None:
            self.load_mapping(Path(mapping_file))

    def load_mapping(self, mapping_file: Path) -> None:
        mapping_file = Path(mapping_file)
        if not mapping_file.exists():
            raise FileNotFoundError(f"Mapping file not found: {mapping_file}")
        sep = "," if mapping_file.suffix == ".csv" else "\t"
        df = pl.read_csv(mapping_file, separator=sep, infer_schema=False, has_header=True)
        if df.width < 2:
            raise ValueError("Mapping file needs >= 2 columns (id, symbol)")
        self.mapping_df = df.select(
            pl.col(df.columns[0]).cast(pl.Utf8).alias("entrez_id"),
            pl.col(df.columns[1]).cast(pl.Utf8).alias("symbol"),
        ).unique(subset=["entrez_id"], keep="first")

    def _require(self) -> None:
        if self.mapping_df is None:
            raise ValueError(
                "No mapping loaded. Pass mapping_file=... (e.g. GAGE's egSymb for "
                "human, or build one with build_kegg_symbol_map for any species)."
            )

    def eg2sym(
        self, entrez_ids: Union[List, pl.Series], as_frame: bool = False
    ) -> Union[List[Optional[str]], pl.DataFrame]:
        self._require()
        ids = entrez_ids.to_list() if isinstance(entrez_ids, pl.Series) else list(entrez_ids)
        lookup = pl.DataFrame({"entrez_id": [str(x) for x in ids]})
        joined = lookup.join(self.mapping_df, on="entrez_id", how="left")
        if as_frame:
            return joined.rename({"entrez_id": "input", "symbol": "output"})
        return joined["symbol"].to_list()

    def sym2eg(
        self, symbols: Union[List, pl.Series], as_frame: bool = False
    ) -> Union[List[Optional[str]], pl.DataFrame]:
        self._require()
        syms = symbols.to_list() if isinstance(symbols, pl.Series) else list(symbols)
        lookup = pl.DataFrame({"symbol": [str(x) for x in syms]})
        rev = self.mapping_df.unique(subset=["symbol"], keep="first")
        joined = lookup.join(rev, on="symbol", how="left")
        if as_frame:
            return joined.rename({"symbol": "input", "entrez_id": "output"})
        return joined["entrez_id"].to_list()


def build_kegg_symbol_map(species: str = "hsa", output_file: Optional[Path] = None) -> pl.DataFrame:
    """Build an Entrez<->symbol map for any KEGG species and optionally save it.

    Uses ``/conv/<org>/ncbi-geneid`` (Entrez) plus ``/list/<org>`` (symbols),
    so the converter is not restricted to human.
    """
    from .pathway_database_utils import KEGGPathwayRetriever

    kg = KEGGPathwayRetriever()
    conv = kg._get(f"conv/{species}/ncbi-geneid")          # entrez <-> kegg gene
    kegg_to_entrez = {}
    for line in conv.strip().split("\n"):
        if not line:
            continue
        a, b = (line.split("\t") + [""])[:2]
        if a.startswith("ncbi-geneid:"):
            entrez, kegg = a.split(":", 1)[1], b.split(":", 1)[1]
        else:
            kegg, entrez = a.split(":", 1)[1], b.split(":", 1)[1]
        kegg_to_entrez[kegg] = entrez

    listing = kg._get(f"list/{species}")                   # kegg gene -> symbol;desc
    rows = []
    for line in listing.strip().split("\n"):
        p = line.split("\t")
        if len(p) < 2:
            continue
        kegg = p[0].split(":", 1)[1] if ":" in p[0] else p[0]
        symbol = p[-1].split(";")[0].split(",")[0].strip()
        entrez = kegg_to_entrez.get(kegg)
        if entrez:
            rows.append({"entrez_id": entrez, "symbol": symbol})

    df = pl.DataFrame(rows).unique(subset=["entrez_id"], keep="first")
    if output_file is not None:
        df.write_csv(output_file, separator="\t")
    return df
