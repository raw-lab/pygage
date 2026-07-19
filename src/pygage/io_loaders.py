"""Input adapters: raw matrices, DE tables, pre-ranked vectors, pandas/AnnData.

Most users arrive with a DESeq2/edgeR/limma result table or a single ranked
statistic, not a raw expression matrix.  These helpers turn any of those into a
gene_col + fold-change/statistic frame that ``GAGEAnalysis.run_gage`` consumes.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Sequence, Union

if TYPE_CHECKING:
    from ._types import FrameLike, GeneSetsLike, PathLike

import numpy as np
import polars as pl

logger = logging.getLogger("pygage.io")

# common column aliases seen in DE tool outputs
_GENE_ALIASES = ("gene", "gene_id", "geneid", "id", "symbol", "row", "", "feature", "ensembl")
_LFC_ALIASES = ("log2foldchange", "log2fc", "logfc", "lfc", "coef", "fold_change", "l2fc")
_STAT_ALIASES = ("stat", "t", "statistic", "wald", "z", "lr")


def _read_csv_str_gene(path: Path, gene_aliases=_GENE_ALIASES) -> pl.DataFrame:
    """Read a CSV/TSV forcing the detected gene column to Utf8 (mixed Entrez/AFFX safe)."""
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    with open(path) as fh:
        header_line = fh.readline().rstrip("\n")
    header = [h.strip().strip('"') for h in header_line.split(sep)]
    low = {c.lower().replace(" ", "").replace("_", ""): c for c in header}
    gc = _find(low, [a.replace("_", "") for a in gene_aliases]) or header[0]
    return pl.read_csv(path, separator=sep, schema_overrides={gc: pl.Utf8},
                       infer_schema_length=None)


def _to_polars(data) -> pl.DataFrame:
    """Accept polars, pandas, numpy-record, dict, or AnnData-like input."""
    if isinstance(data, pl.DataFrame):
        return data
    # pandas
    try:
        import pandas as pd
        if isinstance(data, pd.DataFrame):
            return pl.from_pandas(data.reset_index() if data.index.name else data)
    except ImportError:
        pass
    # AnnData: use .var/.obs-aware matrix (genes as var)
    try:
        import anndata as ad
        if isinstance(data, ad.AnnData):
            X = np.asarray(data.X.todense()) if hasattr(data.X, "todense") else np.asarray(data.X)
            # AnnData is obs x var (samples x genes); transpose to genes x samples
            mat = X.T
            genes = list(data.var_names)
            samples = list(data.obs_names)
            out = {"gene_id": genes}
            for j, sname in enumerate(samples):
                out[str(sname)] = mat[:, j]
            return pl.DataFrame(out)
    except ImportError:
        pass
    if isinstance(data, dict):
        return pl.DataFrame(data)
    raise TypeError(f"Unsupported input type {type(data)!r}; pass polars/pandas/AnnData/dict")


def _find(cols_lower: Dict[str, str], aliases: Sequence[str]) -> Optional[str]:
    for a in aliases:
        if a in cols_lower:
            return cols_lower[a]
    # loose contains-match
    for low, orig in cols_lower.items():
        if any(a and a in low for a in aliases):
            return orig
    return None


def read_matrix(
    source: Union[str, Path, "pl.DataFrame"],
    gene_col: str = "gene_id",
) -> pl.DataFrame:
    """Load a raw expression matrix (genes x samples) from a file or frame."""
    if isinstance(source, (str, Path)):
        df = _read_csv_str_gene(Path(source))
    else:
        df = _to_polars(source)
    if gene_col not in df.columns:
        df = df.rename({df.columns[0]: gene_col})
    return df.with_columns(pl.col(gene_col).cast(pl.Utf8))


def read_de_table(
    source: Union[str, Path, "pl.DataFrame"],
    gene_col: Optional[str] = None,
    value: str = "log2FC",
    stat_col: Optional[str] = None,
    lfc_col: Optional[str] = None,
) -> pl.DataFrame:
    """Load a DESeq2/edgeR/limma-style DE table into a single-column ranked frame.

    Auto-detects gene, log2FoldChange, and stat columns by common aliases.
    ``value`` selects which per-gene statistic to rank on: ``'log2FC'`` (default,
    matches gage use.fold=TRUE) or ``'stat'`` (the Wald/t/LR statistic).
    Returns a frame with ``gene_id`` + one value column ready for run_gage.
    """
    if isinstance(source, (str, Path)):
        df = _read_csv_str_gene(Path(source))
    else:
        df = _to_polars(source)

    cols_lower = {c.lower().replace(" ", "").replace("_", ""): c for c in df.columns}
    gc = gene_col or _find(cols_lower, [a.replace("_", "") for a in _GENE_ALIASES]) or df.columns[0]
    lc = lfc_col or _find(cols_lower, _LFC_ALIASES)
    sc = stat_col or _find(cols_lower, _STAT_ALIASES)

    pick = lc if value == "log2FC" else sc
    if pick is None:
        pick = lc or sc
        if pick is None:
            raise ValueError(
                f"Could not find a log2FC or stat column in {df.columns}. "
                "Pass lfc_col=/stat_col= explicitly."
            )
        logger.warning("value=%r column not found; using %r", value, pick)

    out = df.select(
        pl.col(gc).cast(pl.Utf8).alias("gene_id"),
        pl.col(pick).cast(pl.Float64).alias(value if value in ("log2FC", "stat") else "stat"),
    ).drop_nulls()
    logger.info("Loaded DE table: %d genes, ranking on %r", out.height, pick)
    return out


def read_preranked(
    source: Union[str, Path, Dict[str, float], "pl.DataFrame"],
    gene_col: str = "gene_id",
    score_col: str = "score",
) -> pl.DataFrame:
    """Load a pre-ranked gene->score vector (fgsea-style) into a run_gage frame."""
    if isinstance(source, dict):
        df = pl.DataFrame({"gene_id": list(source.keys()),
                           "score": list(source.values())})
        return df.with_columns(pl.col("gene_id").cast(pl.Utf8))
    if isinstance(source, (str, Path)):
        df = _read_csv_str_gene(Path(source))
    else:
        df = _to_polars(source)
    g = gene_col if gene_col in df.columns else df.columns[0]
    sc = score_col if score_col in df.columns else df.columns[1]
    return df.select(pl.col(g).cast(pl.Utf8).alias("gene_id"),
                     pl.col(sc).cast(pl.Float64).alias("score"))


def gage(
    data: "FrameLike | PathLike",
    gene_sets: "GeneSetsLike",
    ref_indices: Optional[Sequence[int]] = None,
    samp_indices: Optional[Sequence[int]] = None,
    gene_col: str = "gene_id",
    comparison: str = "paired",
    same_dir: bool = True,
    test_method: str = "t-test",
    meta_method: str = "stouffer",
    set_size_range=(10, 500),
    input_logged: bool = True,
    prepared: bool = False,
    tidy: bool = True,
    **run_kwargs,
):
    """One-call GAGE for a raw matrix, a prepared frame, or a DE/pre-ranked table.

    - raw matrix + ref_indices/samp_indices -> gagePrep then run_gage
    - a single-column frame (DE/pre-ranked, or ``prepared=True``) -> run_gage directly
    - ``data`` may be polars, pandas, AnnData, dict, or a file path.

    Returns a tidy polars frame (direction-labelled) by default, else the
    greater/less/stats dict.
    """
    from .core import GAGEAnalysis, GAGEPreparation

    if isinstance(data, (str, Path)):
        frame = read_matrix(data, gene_col=gene_col)
    else:
        frame = _to_polars(data)
        if gene_col not in frame.columns:
            frame = frame.rename({frame.columns[0]: gene_col})
        frame = frame.with_columns(pl.col(gene_col).cast(pl.Utf8))

    if ref_indices is not None and not prepared:
        frame = GAGEPreparation.prepare_expression(
            frame, ref_indices=ref_indices, samp_indices=samp_indices,
            gene_col=gene_col, comparison=comparison, same_dir=same_dir,
            input_logged=input_logged,
        )
    # else: already prepared / single-column DE or pre-ranked -> use as-is

    ga = GAGEAnalysis()
    res = ga.run_gage(
        frame, gene_sets, gene_col=gene_col, set_size_range=set_size_range,
        same_dir=same_dir, test_method=test_method, meta_method=meta_method,
        **run_kwargs,
    )
    if not tidy:
        return res
    parts = []
    for direction in ("greater", "less"):
        if direction in res:
            parts.append(res[direction].with_columns(pl.lit(direction).alias("direction")))
    return pl.concat(parts, how="diagonal").sort(["direction", "p_val"]) if parts else res["greater"]
