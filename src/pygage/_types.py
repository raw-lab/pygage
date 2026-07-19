"""Shared type aliases for PyGAGE's public APIs (PEP 484).

Centralising these keeps loader, engine, and plotting signatures consistent and
gives downstream pipelines precise, importable types to build against.
"""

from __future__ import annotations

from os import PathLike as _OsPathLike
from typing import TYPE_CHECKING, Mapping, Sequence, TypeVar, Union

import polars as pl

if TYPE_CHECKING:  # imported only for type checkers, never at runtime
    import pandas as pd
    from anndata import AnnData

    from .gene_sets import GeneSetCollection

#: A filesystem path accepted by the loaders.
PathLike = Union[str, "_OsPathLike[str]"]

#: Any table-like input accepted by :func:`pygage.io_loaders.gage` and the
#: ``read_*`` loaders: a Polars/pandas frame, an AnnData object, or a mapping.
FrameLike = Union[pl.DataFrame, "pd.DataFrame", "AnnData", Mapping[str, Sequence]]

#: A gene-set mapping ``{set_name: [gene, ...]}``.
GeneSetMapping = Mapping[str, Sequence[str]]

#: Gene sets accepted by the engine: a plain mapping or a
#: :class:`pygage.gene_sets.GeneSetCollection` (both normalise to the same
#: canonical dict internally).
GeneSetsLike = Union[GeneSetMapping, "GeneSetCollection"]

#: Result frame type variable (preserves the concrete frame type through helpers).
FrameT = TypeVar("FrameT", bound=pl.DataFrame)

__all__ = [
    "PathLike", "FrameLike", "GeneSetMapping", "GeneSetsLike", "FrameT",
]
