"""Runtime configuration and data resolution for PyGAGE.

The package ships only a *minimal* set of static assets (the Entrez↔symbol map
and small regression fixtures). Everything that changes over time — KEGG
pathways, GO terms, Reactome, MSigDB — is fetched live (see
:mod:`pygage.pathway_database_utils` and :mod:`pygage.gene_sets`) or supplied by
the user, so gene-set updates never require a package version bump.

Data resolution order for a named asset:

1. an explicit path passed to the calling function/constructor;
2. the ``PYGAGE_DATA_DIR`` environment variable (a directory of overrides);
3. the asset bundled inside the installed package.

This lets deployments point PyGAGE at a curated/updated data directory without
touching the code or reinstalling.
"""

from __future__ import annotations

import logging
import os
from importlib import resources
from pathlib import Path
from typing import Optional

logger = logging.getLogger("pygage.config")

ENV_DATA_DIR = "PYGAGE_DATA_DIR"
ENV_NUM_THREADS = "PYGAGE_NUM_THREADS"


def packaged_data_dir() -> Path:
    """Return the data directory bundled inside the installed package."""
    return Path(str(resources.files("pygage") / "data"))


def data_dir() -> Path:
    """Return the active data directory (``PYGAGE_DATA_DIR`` override or packaged)."""
    override = os.environ.get(ENV_DATA_DIR)
    if override:
        p = Path(override).expanduser()
        if p.is_dir():
            return p
    return packaged_data_dir()


def resolve(name: str, explicit: Optional[str] = None) -> Path:
    """Resolve a data asset by name, honouring explicit > env override > packaged.

    Parameters
    ----------
    name:
        File name relative to the data directory, e.g. ``"egSymb.tsv"``.
    explicit:
        An explicit path supplied by the caller; used verbatim if given.
    """
    if explicit:
        return Path(explicit).expanduser()
    override_dir = os.environ.get(ENV_DATA_DIR)
    if override_dir:
        candidate = Path(override_dir).expanduser() / name
        if candidate.exists():
            return candidate
    return packaged_data_dir() / name


def egsymb_path(explicit: Optional[str] = None) -> Path:
    """Path to the Entrez↔symbol mapping (override-aware)."""
    return resolve("egSymb.tsv", explicit)


# --------------------------------------------------------------------------- #
# Compute / threading controls (for shared & HPC environments)
# --------------------------------------------------------------------------- #
def default_n_jobs() -> int:
    """Default worker count for gene-set parallelism.

    Resolved from the ``PYGAGE_NUM_THREADS`` environment variable, else ``1``
    (serial). Pass ``n_jobs=`` explicitly to override per call.
    """
    raw = os.environ.get(ENV_NUM_THREADS)
    if raw:
        try:
            return max(1, int(raw))
        except ValueError:
            logger.warning("Ignoring non-integer %s=%r", ENV_NUM_THREADS, raw)
    return 1


def set_thread_limits(n_threads: int) -> None:
    """Cap the Polars / BLAS thread pools for shared or HPC environments.

    Sets ``POLARS_MAX_THREADS`` and the common BLAS variables (``OMP``,
    ``OPENBLAS``, ``MKL``, ``NUMEXPR``) plus ``PYGAGE_NUM_THREADS``. To size the
    **Polars** pool this must run *before* Polars is first used (ideally before
    import); a warning is logged if the pool is already initialised at a
    different size. BLAS limits and :func:`default_n_jobs` take effect
    immediately.
    """
    n = max(1, int(n_threads))
    for var in (
        "POLARS_MAX_THREADS", ENV_NUM_THREADS, "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        os.environ[var] = str(n)
    try:
        import polars as pl

        if pl.thread_pool_size() != n:
            logger.warning(
                "Polars thread pool already initialised at %d; set "
                "POLARS_MAX_THREADS/PYGAGE_NUM_THREADS before importing polars "
                "to resize it (BLAS limits and n_jobs still applied).",
                pl.thread_pool_size(),
            )
    except Exception:  # pragma: no cover - polars always present in practice
        pass


def thread_config() -> "dict":
    """Report the active thread configuration (Polars pool + relevant env vars)."""
    cfg = {
        ENV_NUM_THREADS: os.environ.get(ENV_NUM_THREADS),
        "POLARS_MAX_THREADS": os.environ.get("POLARS_MAX_THREADS"),
        "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
    }
    try:
        import polars as pl

        cfg["polars_pool_size"] = pl.thread_pool_size()
    except Exception:  # pragma: no cover
        pass
    return cfg


# --------------------------------------------------------------------------- #
# Live-fetchable reference maps (packaged file is a cached fallback)
# --------------------------------------------------------------------------- #
def ensure_egsymb(
    explicit: Optional[str] = None,
    fetch: bool = False,
    species: str = "hsa",
) -> Path:
    """Return a usable Entrez↔symbol map, optionally building it live.

    Resolution order is the usual explicit → ``PYGAGE_DATA_DIR`` → packaged
    asset. The packaged ``egSymb.tsv`` is a *cached convenience*: with
    ``fetch=True`` (and no local copy found) the map is rebuilt from the KEGG
    REST endpoint via :func:`pygage.gene_id_utils.build_kegg_symbol_map` and
    cached under the active data directory, so the identifier map can be
    refreshed without a package release.
    """
    path = egsymb_path(explicit)
    if path.exists() and not fetch:
        return path
    if fetch:
        from .gene_id_utils import build_kegg_symbol_map

        target = data_dir() / "egSymb.tsv"
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            build_kegg_symbol_map(species=species, output_file=target)
            logger.info("Fetched Entrez↔symbol map from KEGG (%s) -> %s", species, target)
            return target
        except Exception as exc:  # fall back to packaged copy
            logger.warning("Live egSymb fetch failed (%s); using packaged copy", exc)
    return path
