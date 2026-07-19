"""Gene-set sourcing: GMT (MSigDB), Reactome, GO/OBO + a versioned local cache.

Every loader returns a ``GeneSetCollection`` carrying the sets plus provenance
metadata (source, release, retrieval date, checksum) so a run is reproducible
and the gene-set release can be stamped into results.
"""

from __future__ import annotations

import datetime as _dt
import gzip
import hashlib
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional

logger = logging.getLogger("pygage.gene_sets")


@dataclass
class GeneSetCollection:
    gene_sets: Dict[str, List[str]]
    source: str = "unknown"
    release: str = "unknown"
    retrieved: str = field(default_factory=lambda: _dt.date.today().isoformat())
    n_sets: int = 0
    checksum: str = ""
    extra: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        self.n_sets = len(self.gene_sets)
        if not self.checksum:
            blob = json.dumps({k: sorted(v) for k, v in sorted(self.gene_sets.items())}).encode()
            self.checksum = hashlib.sha256(blob).hexdigest()[:16]

    def metadata(self) -> Dict[str, object]:
        d = asdict(self); d.pop("gene_sets"); return d

    def as_dict(self) -> Dict[str, List[str]]:
        """Return the canonical ``{name: [gene, ...]}`` mapping (engine schema)."""
        return {k: list(v) for k, v in self.gene_sets.items()}

    def filter_size(self, lo: int = 5, hi: int = 1000) -> "GeneSetCollection":
        gs = {k: v for k, v in self.gene_sets.items() if lo <= len(v) <= hi}
        return GeneSetCollection(gs, self.source, self.release, self.retrieved, extra=self.extra)


def _open(path: Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def normalize_gene_sets(gene_sets: "object") -> Dict[str, List[str]]:
    """Coerce any accepted gene-set representation to the canonical schema.

    Accepts a :class:`GeneSetCollection` or a plain ``{name: iterable}`` mapping
    (e.g. parsed ``kegg_gs.json``) and returns ``{name: [gene, ...]}`` with gene
    identifiers stringified and de-duplicated **in order**. This is the single
    internal registry format the analysis engine consumes, so JSON dicts, GMT /
    Reactome / GO collections, and ad-hoc mappings all interoperate without the
    caller thinking about representation.
    """
    if isinstance(gene_sets, GeneSetCollection):
        items = gene_sets.gene_sets.items()
    elif isinstance(gene_sets, Mapping):
        items = gene_sets.items()
    else:
        raise TypeError(
            "gene_sets must be a mapping {name: [genes]} or a GeneSetCollection; "
            f"got {type(gene_sets).__name__}"
        )
    normalized: Dict[str, List[str]] = {}
    for name, genes in items:
        if isinstance(genes, (str, bytes)):
            raise TypeError(f"gene set {name!r} must map to an iterable of genes, not a string")
        normalized[str(name)] = list(dict.fromkeys(str(g) for g in genes))
    return normalized


# --------------------------------------------------------------------------- #
# GMT (MSigDB Hallmark / C2 / C5, and any GMT)
# --------------------------------------------------------------------------- #
def load_gmt(path, source: str = "GMT", release: str = "unknown") -> GeneSetCollection:
    """Parse a GMT file: each line = name<TAB>description<TAB>gene1<TAB>gene2...

    Works directly with MSigDB downloads (hallmark.gmt, c2.cp.reactome.gmt,
    c5.go.bp.gmt, etc.).  ``.gz`` supported.
    """
    gene_sets: Dict[str, List[str]] = {}
    descriptions: Dict[str, str] = {}
    with _open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, desc, *genes = parts
            gene_sets[name] = [g for g in genes if g]
            descriptions[name] = desc
    logger.info("Loaded %d gene sets from GMT %s", len(gene_sets), path)
    return GeneSetCollection(gene_sets, source=source, release=release,
                             extra={"file": str(path)})


def load_msigdb(path, collection: str = "H") -> GeneSetCollection:
    """Load an MSigDB collection GMT and tag the collection code (H/C2/C5...)."""
    c = load_gmt(path, source="MSigDB")
    c.extra["collection"] = collection
    return c


# --------------------------------------------------------------------------- #
# Reactome (ReactomePathways.gmt or NCBI2Reactome mapping)
# --------------------------------------------------------------------------- #
def load_reactome(path, id_type: str = "gmt", species: str = "Homo sapiens") -> GeneSetCollection:
    """Load Reactome gene sets.

    ``id_type='gmt'``: a ReactomePathways.gmt export (name<TAB>id<TAB>genes...).
    ``id_type='ncbi2reactome'``: the NCBI2Reactome_All_Levels.txt mapping
    (gene<TAB>pathwayId<TAB>url<TAB>pathwayName<TAB>evidence<TAB>species);
    rows are grouped into pathwayId -> genes for the requested species.
    """
    path = Path(path)
    if id_type == "gmt":
        c = load_gmt(path, source="Reactome")
        c.extra["species"] = species
        return c
    gene_sets: Dict[str, List[str]] = {}
    names: Dict[str, str] = {}
    with _open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            gene, pid, _url, pname, _evi, sp = p[0], p[1], p[2], p[3], p[4], p[5]
            if species and sp != species:
                continue
            gene_sets.setdefault(pid, []).append(gene)
            names[pid] = pname
    gene_sets = {k: sorted(set(v)) for k, v in gene_sets.items()}
    logger.info("Loaded %d Reactome pathways for %s", len(gene_sets), species)
    return GeneSetCollection(gene_sets, source="Reactome", extra={"species": species,
                             "names": json.dumps(names)[:1] and "attached"})


# --------------------------------------------------------------------------- #
# GO via GAF + OBO (goatools-compatible, but dependency-free parser)
# --------------------------------------------------------------------------- #
def load_go(
    gaf_path,
    obo_path: Optional[str] = None,
    aspect: Optional[str] = None,
    id_field: str = "symbol",
    include_iea: bool = True,
    propagate: bool = False,
) -> GeneSetCollection:
    """Build GO gene sets from a GAF (+ optional OBO for names / propagation).

    ``aspect`` filters to a single domain: 'BP'/'MF'/'CC' (None = all).
    ``propagate`` (needs OBO) rolls annotations up the is_a/part_of DAG so a
    gene annotated to a child term also counts for its ancestors (as most GO
    enrichment tools do).  Without an OBO, only direct annotations are used.
    """
    from .pathway_database_utils import GOGeneSetRetriever
    r = GOGeneSetRetriever()
    res = r.get_go_gene_sets(gaf_path, id_field=id_field, include_iea=include_iea,
                             obo_file=obo_path)
    gene_sets = {k: list(v) for k, v in res["gene_sets"].items()}
    cats = res["go_categories"]

    if aspect:
        keep = set(cats.get(aspect, []))
        gene_sets = {k: v for k, v in gene_sets.items() if k in keep}

    if propagate and obo_path:
        parents = _parse_obo_isa(obo_path)
        gene_sets = _propagate(gene_sets, parents)

    c = GeneSetCollection(gene_sets, source="GO",
                          extra={"aspect": aspect or "all",
                                 "propagated": str(bool(propagate and obo_path))})
    c.extra["names_available"] = str(bool(obo_path))
    return c


def _parse_obo_isa(obo_path) -> Dict[str, List[str]]:
    parents: Dict[str, List[str]] = {}
    cur = None
    with _open(obo_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = None
            elif line.startswith("id: GO:"):
                cur = line[4:].strip()
            elif cur and line.startswith("is_a:"):
                parents.setdefault(cur, []).append(line.split("!")[0].replace("is_a:", "").strip())
            elif cur and line.startswith("relationship: part_of"):
                parents.setdefault(cur, []).append(line.split("part_of")[1].split("!")[0].strip())
    return parents


def _propagate(gene_sets: Dict[str, List[str]], parents: Dict[str, List[str]]) -> Dict[str, List[str]]:
    def ancestors(go):
        seen, stack = set(), list(parents.get(go, []))
        while stack:
            p = stack.pop()
            if p in seen:
                continue
            seen.add(p)
            stack.extend(parents.get(p, []))
        return seen
    out: Dict[str, set] = {k: set(v) for k, v in gene_sets.items()}
    for go, genes in list(gene_sets.items()):
        for anc in ancestors(go):
            out.setdefault(anc, set()).update(genes)
    return {k: sorted(v) for k, v in out.items()}


# --------------------------------------------------------------------------- #
# Versioned local cache
# --------------------------------------------------------------------------- #
class GeneSetCache:
    """On-disk cache of GeneSetCollections keyed by name; offline-reproducible."""

    def __init__(self, cache_dir="~/.cache/pygage/gene_sets"):
        self.dir = Path(cache_dir).expanduser()
        self.dir.mkdir(parents=True, exist_ok=True)

    def _p(self, key: str) -> Path:
        safe = key.replace("/", "_")
        return self.dir / f"{safe}.json.gz"

    def save(self, key: str, coll: GeneSetCollection) -> Path:
        payload = {"metadata": coll.metadata(), "gene_sets": coll.gene_sets}
        p = self._p(key)
        with gzip.open(p, "wt") as fh:
            json.dump(payload, fh)
        logger.info("Cached %s (%d sets, checksum %s) -> %s", key, coll.n_sets, coll.checksum, p)
        return p

    def load(self, key: str) -> Optional[GeneSetCollection]:
        p = self._p(key)
        if not p.exists():
            return None
        with gzip.open(p, "rt") as fh:
            payload = json.load(fh)
        md = payload["metadata"]
        return GeneSetCollection(payload["gene_sets"], source=md.get("source", "cache"),
                                 release=md.get("release", "unknown"),
                                 retrieved=md.get("retrieved", ""), extra=md.get("extra", {}))

    def list_keys(self) -> List[str]:
        return sorted(p.name.replace(".json.gz", "") for p in self.dir.glob("*.json.gz"))
