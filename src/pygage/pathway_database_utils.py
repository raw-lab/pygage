#!/usr/bin/env python3
"""
Pathway and Gene Set Database Utilities (corrected + extended).

Fixes vs. the original draft
----------------------------
* KEGG links are fetched with a single bulk ``/link/<org>/pathway`` call
  instead of one HTTP request per pathway (the old ``pathway:<org>...`` URL
  used the wrong entry prefix -- it must be ``path:`` -- so it returned
  nothing, and the N+1 request pattern hit KEGG's rate limits).
* Pathway-list and link parsing strip the ``path:``/``<org>:`` prefixes
  robustly, tolerating both the modern (``hsa00010``) and legacy
  (``path:hsa00010``) formats.
* Entrez conversion tolerates either column order and is optional/robust.
* GO/GAF parsing uses the Aspect column (P/F/C -> BP/MF/CC) for real
  categorisation and honours the ``NOT`` qualifier; it no longer mislabels the
  gene's name as the GO-term name.

Extensions for non-human / non-model use
----------------------------------------
* :meth:`KEGGPathwayRetriever.get_pathway_genes` works for any KEGG organism
  code (mmu, eco, ath, dme, sce, ...), not just ``hsa``.
* :meth:`KEGGPathwayRetriever.get_ko_gene_sets` builds *species-agnostic*
  gene sets keyed by KEGG Orthology (KO) IDs -- the right tool for
  metagenomes, viromes, phages and non-model organisms that lack a dedicated
  KEGG genome.  Annotate your genes with KO IDs (KofamScan / MetaCerberus /
  eggNOG-mapper) and feed those to GAGE.
* :meth:`KEGGPathwayRetriever.get_module_gene_sets` exposes KEGG MODULE gene
  sets (organism-specific or KO-based).
"""

from __future__ import annotations

import datetime as _dt
import json
import logging
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import requests

logger = logging.getLogger("pygage.pathway")


class KEGGPathwayRetriever:
    """Retrieve KEGG pathway / module / KO gene sets for any organism."""

    KEGG_REST_BASE = "https://rest.kegg.jp"

    def __init__(self, cache_dir: Optional[Path] = None, retries: int = 3, pause: float = 0.34):
        self.session = requests.Session()
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.retries = retries
        self.pause = pause          # be polite: KEGG asks for < ~3 req/s
        self.species_info: Optional[Dict[str, str]] = None

    # ------------------------------------------------------------------ #
    # low-level GET with retry + optional on-disk cache                  #
    # ------------------------------------------------------------------ #
    def _get(self, path: str) -> str:
        if self.cache_dir:
            key = path.strip("/").replace("/", "__").replace(":", "_") + ".txt"
            cached = self.cache_dir / key
            if cached.exists():
                return cached.read_text()

        url = f"{self.KEGG_REST_BASE}/{path.lstrip('/')}"
        last_exc: Optional[Exception] = None
        for attempt in range(self.retries):
            try:
                r = self.session.get(url, timeout=60)
                if r.status_code == 200:
                    if self.cache_dir:
                        (self.cache_dir / key).write_text(r.text)
                    time.sleep(self.pause)
                    return r.text
                if r.status_code == 404:
                    return ""       # empty (e.g. no such link) -> caller handles
                last_exc = RuntimeError(f"KEGG {url} -> HTTP {r.status_code}")
            except requests.RequestException as exc:      # network hiccup
                last_exc = exc
            time.sleep(self.pause * (attempt + 1))
        raise RuntimeError(f"KEGG request failed for {url}: {last_exc}")

    # ------------------------------------------------------------------ #
    # helpers                                                            #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _strip(entry: str) -> str:
        """Remove leading db/org prefixes: 'path:hsa00010' -> 'hsa00010'."""
        return entry.split(":", 1)[1] if ":" in entry else entry

    def list_organisms(self) -> Dict[str, str]:
        """Return {org_code: definition} for every KEGG organism."""
        text = self._get("list/organism")
        out: Dict[str, str] = {}
        for line in text.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) >= 3:
                out[p[1]] = p[2]      # p[0]=T-number, p[1]=code, p[2]=name
        return out

    def get_species_code(self, species: str = "hsa") -> Dict[str, str]:
        """Validate a KEGG organism code and return its record."""
        text = self._get("list/organism")
        for line in text.strip().split("\n"):
            p = line.split("\t")
            if len(p) >= 2 and p[1] == species:
                self.species_info = {
                    "t_number": p[0],
                    "code": p[1],
                    "name": p[2] if len(p) > 2 else "",
                    "taxonomy": p[3] if len(p) > 3 else "",
                }
                return self.species_info
        raise ValueError(
            f"Unknown KEGG organism code: {species!r}. "
            f"Use list_organisms() to browse valid codes, or get_ko_gene_sets() "
            f"for organisms without a KEGG genome."
        )

    # ------------------------------------------------------------------ #
    # organism-specific pathway gene sets                                #
    # ------------------------------------------------------------------ #
    def get_pathway_names(self, species: str = "hsa") -> Dict[str, str]:
        text = self._get(f"list/pathway/{species}")
        names: Dict[str, str] = {}
        for line in text.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) >= 2:
                names[self._strip(p[0])] = p[1]
        return names

    def get_pathway_genes(
        self, species: str = "hsa", id_type: str = "kegg"
    ) -> Dict[str, object]:
        """Return pathway gene sets for a KEGG organism (any species).

        Uses one bulk link call for all gene<->pathway edges.

        id_type: ``kegg`` (native ``<org>:gene`` IDs, prefix stripped) or
        ``entrez`` (NCBI Gene IDs via ``/conv``).
        """
        self.get_species_code(species)
        pathway_names = self.get_pathway_names(species)

        link_text = self._get(f"link/{species}/pathway")   # bulk: all edges
        gene_sets: Dict[str, List[str]] = defaultdict(list)
        for line in link_text.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 2:
                continue
            # order is path:<org>#####  <TAB>  <org>:gene  (bulk link direction)
            a, b = p[0], p[1]
            if a.startswith("path:"):
                pid, gene = self._strip(a), self._strip(b)
            else:
                pid, gene = self._strip(b), self._strip(a)
            gene_sets[pid].append(gene)

        gene_sets = {k: sorted(set(v)) for k, v in gene_sets.items()}

        if id_type == "entrez":
            gene_sets = self._convert_to_entrez(species, gene_sets)

        return {
            "gene_sets": gene_sets,
            "pathway_names": pathway_names,
            "categories": self.categorize_pathways(list(pathway_names.keys())),
        }

    # ------------------------------------------------------------------ #
    # KO-based (species-agnostic) gene sets -- key to "beyond human"     #
    # ------------------------------------------------------------------ #
    def get_ko_gene_sets(self, reference: str = "pathway") -> Dict[str, object]:
        """Build KO-keyed gene sets from KEGG reference pathways or modules.

        These sets are organism-independent: each set lists the KO IDs
        (e.g. ``K00844``) assigned to a reference map (``map#####``) or module.
        Annotate query genes with KO IDs (KofamScan / MetaCerberus / eggNOG)
        and run GAGE directly against these -- ideal for metagenomes, viromes,
        phages and non-model organisms.

        reference: ``pathway`` or ``module``.
        """
        target = "pathway" if reference == "pathway" else "module"
        link_text = self._get(f"link/ko/{target}")
        gene_sets: Dict[str, List[str]] = defaultdict(list)
        for line in link_text.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 2:
                continue
            a, b = p[0], p[1]
            if a.startswith("ko:") or a.startswith("K"):
                ko, mid = self._strip(a), self._strip(b)
            else:
                ko, mid = self._strip(b), self._strip(a)
            gene_sets[mid].append(ko)
        gene_sets = {k: sorted(set(v)) for k, v in gene_sets.items()}

        names_text = self._get(f"list/{target}")
        names: Dict[str, str] = {}
        for line in names_text.strip().split("\n"):
            p = line.split("\t")
            if len(p) >= 2:
                names[self._strip(p[0])] = p[1]

        cats = (
            self.categorize_pathways(list(names.keys()))
            if reference == "pathway"
            else {}
        )
        logger.info("Built %d KO gene sets from KEGG %s reference", len(gene_sets), target)
        return {"gene_sets": gene_sets, "pathway_names": names, "categories": cats}

    def list_all_kos(self) -> Dict[str, str]:
        """Return **every** KEGG Orthology entry: {KO_id: description}.

        e.g. ``{'K00001': 'E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]', ...}``.
        This is the full KO namespace (~26k entries) used to annotate genes for
        species-agnostic analysis.
        """
        text = self._get("list/ko")
        out: Dict[str, str] = {}
        for line in text.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t", 1)
            if len(p) >= 2:
                out[self._strip(p[0])] = p[1]
        logger.info("Listed %d KEGG Orthology (KO) entries", len(out))
        return out

    def download_ko_gene_sets(
        self,
        output_file: Path,
        reference: str = "pathway",
        id_field: str = "ko",
    ) -> Dict[str, object]:
        """Download KO gene sets and write them to ``output_file`` (JSON).

        Bundles the KO gene sets, their names, categories, and the full KO list
        with provenance (source + retrieval date) so the download is reproducible
        and can be re-loaded offline.  Progress is logged.
        """
        logger.info("Downloading KO gene sets (reference=%s) ...", reference)
        res = self.get_ko_gene_sets(reference=reference)
        logger.info("Fetching full KO namespace for annotation ...")
        ko_list = self.list_all_kos()
        payload = {
            "gene_sets": res["gene_sets"],
            "set_names": res.get("pathway_names", {}),
            "categories": res.get("categories", {}),
            "ko_catalog": ko_list,
            "provenance": {
                "source": "KEGG REST (rest.kegg.jp)",
                "reference": reference,
                "n_sets": len(res["gene_sets"]),
                "n_kos_in_sets": len({k for v in res["gene_sets"].values() for k in v}),
                "n_kos_total": len(ko_list),
                "retrieved": _dt.date.today().isoformat(),
                "note": "KEGG is academic-use only; cite Kanehisa et al.",
            },
        }
        Path(output_file).write_text(json.dumps(payload))
        logger.info("Wrote %d KO gene sets (%d KOs used of %d total) -> %s",
                    payload["provenance"]["n_sets"],
                    payload["provenance"]["n_kos_in_sets"],
                    payload["provenance"]["n_kos_total"], output_file)
        return payload

    def get_module_gene_sets(self, species: str = "hsa") -> Dict[str, object]:
        """KEGG MODULE gene sets for a given organism."""
        self.get_species_code(species)
        link_text = self._get(f"link/{species}/module")
        gene_sets: Dict[str, List[str]] = defaultdict(list)
        for line in link_text.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 2:
                continue
            a, b = p[0], p[1]
            if a.startswith("md:") or "_M" in a:
                mid, gene = self._strip(a), self._strip(b)
            else:
                mid, gene = self._strip(b), self._strip(a)
            gene_sets[mid].append(gene)
        gene_sets = {k: sorted(set(v)) for k, v in gene_sets.items()}
        names_text = self._get("list/module")
        names = {}
        for line in names_text.strip().split("\n"):
            p = line.split("\t")
            if len(p) >= 2:
                names[self._strip(p[0])] = p[1]
        return {"gene_sets": gene_sets, "module_names": names}

    # ------------------------------------------------------------------ #
    # ID conversion                                                      #
    # ------------------------------------------------------------------ #
    def _convert_to_entrez(
        self, species: str, gene_sets: Dict[str, List[str]]
    ) -> Dict[str, List[str]]:
        """Map native KEGG gene IDs to NCBI Entrez Gene IDs (order-tolerant)."""
        try:
            text = self._get(f"conv/{species}/ncbi-geneid")
        except RuntimeError:
            logger.warning("could not fetch KEGG<->Entrez conversion; keeping KEGG IDs")
            return gene_sets

        conv: Dict[str, str] = {}
        for line in text.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 2:
                continue
            # one side is 'ncbi-geneid:1234', the other '<org>:gene' -- detect it
            left, right = p[0], p[1]
            if left.startswith("ncbi-geneid:"):
                entrez, kegg = self._strip(left), self._strip(right)
            else:
                kegg, entrez = self._strip(left), self._strip(right)
            conv[kegg] = entrez

        return {
            pid: [conv.get(g, g) for g in genes] for pid, genes in gene_sets.items()
        }

    # ------------------------------------------------------------------ #
    # BRITE-based pathway categorisation (GAGE-style), numeric fallback  #
    # ------------------------------------------------------------------ #
    def categorize_pathways(self, pathway_ids: List[str]) -> Dict[str, List[str]]:
        """Split pathways into GAGE-style groups.

        Preferred: parse the KEGG pathway BRITE hierarchy (br08901) so the
        grouping matches KEGG's own top-level categories.  If that fetch fails,
        fall back to the pathway-number ranges KEGG uses to lay out the BRITE
        (Metabolism 00/01, Genetic info 03, Env/Cellular/Organismal signalling
        02/04, Human Diseases 05, Drug development 07).
        """
        num = {pid: "".join(ch for ch in pid if ch.isdigit())[-5:] for pid in pathway_ids}

        cats: Dict[str, List[str]] = {
            "metabolism": [],
            "genetic_info": [],
            "signaling": [],       # environmental/cellular/organismal processing
            "disease": [],
            "drug": [],
            "other": [],
        }
        try:
            brite = self._get("get/br:br08901")
            mapping = self._parse_brite08901(brite)
        except RuntimeError:
            mapping = {}

        for pid in pathway_ids:
            top = mapping.get(num.get(pid, ""))
            if top is None:
                n = num.get(pid, "")
                first2 = n[:2] if len(n) >= 2 else ""
                if first2 in ("00", "01"):
                    top = "metabolism"
                elif first2 == "03":
                    top = "genetic_info"
                elif first2 in ("02", "04"):
                    top = "signaling"
                elif first2 == "05":
                    top = "disease"
                elif first2 == "07":
                    top = "drug"
                else:
                    top = "other"
            cats[top].append(pid)

        # GAGE convenience aliases
        cats["sigmet"] = cats["metabolism"] + cats["signaling"]
        cats["sig"] = cats["signaling"]
        cats["met"] = cats["metabolism"]
        cats["dise"] = cats["disease"]
        return cats

    @staticmethod
    def _parse_brite08901(brite_text: str) -> Dict[str, str]:
        """Map 5-digit pathway number -> top-level BRITE category name."""
        top_map = {
            "Metabolism": "metabolism",
            "Genetic Information Processing": "genetic_info",
            "Environmental Information Processing": "signaling",
            "Cellular Processes": "signaling",
            "Organismal Systems": "signaling",
            "Human Diseases": "disease",
            "Drug Development": "drug",
        }
        current = "other"
        out: Dict[str, str] = {}
        for line in brite_text.split("\n"):
            if not line or line[0] not in "ABCDEF":
                continue
            level, rest = line[0], line[1:].strip()
            rest = rest.lstrip("+-").strip()
            if level == "A":
                name = rest.split("  ", 1)[-1].strip()
                # strip leading numeric code if present
                name = name.split(" ", 1)[-1].strip() if name[:5].isdigit() else name
                current = next((v for k, v in top_map.items() if k in rest), "other")
            else:
                tok = rest.split(" ", 1)[0]
                digits = "".join(ch for ch in tok if ch.isdigit())
                if len(digits) == 5:
                    out[digits] = current
        return out


class GOGeneSetRetriever:
    """Retrieve Gene Ontology gene sets from a GAF annotation file."""

    ASPECT_TO_DOMAIN = {"P": "BP", "F": "MF", "C": "CC"}

    def get_go_gene_sets(
        self,
        annotation_file: Path,
        id_field: str = "symbol",
        include_iea: bool = True,
        obo_file: Optional[Path] = None,
    ) -> Dict[str, object]:
        """Parse a GAF file into GO gene sets, correctly categorised by domain.

        Args:
            annotation_file: GAF 2.x file (any species; download the matching
                organism's GAF from the GO or an organism annotation DB).
            id_field: ``symbol`` (col 3, default) or ``object_id`` (col 2).
            include_iea: keep electronic (IEA) annotations if ``True``.
            obo_file: optional go-basic.obo to attach real GO-term names
                (the GAF itself does *not* contain GO-term names).

        Returns:
            ``gene_sets`` (GO ID -> genes), ``go_categories`` (BP/MF/CC ->
            GO IDs) and ``go_names`` (GO ID -> term name, empty unless
            ``obo_file`` is given).
        """
        annotation_file = Path(annotation_file)
        gene_sets: Dict[str, set] = defaultdict(set)
        go_domain: Dict[str, str] = {}

        with open(annotation_file) as fh:
            for line in fh:
                if line.startswith("!") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 15:
                    continue
                qualifier = p[3]
                if "NOT" in qualifier.split("|"):
                    continue                       # negative annotation
                evidence = p[6]
                if not include_iea and evidence == "IEA":
                    continue
                gene = p[2] if id_field == "symbol" else p[1]
                go_id = p[4]
                aspect = p[8]                      # P / F / C  <-- the real fix
                gene_sets[go_id].add(gene)
                go_domain[go_id] = self.ASPECT_TO_DOMAIN.get(aspect, "BP")

        categories: Dict[str, List[str]] = {"BP": [], "MF": [], "CC": []}
        for go_id, dom in go_domain.items():
            categories[dom].append(go_id)

        go_names = self._parse_obo_names(obo_file) if obo_file else {}
        return {
            "gene_sets": {k: sorted(v) for k, v in gene_sets.items()},
            "go_categories": categories,
            "go_names": {k: go_names.get(k, "") for k in gene_sets},
        }

    @staticmethod
    def _parse_obo_names(obo_file: Path) -> Dict[str, str]:
        names: Dict[str, str] = {}
        cur: Optional[str] = None
        with open(obo_file) as fh:
            for line in fh:
                line = line.strip()
                if line == "[Term]":
                    cur = None
                elif line.startswith("id: GO:"):
                    cur = line.split("id:", 1)[1].strip()
                elif line.startswith("name:") and cur:
                    names[cur] = line.split("name:", 1)[1].strip()
        return names
