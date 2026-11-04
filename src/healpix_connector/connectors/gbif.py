"""GBIF occurrence connector.

Rules (see CLAUDE.md): always send ``checklistKey`` and record the taxonomy
used; identify records by ``gbifID``; up to 100,000 records via the search API,
beyond that a GBIF download (record-level, with a DOI). Nothing is mirrored.
"""

from __future__ import annotations

import time
from datetime import datetime, timezone

from healpix_connector import __version__
from healpix_connector.connectors._common import (  # noqa: F401  (re-exported)
    USER_AGENT, Result, TooManyRecords, keep_in_region, to_cells)
from healpix_connector.connectors._common import get_json as _get
from healpix_connector.region import Region

API = "https://api.gbif.org/v1"
COL_XR = "7ddf754f-d193-4cc9-b351-99906754a03b"  # Catalogue of Life eXtended Release
LEGACY_BACKBONE = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"  # GBIF Backbone Taxonomy (last built 2023)
SEARCH_LIMIT = 100_000  # GBIF refuses offset + limit beyond this (HTTP 400)
# Deep paging collapses, with or without a geometry filter. Measured 2026-09-22:
# one page took 0.5 s at offset 0 and 0.4 s at offset 6,000, but 361 s at offset
# 12,000 (394 s for the same query without geometry). Past a few thousand
# records, use a download.
PRACTICAL_PAGING_LIMIT = 5_000
PAGE = 300  # maximum page size of the occurrence search API
PAGE_PAUSE_S = 0.25  # between successive pages, to stay under GBIF's rate limit
#: Names of the checklists a record's taxonomy can come from; any other key is reported as itself.
CHECKLISTS = {COL_XR: "Catalogue of Life eXtended Release", LEGACY_BACKBONE: "GBIF Backbone Taxonomy"}


def normalise(rec: dict, checklist_key: str) -> dict:
    """Flatten one search-API record, keeping identity, position, time, taxonomy and the
    fields filters need: country, GBIF's own issue flags (passed through, not acted on)
    and ``occurrenceStatus``, so that absence records are never read as presences.

    Returns the core record (``_common.CORE_FIELDS``) plus GBIF's native fields."""
    cls = (rec.get("classifications") or {}).get(checklist_key) or {}
    usage, accepted = cls.get("usage") or {}, cls.get("acceptedUsage") or {}
    gbif_id = str(rec.get("gbifID") or rec.get("key"))
    issues = list(rec.get("issues") or [])
    return {
        "source": "GBIF",
        "record_id": gbif_id,
        "dataset_id": rec.get("datasetKey"),
        "taxon_id": usage.get("key"),
        "taxonomy": CHECKLISTS.get(checklist_key, checklist_key),
        "taxonomy_key": checklist_key,
        "flags": issues,
        "gbifID": gbif_id,
        "occurrenceID": rec.get("occurrenceID"),
        "datasetKey": rec.get("datasetKey"),
        "eventDate": rec.get("eventDate"),
        "decimalLongitude": rec.get("decimalLongitude"),
        "decimalLatitude": rec.get("decimalLatitude"),
        "coordinateUncertaintyInMeters": rec.get("coordinateUncertaintyInMeters"),
        "countryCode": rec.get("countryCode"),
        "year": rec.get("year"),
        "occurrenceStatus": rec.get("occurrenceStatus"),
        "issues": issues,
        "basisOfRecord": rec.get("basisOfRecord"),
        "license": rec.get("license"),
        "taxon_key": usage.get("key"),
        "taxon_name": usage.get("name"),
        "accepted_key": accepted.get("key"),
        "taxonomic_status": cls.get("taxonomicStatus"),
        "checklist_key": checklist_key,
    }


def match_name(name: str, *, checklist_key: str = COL_XR, session=None) -> dict:
    """Resolve a scientific name to a taxon key under ``checklist_key`` (GBIF v2 match).

    Returns the key, the matched name, its rank and GBIF's match type. Callers
    should treat anything other than an EXACT match as needing review, as GBIF
    Alert does, rather than guess.
    """
    import requests

    session = session or requests.Session()
    d = _get(session, "https://api.gbif.org/v2/species/match",
             {"checklistKey": checklist_key, "scientificName": name})
    usage = d.get("usage") or {}
    return {"query": name, "checklist_key": checklist_key, "taxon_key": usage.get("key"),
            "name": usage.get("name"), "rank": usage.get("rank"),
            "match_type": (d.get("diagnostics") or {}).get("matchType")}


def search(region: Region, *, checklist_key: str = COL_XR, taxon_key: str | None = None,
           filters: dict | None = None, pad_deg: float = 0.05, session=None,
           max_records: int = PRACTICAL_PAGING_LIMIT) -> Result:
    """Records in ``region`` via the search API, filtered exactly to the region.

    ``max_records`` bounds the work. It defaults to ``PRACTICAL_PAGING_LIMIT``,
    well below the API's own ceiling, because deep paging becomes pathologically
    slow (see that constant). Raise it deliberately, or use a download.
    """
    import requests

    session = session or requests.Session()
    params = {"geometry": region.query_wkt(pad_deg), "hasCoordinate": "true",
              "checklistKey": checklist_key, "limit": PAGE}
    if taxon_key is not None:
        params["taxonKey"] = taxon_key
    params.update(filters or {})

    first = _get(session, f"{API}/occurrence/search", {**params, "limit": 0})
    count = int(first["count"])
    cap = min(int(max_records), SEARCH_LIMIT)
    if count > cap:
        raise TooManyRecords(
            f"{count:,} records match; this call allows {cap:,}. Deep paging is pathologically "
            f"slow (~361 s for one page at offset 12,000) and the API stops at {SEARCH_LIMIT:,} "
            "anyway. Narrow the region or filters, raise max_records deliberately, or use a "
            "GBIF download (record-level, with a DOI).")
    raw = []
    for offset in range(0, count, PAGE):
        if offset:
            time.sleep(PAGE_PAUSE_S)
        page = _get(session, f"{API}/occurrence/search", {**params, "offset": offset})
        raw.extend(page["results"])
        if page.get("endOfRecords"):
            break
    recs = keep_in_region([normalise(r, checklist_key) for r in raw], region)
    return Result(recs, {
        "source": "GBIF occurrence search API",
        "api": f"{API}/occurrence/search",
        "params": {k: v for k, v in params.items() if k != "limit"},
        "checklist_key": checklist_key,
        "matched_in_query_polygon": count,
        "kept_in_region": len(recs),
        "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "producer": f"healpix-connector {__version__}",
    })


def download_metadata(key: str, session=None) -> dict:
    """Metadata of an existing GBIF download: DOI, request, record count, taxonomy."""
    import requests

    session = session or requests.Session()
    d = _get(session, f"{API}/occurrence/download/{key}")
    req = d.get("request") or {}
    return {
        "key": d["key"], "doi": d.get("doi"), "status": d.get("status"),
        "created": d.get("created"), "total_records": d.get("totalRecords"),
        "format": req.get("format"), "predicate": req.get("predicate"),
        # Downloads made before GBIF's 2026 switch to COL XR usually carry none.
        "checklist_key": req.get("checklistKey") or "not specified in the download request",
        "download_link": d.get("downloadLink"),
    }
