"""OBIS occurrence connector.

Same contract as the GBIF connector: records are read from OBIS on request and
never mirrored, each one keeps its own identifier, and the taxonomy it was
grouped by is part of the result.

Three things differ from GBIF and are visible in the API here:

* **Taxonomy is WoRMS**, not the Catalogue of Life, so a name resolves to an
  AphiaID. Reconciling the two is the caller's job (ChecklistBank), and
  :func:`match_name` states which backbone it used.
* **Paging is a cursor**, not an offset: ``after`` takes the last record id
  seen, so deep paging does not collapse the way GBIF's does.
* **There is no per-query DOI.** GBIF mints one for a download; OBIS does not,
  so a published analysis pins its records by listing the dataset ids and the
  retrieval date, which :func:`search` records in its provenance.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np

from healpix_connector import __version__
from healpix_connector.cells import assign_cells, footprint_cells
from healpix_connector.region import Region

API = "https://api.obis.org/v3"
BACKBONE = "WoRMS"  # OBIS groups occurrences by the World Register of Marine Species
PAGE = 5_000  # the API refuses size > 10,000 (HTTP 400); this leaves headroom
PAGE_PAUSE_S = 0.25
#: Fields kept from each record; asking for them cuts the payload from 68 fields.
FIELDS = ("id", "occurrenceID", "dataset_id", "scientificName", "aphiaID", "eventDate",
          "date_year", "decimalLongitude", "decimalLatitude",
          "coordinateUncertaintyInMeters", "basisOfRecord", "flags")
USER_AGENT = f"healpix-connector/{__version__} (+https://github.com/annefou/healpix-connector)"


class TooManyRecords(RuntimeError):
    """The query matches more records than this call was allowed to fetch."""


@dataclass
class Result:
    records: list[dict]
    provenance: dict = field(default_factory=dict)


def _get(session, url, params=None, retries=6, sleep=time.sleep):
    """GET with retries on rate limiting (429) and server errors (5xx)."""
    for attempt in range(retries):
        r = session.get(url, params=params, timeout=60, headers={"User-Agent": USER_AGENT})
        if r.status_code != 429 and r.status_code < 500:
            r.raise_for_status()
            return r.json()
        wait = min(2 ** attempt, 60)
        retry_after = (getattr(r, "headers", None) or {}).get("Retry-After")
        if retry_after and str(retry_after).isdigit():
            wait = max(wait, int(retry_after))
        sleep(wait)
    r.raise_for_status()


def normalise(rec: dict) -> dict:
    """Keep identity, position, time, taxonomy and OBIS's own quality flags.

    ``flags`` is passed through rather than acted on: ``ON_LAND`` on a marine
    record, or ``NO_ACCEPTED_NAME``, is the repository's own judgement and
    belongs in the result, not in a filter we invent.
    """
    return {
        "obisID": rec.get("id"),
        "occurrenceID": rec.get("occurrenceID"),
        "dataset_id": rec.get("dataset_id"),
        "eventDate": rec.get("eventDate"),
        "year": rec.get("date_year"),
        "decimalLongitude": rec.get("decimalLongitude"),
        "decimalLatitude": rec.get("decimalLatitude"),
        "coordinateUncertaintyInMeters": rec.get("coordinateUncertaintyInMeters"),
        "basisOfRecord": rec.get("basisOfRecord"),
        "taxon_name": rec.get("scientificName"),
        "aphia_id": rec.get("aphiaID"),
        "flags": rec.get("flags") or [],
        "backbone": BACKBONE,
    }


def match_name(name: str, session=None) -> dict:
    """Resolve a scientific name to an AphiaID through OBIS's taxon endpoint."""
    import requests

    session = session or requests.Session()
    d = _get(session, f"{API}/taxon/{name}")
    hits = d.get("results") or []
    top = hits[0] if hits else {}
    return {"query": name, "backbone": BACKBONE,
            # OBIS returns the AphiaID as taxonID; its aphiaID field is often null.
            "aphia_id": top.get("taxonID") or top.get("aphiaID"),
            "name": top.get("scientificName"), "rank": top.get("taxonRank"),
            "accepted_name": top.get("acceptedNameUsage"), "matches": len(hits)}


def statistics(region: Region, *, aphia_id: int | None = None, filters: dict | None = None,
               pad_deg: float = 0.05, session=None) -> dict:
    """Counts for a region without paging through it: records, species, datasets, years.

    Takes the same taxon and filters as :func:`search`, so the count describes
    the query that will actually be paged.
    """
    import requests

    session = session or requests.Session()
    params = {"geometry": region.query_wkt(pad_deg), **(filters or {})}
    if aphia_id is not None:
        params["taxonid"] = int(aphia_id)
    d = _get(session, f"{API}/statistics", params)
    return {"records": d.get("records"), "species": d.get("species"), "taxa": d.get("taxa"),
            "datasets": d.get("datasets"), "year_range": d.get("yearrange")}


def search(region: Region, *, aphia_id: int | None = None, filters: dict | None = None,
           pad_deg: float = 0.05, session=None, max_records: int = 50_000) -> Result:
    """Records in ``region``, filtered exactly to it after the repository query.

    The query polygon is padded, because HEALPix cell edges are not straight in
    longitude and latitude; the padding is then removed by an exact test, so it
    never leaks into the result.
    """
    import requests

    session = session or requests.Session()
    params = {"geometry": region.query_wkt(pad_deg), "size": PAGE,
              "fields": ",".join(FIELDS)}
    if aphia_id is not None:
        params["taxonid"] = int(aphia_id)
    params.update(filters or {})

    counts = statistics(region, aphia_id=aphia_id, filters=filters,
                        pad_deg=pad_deg, session=session)
    count = int(counts["records"] or 0)
    if count > int(max_records):
        raise TooManyRecords(
            f"{count:,} records match; this call allows {int(max_records):,}. Narrow the "
            "region, the taxon or the dates, or raise max_records deliberately.")

    raw, after = [], None
    while len(raw) < count:
        if after is not None:
            time.sleep(PAGE_PAUSE_S)
        page = _get(session, f"{API}/occurrence", {**params, **({"after": after} if after else {})})
        results = page.get("results") or []
        if not results:
            break
        raw.extend(results)
        after = results[-1].get("id")

    recs = [normalise(r) for r in raw]
    recs = [r for r in recs
            if r["decimalLongitude"] is not None and r["decimalLatitude"] is not None]
    lon = np.array([r["decimalLongitude"] for r in recs], dtype=float)
    lat = np.array([r["decimalLatitude"] for r in recs], dtype=float)
    keep = region.contains(lon, lat) if recs else np.zeros(0, bool)
    recs = [r for r, k in zip(recs, keep) if k]
    return Result(recs, {
        "source": "OBIS occurrence API",
        "api": f"{API}/occurrence",
        "params": {k: v for k, v in params.items() if k not in ("size", "fields")},
        "backbone": BACKBONE,
        "matched_in_query_polygon": count,
        "kept_in_region": len(recs),
        "datasets": sorted({r["dataset_id"] for r in recs if r["dataset_id"]}),
        # OBIS mints no DOI for a query, so the dataset ids above and this
        # timestamp are what a published result pins itself to.
        "query_doi": None,
        "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "producer": f"healpix-connector {__version__}",
    })


def dataset_metadata(dataset_id: str, session=None) -> dict:
    """Title, citation, publisher and DOI of one OBIS dataset.

    Many datasets carry no DOI; the field is then ``None`` rather than a guess.
    """
    import requests

    session = session or requests.Session()
    d = _get(session, f"{API}/dataset/{dataset_id}")
    hits = d.get("results") or []
    top = hits[0] if hits else {}
    return {"dataset_id": dataset_id, "title": top.get("title"), "doi": top.get("doi"),
            "citation": top.get("citation"), "url": top.get("url"),
            "published": top.get("published"),
            "institutes": [i.get("name") for i in (top.get("institutes") or [])]}


def to_cells(records: list[dict], depth: int, *, footprints: bool = True) -> list[dict]:
    """Add each record's cell and, optionally, the cells its uncertainty disc covers."""
    if not records:
        return []
    lon = np.array([r["decimalLongitude"] for r in records], dtype=float)
    lat = np.array([r["decimalLatitude"] for r in records], dtype=float)
    cells = assign_cells(lon, lat, depth)
    out = []
    for r, c in zip(records, cells):
        row = {**r, "depth": depth, "cell": int(c)}
        if footprints:
            unc = r.get("coordinateUncertaintyInMeters")
            fp, known = footprint_cells(r["decimalLongitude"], r["decimalLatitude"],
                                        None if unc is None else float(unc), depth)
            row.update(footprint_cells=[int(x) for x in fp], footprint_known=known)
        out.append(row)
    return out
