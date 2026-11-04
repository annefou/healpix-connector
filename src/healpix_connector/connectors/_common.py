"""What every occurrence connector shares: HTTP retries, the result type, the exact
region filter, cell assignment, and the core record every connector returns.

The core record is Darwin Core position, time and identity, plus the fields that
make records from different repositories comparable: which repository it came
from (``source``), its identifier there (``record_id``), the dataset, the taxon and
the taxonomy that taxon belongs to, and the repository's own quality flags. Each
connector keeps its repository's native fields beside these (``gbifID``,
``aphia_id``, ...), so nothing a caller already reads goes away.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from healpix_connector import __version__
from healpix_connector.cells import assign_cells, footprint_cells
from healpix_connector.region import Region

USER_AGENT = f"healpix-connector/{__version__} (+https://github.com/annefou/healpix-connector)"

#: Fields every connector's ``normalise`` returns, whatever the repository.
CORE_FIELDS = (
    "source",                          # repository name, e.g. "GBIF", "OBIS"
    "record_id",                       # the repository's own identifier, as a string
    "occurrenceID",                    # Darwin Core; the publisher's id, shared across repositories
    "dataset_id",                      # the dataset the record was published in
    "eventDate", "year",               # time matching keys on eventDate, never publication date
    "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters",
    "basisOfRecord",
    "taxon_id", "taxon_name",          # in the taxonomy named next
    "taxonomy", "taxonomy_key",        # e.g. "WoRMS"; key = checklist id where there is one
    "flags",                           # the repository's own quality flags, passed through
)


class TooManyRecords(RuntimeError):
    """The query matches more records than this call was allowed to fetch."""


@dataclass
class Result:
    records: list[dict]
    provenance: dict = field(default_factory=dict)


def get_json(session, url, params=None, retries=6, sleep=time.sleep):
    """GET with retries on rate limiting (429) and server errors (5xx).

    Repositories answer 429 when requests come too fast and 503 during backend
    hiccups; both are retried with exponential backoff, honouring ``Retry-After``
    when given. Any other error status is raised at once.
    """
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


def keep_in_region(records: list[dict], region: Region) -> list[dict]:
    """Records with coordinates that lie inside ``region`` (the exact test after a padded query)."""
    recs = [r for r in records
            if r["decimalLongitude"] is not None and r["decimalLatitude"] is not None]
    if not recs:
        return []
    lon = np.array([r["decimalLongitude"] for r in recs], dtype=float)
    lat = np.array([r["decimalLatitude"] for r in recs], dtype=float)
    return [r for r, k in zip(recs, region.contains(lon, lat)) if k]


def missing_core_fields(record: dict) -> list[str]:
    """Core fields absent from ``record`` (an empty list means it conforms)."""
    return [f for f in CORE_FIELDS if f not in record]


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
