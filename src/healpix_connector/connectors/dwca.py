"""Darwin Core Archives: occurrence records from a download, as core records.

A Darwin Core Archive (DwC-A) is how repositories hand over record-level data in
bulk: a GBIF download (with its DOI), an IPT resource, an OBIS export. It is a
zip holding data tables, a ``meta.xml`` that says which column is which Darwin
Core term, and usually an ``eml.xml`` with the dataset's metadata. GBIF's
``SIMPLE_CSV`` download is the same idea without ``meta.xml``: one tab-separated
table whose header row holds the term names; it is read too.

Two layouts are handled:

* **Occurrence core** (GBIF downloads, most IPT occurrence datasets): one row per
  record.
* **Event core with an Occurrence extension** (sampling-event datasets, common in
  OBIS): an occurrence often states no position or date of its own. Each one
  inherits the fields it lacks from its event, and an event from its parent
  events (``parentEventID``), the way the archive defines them.

The archive is read in chunks and filtered to the region chunk by chunk, so a
download of millions of records never has to fit in memory. Nothing is mirrored:
:func:`healpix_connector.connectors.gbif.download_records` deletes the archive
after reading it unless asked to keep it.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

import numpy as np

from healpix_connector import __version__
from healpix_connector.connectors._common import Result
from healpix_connector.region import Region

CHUNK_ROWS = 200_000
#: Fields an occurrence may take from its event when it states none itself.
EVENT_FIELDS = ("eventDate", "year", "decimalLatitude", "decimalLongitude",
                "coordinateUncertaintyInMeters", "datasetID", "countryCode")
FLOATS = ("decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters")
#: Columns read for every record; anything else only when asked for (``keep_terms``).
TERMS = ("id", "coreid", "gbifID", "occurrenceID", "datasetKey", "datasetID", "eventDate", "year",
         "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "basisOfRecord",
         "occurrenceStatus", "countryCode", "license", "scientificName", "scientificNameID",
         "taxonID", "taxonKey", "acceptedTaxonKey", "taxonomicStatus", "issue", "flags")
_DOI = re.compile(r"(10\.\d{4,9}/\S+)")


def _local(term: str) -> str:
    """``http://rs.tdwg.org/dwc/terms/decimalLatitude`` -> ``decimalLatitude``."""
    return re.split(r"[/#]", term.rstrip("/"))[-1]


def _unescape(s: str | None, default: str) -> str:
    if s is None:
        return default
    return s.replace("\\t", "\t").replace("\\n", "\n").replace("\\r", "\r")


# ------------------------------------------------------------- archive layout
def _table(el, ns) -> dict:
    """One ``<core>`` or ``<extension>`` element as a reading recipe."""
    cols, defaults = {}, {}
    for f in el.findall(f"{ns}field"):
        name = _local(f.get("term"))
        if f.get("index") is not None:
            cols[int(f.get("index"))] = name
        elif f.get("default") is not None:
            defaults[name] = f.get("default")
    id_el, coreid_el = el.find(f"{ns}id"), el.find(f"{ns}coreid")
    key = id_el if id_el is not None else coreid_el
    key_name = None if key is None else ("id" if key is id_el else "coreid")
    return {
        "row_type": _local(el.get("rowType", "")),
        "files": [loc.text.strip() for loc in el.iter(f"{ns}location")],
        "sep": _unescape(el.get("fieldsTerminatedBy"), ","),
        "quote": _unescape(el.get("fieldsEnclosedBy"), '"'),
        "skip": int(el.get("ignoreHeaderLines", "0")),
        "encoding": el.get("encoding", "UTF-8"),
        "columns": cols, "defaults": defaults, "key": key_name,
        "key_index": None if key is None else int(key.get("index")),
    }


def layout(z: zipfile.ZipFile) -> dict:
    """The archive's core and extensions as reading recipes (from ``meta.xml``, or guessed
    from a single header-row table such as a GBIF ``SIMPLE_CSV`` download)."""
    names = z.namelist()
    meta = next((n for n in names if n.rsplit("/", 1)[-1] == "meta.xml"), None)
    if meta is None:
        data = [n for n in names if n.lower().endswith((".csv", ".txt", ".tsv")) and not n.endswith("/")]
        if len(data) != 1:
            raise ValueError(f"no meta.xml and not one data table: {data}")
        with z.open(data[0]) as f:
            header = io.TextIOWrapper(f, encoding="utf-8").readline().rstrip("\r\n")
        sep = "\t" if header.count("\t") >= header.count(",") else ","
        cols = {i: _local(t) for i, t in enumerate(next(csv.reader([header], delimiter=sep)))}
        core = {"row_type": "Occurrence", "files": data, "sep": sep, "quote": "" if sep == "\t" else '"',
                "skip": 1, "encoding": "UTF-8", "columns": cols, "defaults": {}, "key": None,
                "key_index": None}
        return {"core": core, "extensions": [], "metadata": None, "prefix": ""}
    root = ElementTree.fromstring(z.read(meta))
    ns = root.tag[: root.tag.index("}") + 1] if root.tag.startswith("{") else ""
    prefix = meta[: -len("meta.xml")]
    return {"core": _table(root.find(f"{ns}core"), ns),
            "extensions": [_table(e, ns) for e in root.findall(f"{ns}extension")],
            "metadata": root.get("metadata"), "prefix": prefix}


def eml_summary(z: zipfile.ZipFile, lay: dict | None = None) -> dict:
    """Title, DOI and package id from the archive's EML, where it has one."""
    lay = lay or layout(z)
    name = lay["prefix"] + (lay["metadata"] or "eml.xml")
    if name not in z.namelist():
        return {}
    root = ElementTree.fromstring(z.read(name))
    texts = lambda tag: [e.text.strip() for e in root.iter() if e.tag.rsplit("}", 1)[-1] == tag and e.text]
    doi = next((m.group(1) for t in texts("alternateIdentifier") + texts("identifier")
                for m in [_DOI.search(t)] if m), None)
    title = next(iter(texts("title")), None)
    return {"title": title, "doi": doi, "package_id": root.get("packageId")}


# --------------------------------------------------------------- reading
def _read(z, prefix, t, wanted, chunk_rows):
    """Chunks of table ``t`` as DataFrames holding the ``wanted`` columns present (strings)."""
    import pandas as pd

    use = {i: n for i, n in t["columns"].items() if n in wanted}
    key, ki = t["key"], t["key_index"]
    # The key column is often also a term's column (GBIF: id = gbifID); it then gets both names.
    key_alias = ki is not None and ki in use
    if ki is not None and not key_alias:
        use[ki] = key
    for name in t["files"]:
        with z.open(prefix + name) as f:
            reader = pd.read_csv(
                f, sep=t["sep"], header=None, skiprows=t["skip"], usecols=sorted(use),
                dtype=str, keep_default_na=False, na_values=[""], encoding=t["encoding"],
                quoting=csv.QUOTE_NONE if not t["quote"] else csv.QUOTE_MINIMAL,
                quotechar=t["quote"] or '"', chunksize=chunk_rows, on_bad_lines="error")
            for chunk in reader:
                chunk = chunk.rename(columns=use)
                if key_alias:
                    chunk[key] = chunk[use[ki]]
                for k, v in t["defaults"].items():
                    if k in wanted:
                        chunk[k] = chunk[k].fillna(v) if k in chunk else v
                yield chunk


def _events(z, prefix, core) -> "pd.DataFrame":
    """Event table indexed by event id, each event carrying what its parents state."""
    import pandas as pd

    ev = pd.concat(list(_read(z, prefix, core, set(EVENT_FIELDS) | {"parentEventID"}, CHUNK_ROWS)))
    ev = ev.set_index("id")
    if "parentEventID" in ev:
        for _ in range(20):  # sampling-event hierarchies are a few levels deep
            parent = ev["parentEventID"]
            changed = False
            for col in (c for c in EVENT_FIELDS if c in ev):
                inherited = ev[col].isna() & parent.notna()
                if inherited.any():
                    fill = ev[col].reindex(parent[inherited]).to_numpy()
                    got = pd.notna(fill)
                    if got.any():
                        idx = inherited[inherited].index[got]
                        ev.loc[idx, col] = fill[got]
                        changed = True
            if not changed:
                break
    return ev


def _text(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)


def _number(v):
    """A float, or None for an empty or unreadable value (the record is kept, the field is not)."""
    try:
        x = float(_text(v))
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _core_record(t: dict, profile: str, taxonomy: str, taxonomy_key, source: str, package_id,
                 keep_terms) -> dict:
    g = t.get
    year = _number(g("year"))
    rec = {
        "source": source,
        "occurrenceID": _text(g("occurrenceID")),
        "eventDate": _text(g("eventDate")),
        "year": None if year is None else int(year),
        **{k: _number(g(k)) for k in FLOATS},
        "basisOfRecord": _text(g("basisOfRecord")),
        "taxon_name": _text(g("scientificName")),
        "taxonomy": taxonomy, "taxonomy_key": taxonomy_key,
        "occurrenceStatus": _text(g("occurrenceStatus")),
        "countryCode": _text(g("countryCode")),
        "license": _text(g("license")),
    }
    if profile == "gbif":
        issues = [s for s in (_text(g("issue")) or "").split(";") if s]
        rec.update(record_id=_text(g("gbifID")) or _text(g("id")), dataset_id=_text(g("datasetKey")),
                   taxon_id=_text(g("taxonKey")), flags=issues,
                   # the same native fields as a search-API record
                   gbifID=_text(g("gbifID")) or _text(g("id")), datasetKey=_text(g("datasetKey")),
                   issues=issues, taxon_key=_text(g("taxonKey")), accepted_key=_text(g("acceptedTaxonKey")),
                   taxonomic_status=_text(g("taxonomicStatus")), checklist_key=taxonomy_key)
    else:
        rec.update(record_id=_text(g("occurrenceID")) or _text(g("id")),
                   dataset_id=_text(g("datasetID")) or package_id,
                   taxon_id=_text(g("taxonID")) or _text(g("scientificNameID")),
                   flags=[s for s in re.split(r"[;,]", _text(g("flags")) or "") if s],
                   scientificNameID=_text(g("scientificNameID")))
    for k in keep_terms:
        rec.setdefault(k, _text(g(k)))
    return rec


def read(path, region: Region | None = None, *, source: str | None = None,
         taxonomy: str | None = None, taxonomy_key: str | None = None,
         keep_terms: tuple[str, ...] = (), chunk_rows: int = CHUNK_ROWS) -> Result:
    """Occurrence records of a Darwin Core Archive, as core records.

    ``region`` keeps only the records inside it (exact test); records without
    coordinates are dropped and counted either way, since they have no cell.
    A GBIF archive (it has ``gbifID``) gets the same native fields as a search-API
    record; any other gets ``record_id`` = ``occurrenceID``. ``taxonomy`` names
    the taxonomy its taxon ids belong to: an archive does not say, so say it when
    you know (a GBIF download's ``checklistKey``, WoRMS for OBIS); otherwise it is
    reported as ``"as published"``. ``keep_terms`` keeps further Darwin Core terms.
    """
    import pandas as pd

    path = Path(path)
    sha = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            sha.update(block)

    with zipfile.ZipFile(path) as z:
        lay = layout(z)
        eml = eml_summary(z, lay)
        core = lay["core"]
        if core["row_type"] == "Occurrence":
            occ, events = core, None
        elif core["row_type"] == "Event":
            occ = next((e for e in lay["extensions"] if e["row_type"] == "Occurrence"), None)
            if occ is None:
                raise ValueError(f"{path.name}: an Event core without an Occurrence extension holds no records")
            events = _events(z, lay["prefix"], core)
        else:
            raise ValueError(f"{path.name}: core {core['row_type']!r} holds no occurrence records")

        names = set(occ["columns"].values())
        profile = "gbif" if "gbifID" in names else "dwc"
        source = source or ("GBIF" if profile == "gbif" else (eml.get("title") or path.stem))
        taxonomy = taxonomy or ("not specified in the archive" if profile == "gbif" else "as published")
        wanted = set(TERMS) | set(keep_terms)

        records, rows, with_coords = [], 0, 0
        for chunk in _read(z, lay["prefix"], occ, wanted, chunk_rows):
            rows += len(chunk)
            if events is not None and "coreid" in chunk:
                ev = events.reindex(chunk["coreid"].to_numpy())
                for col in (c for c in EVENT_FIELDS if c in ev):
                    own = chunk[col] if col in chunk else pd.Series(pd.NA, index=chunk.index, dtype=object)
                    chunk[col] = own.where(own.notna(), ev[col].to_numpy())
            if "decimalLongitude" not in chunk or "decimalLatitude" not in chunk:
                continue
            lon = pd.to_numeric(chunk["decimalLongitude"], errors="coerce")
            lat = pd.to_numeric(chunk["decimalLatitude"], errors="coerce")
            ok = np.array(lon.notna() & lat.notna(), dtype=bool)  # a writable copy
            with_coords += int(ok.sum())
            if region is not None:
                ok[ok] = region.contains(lon.to_numpy()[ok], lat.to_numpy()[ok])
            for row in chunk[ok].to_dict("records"):
                records.append(_core_record(row, profile, taxonomy, taxonomy_key, source,
                                            eml.get("package_id"), keep_terms))

    return Result(records, {
        "source": f"{source} (Darwin Core Archive)",
        "api": None,
        "params": {"region": None if region is None else region.query_wkt()},
        "archive": path.name, "sha256": sha.hexdigest(),
        "core": core["row_type"], "occurrences_from": occ["row_type"] + ("" if events is None else " extension"),
        "archive_title": eml.get("title"), "archive_doi": eml.get("doi"),
        "taxonomy": taxonomy, "taxonomy_key": taxonomy_key,
        "rows_read": rows, "with_coordinates": with_coords, "kept_in_region": len(records),
        "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "producer": f"healpix-connector {__version__}",
    })
