"""Darwin Core Archives: the layouts repositories actually ship, built here in miniature."""

import zipfile

import pytest

from healpix_connector.connectors import _common, dwca, gbif
from healpix_connector.region import Region

DWC, GBIF_T, DC = "http://rs.tdwg.org/dwc/terms/", "http://rs.gbif.org/terms/1.0/", "http://purl.org/dc/terms/"
IBERIA = Region.from_bbox(-10, 35, 4, 44)


def _zip(path, files: dict):
    with zipfile.ZipFile(path, "w") as z:
        for name, text in files.items():
            z.writestr(name, text)
    return path


def _fields(terms):
    return "\n".join(f'    <field index="{i}" term="{t}"/>' for i, t in enumerate(terms))


# A GBIF DWCA download: tab-separated, no quoting, gbifID doubles as the row id, issues joined by ';',
# the download's EML is metadata.xml.
GBIF_TERMS = [GBIF_T + "gbifID", GBIF_T + "datasetKey", DWC + "occurrenceID", DWC + "eventDate", DWC + "year",
              DWC + "decimalLatitude", DWC + "decimalLongitude", DWC + "coordinateUncertaintyInMeters",
              DWC + "scientificName", GBIF_T + "taxonKey", DWC + "occurrenceStatus", GBIF_T + "issue",
              DWC + "countryCode", DWC + "recordedBy"]
GBIF_META = f"""<archive xmlns="http://rs.tdwg.org/dwc/text/" metadata="metadata.xml">
  <core encoding="UTF-8" fieldsTerminatedBy="\\t" linesTerminatedBy="\\n" fieldsEnclosedBy=""
        ignoreHeaderLines="1" rowType="{DWC}Occurrence">
    <files><location>occurrence.txt</location></files>
    <id index="0"/>
{_fields(GBIF_TERMS)}
    <field term="{DC}license" default="CC_BY_4_0"/>
  </core>
</archive>"""
GBIF_ROWS = [
    ["1001", "ds-1", "occ-a", "2019-05-01", "2019", "40.4", "-3.7", "250", "Calopteryx virgo", "Q2M4",
     "PRESENT", "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES", "ES", 'Ann "the" Observer'],
    ["1002", "ds-1", "occ-b", "2020-06", "2020", "41.0", "-4.0", "", "Calopteryx virgo", "Q2M4",
     "ABSENT", "", "ES", ""],
    ["1003", "ds-2", "occ-c", "", "", "", "", "", "Calopteryx virgo", "Q2M4", "PRESENT", "", "", ""],  # no position
    ["1004", "ds-2", "occ-d", "2021-01-01", "2021", "52.0", "5.0", "10", "Calopteryx virgo", "Q2M4",
     "PRESENT", "", "NL", ""],  # outside Iberia
]
EML = """<eml:eml xmlns:eml="https://eml.ecoinformatics.org/eml-2.2.0" packageId="pkg-42">
  <dataset><alternateIdentifier>https://doi.org/10.15468/dl.abc123</alternateIdentifier>
  <title>GBIF Occurrence Download 0001234-261104</title></dataset></eml:eml>"""


@pytest.fixture
def gbif_archive(tmp_path):
    table = "\t".join(t.rsplit("/", 1)[-1] for t in GBIF_TERMS) + "\n"
    table += "".join("\t".join(r) + "\n" for r in GBIF_ROWS)
    return _zip(tmp_path / "0001234-261104.zip",
                {"meta.xml": GBIF_META, "occurrence.txt": table, "metadata.xml": EML})


def test_gbif_archive_gives_core_records_like_the_search_api(gbif_archive):
    res = dwca.read(gbif_archive, IBERIA, taxonomy="Catalogue of Life eXtended Release",
                    taxonomy_key=gbif.COL_XR)
    assert [r["gbifID"] for r in res.records] == ["1001", "1002"]
    a, b = res.records
    assert _common.missing_core_fields(a) == []
    assert a["source"] == "GBIF" and a["record_id"] == "1001" and a["dataset_id"] == "ds-1"
    assert a["flags"] == a["issues"] == ["COORDINATE_ROUNDED", "COUNTRY_DERIVED_FROM_COORDINATES"]
    assert a["decimalLatitude"] == 40.4 and a["coordinateUncertaintyInMeters"] == 250.0 and a["year"] == 2019
    assert a["taxon_id"] == a["taxon_key"] == "Q2M4" and a["checklist_key"] == gbif.COL_XR
    assert a["license"] == "CC_BY_4_0"  # a meta.xml default, applied to every row
    assert b["occurrenceStatus"] == "ABSENT" and b["coordinateUncertaintyInMeters"] is None
    p = res.provenance
    assert (p["rows_read"], p["with_coordinates"], p["kept_in_region"]) == (4, 3, 2)
    assert p["archive_doi"] == "10.15468/dl.abc123" and len(p["sha256"]) == 64


def test_extra_terms_are_kept_only_when_asked_for(gbif_archive):
    plain = dwca.read(gbif_archive, IBERIA).records[0]
    extra = dwca.read(gbif_archive, IBERIA, keep_terms=("recordedBy",)).records[0]
    assert "recordedBy" not in plain and extra["recordedBy"] == 'Ann "the" Observer'  # unquoted TSV, verbatim


def test_without_a_region_every_record_with_coordinates_is_kept(gbif_archive):
    res = dwca.read(gbif_archive)
    assert len(res.records) == 3 and res.provenance["taxonomy"] == "not specified in the archive"


def test_reading_in_small_chunks_gives_the_same_records(gbif_archive):
    assert dwca.read(gbif_archive, IBERIA, chunk_rows=1).records == dwca.read(gbif_archive, IBERIA).records


# An IPT sampling-event archive: quoted CSV, an Event core whose child events take their position from a
# parent, and an Occurrence extension that states neither position nor date.
EVENT_META = f"""<archive xmlns="http://rs.tdwg.org/dwc/text/" metadata="eml.xml">
  <core encoding="UTF-8" fieldsTerminatedBy="," fieldsEnclosedBy='"' ignoreHeaderLines="1" rowType="{DWC}Event">
    <files><location>event.txt</location></files>
    <id index="0"/>
{_fields([DWC + "eventID", DWC + "parentEventID", DWC + "eventDate", DWC + "decimalLatitude",
          DWC + "decimalLongitude", DWC + "datasetID"])}
  </core>
  <extension encoding="UTF-8" fieldsTerminatedBy="," fieldsEnclosedBy='"' ignoreHeaderLines="1"
             rowType="{DWC}Occurrence">
    <files><location>occurrence.txt</location></files>
    <coreid index="0"/>
{_fields([DWC + "eventID", DWC + "occurrenceID", DWC + "scientificName", DWC + "scientificNameID",
          DWC + "decimalLatitude", DWC + "decimalLongitude"])}
  </extension>
</archive>"""
EVENTS = """id,parentEventID,eventDate,decimalLatitude,decimalLongitude,datasetID
cruise,,2018-07,42.5,-9.2,ds-cruise
station-1,cruise,2018-07-03,,,
haul-1,station-1,2018-07-03/04,,,
station-2,cruise,2018-07-09,43.1,-8.5,
"""
OCCURRENCES = """coreid,occurrenceID,scientificName,scientificNameID,decimalLatitude,decimalLongitude
haul-1,o1,"Delphinus delphis, Linnaeus 1758",urn:lsid:marinespecies.org:taxname:137094,,
station-2,o2,Tursiops truncatus,urn:lsid:marinespecies.org:taxname:137111,,
station-2,o3,Phocoena phocoena,urn:lsid:marinespecies.org:taxname:137117,43.2,-8.4
"""


@pytest.fixture
def event_archive(tmp_path):
    return _zip(tmp_path / "ipt.zip", {"meta.xml": EVENT_META, "event.txt": EVENTS,
                                       "occurrence.txt": OCCURRENCES, "eml.xml": EML})


def test_occurrences_inherit_position_and_date_through_the_event_hierarchy(event_archive):
    res = dwca.read(event_archive, IBERIA, taxonomy="WoRMS")
    by_id = {r["record_id"]: r for r in res.records}
    o1, o2, o3 = by_id["o1"], by_id["o2"], by_id["o3"]
    # o1 -> haul-1 (its date) -> station-1 -> cruise (its position and dataset)
    assert (o1["decimalLatitude"], o1["decimalLongitude"], o1["eventDate"]) == (42.5, -9.2, "2018-07-03/04")
    assert o1["dataset_id"] == "ds-cruise"
    assert (o2["decimalLatitude"], o2["eventDate"]) == (43.1, "2018-07-09")
    assert o3["decimalLatitude"] == 43.2  # an occurrence's own position wins over its event's
    assert o1["taxon_name"] == "Delphinus delphis, Linnaeus 1758"  # quoted comma survives
    assert o1["taxon_id"] == "urn:lsid:marinespecies.org:taxname:137094" and o1["taxonomy"] == "WoRMS"
    assert all(_common.missing_core_fields(r) == [] for r in res.records)
    assert res.provenance["occurrences_from"] == "Occurrence extension"


def test_an_event_archive_without_occurrences_is_refused(tmp_path):
    meta = EVENT_META[: EVENT_META.index("<extension")] + "</archive>"
    path = _zip(tmp_path / "events.zip", {"meta.xml": meta, "event.txt": EVENTS})
    with pytest.raises(ValueError, match="without an Occurrence extension"):
        dwca.read(path)


def test_a_simple_csv_download_is_read_from_its_header(tmp_path):
    header = ["gbifID", "datasetKey", "occurrenceID", "eventDate", "year", "decimalLatitude",
              "decimalLongitude", "coordinateUncertaintyInMeters", "taxonKey", "issue"]
    rows = [["7", "ds", "x", "2022-03-04", "2022", "40.0", "-3.0", "30", "Q2M4", "RECORDED_DATE_INVALID"]]
    text = "\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in rows)
    res = dwca.read(_zip(tmp_path / "simple.zip", {"0009-261104.csv": text}), IBERIA)
    (r,) = res.records
    assert r["gbifID"] == "7" and r["flags"] == ["RECORDED_DATE_INVALID"] and r["decimalLongitude"] == -3.0


def test_records_from_an_archive_go_straight_to_cells(event_archive):
    rows = _common.to_cells(dwca.read(event_archive, IBERIA).records, depth=8)
    assert all(r["footprint_known"] is False for r in rows)  # no uncertainty in the archive


# ------------------------------------------------------- gbif.download_records
class _Resp:
    def __init__(self, payload=None, body=b""):
        self.payload, self.body, self.status_code, self.headers = payload, body, 200, {}

    def json(self):
        return self.payload

    def raise_for_status(self):
        pass

    def iter_content(self, n):
        for i in range(0, len(self.body), n):
            yield self.body[i:i + n]

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _Session:
    def __init__(self, archive_bytes, status="SUCCEEDED", fmt="DWCA", checklist=gbif.COL_XR):
        self.archive, self.calls = archive_bytes, []
        self.meta = {"key": "0001234-261104", "doi": "10.15468/dl.abc123", "status": status,
                     "created": "2026-11-04", "totalRecords": 4, "downloadLink": "https://x/0001234.zip",
                     "request": {"format": fmt, "predicate": {"type": "within"},
                                 **({"checklistKey": checklist} if checklist else {})}}

    def get(self, url, params=None, timeout=None, headers=None, stream=False):
        self.calls.append(url)
        return _Resp(body=self.archive) if stream else _Resp(self.meta)


def test_download_records_reads_the_archive_and_keeps_nothing(gbif_archive, tmp_path):
    s = _Session(gbif_archive.read_bytes())
    res = gbif.download_records("0001234-261104", IBERIA, session=s)
    assert [r["record_id"] for r in res.records] == ["1001", "1002"]
    assert res.records[0]["taxonomy"] == "Catalogue of Life eXtended Release"
    p = res.provenance
    assert p["doi"] == "10.15468/dl.abc123" and p["download_key"] == "0001234-261104"
    assert p["checklist_key"] == gbif.COL_XR and p["kept_in_region"] == 2


def test_a_kept_archive_is_reused_rather_than_fetched_again(gbif_archive, tmp_path):
    keep = tmp_path / "kept" / "dl.zip"
    s = _Session(gbif_archive.read_bytes())
    gbif.download_records("k", IBERIA, session=s, keep_archive=keep)
    gbif.download_records("k", IBERIA, session=s, keep_archive=keep)
    assert keep.exists() and s.calls.count("https://x/0001234.zip") == 1


def test_a_download_without_a_checklist_says_so(gbif_archive):
    res = gbif.download_records("k", session=_Session(gbif_archive.read_bytes(), checklist=None))
    assert res.records[0]["taxonomy"] == "not specified in the download request"
    assert res.records[0]["taxonomy_key"] is None


@pytest.mark.parametrize("kw, err", [({"status": "RUNNING"}, RuntimeError), ({"fmt": "SPECIES_LIST"}, ValueError)])
def test_unfinished_or_recordless_downloads_are_refused(kw, err):
    with pytest.raises(err):
        gbif.download_records("k", session=_Session(b"", **kw))

