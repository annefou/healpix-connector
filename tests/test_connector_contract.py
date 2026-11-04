"""The contract every occurrence connector meets, run against each one.

A new connector joins ``CONNECTORS`` with a recorded fixture and a fake session;
if it passes these tests, records from it can be mixed with any other
connector's and passed to ``sample.sample_records`` unchanged.
"""

import pytest

from healpix_connector.connectors import _common, gbif, obis
from healpix_connector.region import Region
from test_gbif import FIX as GBIF_FIX, FakeSession as GbifSession
from test_obis import FIX as OBIS_FIX, FakeSession as ObisSession

CONNECTORS = {
    "GBIF": dict(module=gbif, fixture=GBIF_FIX, session=GbifSession,
                 normalise=lambda rec: gbif.normalise(rec, gbif.COL_XR)),
    "OBIS": dict(module=obis, fixture=OBIS_FIX, session=ObisSession, normalise=obis.normalise),
}
#: What a result must carry to be cited and reproduced, whatever the repository.
PROVENANCE = ("source", "api", "params", "matched_in_query_polygon", "kept_in_region",
              "retrieved_at", "producer")

each = pytest.mark.parametrize("name", sorted(CONNECTORS))


def _around(fix, d):
    lon, lat = fix["decimalLongitude"], fix["decimalLatitude"]
    return Region.from_bbox(lon - d, lat - d, lon + d, lat + d)


@each
def test_normalise_returns_the_core_record(name):
    c = CONNECTORS[name]
    rec = c["normalise"](c["fixture"])
    assert _common.missing_core_fields(rec) == []
    assert rec["source"] == name
    assert isinstance(rec["record_id"], str) and rec["record_id"]
    assert isinstance(rec["taxonomy"], str) and rec["taxonomy"]
    assert rec["taxon_id"] is not None
    assert isinstance(rec["flags"], list)


@each
def test_search_returns_core_records_with_provenance_filtered_exactly(name):
    c = CONNECTORS[name]
    res = c["module"].search(_around(c["fixture"], 0.1), session=c["session"](count=2))
    assert isinstance(res, _common.Result)
    assert len(res.records) == 2
    assert all(_common.missing_core_fields(r) == [] for r in res.records)
    assert [k for k in PROVENANCE if k not in res.provenance] == []
    lon, lat = c["fixture"]["decimalLongitude"], c["fixture"]["decimalLatitude"]
    far = Region.from_bbox(lon + 1, lat + 1, lon + 2, lat + 2)
    assert c["module"].search(far, session=c["session"](count=1)).records == []


@each
def test_one_exception_type_for_too_many_records(name):
    c = CONNECTORS[name]
    with pytest.raises(_common.TooManyRecords):
        c["module"].search(_around(c["fixture"], 1), session=c["session"](count=10_000), max_records=10)


@each
def test_cells_come_from_the_shared_implementation(name):
    c = CONNECTORS[name]
    assert c["module"].to_cells is _common.to_cells
    row = c["module"].to_cells([c["normalise"](c["fixture"])], depth=8)[0]
    assert row["cell"] in row["footprint_cells"]


def test_records_from_different_repositories_mix():
    recs = [CONNECTORS[n]["normalise"](CONNECTORS[n]["fixture"]) for n in sorted(CONNECTORS)]
    rows = _common.to_cells(recs, depth=6)
    assert {r["source"] for r in rows} == set(CONNECTORS)
    assert len({(r["source"], r["record_id"]) for r in rows}) == len(rows)


def test_keep_in_region_drops_records_without_coordinates():
    rec = {"decimalLongitude": None, "decimalLatitude": 40.0}
    assert _common.keep_in_region([rec], Region.from_bbox(-180, -90, 180, 90)) == []


def test_missing_core_fields_names_what_is_absent():
    assert _common.missing_core_fields({"source": "X"})[0] == "record_id"
