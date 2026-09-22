import json
from pathlib import Path

import pytest

from healpix_connector.connectors import obis
from healpix_connector.region import Region

FIX = json.loads((Path(__file__).parent / "fixtures" / "obis_record_delphinus.json").read_text())


def test_normalise_keeps_identity_taxonomy_and_the_repositorys_own_flags():
    r = obis.normalise(FIX)
    assert r["obisID"] == FIX["id"] and r["occurrenceID"] == FIX["occurrenceID"]
    assert r["aphia_id"] == 137094 and r["backbone"] == "WoRMS"
    # OBIS's quality judgement travels with the record rather than being applied.
    assert r["flags"] == ["NO_DEPTH", "ON_LAND"]
    assert r["coordinateUncertaintyInMeters"] is None


class FakeResponse:
    def __init__(self, payload, status=200, headers=None):
        self.payload, self.status_code, self.headers = payload, status, headers or {}

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(self.status_code)


class FakeSession:
    """Answers /statistics with a count, then pages records by the ``after`` cursor."""

    def __init__(self, count, record=FIX):
        self.count, self.record, self.calls, self.served = count, record, [], 0

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append({"url": url, **(params or {})})
        if url.endswith("/statistics"):
            return FakeResponse({"records": self.count, "species": 1, "taxa": 1,
                                 "datasets": 1, "yearrange": [2000, 2026]})
        n = min(params["size"], self.count - self.served)
        self.served += n
        # each record needs a distinct id, so the cursor advances
        return FakeResponse({"results": [{**self.record, "id": f"{i}"} for i in range(n)]})


def test_search_counts_before_paging_and_filters_exactly():
    lon, lat = FIX["decimalLongitude"], FIX["decimalLatitude"]
    inside = Region.from_bbox(lon - 0.1, lat - 0.1, lon + 0.1, lat + 0.1)
    s = FakeSession(count=3)
    res = obis.search(inside, aphia_id=137094, session=s)
    assert s.calls[0]["url"].endswith("/statistics")  # counted, not paged blindly
    # the count must describe the query that gets paged, taxon filter included
    assert s.calls[0]["taxonid"] == 137094
    assert len(res.records) == 3
    assert res.provenance["params"]["taxonid"] == 137094
    assert res.provenance["backbone"] == "WoRMS"
    # No DOI exists for an OBIS query: the dataset ids and the timestamp stand in.
    assert res.provenance["query_doi"] is None
    assert res.provenance["datasets"] == [FIX["dataset_id"]]

    far = Region.from_bbox(lon + 1, lat + 1, lon + 2, lat + 2)
    assert obis.search(far, session=FakeSession(count=1)).records == []


def test_search_refuses_more_than_it_was_allowed():
    inside = Region.from_bbox(-9, 41, -8, 43)
    with pytest.raises(obis.TooManyRecords, match="100,000"):
        obis.search(inside, session=FakeSession(count=100_000), max_records=1_000)


def test_pages_follow_the_cursor_not_an_offset():
    inside = Region.from_bbox(-9, 41, -8, 43)
    s = FakeSession(count=obis.PAGE + 2)
    obis.search(inside, session=s, max_records=obis.PAGE * 2)
    pages = [c for c in s.calls if not c["url"].endswith("/statistics")]
    assert len(pages) == 2 and "offset" not in pages[1]
    assert pages[1]["after"] == str(obis.PAGE - 1)


def test_to_cells_marks_records_whose_uncertainty_is_unknown():
    rows = obis.to_cells([obis.normalise(FIX)], depth=8)
    assert rows[0]["cell"] >= 0 and rows[0]["footprint_known"] is False
    assert rows[0]["footprint_cells"] == [rows[0]["cell"]]
