import json
from pathlib import Path

import pytest

from healpix_connector.connectors import gbif
from healpix_connector.region import Region

FIX = json.loads((Path(__file__).parent / "fixtures" / "gbif_record_calopteryx.json").read_text())


def test_normalise_keeps_identity_and_the_requested_taxonomy():
    col = gbif.normalise(FIX, gbif.COL_XR)
    legacy = gbif.normalise(FIX, gbif.LEGACY_BACKBONE)
    assert col["gbifID"] == str(FIX["gbifID"])
    assert col["taxon_key"] == "Q2M4" and col["checklist_key"] == gbif.COL_XR
    assert legacy["taxon_key"] == "1427067"


class FakeResponse:
    def __init__(self, payload, status=200, headers=None):
        self.payload, self.status_code, self.headers = payload, status, headers or {}

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(self.status_code)


class FakeSession:
    """Returns ``count`` then pages of the fixture record; records every call."""

    def __init__(self, count, record=FIX):
        self.count, self.record, self.calls = count, record, []

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append(dict(params or {}))
        if params.get("limit") == 0:
            return FakeResponse({"count": self.count})
        n = min(params["limit"], self.count - params.get("offset", 0))
        return FakeResponse({"results": [self.record] * n, "endOfRecords": params.get("offset", 0) + n >= self.count})


def test_search_always_sends_checklist_key_and_filters_exactly():
    lon, lat = FIX["decimalLongitude"], FIX["decimalLatitude"]
    inside = Region.from_bbox(lon - 0.1, lat - 0.1, lon + 0.1, lat + 0.1)
    s = FakeSession(count=2)
    res = gbif.search(inside, taxon_key="Q2M4", session=s)
    assert all(c["checklistKey"] == gbif.COL_XR for c in s.calls)
    assert len(res.records) == 2
    assert res.provenance["checklist_key"] == gbif.COL_XR
    assert res.provenance["params"]["taxonKey"] == "Q2M4"
    # A region the record is outside of: returned by the query polygon, dropped by the exact test.
    far = Region.from_bbox(lon + 1, lat + 1, lon + 2, lat + 2)
    assert gbif.search(far, session=FakeSession(count=1)).records == []


def test_search_refuses_more_than_the_api_can_page():
    with pytest.raises(gbif.TooManyRecords, match="download"):
        gbif.search(Region.from_bbox(-10, 35, 4, 44), session=FakeSession(count=100_001),
                    max_records=gbif.SEARCH_LIMIT)


def test_search_defaults_to_the_practical_paging_limit():
    region = Region.from_bbox(-10, 35, 4, 44)
    with pytest.raises(gbif.TooManyRecords, match="Deep paging"):
        gbif.search(region, session=FakeSession(count=gbif.PRACTICAL_PAGING_LIMIT + 1))


def test_search_respects_a_lower_max_records_cap():
    region = Region.from_bbox(-10, 35, 4, 44)
    with pytest.raises(gbif.TooManyRecords, match="this call allows 500"):
        gbif.search(region, session=FakeSession(count=5_000), max_records=500)


def test_to_cells_adds_cell_and_footprint():
    rec = gbif.normalise(FIX, gbif.COL_XR)
    out = gbif.to_cells([rec], depth=8)[0]
    assert out["depth"] == 8 and out["cell"] in out["footprint_cells"]
    assert out["footprint_known"] == (rec["coordinateUncertaintyInMeters"] is not None)


def test_match_name_uses_v2_with_checklist_and_reports_match_type():
    class S:
        def __init__(self):
            self.calls = []

        def get(self, url, params=None, timeout=None, headers=None):
            self.calls.append((url, dict(params)))
            return FakeResponse({"usage": {"key": "5Z5T3", "name": "Ciconia ciconia (Linnaeus, 1758)", "rank": "SPECIES"},
                                 "diagnostics": {"matchType": "EXACT"}})

    s = S()
    m = gbif.match_name("Ciconia ciconia", session=s)
    assert m["taxon_key"] == "5Z5T3" and m["match_type"] == "EXACT"
    url, params = s.calls[0]
    assert "/v2/species/match" in url and params["checklistKey"] == gbif.COL_XR


def test_rate_limit_and_server_errors_are_retried_honouring_retry_after():
    responses = [FakeResponse({}, 429, {"Retry-After": "7"}), FakeResponse({}, 503), FakeResponse({"ok": 1})]

    class S:
        def get(self, url, params=None, timeout=None, headers=None):
            return responses.pop(0)

    waits = []
    assert gbif._get(S(), "u", sleep=waits.append) == {"ok": 1}
    assert waits == [7, 2]  # Retry-After wins over the 1 s backoff; then 2 s
