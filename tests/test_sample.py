import numpy as np
import pytest

from healpix_connector import sample as S
from healpix_connector.binning import bin_to_cells
from healpix_connector.cells import assign_cells
from healpix_connector.sources import chelsa

RES = 1 / 120


@pytest.fixture(scope="module")
def source(tmp_path_factory):
    """A small CHELSA-like source over Madrid with a latitudinal gradient."""
    lon = np.arange(-4.2958, -3.0, RES)
    lat = np.arange(41.2958, 40.0, -RES)
    vals = 20.0 - 0.5 * np.broadcast_to(lat[:, None], (lat.size, lon.size))
    ds = chelsa.to_dataset(bin_to_cells(vals, lon, lat, 8, RES), 1)
    path = chelsa.write(ds, tmp_path_factory.mktemp("src"), (-4.3, 40.0, -3.0, 41.3))
    return S.open_source(path)


def test_open_source_reads_grid_and_period(source):
    assert source.depth == 8 and source.variable == "bio1"
    assert source.valid_time == "climatology"
    assert source.period == ("1981-01-01", "2010-12-31")
    assert source.doi == "10.16904/envidat.228"


def test_direct_sampling_returns_the_cells_own_value(source):
    cells = assign_cells([-3.7, -3.5], [40.4, 40.9], 8)
    got = S.sample_cells(source, cells, 8)
    assert got["support"]["how"] == "direct"
    # A latitudinal gradient: the northern cell must be colder.
    assert got["value"][0] > got["value"][1]
    assert np.isfinite(got["value"]).all()


def test_finer_request_inherits_the_coarser_value_and_says_so(source):
    cell8 = assign_cells([-3.7], [40.4], 8)
    cells9 = assign_cells([-3.7], [40.4], 9)
    at8 = S.sample_cells(source, cell8, 8)["value"][0]
    at9 = S.sample_cells(source, cells9, 9)
    assert at9["value"][0] == pytest.approx(at8)
    assert at9["support"]["how"] == "inherited from depth 8"


def test_coarser_request_aggregates_children_and_says_so(source):
    cell7 = assign_cells([-3.7], [40.4], 7)
    got = S.sample_cells(source, cell7, 7)
    assert got["support"]["how"] == "aggregated from depth 8"
    kids = np.asarray(cell7, dtype=np.uint64)[0] * 4 + np.arange(4, dtype=np.uint64)
    child_vals = S.sample_cells(source, kids, 8)["value"]
    assert got["value"][0] == pytest.approx(np.nanmean(child_vals))


def test_cells_outside_the_source_are_nan(source):
    far = assign_cells([10.0], [50.0], 8)  # not in the Madrid window
    assert np.isnan(S.sample_cells(source, far, 8)["value"][0])


@pytest.mark.parametrize(
    "date, expected",
    [("1995-06-01", "inside the 1981-2010 climatology"),
     ("2021-06-01", "outside the 1981-2010 climatology"),
     (2021, "outside the 1981-2010 climatology"),
     (None, "record date unknown")],
)
def test_time_status_compares_the_record_date_with_the_period(source, date, expected):
    assert S.time_status(date, source) == expected


def test_record_sampling_reports_positional_spread_only_across_its_footprint(source):
    lon, lat = -3.7, 40.4
    own = int(assign_cells([lon], [lat], 8)[0])
    north = int(assign_cells([lon], [lat + 0.5], 8)[0])
    records = [
        {"cell": own, "footprint_cells": [own], "footprint_known": True, "eventDate": "2021-06-01"},
        {"cell": own, "footprint_cells": [own, north], "footprint_known": True, "eventDate": "1995-06-01"},
        {"cell": own, "footprint_cells": [own], "footprint_known": False, "eventDate": None},
    ]
    out = S.sample_records(records, source, 8)
    assert out[0]["positional_spread"] == 0.0
    # Two cells: the spread is half the gap between them, and it is not zero.
    gap = out[1]["positional_max"] - out[1]["positional_min"]
    assert gap > 0
    assert out[1]["positional_spread"] == pytest.approx(gap / 2)
    assert out[0]["time_status"].startswith("outside") and out[1]["time_status"].startswith("inside")
    assert out[2]["positional_known"] is False and out[2]["time_status"] == "record date unknown"
    assert out[0]["source_std"] >= 0 and out[0]["support"]["depth"] == 8


def test_provenance_carries_what_a_citation_needs(source):
    p = S.provenance(source, 8)
    assert p["dataset_doi"] == "10.16904/envidat.228" and p["license"] == "CC0-1.0"
    assert p["period"] == ["1981-01-01", "2010-12-31"] and p["sampled_at_depth"] == 8
    assert "binning" in p["resampling_method"]
