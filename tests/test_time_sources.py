"""Sources with a time axis, and the live sources that build DGGS datasets in memory."""

import numpy as np
import pytest
import xarray as xr

from healpix_connector import sample as S
from healpix_connector.cells import assign_cells
from healpix_connector.dggs_zarr import cell_dataset, concat_time, month_bounds
from healpix_connector.sources.climatedt_monthly import YearOfMonths
from healpix_connector.sources.gridded import sample_cells_from_grid

DEPTH = 8
LON, LAT = np.array([-3.7, -3.5, -3.3]), np.array([40.4, 40.6, 40.8])
CELLS = np.asarray(assign_cells(LON, LAT, DEPTH), dtype=np.uint64)


def monthly(year=2020, depth=DEPTH, cells=CELLS):
    """value = 100 * month + position of the cell (0, 1, 2), so every lookup is checkable."""
    vals = 100.0 * np.arange(1, 13)[:, None] + np.arange(cells.size)[None, :]
    return cell_dataset(cells, depth, {"t": (vals, {"units": "degC"})},
                        {"title": f"test monthly {year}", "valid_time": "monthly mean"},
                        time_bounds=month_bounds(year, range(1, 13)))


@pytest.fixture(scope="module")
def src():
    return S.from_dataset(monthly())


@pytest.mark.parametrize("date, start, end", [
    ("2020", "2020-01-01", "2021-01-01"),
    (2020, "2020-01-01", "2021-01-01"),
    ("2020-02", "2020-02-01", "2020-03-01"),
    ("2020-02-29", "2020-02-29", "2020-03-01"),
    ("2020-02-29T23:10:00Z", "2020-02-29", "2020-03-01"),
    ("2020-05-01/2020-05-09", "2020-05-01", "2020-05-10"),
    ("2020-05-01/09", "2020-05-01", "2020-05-10"),
])
def test_event_dates_are_read_at_their_own_precision(date, start, end):
    s, e = S.date_interval(date)
    assert (str(s)[:10], str(e)[:10]) == (start, end)


@pytest.mark.parametrize("bad", ["20", "May 2020", "2020-13-01"])
def test_dates_that_cannot_be_read_are_refused(bad):
    with pytest.raises(ValueError):
        S.date_interval(bad)


def test_time_bounds_and_period_come_from_cf(src):
    assert src.has_time and src.time_bounds.shape == (12, 2)
    assert src.period == ("2020-01-01", "2020-12-31")


@pytest.mark.parametrize("date, step, status", [
    ("2020-05-17", 4, "inside time step 2020-05-01..2020-05-31"),
    ("2020-05-01/09", 4, "inside time step 2020-05-01..2020-05-31"),
    ("2020-05-30/2020-06-02", -1, "record date spans 2 time steps"),
    ("2020", -1, "record date spans 12 time steps"),
    ("2019-12-31", -1, "outside 2020-01-01..2020-12-31"),
    (None, -1, "record date unknown"),
    ("someday", -1, "record date not understood: 'someday'"),
])
def test_a_record_gets_the_one_step_that_contains_its_whole_date(src, date, step, status):
    assert S.time_step(date, src) == (step, status)
    assert S.time_status(date, src) == status


def test_a_time_source_will_not_guess_the_step(src):
    with pytest.raises(ValueError, match="time steps"):
        S.sample_cells(src, CELLS, DEPTH)
    assert S.sample_cells(src, CELLS, DEPTH, steps=2)["value"].tolist() == [300.0, 301.0, 302.0]


def test_records_are_sampled_at_their_own_month(src):
    recs = [{"cell": int(CELLS[1]), "eventDate": "2020-07-04"},
            {"cell": int(CELLS[2]), "eventDate": "2020", "footprint_cells": [int(c) for c in CELLS]},
            {"cell": int(CELLS[0]), "eventDate": None, "year": 2020},
            {"cell": int(CELLS[0]), "eventDate": "2020-03-02", "footprint_cells": [int(c) for c in CELLS]}]
    out = S.sample_records(recs, src, DEPTH)
    assert out[0]["value"] == 701.0 and out[0]["time_step"] == ["2020-07-01", "2020-07-31"]
    assert np.isnan(out[1]["value"]) and out[1]["time_step"] is None  # a year spans 12 monthly steps
    assert out[2]["time_status"] == "record date spans 12 time steps"     # falls back to `year`
    # the footprint is sampled at the record's own month too
    assert (out[3]["positional_min"], out[3]["positional_max"]) == (300.0, 302.0)


def test_coarser_and_finer_requests_keep_the_step(src):
    parent = np.asarray(assign_cells(LON[:1], LAT[:1], DEPTH - 1), dtype=np.uint64)
    got = S.sample_cells(src, parent, DEPTH - 1, steps=0)
    assert got["support"]["how"] == f"aggregated from depth {DEPTH}" and got["value"][0] == 100.0
    child = np.asarray(assign_cells(LON[:1], LAT[:1], DEPTH + 1), dtype=np.uint64)
    assert S.sample_cells(src, child, DEPTH + 1, steps=11)["value"][0] == 1200.0


def test_zarr_round_trip_keeps_the_time_axis(tmp_path):
    path = tmp_path / "m.zarr"
    monthly().to_zarr(path, zarr_format=3, consolidated=False)
    back = S.open_source(path)
    assert back.has_time and back.period == ("2020-01-01", "2020-12-31")
    assert S.time_step("2020-12-25", back)[0] == 11


def test_years_join_along_time():
    both = S.from_dataset(concat_time([monthly(2021), monthly(2020)]))
    assert both.time_bounds.shape == (24, 2) and both.period == ("2020-01-01", "2021-12-31")
    assert S.time_step("2021-01-15", both)[0] == 12
    with pytest.raises(ValueError, match="twice"):
        concat_time([monthly(2020), monthly(2020)])


def test_a_time_axis_without_bounds_is_refused():
    ds = monthly().drop_vars("time_bnds")
    with pytest.raises(ValueError, match="bounds"):
        S.from_dataset(ds)


def test_unsorted_cells_in_a_dataset_are_looked_up_correctly():
    ds = cell_dataset(CELLS, DEPTH, {"v": (np.array([1.0, 2.0, 3.0]), {"units": "1"})}, {"title": "x"})
    rev = ds.isel(cells=slice(None, None, -1))  # cells now in descending order on disk
    got = S.sample_cells(S.from_dataset(rev), CELLS, DEPTH)["value"]
    assert got.tolist() == [1.0, 2.0, 3.0]


def test_a_coarse_grid_sample_carries_its_support_into_every_value():
    lon = np.arange(-5.0, -2.0, 0.25)
    lat = np.arange(42.0, 39.0, -0.25)
    vals = np.broadcast_to(lat[:, None], (lat.size, lon.size)).copy()
    cs = sample_cells_from_grid(vals, lon, lat, CELLS, DEPTH, 0.25)
    src = S.from_dataset(cs.to_dataset(DEPTH, "lat", units="degrees_north", period=("1981", "2010"),
                                       valid_time="climatology"))
    got = S.sample_cells(src, CELLS, DEPTH)
    assert got["support"]["support_m"] == pytest.approx(cs.support_m)
    assert np.array_equal(got["value"], cs.value)
    assert S.provenance(src, DEPTH)["support_m"] == pytest.approx(cs.support_m)
    assert S.time_status(1995, src) == "inside the 1981-2010 climatology"


def test_climate_dt_generation_2_year_becomes_a_monthly_source():
    n = CELLS.size
    year = YearOfMonths(year=2030, cell_ids=CELLS, values={"avg_2t": np.full((12, n), 15.0)},
                        units={"avg_2t": "degC"}, valid_share={"avg_2t": np.ones((12, n))},
                        method="bilinear_conservative", depth=DEPTH,
                        requests=[{"model": "IFS-FESOM", "experiment": "SSP3-7.0"}])
    ds = year.to_dataset()
    assert set(ds.data_vars) >= {"avg_2t", "avg_2t_valid_share", "time_bnds", "crs"}
    assert ds.attrs["dggs"]["refinement_level"] == DEPTH
    src = S.from_dataset(ds)  # the valid share is a companion, never the sampled variable
    assert src.variable == "avg_2t" and src.time_bounds.shape == (12, 2)
    assert "IFS-FESOM SSP3-7.0" in src.title and "bilinear_conservative" in src.resampling_method
    assert isinstance(ds, xr.Dataset)
