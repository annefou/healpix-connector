"""Climate DT source: the parts that need no Polytope call.

The converter itself belongs to healpix-convert and is tested there; what is
ours is the month-long request, the cache name that survives a date range, and
the parameter table.
"""

import numpy as np
import pytest

pytest.importorskip("healpix_convert")

from healpix_connector.sources import climatedt  # noqa: E402


def test_source_is_read_from_healpix_convert():
    from healpix_convert.settings import climatedt as cfg

    assert climatedt.SOURCE["experiment"] == cfg.CDT_EXPERIMENT
    assert climatedt.SOURCE["model"] == cfg.CDT_MODEL
    assert cfg.POLYTOPE_ADDRESS in climatedt.SOURCE["access"]


def test_param_ids_cover_the_converters_variables_and_precipitation():
    from healpix_convert.settings import climatedt as cfg

    for name, (param, *_ ) in cfg.CDT_SFC_VARIABLE_META.items():
        assert climatedt.PARAM_ID[name] == param
    assert climatedt.PARAM_ID["tprate"] == "260048"


def test_a_date_range_gives_a_usable_cache_name(tmp_path):
    conv = climatedt.MonthConverter(date="20300101/to/20300131", time="0000/1200",
                                    params="167", local_dir=tmp_path)
    path = conv._grib_path()
    assert "/" not in path.name and path.parent == tmp_path
    assert path.name == "cdt_20300101-to-20300131T0000-1200_sfc.grib"


def test_the_request_is_the_converters_own(tmp_path):
    from healpix_convert.settings import climatedt as cfg

    conv = climatedt.MonthConverter(date="20300101/to/20300131", time="0000/1200",
                                    params="167/260048", local_dir=tmp_path)
    assert (conv.experiment, conv.expver, conv.model, conv.resolution) == (
        cfg.CDT_EXPERIMENT, cfg.CDT_EXPVER, cfg.CDT_MODEL, cfg.CDT_RESOLUTION)
    assert conv.levtype == "sfc"


def test_sampler_declares_the_source_cell_as_support():
    from healpix_connector.conventions import cell_size_m

    sampler = climatedt.MonthlySampler([10786, 10792], ["t2m"])
    assert sampler._resampler is None  # built on the first month, not before
    assert round(cell_size_m(climatedt.DEPTH) / 1000) == 51


def _sphere_grid_and_wgs84_cells():
    """Climate DT's native grid near Iberia (depth 7 on a sphere), and WGS84 cells inside it."""
    from types import SimpleNamespace

    from healpix_geo import nested

    from healpix_connector.conventions import ELLIPSOID

    lon, lat = np.meshgrid(np.arange(-14.0, 8.0, 0.25), np.arange(32.0, 48.0, 0.25))
    native = np.unique(nested.lonlat_to_healpix(lon.ravel(), lat.ravel(), climatedt.DEPTH,
                                                ellipsoid="sphere"))
    nlon, nlat = nested.healpix_to_lonlat(native, climatedt.DEPTH, ellipsoid="sphere")
    nlon = (np.asarray(nlon) + 180.0) % 360.0 - 180.0  # healpix-geo returns 0..360
    lon, lat = np.meshgrid(np.arange(-9.0, 3.0, 0.25), np.arange(36.0, 44.0, 0.25))
    cells = np.unique(nested.lonlat_to_healpix(lon.ravel(), lat.ravel(), climatedt.DEPTH,
                                               ellipsoid=ELLIPSOID)).astype(np.int64)
    grid = SimpleNamespace(longitude=SimpleNamespace(values=np.asarray(nlon)),
                           latitude=SimpleNamespace(values=np.asarray(nlat)))
    return grid, np.asarray(nlon), np.asarray(nlat), cells


def test_psf_on_precipitation_is_a_bounded_fit():
    """A sharp edge makes the unbounded PSF ring below zero; ``tprate`` must not."""
    pytest.importorskip("healpix_resample")
    grid, lon, lat, cells = _sphere_grid_and_wgs84_cells()
    # a rain band with a sharp edge, zero elsewhere: the worst case for a deconvolution
    rain = np.where((lon > -4.0) & (lon < -1.0), 1e-4, 0.0)

    sampler = climatedt.MonthlySampler(cells, ["tprate", "u10"], method="psf")
    sampler._build_resampler(grid)
    free, n_free = sampler._resample("u10", rain)  # a signed variable: no bound
    bounded, n_held = sampler._resample("tprate", rain)

    assert n_free is None
    # the ringing this guards against: measured -0.67 x the band, 148 of 381 cells
    assert np.nanmin(free) < -0.3e-4
    assert np.nanmin(bounded) >= 0
    assert n_held >= np.sum(free < 0)
    # the bound reshapes the field near the edges only; 3 degrees away it agrees
    from healpix_geo import nested

    clon, _ = nested.healpix_to_lonlat(cells.astype(np.uint64), climatedt.DEPTH,
                                       ellipsoid="WGS84")
    clon = (np.asarray(clon) + 180.0) % 360.0 - 180.0
    far = (np.abs(clon + 4.0) > 3.0) & (np.abs(clon + 1.0) > 3.0)
    assert far.sum() > 50
    np.testing.assert_allclose(bounded[far], free[far], atol=1e-6)


def test_bilinear_conservative_reports_no_clipping():
    grid, lon, lat, cells = _sphere_grid_and_wgs84_cells()
    sampler = climatedt.MonthlySampler(cells, ["tprate"], method="bilinear_conservative")
    sampler._build_resampler(grid)
    values, n = sampler._resample("tprate", np.full(lon.size, 2e-5))
    assert n is None
    np.testing.assert_allclose(values, 2e-5, rtol=1e-6)
