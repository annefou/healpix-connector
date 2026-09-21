import numpy as np
import pytest

from healpix_connector.binning import bin_to_cells, check_binning_valid
from healpix_connector.conventions import ELLIPSOID

RES = 1 / 120  # 30 arc-seconds, like CHELSA


def grid(lon0, lon1, lat0, lat1, res=RES):
    lon = np.arange(lon0 + res / 2, lon1, res)
    lat = np.arange(lat1 - res / 2, lat0, -res)  # north first, like GeoTIFF rows
    return lon, lat


def interior(stats, lon0, lon1, lat0, lat1, margin=0.5):
    """Cells whose centre is well inside the grid, so they are fully covered."""
    from healpix_geo import nested

    clon, clat = nested.healpix_to_lonlat(stats.cell_ids, np.uint8(stats.depth), ellipsoid=ELLIPSOID)
    clon = np.where(clon > 180, clon - 360, clon)
    return ((clon > lon0 + margin) & (clon < lon1 - margin)
            & (clat > lat0 + margin) & (clat < lat1 - margin)), clat


def test_constant_field_has_no_spread_and_full_coverage():
    lon, lat = grid(-4, -2, 40, 42)
    vals = np.full((lat.size, lon.size), 14.0)
    s = bin_to_cells(vals, lon, lat, depth=8, res_deg=RES)
    inside, _ = interior(s, -4, -2, 40, 42)
    assert inside.any()
    assert np.allclose(s.mean, 14.0)
    assert np.allclose(s.std, 0.0)
    assert np.allclose(s.coverage[inside], 1.0, atol=0.02)
    assert s.pixel_count.sum() == vals.size


def test_linear_field_mean_matches_cell_centre_and_spread_is_nonzero():
    lon, lat = grid(-4, -2, 40, 42)
    LAT = np.broadcast_to(lat[:, None], (lat.size, lon.size))
    vals = 20.0 - 0.5 * LAT  # 0.5 degC per degree of latitude
    s = bin_to_cells(vals, lon, lat, depth=8, res_deg=RES)
    inside, clat = interior(s, -4, -2, 40, 42)
    # The area-weighted mean of a linear field is its value at the cell centroid.
    assert np.allclose(s.mean[inside], 20.0 - 0.5 * clat[inside], atol=0.01)
    assert (s.std[inside] > 0.01).all()
    assert (s.minimum <= s.mean).all() and (s.mean <= s.maximum).all()


def test_missing_pixels_reduce_coverage():
    lon, lat = grid(-4, -2, 40, 42)
    vals = np.full((lat.size, lon.size), 1.0)
    vals[:, : lon.size // 2] = np.nan
    s = bin_to_cells(vals, lon, lat, depth=8, res_deg=RES)
    assert s.pixel_count.sum() == np.isfinite(vals).sum()
    assert s.coverage.max() <= 1.02


def test_depth_too_fine_for_binning_is_refused():
    with pytest.raises(ValueError, match="overlap-conservative"):
        check_binning_valid(depth=12, res_deg=RES, lat=40.0)
    check_binning_valid(depth=9, res_deg=RES, lat=40.0)  # ~13 km cells: fine
