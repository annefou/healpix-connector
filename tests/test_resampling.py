"""Conservative bilinear onto WGS84 cells, on synthetic fields (no network)."""

import numpy as np
import pytest

pytest.importorskip("healpix_resample")
hn = pytest.importorskip("healpix_geo.nested")

from healpix_connector.resampling import BilinearConservative  # noqa: E402

DEPTH = 10


def _cells(lon0=33.0, lat0=38.0, half=0.5):
    ring = np.array([[lon0 - half, lat0 - half], [lon0 + half, lat0 - half],
                     [lon0 + half, lat0 + half], [lon0 - half, lat0 + half]])
    ids, _, _ = hn.polygon_coverage(ring, DEPTH, ellipsoid="WGS84", flat=True)
    return np.unique(np.asarray(ids, np.uint64))


def _sphere_source(cells):
    clon, clat = hn.healpix_to_lonlat(cells, DEPTH, ellipsoid="WGS84")
    seed = np.asarray(hn.lonlat_to_healpix(np.asarray(clon), np.asarray(clat), np.uint8(DEPTH),
                                           ellipsoid="sphere"), np.uint64)
    src = np.unique(np.asarray(hn.kth_neighbourhood(seed, DEPTH, 2), np.uint64).ravel())
    slon, slat = hn.healpix_to_lonlat(src, DEPTH, ellipsoid="sphere")
    return np.asarray(slon), np.asarray(slat)


def test_constant_field_stays_constant_and_every_cell_is_filled():
    cells = _cells()
    slon, slat = _sphere_source(cells)
    op = BilinearConservative(slon, slat, cells, DEPTH)
    mean, share = op(np.full(slon.size, 12.5))
    assert np.isfinite(mean).all()
    np.testing.assert_allclose(mean, 12.5, atol=1e-9)
    assert (share > 0.99).all()


def test_missing_samples_are_excluded_not_counted_as_zero():
    cells = _cells()
    slon, slat = _sphere_source(cells)
    op = BilinearConservative(slon, slat, cells, DEPTH)
    v = np.full(slon.size, 3.0)
    v[::3] = np.nan
    mean, share = op(v)
    ok = np.isfinite(mean)
    np.testing.assert_allclose(mean[ok], 3.0, atol=1e-9)
    assert share[ok].min() < 1.0                  # cells next to a gap report a partial share


def test_sphere_to_wgs84_shift_is_corrected():
    """A field that varies with geodetic latitude: matching ids would misplace it by ~0.12 deg."""
    cells = _cells()
    slon, slat = _sphere_source(cells)
    f = lambda lat: 10.0 * (lat - 38.0)          # 10 units per degree of latitude
    op = BilinearConservative(slon, slat, cells, DEPTH)
    mean, _ = op(f(slat))
    clon, clat = hn.healpix_to_lonlat(cells, DEPTH, ellipsoid="WGS84")
    truth = f(np.asarray(clat))
    _, same_id_lat = hn.healpix_to_lonlat(cells, DEPTH, ellipsoid="sphere")
    same_id_error = np.abs(f(np.asarray(same_id_lat)) - truth).mean()
    assert same_id_error > 1.0                   # ~1.2 units: the shift is real
    assert np.abs(mean - truth).mean() < 0.2 * same_id_error


def test_unsorted_cells_are_refused():
    cells = _cells()
    slon, slat = _sphere_source(cells)
    with pytest.raises(ValueError):
        BilinearConservative(slon, slat, cells[::-1], DEPTH)
