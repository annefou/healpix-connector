"""Records to HEALPix cells, with positional-uncertainty footprints."""

from __future__ import annotations

import numpy as np

from healpix_connector.conventions import ELLIPSOID

# Metres per degree of great-circle arc on the mean Earth sphere.
_M_PER_DEG = 6_371_008.8 * np.pi / 180.0


def assign_cells(lon, lat, depth: int) -> np.ndarray:
    """WGS84 NESTED cell of each point (NaN coordinates are not allowed)."""
    from healpix_geo import nested

    return nested.lonlat_to_healpix(np.asarray(lon, dtype=np.float64),
                                    np.asarray(lat, dtype=np.float64),
                                    np.uint8(depth), ellipsoid=ELLIPSOID)


def footprint_cells(lon: float, lat: float, radius_m: float | None, depth: int) -> tuple[np.ndarray, bool]:
    """Cells covered by a record's uncertainty disc.

    Returns ``(cells, known)``. With no stated uncertainty the footprint is the
    record's own cell and ``known`` is False, so callers can flag it rather than
    silently treating the point as exact.
    """
    from healpix_geo import nested

    own = assign_cells([lon], [lat], depth)
    if radius_m is None or not np.isfinite(radius_m) or radius_m <= 0:
        return own, False
    cov = nested.cone_coverage((float(lon), float(lat)), float(radius_m) / _M_PER_DEG,
                               np.uint8(depth), ellipsoid=ELLIPSOID, flat=True)
    cells = np.asarray(cov[0], dtype=np.uint64)
    return np.union1d(cells, own), True
