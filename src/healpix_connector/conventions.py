"""GRID4EARTH conventions, stated once so no module can drift from them.

- NESTED ordering only: the multi-resolution join relies on ``parent = cell >> 2``.
- WGS84 always, passed explicitly: healpix-geo defaults to ``"sphere"``.
- Resolution is called *depth* here. healpix-resample calls it ``level`` and the
  zarr-conventions/dggs metadata calls it ``refinement_level``; all three are
  the same number (nside = 2**depth).
"""

from __future__ import annotations

import math

ELLIPSOID = "WGS84"
INDEXING_SCHEME = "nested"
MAX_DEPTH = 29

# Mean Earth radius (m), as healpix-resample's ``cell_size_m`` uses.
_EARTH_RADIUS_M = 6_371_000.0


def cell_size_m(depth: int) -> float:
    """Equal-area HEALPix cell width at ``depth``, in metres.

    Same formula as ``healpix_resample.cell_size_m``: the square root of one
    cell's solid angle, times the Earth radius.
    """
    if not 0 <= depth <= MAX_DEPTH:
        raise ValueError(f"depth must be in [0, {MAX_DEPTH}], got {depth}")
    return math.sqrt(4.0 * math.pi / (12.0 * 4.0**depth)) * _EARTH_RADIUS_M
