"""Area-weighted binning of a fine lon/lat grid into HEALPix cells.

Each source pixel is assigned to the WGS84 HEALPix cell that contains its
centre, and per-cell statistics are computed with the pixel's area as weight.
This is accurate when source pixels are much smaller than the target cells:
the only error is at cell edges, where a pixel is counted wholly in one cell.
``MIN_CELL_TO_PIXEL_RATIO`` enforces that; below it, use an exact overlap
method (healpix-resample's OverlapConservativeResampler, with authalic-latitude
bounds so its cells match healpix-geo's WGS84 cells).

Unlike a conservative mean alone, binning also yields the spread of the source
field inside each cell: the part of the signal that one number per cell hides.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from healpix_connector.conventions import ELLIPSOID, cell_size_m

# A cell must be at least this many source pixels wide for centre-binning.
MIN_CELL_TO_PIXEL_RATIO = 10.0

# Mean Earth radius (m) used for pixel areas; relative weights only.
_R = 6_371_007.2


@dataclass(frozen=True)
class CellStats:
    """Per-cell statistics, all arrays aligned on ``cell_ids`` (sorted)."""

    depth: int
    cell_ids: np.ndarray  # uint64
    mean: np.ndarray
    std: np.ndarray  # area-weighted population standard deviation
    minimum: np.ndarray
    maximum: np.ndarray
    pixel_count: np.ndarray  # valid source pixels in the cell
    coverage: np.ndarray  # valid pixel area / cell area (1.0 = fully covered)


def pixel_size_m(res_deg: float, lat: float = 0.0) -> float:
    """Equivalent square width (sqrt of area) of a ``res_deg`` pixel at ``lat``, in metres."""
    return np.radians(res_deg) * _R * np.cos(np.radians(lat)) ** 0.5


def check_binning_valid(depth: int, res_deg: float, lat: float = 0.0) -> None:
    ratio = cell_size_m(depth) / pixel_size_m(res_deg, lat)
    if ratio < MIN_CELL_TO_PIXEL_RATIO:
        raise ValueError(
            f"depth {depth} cells are only {ratio:.1f} source pixels wide "
            f"(< {MIN_CELL_TO_PIXEL_RATIO:.0f}); centre-binning would be inaccurate. "
            "Use an exact overlap-conservative remap instead."
        )


def bin_to_cells(
    values: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    depth: int,
    res_deg: float,
) -> CellStats:
    """Bin a regular lon/lat grid into HEALPix cells at ``depth``.

    ``values`` has shape (len(lat), len(lon)); ``lon``/``lat`` are pixel-centre
    coordinates in degrees; ``res_deg`` is the pixel size. NaN pixels are
    ignored and reduce ``coverage``.
    """
    from healpix_geo import nested

    values = np.asarray(values, dtype=np.float64)
    if values.shape != (lat.size, lon.size):
        raise ValueError(f"values shape {values.shape} != (lat, lon) {(lat.size, lon.size)}")
    check_binning_valid(depth, res_deg, float(np.mean(lat)))

    LON, LAT = np.meshgrid(lon.astype(np.float64), lat.astype(np.float64))
    ok = np.isfinite(values)
    v, lo, la = values[ok], LON[ok], LAT[ok]
    # Pixel area on the sphere of radius _R: res^2 * cos(lat).
    area = (np.radians(res_deg) * _R) ** 2 * np.cos(np.radians(la))

    cells = nested.lonlat_to_healpix(lo, la, np.uint8(depth), ellipsoid=ELLIPSOID)
    ids, inv = np.unique(cells, return_inverse=True)

    wsum = np.bincount(inv, weights=area)
    mean = np.bincount(inv, weights=area * v) / wsum
    var = np.bincount(inv, weights=area * (v - mean[inv]) ** 2) / wsum
    count = np.bincount(inv).astype(np.int64)
    vmin = np.full(ids.size, np.inf)
    vmax = np.full(ids.size, -np.inf)
    np.minimum.at(vmin, inv, v)
    np.maximum.at(vmax, inv, v)
    cell_area = 4.0 * np.pi * _R**2 / (12.0 * 4.0**depth)

    return CellStats(
        depth=depth,
        cell_ids=ids.astype(np.uint64),
        mean=mean,
        std=np.sqrt(np.maximum(var, 0.0)),
        minimum=vmin,
        maximum=vmax,
        pixel_count=count,
        coverage=wsum / cell_area,
    )
