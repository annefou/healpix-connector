"""A coarse regular grid sampled onto HEALPix cells.

When the source grid box is comparable to or larger than the target cell
(ERA5's 0.25 degrees against a depth-8 cell, say), binning is meaningless: a
cell contains no source pixels at all. The honest operation is to take the
value of the source box containing the cell's centre and to declare the
support as *the source box*, not the cell. ``support_m`` says so, and callers
carry it into the result rather than implying the value describes the cell.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from healpix_connector.binning import pixel_size_m
from healpix_connector.conventions import ELLIPSOID, cell_size_m


@dataclass(frozen=True)
class CoarseSample:
    cell_ids: np.ndarray
    value: np.ndarray
    support_m: float  # width of the source box, not of the cell
    cell_width_m: float
    method: str


def sample_cells_from_grid(values: np.ndarray, lon: np.ndarray, lat: np.ndarray,
                           cell_ids: np.ndarray, depth: int, res_deg: float) -> CoarseSample:
    """Value of the source box containing each cell's centre.

    ``values`` has shape (len(lat), len(lon)) on a regular lon/lat grid.
    """
    from healpix_geo import nested

    values = np.asarray(values, dtype=np.float64)
    if values.shape != (lat.size, lon.size):
        raise ValueError(f"values shape {values.shape} != (lat, lon) {(lat.size, lon.size)}")
    cell_ids = np.asarray(cell_ids, dtype=np.uint64)
    clon, clat = nested.healpix_to_lonlat(cell_ids, np.uint8(depth), ellipsoid=ELLIPSOID)
    clon = np.where(clon > 180, clon - 360, clon)
    i = np.abs(lat[:, None] - clat[None, :]).argmin(axis=0)
    j = np.abs(lon[:, None] - clon[None, :]).argmin(axis=0)
    support = pixel_size_m(res_deg, float(np.mean(lat)))
    width = cell_size_m(depth)
    return CoarseSample(
        cell_ids=cell_ids, value=values[i, j], support_m=support, cell_width_m=width,
        method=(f"value of the source grid box containing the cell centre; support is the "
                f"source box (~{support/1000:.0f} km), wider than the cell (~{width/1000:.0f} km)"
                if support >= width else
                f"value of the source grid box containing the cell centre (~{support/1000:.0f} km)"),
    )
