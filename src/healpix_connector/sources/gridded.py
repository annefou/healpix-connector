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

    def to_dataset(self, depth: int, name: str, *, units: str, long_name: str | None = None,
                   title: str | None = None, source_doi: str | None = None, license: str | None = None,
                   period: tuple[str, str] | None = None, valid_time: str | None = None):
        """This sample as a DGGS dataset, ready for ``sample.from_dataset``.

        ``support_m`` travels with it, so every value sampled from it says that it
        describes the source box rather than the cell.
        """
        from healpix_connector import __version__
        from healpix_connector.dggs_zarr import cell_dataset

        attrs = {"title": title or f"{name} on HEALPix depth {depth}", "support_m": float(self.support_m),
                 "resampling_method": self.method, "producer": f"healpix-connector {__version__}"}
        if source_doi:
            attrs["source_doi"] = source_doi
        if license:
            attrs["license"] = license
        if period:
            attrs["climatology_bounds" if valid_time == "climatology" else "time_bounds"] = list(period)
        if valid_time:
            attrs["valid_time"] = valid_time
        var_attrs = {"units": units, **({"long_name": long_name} if long_name else {})}
        return cell_dataset(self.cell_ids, depth, {name: (self.value, var_attrs)}, attrs)


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
    # Every cell centre must fall inside a source box: the nearest-point lookup below would otherwise return
    # the edge value for a cell outside the grid (or for a grid in 0-360 longitudes and a cell west of 0).
    half = res_deg / 2
    outside = ((clat < lat.min() - half) | (clat > lat.max() + half)
               | (clon < lon.min() - half) | (clon > lon.max() + half))
    if outside.any():
        raise ValueError(f"{int(outside.sum())} of {cell_ids.size} cell centres lie outside the source grid "
                         f"(lon {lon.min()}..{lon.max()}, lat {lat.min()}..{lat.max()}, longitudes -180..180 "
                         "expected); request a larger area or drop those cells")
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
