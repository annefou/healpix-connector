"""Regions: a bounding box, a polygon, or a set of HEALPix cells.

Repositories such as GBIF and OBIS know nothing about HEALPix. A region is
therefore turned into a (slightly padded) query polygon for the repository, and
the returned records are then filtered exactly with ``contains``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from healpix_connector.conventions import ELLIPSOID


@dataclass(frozen=True)
class Region:
    kind: str  # "bbox" | "polygon" | "cells"
    bbox: tuple[float, float, float, float] | None = None  # lon_min, lat_min, lon_max, lat_max
    polygon: tuple[tuple[float, float], ...] | None = None  # (lon, lat) vertices, not closed
    depth: int | None = None
    cells: np.ndarray | None = field(default=None, compare=False)

    @classmethod
    def from_bbox(cls, lon_min, lat_min, lon_max, lat_max) -> "Region":
        if not (lon_min < lon_max and lat_min < lat_max):
            raise ValueError("bbox must be (lon_min, lat_min, lon_max, lat_max) with min < max")
        return cls("bbox", bbox=(float(lon_min), float(lat_min), float(lon_max), float(lat_max)))

    @classmethod
    def from_polygon(cls, vertices) -> "Region":
        v = [(float(x), float(y)) for x, y in vertices]
        if len(v) > 1 and v[0] == v[-1]:
            v = v[:-1]
        if len(v) < 3:
            raise ValueError("a polygon needs at least 3 distinct vertices")
        return cls("polygon", polygon=tuple(v))

    @classmethod
    def from_cells(cls, cell_ids, depth: int) -> "Region":
        ids = np.unique(np.asarray(cell_ids, dtype=np.uint64))
        if ids.size == 0:
            raise ValueError("empty cell set")
        return cls("cells", depth=int(depth), cells=ids)

    # ------------------------------------------------------------------ query
    def query_bounds(self, pad_deg: float = 0.0) -> tuple[float, float, float, float]:
        """Lon/lat bounds enclosing the region, padded by ``pad_deg``."""
        if self.kind == "bbox":
            x0, y0, x1, y1 = self.bbox
        elif self.kind == "polygon":
            xs, ys = zip(*self.polygon)
            x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
        else:
            from healpix_geo import nested

            lon, lat = nested.vertices(self.cells, np.uint8(self.depth), ellipsoid=ELLIPSOID)
            lon = np.where(lon > 180, lon - 360, lon)
            if lon.max() - lon.min() > 180:
                raise ValueError("cell region crosses the antimeridian; split it")
            # HEALPix edges are not straight in lon/lat: pad by a fraction of a cell.
            x0, y0, x1, y1 = lon.min(), lat.min(), lon.max(), lat.max()
        return (max(x0 - pad_deg, -180.0), max(y0 - pad_deg, -90.0),
                min(x1 + pad_deg, 180.0), min(y1 + pad_deg, 90.0))

    def query_wkt(self, pad_deg: float = 0.0) -> str:
        """Counter-clockwise WKT polygon for repository queries (GBIF requires CCW)."""
        if self.kind == "polygon" and pad_deg == 0.0:
            ring = list(self.polygon)
            if _signed_area(ring) < 0:
                ring = ring[::-1]
        else:
            x0, y0, x1, y1 = self.query_bounds(pad_deg)
            ring = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        ring = ring + [ring[0]]
        return "POLYGON((" + ",".join(f"{x:.6f} {y:.6f}" for x, y in ring) + "))"

    # ------------------------------------------------------------ exact test
    def contains(self, lon, lat) -> np.ndarray:
        lon = np.asarray(lon, dtype=np.float64)
        lat = np.asarray(lat, dtype=np.float64)
        if self.kind == "bbox":
            x0, y0, x1, y1 = self.bbox
            return (lon >= x0) & (lon <= x1) & (lat >= y0) & (lat <= y1)
        if self.kind == "polygon":
            return _points_in_polygon(lon, lat, np.asarray(self.polygon))
        from healpix_geo import nested

        ok = np.isfinite(lon) & np.isfinite(lat)
        out = np.zeros(lon.shape, dtype=bool)
        cells = nested.lonlat_to_healpix(lon[ok], lat[ok], np.uint8(self.depth), ellipsoid=ELLIPSOID)
        out[ok] = np.isin(cells, self.cells)
        return out


def _signed_area(ring) -> float:
    x = np.array([p[0] for p in ring])
    y = np.array([p[1] for p in ring])
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _points_in_polygon(x, y, poly) -> np.ndarray:
    """Even-odd ray casting; vectorised over points."""
    inside = np.zeros(x.shape, dtype=bool)
    px, py = poly[:, 0], poly[:, 1]
    j = len(poly) - 1
    for i in range(len(poly)):
        crosses = (py[i] > y) != (py[j] > y)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i]) + px[i]
        inside ^= crosses & (x < xint)
        j = i
    return inside
