"""Conservative bilinear resampling onto a fixed set of WGS84 cells.

healpix-resample's ``BilinearResampler`` in ``conservative=True`` mode spreads each
source sample over its four nearest target cells with weights that sum to one per
sample, so a field's integral is preserved exactly. For an *intensive* field
(temperature, a rate, a reflectance) the conservative result is a sum of
``value x area``; the cell mean is that sum divided by the same operator applied to
the valid samples' areas. Missing samples (NaN: land, cloud, a gap) are excluded
from both, so the mean is over the valid samples only, and the valid share is
reported per cell.

Why this one: when source and target cells are about the same size (Climate DT at
depth 7 or 10 onto WGS84 cells of the same depth, a 0.05 deg grid onto depth 10),
binning by sample centre (``ConservativeResampler``, ``GroupByResampler``) leaves
some cells with two samples and others with none, and ``nearest`` copies one
sample. PSF is a deconvolution and rings below zero on non-negative fields. Tested
against an exact polygon overlap (2026-10-23, Anatolia, depth 10): mean absolute
difference 0.07 degC for Climate DT 2 m temperature (sphere -> WGS84), 0.09 degC and
0.05 mm for WorldClim bio9 and bio14, regional means preserved, every cell filled;
taking the Climate DT value of the cell with the same id instead (no correction)
gives 0.80 degC on average and 6.7 degC at the 99th percentile.
"""

from __future__ import annotations

import numpy as np

from healpix_connector.conventions import ELLIPSOID


class BilinearConservative:
    """Operator from fixed source samples (geodetic lon/lat) onto fixed WGS84 cells.

    ``area`` is each sample's footprint (1 for equal-area cells such as a HEALPix
    source; ``cos(lat)`` for a regular lon/lat grid). Build once, apply to many
    fields: ``mean, valid_share = op(values)``; ``values`` may be ``(N,)`` or
    ``(B, N)``.
    """

    def __init__(self, lon, lat, cell_ids, depth: int, *, area=None):
        from healpix_resample import BilinearResampler

        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        self.cell_ids = np.asarray(cell_ids, dtype=np.uint64)
        if np.any(np.diff(self.cell_ids.astype(np.int64)) <= 0):
            raise ValueError("cell_ids must be sorted and unique")
        self.depth = int(depth)
        self.area = np.ones(lon.size) if area is None else np.asarray(area, dtype=float)
        self._op = BilinearResampler(lon_deg=lon, lat_deg=lat, level=self.depth, nest=True,
                                     ellipsoid=ELLIPSOID, verbose=False, area=self.area)
        self._full = self._place(self._op.resample(np.ones(lon.size), conservative=True))

    def _place(self, res) -> np.ndarray:
        cid = np.asarray(res.cell_ids, dtype=np.uint64)
        data = np.asarray(res.cell_data, dtype=float)
        out = np.full(data.shape[:-1] + (self.cell_ids.size,), np.nan)
        pos = np.clip(np.searchsorted(self.cell_ids, cid), 0, self.cell_ids.size - 1)
        hit = self.cell_ids[pos] == cid
        out[..., pos[hit]] = data[..., hit]
        return out

    def __call__(self, values) -> tuple[np.ndarray, np.ndarray]:
        v = np.asarray(values, dtype=float)
        ok = np.isfinite(v)
        num = self._place(self._op.resample(np.where(ok, v, 0.0), conservative=True))
        den = self._place(self._op.resample(ok.astype(float), conservative=True))
        with np.errstate(invalid="ignore", divide="ignore"):
            return num / den, den / self._full
