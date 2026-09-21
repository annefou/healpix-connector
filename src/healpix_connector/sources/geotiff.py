"""Windowed reads of cloud-hosted GeoTIFFs (EPSG:4326), with scale/offset applied."""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Window:
    values: np.ndarray  # (nlat, nlon), physical units, NaN = missing
    lon: np.ndarray  # pixel-centre longitudes, increasing
    lat: np.ndarray  # pixel-centre latitudes, decreasing (north first)
    res_deg: float
    url: str


def read_window(url: str, bbox: tuple[float, float, float, float], pad_deg: float = 0.0) -> Window:
    """Read ``bbox`` = (lon_min, lat_min, lon_max, lat_max), padded by ``pad_deg``.

    Only the byte ranges covering the window are fetched (GDAL ``/vsicurl/``).
    The file's scale and offset are applied, and its nodata value (if any)
    becomes NaN.
    """
    import rasterio
    from rasterio.windows import from_bounds

    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
    lon_min, lat_min, lon_max, lat_max = bbox
    path = url if url.startswith("/vsi") or os.path.exists(url) else "/vsicurl/" + url
    with rasterio.open(path) as src:
        if src.crs is None or src.crs.to_epsg() != 4326:
            raise ValueError(f"expected EPSG:4326, got {src.crs}")
        win = from_bounds(lon_min - pad_deg, lat_min - pad_deg,
                          lon_max + pad_deg, lat_max + pad_deg, src.transform)
        win = win.round_offsets().round_lengths()
        raw = src.read(1, window=win)
        wt = src.window_transform(win)
        scale, offset, nodata = src.scales[0], src.offsets[0], src.nodata
    vals = raw.astype(np.float64)
    if nodata is not None:
        vals[raw == nodata] = np.nan
    vals = vals * scale + offset
    h, w = vals.shape
    lon = wt.c + (np.arange(w) + 0.5) * wt.a
    lat = wt.f + (np.arange(h) + 0.5) * wt.e
    return Window(values=vals, lon=lon, lat=lat, res_deg=abs(wt.a), url=url)
