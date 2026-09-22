"""CHELSA v2.1 bioclimatic climatologies (1981-2010, ~1 km) on HEALPix.

Source: Karger et al., CHELSA v2.1, EnviDat, doi:10.16904/envidat.228, CC0-1.0.
Only variables whose physical meaning is verified are listed in ``BIO``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from healpix_connector import __version__
from healpix_connector.binning import CellStats, bin_to_cells
from healpix_connector.dggs_zarr import cf_grid_mapping_attrs, dggs_attrs
from healpix_connector.sources.geotiff import read_window

SOURCE = {
    "name": "CHELSA v2.1 climatologies",
    "doi": "10.16904/envidat.228",
    "license": "CC0-1.0",
    "version": "2.1",
    "reference_period": ["1981-01-01", "2010-12-31"],
    "url_template": (
        "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/"
        "CHELSA_bio{n}_1981-2010_V.2.1.tif"
    ),
}

# bio number -> (CF standard_name, units, long name).
# Units are those documented in CHELSA's own file specification (all six are
# stored as uint16 with scale 0.1; temperatures also carry offset -273.15,
# which read_window applies).
BIO = {
    1: ("air_temperature", "degC", "Mean annual air temperature (BIO1)"),
    4: ("air_temperature", "degC", "Temperature seasonality: standard deviation of monthly mean temperatures (BIO4)"),
    5: ("air_temperature", "degC", "Mean daily maximum air temperature of the warmest month (BIO5)"),
    6: ("air_temperature", "degC", "Mean daily minimum air temperature of the coldest month (BIO6)"),
    12: ("precipitation_amount", "kg m-2", "Annual precipitation amount (BIO12)"),
    # CHELSA documents bio15 as kg m-2, though it is a coefficient of variation.
    15: ("precipitation_amount", "kg m-2", "Precipitation seasonality: coefficient of variation of monthly precipitation (BIO15)"),
}

RESAMPLING = "area-weighted binning of source pixel centres into WGS84 HEALPix cells"


def url(n: int) -> str:
    return SOURCE["url_template"].format(n=n)


def to_dataset(stats: CellStats, n: int):
    """Build a DGGS-Zarr-conformant xarray Dataset from per-cell statistics."""
    import xarray as xr

    if n not in BIO:
        raise KeyError(f"bio{n} metadata not verified yet; known: {sorted(BIO)}")
    std_name, units, long_name = BIO[n]
    var = f"bio{n}"

    def field(a, what):
        return ("cells", a, {"units": units, "long_name": f"{long_name}: {what} over the cell",
                             "grid_mapping": "crs"})

    ds = xr.Dataset(
        {
            var: ("cells", stats.mean, {"standard_name": std_name, "units": units,
                                        "long_name": long_name, "cell_methods": "cells: mean (area-weighted)",
                                        "grid_mapping": "crs"}),
            f"{var}_std": field(stats.std, "area-weighted standard deviation"),
            f"{var}_min": field(stats.minimum, "minimum"),
            f"{var}_max": field(stats.maximum, "maximum"),
            "pixel_count": ("cells", stats.pixel_count, {"long_name": "valid source pixels in the cell", "units": "1"}),
            "coverage": ("cells", stats.coverage, {"long_name": "valid source area / cell area", "units": "1"}),
            "crs": ((), np.int8(0), cf_grid_mapping_attrs(stats.depth)),
        },
        coords={"cell_ids": ("cells", stats.cell_ids, {"standard_name": "healpix_index", "units": "1"})},
    )
    ds.attrs.update(dggs_attrs(stats.depth))
    # Layer 3: what the numbers are.
    ds.attrs.update({
        "title": f"CHELSA v2.1 {var} on HEALPix depth {stats.depth}",
        "source": SOURCE["name"],
        "source_doi": SOURCE["doi"],
        "source_url": url(n),
        "license": SOURCE["license"],
        "support_depth": int(stats.depth),
        "valid_time": "climatology",
        "climatology_bounds": SOURCE["reference_period"],
        "resampling_method": RESAMPLING,
        "uncertainty_kind": "within-cell spread (std, min, max) of the ~1 km source field",
        "producer": f"healpix-connector {__version__}",
    })
    return ds


def convert(n: int, bbox: tuple[float, float, float, float], depth: int, pad_deg: float = 1.0):
    """Read bio``n`` over ``bbox`` (+pad) and bin it to HEALPix ``depth``."""
    win = read_window(url(n), bbox, pad_deg=pad_deg)
    stats = bin_to_cells(win.values, win.lon, win.lat, depth, win.res_deg)
    return to_dataset(stats, n), win


def stac_item(ds, zarr_href: str, bbox: tuple[float, float, float, float]) -> dict:
    """A minimal STAC 1.0 item describing a converted dataset."""
    lon0, lat0, lon1, lat1 = bbox
    a = ds.attrs
    return {
        "type": "Feature",
        "stac_version": "1.0.0",
        "id": a["title"].replace(" ", "-").lower(),
        "bbox": [lon0, lat0, lon1, lat1],
        "geometry": {"type": "Polygon", "coordinates": [[[lon0, lat0], [lon1, lat0], [lon1, lat1], [lon0, lat1], [lon0, lat0]]]},
        "properties": {
            "datetime": None,
            "start_datetime": a["climatology_bounds"][0] + "T00:00:00Z",
            "end_datetime": a["climatology_bounds"][1] + "T23:59:59Z",
            "title": a["title"],
            "license": a["license"],
            "sci:doi": a["source_doi"],
            "dggs": a["dggs"],
            "support_depth": a["support_depth"],
            "valid_time": a["valid_time"],
            "resampling_method": a["resampling_method"],
            "uncertainty_kind": a["uncertainty_kind"],
        },
        "links": [{"rel": "cite-as", "href": f"https://doi.org/{a['source_doi']}"}],
        "assets": {"data": {"href": zarr_href, "type": "application/vnd+zarr", "roles": ["data"]}},
    }


def write(ds, out_dir: Path, bbox) -> Path:
    """Write ``ds`` as Zarr v3 plus a STAC item next to it; return the Zarr path."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = ds.attrs["title"].replace(" ", "_").replace(".", "")
    zpath = out_dir / f"{name}.zarr"
    ds.to_zarr(zpath, mode="w", zarr_format=3, consolidated=False)
    (out_dir / f"{name}.stac.json").write_text(json.dumps(stac_item(ds, zpath.name, bbox), indent=2))
    return zpath
