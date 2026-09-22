"""Test 1: rebuild a published study's covariates as declared-support cell means.

sdm-hotspot-spatial-effort samples CHELSA at each HEALPix cell's centre (one
~1 km pixel per cell, ~650 km2 at depth 8). This script writes the same file in
the same format, with the same cells and the same source data, but with each
value replaced by the area-weighted mean over the whole cell, as produced by
healpix-connector. Nothing else about the study changes, so any difference in
its verdict is attributable to the covariate sampling.

Usage:
    python examples/test1_swap_predictors.py <path-to-sdm-hotspot-spatial-effort> [nside]
"""

import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr

from healpix_connector.binning import bin_to_cells
from healpix_connector.sources import chelsa
from healpix_connector.sources.geotiff import read_window

IBERIA = (-10.0, 35.0, 4.0, 44.0)
PAD_DEG = 1.5
# CHELSA stores uint16 with scale 0.1; temperatures also carry offset -273.15.
# The study reads the raw integers, so convert its values for comparison.
OFFSET = {1: -273.15, 4: 0.0, 5: -273.15, 6: -273.15, 12: 0.0, 15: 0.0}


def main(repo: Path, nside: int = 256):
    depth = int(np.log2(nside))
    src = repo / "data/clean/predictors_cells.nc"
    out = repo / "data/clean/predictors_cells_cellmean.nc"
    theirs = xr.open_dataset(src, group=f"nside_{nside}")
    cells = theirs["cell"].values.astype(np.uint64)
    bio_vars = [v for v in theirs.data_vars]
    print(f"{src.name}: {len(cells)} cells at nside {nside} (depth {depth}), vars {bio_vars}")

    data, report = {}, {}
    for v in bio_vars:
        n = int(v.replace("bio", ""))
        win = read_window(chelsa.url(n), IBERIA, pad_deg=PAD_DEG)
        stats = bin_to_cells(win.values, win.lon, win.lat, depth, win.res_deg)
        pos = np.searchsorted(stats.cell_ids, cells)
        pos = np.clip(pos, 0, stats.cell_ids.size - 1)
        hit = stats.cell_ids[pos] == cells
        mean = np.where(hit, stats.mean[pos], np.nan)
        spread = np.where(hit, stats.std[pos], np.nan)
        data[v] = mean
        centre = theirs[v].values * 0.1 + OFFSET[n]  # their raw integers -> physical units
        d = np.abs(centre - mean)
        report[v] = {
            "cells": int(len(cells)),
            "unit": chelsa.BIO[n][1],
            "median_abs_diff": round(float(np.nanmedian(d)), 3),
            "p95_abs_diff": round(float(np.nanpercentile(d, 95)), 3),
            "max_abs_diff": round(float(np.nanmax(d)), 3),
            "median_within_cell_std": round(float(np.nanmedian(spread)), 3),
        }
        print(f"  {v}: centre vs cell mean, median |diff| {report[v]['median_abs_diff']} "
              f"{report[v]['unit']}, p95 {report[v]['p95_abs_diff']}, max {report[v]['max_abs_diff']}")

    # Same structure as 02_data_clean.py writes, so 03_analysis.py reads it unchanged.
    # Values are physical units here; the study standardises its covariates, so a
    # linear rescaling of the inputs does not change the model.
    ds = xr.Dataset(
        {v: (("cell",), data[v]) for v in bio_vars},
        coords={"cell": ("cell", cells.astype(np.int64)),
                "lon": ("cell", theirs["lon"].values), "lat": ("cell", theirs["lat"].values)},
        attrs={**theirs.attrs, "source": "CHELSA v2.1 via healpix-connector (area-weighted cell means)",
               "sampling": "area-weighted mean over the whole cell", "support_depth": depth},
    )
    if out.exists():
        out.unlink()
    ds.to_netcdf(out, mode="w", group=f"nside_{nside}", engine="netcdf4",
                 encoding={v: {"zlib": True, "complevel": 4} for v in bio_vars})
    print(f"wrote {out}")
    Path(__file__).parent.joinpath("results/test1_covariate_difference.json").write_text(
        json.dumps({"nside": nside, "depth": depth, "per_variable": report}, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), int(sys.argv[2]) if len(sys.argv) > 2 else 256)
