"""Preview of test 1: cell-centre nearest pixel vs area-weighted cell mean.

Reproduces how sdm-hotspot-spatial-effort samples CHELSA (nearest source pixel
at each HEALPix cell centre, cells = centre inside the Iberia bbox) and compares
it with the declared-support value from healpix-connector (area-weighted mean
over the whole cell), for CHELSA v2.1 bio1 at depths 6-9.

Run: pixi run -e test python examples/compare_centre_vs_cell_mean.py
"""

import json
from pathlib import Path

import numpy as np
from healpix_geo import nested

from healpix_connector.binning import bin_to_cells
from healpix_connector.conventions import ELLIPSOID, cell_size_m
from healpix_connector.sources import chelsa
from healpix_connector.sources.geotiff import read_window

BBOX = (-10.0, 35.0, 4.0, 44.0)  # as in sdm-hotspot-spatial-effort
DEPTHS = [6, 7, 8, 9]
OUT = Path(__file__).parent / "results"


def iberian_cells(depth):
    """Test-case definition: all NESTED cells whose centre lies in the bbox."""
    ids = np.arange(12 * 4**depth, dtype=np.uint64)
    lon, lat = nested.healpix_to_lonlat(ids, np.uint8(depth), ellipsoid=ELLIPSOID)
    lon = np.where(lon > 180, lon - 360, lon)
    m = (lon >= BBOX[0]) & (lon <= BBOX[2]) & (lat >= BBOX[1]) & (lat <= BBOX[3])
    return ids[m], lon[m], lat[m]


def nearest(win, lon, lat):
    i = np.abs(win.lat[:, None] - lat[None, :]).argmin(axis=0)
    j = np.abs(win.lon[:, None] - lon[None, :]).argmin(axis=0)
    return win.values[i, j]


def main():
    win = read_window(chelsa.url(1), BBOX, pad_deg=1.5)
    print(f"CHELSA bio1 window {win.values.shape}, {win.res_deg*3600:.0f} arcsec, "
          f"range {np.nanmin(win.values):.2f}..{np.nanmax(win.values):.2f} degC")
    report = {"source": chelsa.SOURCE, "bbox": BBOX, "variable": "bio1 (degC)", "depths": {}}
    for d in DEPTHS:
        stats = bin_to_cells(win.values, win.lon, win.lat, d, win.res_deg)
        ids, clon, clat = iberian_cells(d)
        pos = np.searchsorted(stats.cell_ids, ids)
        assert (stats.cell_ids[pos] == ids).all(), "target cell missing from window"
        cov = stats.coverage[pos]
        centre = nearest(win, clon, clat)
        mean, std = stats.mean[pos], stats.std[pos]
        diff = centre - mean
        a = np.abs(diff)
        r = {
            "cell_width_km": round(cell_size_m(d) / 1000, 2),
            "cells": int(ids.size),
            "min_coverage": round(float(cov.min()), 3),
            "abs_diff_mean_degC": round(float(a.mean()), 3),
            "abs_diff_p95_degC": round(float(np.percentile(a, 95)), 3),
            "abs_diff_max_degC": round(float(a.max()), 3),
            "share_abs_diff_gt_0.5C": round(float((a > 0.5).mean()), 4),
            "share_abs_diff_gt_1C": round(float((a > 1.0).mean()), 4),
            "within_cell_std_median_degC": round(float(np.median(std)), 3),
            "within_cell_std_p95_degC": round(float(np.percentile(std, 95)), 3),
            "within_cell_range_max_degC": round(float((stats.maximum[pos] - stats.minimum[pos]).max()), 2),
        }
        report["depths"][str(d)] = r
        print(f"depth {d} ({r['cell_width_km']} km, {r['cells']} cells, coverage>={r['min_coverage']}): "
              f"|centre-mean| mean {r['abs_diff_mean_degC']} p95 {r['abs_diff_p95_degC']} max {r['abs_diff_max_degC']} degC; "
              f">0.5C {100*r['share_abs_diff_gt_0.5C']:.1f}%, >1C {100*r['share_abs_diff_gt_1C']:.1f}%; "
              f"within-cell std median {r['within_cell_std_median_degC']} p95 {r['within_cell_std_p95_degC']}")
        if d == 8:
            z = chelsa.write(chelsa.to_dataset(stats, 1), Path("data/converted"), BBOX)
            print(f"  wrote {z}")
    OUT.mkdir(exist_ok=True)
    (OUT / "centre_vs_cell_mean_bio1_iberia.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
