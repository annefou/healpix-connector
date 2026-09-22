"""Step 4 on the test case: CHELSA bio1 attached to GBIF records.

Every value carries its support (depth 8), its time status against CHELSA's
1981-2010 climatology, and two uncertainty components (the source's own spread
inside the cell, and the spread across the cells a record's positional
uncertainty covers).

It also answers the record-weighted version of test 1: what share of *records*
sit in cells where the cell-centre pixel (what sdm-hotspot-spatial-effort uses)
differs from the cell mean.

Run: pixi run -e test python examples/sample_records_demo.py
"""

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from healpix_connector import sample as S
from healpix_connector.cells import assign_cells
from healpix_connector.connectors import gbif
from healpix_connector.region import Region
from healpix_connector.sources import chelsa
from healpix_connector.sources.geotiff import read_window

DEPTH = 8
IBERIA = Region.from_bbox(-10.0, 35.0, 4.0, 44.0)
ZARR = Path("data/converted/CHELSA_v21_bio1_on_HEALPix_depth_8.zarr")
ZIP = Path("data/gbif/0008222-260519110011954.zip")
OUT = Path(__file__).parent / "results"


def load_records() -> pd.DataFrame:
    cols = ["gbifID", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
            "eventDate", "year", "species"]
    with zipfile.ZipFile(ZIP) as z:
        name = next(n for n in z.namelist() if n.endswith(".csv"))
        with z.open(name) as f:
            df = pd.read_csv(f, sep="\t", usecols=cols, dtype={"gbifID": str},
                             quoting=3, on_bad_lines="skip", low_memory=False)
    df = df.dropna(subset=["decimalLatitude", "decimalLongitude"])
    df = df[IBERIA.contains(df.decimalLongitude.values, df.decimalLatitude.values)].copy()
    df["cell"] = assign_cells(df.decimalLongitude.values, df.decimalLatitude.values, DEPTH)
    return df


def main():
    src = S.open_source(ZARR)
    df = load_records()

    got = S.sample_cells(src, df.cell.values, DEPTH)
    df["value"], df["source_std"] = got["value"], got["source_std"]
    df["time_status"] = [S.time_status(d, src) for d in df.year.values]

    # Positional component on a reproducible sample (cone coverage per record).
    with_unc = df[df.coordinateUncertaintyInMeters.notna() & (df.coordinateUncertaintyInMeters > 0)]
    sample = with_unc.sample(n=min(5000, len(with_unc)), random_state=20260922)
    recs = gbif.to_cells([{"decimalLongitude": x, "decimalLatitude": y,
                           "coordinateUncertaintyInMeters": u, "eventDate": e}
                          for x, y, u, e in zip(sample.decimalLongitude, sample.decimalLatitude,
                                                sample.coordinateUncertaintyInMeters, sample.eventDate)], DEPTH)
    sampled = S.sample_records(recs, src, DEPTH)
    spread = np.array([r["positional_spread"] for r in sampled])

    # Record-weighted test 1: cell-centre pixel vs cell mean, per record.
    win = read_window(chelsa.url(1), IBERIA.bbox, pad_deg=1.0)
    from healpix_geo import nested
    cells = np.unique(df.cell.values)
    clon, clat = nested.healpix_to_lonlat(cells, np.uint8(DEPTH), ellipsoid="WGS84")
    clon = np.where(clon > 180, clon - 360, clon)
    i = np.abs(win.lat[:, None] - clat[None, :]).argmin(axis=0)
    j = np.abs(win.lon[:, None] - clon[None, :]).argmin(axis=0)
    centre = pd.Series(win.values[i, j], index=cells)
    mean = pd.Series(S.sample_cells(src, cells, DEPTH)["value"], index=cells)
    diff = (centre - mean).abs()
    df["centre_minus_mean"] = df.cell.map(diff)

    report = {
        "provenance": S.provenance(src, DEPTH),
        "records": {
            "in_iberia": int(len(df)),
            "with_a_value": int(df.value.notna().sum()),
            "cells_occupied": int(df.cell.nunique()),
        },
        "time": {
            "outside_the_climatology": round(float(df.time_status.str.startswith("outside").mean()), 4),
            "inside_the_climatology": round(float(df.time_status.str.startswith("inside").mean()), 4),
            "date_unknown": round(float(df.time_status.str.contains("unknown").mean()), 4),
        },
        "uncertainty_source_within_cell_degC": {
            "median": round(float(df.source_std.median()), 3),
            "p95": round(float(df.source_std.quantile(0.95)), 3),
            "max": round(float(df.source_std.max()), 3),
        },
        "uncertainty_positional_degC_sample": {
            "n": len(spread),
            "share_nonzero": round(float((spread > 0).mean()), 4),
            "median_when_nonzero": round(float(np.median(spread[spread > 0])) if (spread > 0).any() else 0.0, 3),
            "max": round(float(spread.max()), 3),
        },
        "record_weighted_test1_centre_vs_mean_degC": {
            "median": round(float(df.centre_minus_mean.median()), 3),
            "p95": round(float(df.centre_minus_mean.quantile(0.95)), 3),
            "share_records_gt_0.5C": round(float((df.centre_minus_mean > 0.5).mean()), 4),
            "share_records_gt_1C": round(float((df.centre_minus_mean > 1.0).mean()), 4),
            "share_cells_gt_1C": round(float((diff > 1.0).mean()), 4),
        },
    }
    print(json.dumps(report, indent=1))
    OUT.mkdir(exist_ok=True)
    (OUT / "sample_records_chelsa_bio1.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
