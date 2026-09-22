"""GBIF connector on the test case (sdm-hotspot-spatial-effort).

1. Search: white stork (Ciconia ciconia) records of 2025 in a region given as
   HEALPix cells, its name resolved under COL XR and the taxonomy stated
   explicitly. The year filter keeps the query inside the search API's
   practical paging limit (see gbif.PRACTICAL_PAGING_LIMIT).
2. Download: the test case's existing museum download (doi:10.15468/dl.r8pcat),
   filtered to Iberia and put on depth-8 cells, with the positional-uncertainty
   and time facts that tests 2 and 3 depend on.

Run: pixi run -e test python examples/gbif_test_case_records.py
"""

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from healpix_connector.cells import assign_cells, footprint_cells
from healpix_connector.connectors import gbif
from healpix_connector.conventions import cell_size_m
from healpix_connector.region import Region

DEPTH = 8
IBERIA = Region.from_bbox(-10.0, 35.0, 4.0, 44.0)
MUSEUM_KEY = "0008222-260519110011954"
DATA = Path("data/gbif")
OUT = Path(__file__).parent / "results"


def part1_search():
    madrid = assign_cells([-3.70], [40.42], 8)  # one ~25 km cell in Madrid
    region = Region.from_cells(madrid, 8)
    match = gbif.match_name("Ciconia ciconia")
    res = gbif.search(region, taxon_key=match["taxon_key"], filters={"year": 2025})
    rows = gbif.to_cells(res.records, DEPTH)
    unc = [r["coordinateUncertaintyInMeters"] for r in rows]
    return {
        "region": f"HEALPix depth-8 cell {int(madrid[0])} (Madrid)",
        "taxon": match,
        "provenance": res.provenance,
        "records_in_region": len(rows),
        "with_uncertainty": int(sum(u is not None for u in unc)),
        "footprint_cells_depth8": {"median": float(np.median([len(r["footprint_cells"]) for r in rows])) if rows else None,
                                   "max": max((len(r["footprint_cells"]) for r in rows), default=None)},
    }


def part2_download():
    meta = gbif.download_metadata(MUSEUM_KEY)
    DATA.mkdir(parents=True, exist_ok=True)
    zpath = DATA / f"{MUSEUM_KEY}.zip"
    if not zpath.exists():
        with requests.get(meta["download_link"], stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(zpath, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
    cols = ["gbifID", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
            "year", "basisOfRecord", "taxonKey", "species"]
    with zipfile.ZipFile(zpath) as z:
        name = next(n for n in z.namelist() if n.endswith(".csv"))
        with z.open(name) as f:
            df = pd.read_csv(f, sep="\t", usecols=cols, dtype={"gbifID": str, "taxonKey": str},
                             quoting=3, on_bad_lines="skip", low_memory=False)
    n_all = len(df)
    df = df.dropna(subset=["decimalLatitude", "decimalLongitude"])
    df = df[IBERIA.contains(df.decimalLongitude.values, df.decimalLatitude.values)].copy()
    df["cell"] = assign_cells(df.decimalLongitude.values, df.decimalLatitude.values, DEPTH)

    unc = df.coordinateUncertaintyInMeters
    half_cell = cell_size_m(DEPTH) / 2
    known = unc.notna() & (unc > 0)
    # Exact footprints on a reproducible sample of records with stated uncertainty.
    sample = df[known].sample(n=min(5000, int(known.sum())), random_state=20260921)
    fp = [footprint_cells(x, y, u, DEPTH)[0].size for x, y, u in
          zip(sample.decimalLongitude, sample.decimalLatitude, sample.coordinateUncertaintyInMeters)]
    yr = df.year
    return {
        "download": {k: meta[k] for k in ("key", "doi", "created", "total_records", "format", "checklist_key")},
        "rows_read": int(n_all),
        "records_in_iberia": int(len(df)),
        "cells_occupied_depth8": int(df.cell.nunique()),
        "species": int(df.species.nunique()),
        "uncertainty": {
            "share_stated": round(float(known.mean()), 4),
            "median_m_when_stated": float(unc[known].median()),
            "p90_m_when_stated": float(unc[known].quantile(0.9)),
            "share_stated_wider_than_half_a_depth8_cell": round(float((unc[known] > half_cell).mean()), 4),
            "sample_share_footprint_spans_more_than_one_cell": round(float(np.mean(np.array(fp) > 1)), 4),
            "sample_size": len(fp),
        },
        "time": {
            "share_year_missing": round(float(yr.isna().mean()), 4),
            "share_in_1981_2010": round(float(yr.between(1981, 2010).mean()), 4),
            "share_before_1981": round(float((yr < 1981).mean()), 4),
            "share_after_2010": round(float((yr > 2010).mean()), 4),
            "median_year": float(yr.median()),
        },
    }


def main():
    report = {"part1_search": part1_search()}
    print(json.dumps(report["part1_search"], indent=1, default=str)[:1500])
    report["part2_download"] = part2_download()
    print(json.dumps(report["part2_download"], indent=1, default=str))
    OUT.mkdir(exist_ok=True)
    (OUT / "gbif_test_case_records.json").write_text(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
