"""Test 4: do the study's hotspots hold under a projected 2030s climate?

Climate DT (DestinE, IFS-NEMO, SSP3-7.0) is a free-running projection, so it
cannot say what a 2021 record experienced (that was test 3, with ERA5). It can
say how the climate itself changes, by comparing two decades of the same
simulation: 2020-2029 and 2030-2039, which is what Polytope serves.

Everything about the source - the MARS request, the endpoint, the HEALPix
conventions and the sphere -> WGS84 correction - comes from GRID4EARTH's
healpix-convert and healpix-resample, through healpix_connector.sources.climatedt. This script only
turns monthly means into bioclimatic variables and swaps them into the study.

bio1 and bio4 come from 2 m temperature; bio12 and bio15 from the total
precipitation rate sampled at 00 and 12 UTC, which estimates the monthly mean
rate rather than accumulating every hour.

Usage: python examples/climatedt_projection.py <repo> [nside]
"""

import calendar
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr

from healpix_connector.sources import climatedt

# The SSP3-7.0 run is served for 2020-2039 (2018 and 2040 are refused).
PERIODS = {"2020-2029": (2020, 2029), "2030-2039": (2030, 2039)}
VARIABLES = ["t2m", "tprate"]
OUT = Path(__file__).parent / "results"


def climatology(sampler: climatedt.MonthlySampler, period: str) -> dict:
    """bio1, bio4, bio12, bio15 over the period, per cell."""
    y0, y1 = PERIODS[period]
    t_by_month = {m: [] for m in range(1, 13)}
    p_by_month = {m: [] for m in range(1, 13)}
    for year in range(y0, y1 + 1):
        for month in range(1, 13):
            means = sampler.month(year, month)
            t_by_month[month].append(means.values["t2m"] - 273.15)
            p_by_month[month].append(means.values["tprate"])
        print(f"  {period}: {year} done", flush=True)

    t_m = np.stack([np.mean(t_by_month[m], axis=0) for m in range(1, 13)])
    days = np.array([calendar.monthrange(2001, m)[1] for m in range(1, 13)])[:, None]
    p_m = np.stack([np.mean(p_by_month[m], axis=0) for m in range(1, 13)]) * days * 86400.0
    return {"bio1": t_m.mean(axis=0), "bio4": t_m.std(axis=0),
            "bio12": p_m.sum(axis=0), "bio15": 100.0 * p_m.std(axis=0) / p_m.mean(axis=0)}


def main(repo: Path, nside: int = 256):
    depth = int(np.log2(nside))
    theirs = xr.open_dataset(repo / "data/clean/predictors_cells.nc", group=f"nside_{nside}")
    cells = theirs["cell"].values.astype(np.uint64)
    # Climate DT's cells are depth 7; the study's are depth 8, so each study
    # cell inherits its parent's value and the support is the ~51 km source cell.
    parents = cells >> np.uint64(2 * (depth - climatedt.DEPTH))
    source_cells = np.unique(parents).astype(np.int64)
    sampler = climatedt.MonthlySampler(source_cells, VARIABLES)
    print(f"{len(cells)} study cells at depth {depth} -> {len(source_cells)} "
          f"Climate DT cells at depth {climatedt.DEPTH}")

    report = {"nside": nside, "source_depth": climatedt.DEPTH, "source": climatedt.SOURCE,
              "sampling": f"inherited from depth {climatedt.DEPTH}", "periods": {}}
    fields = {}
    for period in PERIODS:
        coarse = climatology(sampler, period)
        take = np.searchsorted(source_cells, parents.astype(np.int64))
        fields[period] = {v: coarse[v][take] for v in coarse}  # depth 7 -> depth 8
        out = repo / f"data/clean/predictors_cells_cdt_{period}.nc"
        out.unlink(missing_ok=True)
        xr.Dataset(
            {v: (("cell",), fields[period][v]) for v in fields[period]},
            coords={"cell": ("cell", cells.astype(np.int64)),
                    "lon": ("cell", theirs["lon"].values), "lat": ("cell", theirs["lat"].values)},
            attrs={**theirs.attrs, "period": period, "support_depth": climatedt.DEPTH,
                   "source": f"Climate DT SSP3-7.0 {period} via healpix-connector "
                             f"(healpix-convert ClimateDTConverter, PSF sphere->WGS84)"},
        ).to_netcdf(out, mode="w", group=f"nside_{nside}", engine="netcdf4")
        report["periods"][period] = {v: round(float(np.nanmean(a)), 2) for v, a in coarse.items()}
        print(f"{period}: " + ", ".join(f"{v} {report['periods'][period][v]}" for v in coarse))

    a, b = fields["2020-2029"], fields["2030-2039"]
    report["ellipsoid_correction"] = {
        "method": f"healpix_resample {sampler.method} resampler, sphere -> WGS84",
        "why": ("nearest is healpix-convert's declared resampler for these fields; its "
                "optional PSF correction is a deconvolution and drove 114 of 480 cells "
                "negative on precipitation, which a rate cannot be"),
    }
    report["change_2030_2039_minus_2020_2029"] = {
        v: {"mean": round(float(np.nanmean(b[v] - a[v])), 2),
            "max_abs": round(float(np.nanmax(np.abs(b[v] - a[v]))), 2)} for v in a}
    print("\nprojected change:", json.dumps(report["change_2030_2039_minus_2020_2029"]))
    OUT.mkdir(exist_ok=True)
    (OUT / "test4_climatedt_projection.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), int(sys.argv[2]) if len(sys.argv) > 2 else 256)
