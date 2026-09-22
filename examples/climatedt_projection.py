"""Test 4: do the hotspots hold under 2030s-2040s climate?

Climate DT (DestinE, IFS-NEMO, SSP3-7.0) is a free-running projection, so it
cannot say what a 2021 record experienced (that was test 3, with ERA5). It can
say how the climate itself changes, by comparing two decades of the *same*
simulation: 2020-2029 and 2030-2039 (the years Polytope serves).

Climate DT is delivered on HEALPix nside 128 (depth 7, ~50 km) in NESTED order,
on a sphere. The study's cells are WGS84 at depth 8, so each study cell's centre
is looked up in the *spherical* depth-7 grid rather than assuming ids match; the
support is declared as the depth-7 cell (~51 km), wider than the study's cell.

bio1 and bio4 come from 2 m temperature; bio12 and bio15 from the total
precipitation rate, sampled twice a day (00 and 12 UTC), which estimates the
monthly mean rate rather than accumulating every hour.

Usage: python examples/climatedt_projection.py <repo> [nside]
"""

import calendar
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr

# The SSP3-7.0 run is served for 2020-2039 (2018 and 2040 are refused by Polytope).
PERIODS = {"2020-2029": (2020, 2029), "2030-2039": (2030, 2039)}
PARAMS = "167/260048"  # 2 m temperature, total precipitation rate
TIMES = "0000/1200"
SRC_DEPTH = 7
ADDRESS = "polytope.lumi.apps.dte.destination-earth.eu"
DATA = Path("data/climatedt")
OUT = Path(__file__).parent / "results"


def month_file(year: int, month: int) -> Path:
    import polytope.api as polytope

    DATA.mkdir(parents=True, exist_ok=True)
    path = DATA / f"cdt_{year}{month:02d}.grib"
    if path.exists():
        return path
    last = calendar.monthrange(year, month)[1]
    req = {"activity": "ScenarioMIP", "class": "d1", "dataset": "climate-dt",
           "date": f"{year}{month:02d}01/to/{year}{month:02d}{last:02d}",
           "experiment": "SSP3-7.0", "expver": "0001", "generation": "1",
           "levtype": "sfc", "model": "IFS-NEMO", "param": PARAMS, "realization": "1",
           "resolution": "standard", "stream": "clte", "time": TIMES, "type": "fc"}
    try:
        polytope.Client(address=ADDRESS, quiet=True).retrieve("destination-earth", req, str(path))
    except SystemExit as exc:  # the polytope client exits the process on a refused request
        raise RuntimeError(f"Climate DT refused {year}-{month:02d}") from exc
    return path


def monthly_means(year: int, month: int, idx: np.ndarray) -> tuple[float, float]:
    """(mean 2 m temperature in degC, mean precipitation rate in kg m-2 s-1) per cell."""
    import cfgrib

    path = month_file(year, month)
    t = p = None
    for ds in cfgrib.open_datasets(str(path)):
        if "t2m" in ds:
            t = ds["t2m"].mean("time").values[idx] - 273.15
        if "tprate" in ds:
            p = ds["tprate"].mean("time").values[idx]
    path.unlink()  # processed; keep the download footprint small
    return t, p


def climatology(period: str, idx: np.ndarray) -> dict:
    y0, y1 = PERIODS[period]
    t_by_month, p_by_month = {m: [] for m in range(1, 13)}, {m: [] for m in range(1, 13)}
    for year in range(y0, y1 + 1):
        for month in range(1, 13):
            t, p = monthly_means(year, month, idx)
            t_by_month[month].append(t)
            p_by_month[month].append(p)
        print(f"  {period}: {year} done", flush=True)
    t_m = np.stack([np.mean(t_by_month[m], axis=0) for m in range(1, 13)])
    # monthly precipitation totals (mm) from the mean rate
    days = np.array([calendar.monthrange(2001, m)[1] for m in range(1, 13)])[:, None]
    p_m = np.stack([np.mean(p_by_month[m], axis=0) for m in range(1, 13)]) * days * 86400.0
    return {"bio1": t_m.mean(axis=0), "bio4": t_m.std(axis=0),
            "bio12": p_m.sum(axis=0), "bio15": 100.0 * p_m.std(axis=0) / p_m.mean(axis=0)}


def main(repo: Path, nside: int = 256):
    from healpix_geo import nested

    depth = int(np.log2(nside))
    theirs = xr.open_dataset(repo / "data/clean/predictors_cells.nc", group=f"nside_{nside}")
    cells = theirs["cell"].values.astype(np.uint64)
    # Each WGS84 study cell's centre, looked up in Climate DT's spherical depth-7 grid.
    clon, clat = nested.healpix_to_lonlat(cells, np.uint8(depth), ellipsoid="WGS84")
    idx = nested.lonlat_to_healpix(np.where(clon > 180, clon - 360, clon), clat,
                                   np.uint8(SRC_DEPTH), ellipsoid="sphere").astype(np.int64)
    print(f"{len(cells)} study cells -> {len(np.unique(idx))} Climate DT cells at depth {SRC_DEPTH}")

    report, fields = {"nside": nside, "source_depth": SRC_DEPTH,
                      "support": "Climate DT depth-7 cell (~51 km), wider than the study's ~25 km cell",
                      "periods": {}}, {}
    for period in PERIODS:
        fields[period] = climatology(period, idx)
        out = repo / f"data/clean/predictors_cells_cdt_{period}.nc"
        if out.exists():
            out.unlink()
        xr.Dataset({v: (("cell",), fields[period][v]) for v in fields[period]},
                   coords={"cell": ("cell", cells.astype(np.int64)),
                           "lon": ("cell", theirs["lon"].values), "lat": ("cell", theirs["lat"].values)},
                   attrs={**theirs.attrs, "source": f"Climate DT SSP3-7.0 {period} via healpix-connector",
                          "support_depth": SRC_DEPTH, "period": period}).to_netcdf(
            out, mode="w", group=f"nside_{nside}", engine="netcdf4")
        report["periods"][period] = {v: round(float(np.nanmean(a)), 2) for v, a in fields[period].items()}
        print(f"{period}: " + ", ".join(f"{v} {report['periods'][period][v]}" for v in fields[period]))

    a, b = fields["2020-2029"], fields["2030-2039"]
    report["change_2030_2039_minus_2020_2029"] = {
        v: {"mean": round(float(np.nanmean(b[v] - a[v])), 2),
            "max_abs": round(float(np.nanmax(np.abs(b[v] - a[v]))), 2)} for v in a}
    print("\nprojected change:", json.dumps(report["change_2030_2039_minus_2020_2029"]))
    OUT.mkdir(exist_ok=True)
    (OUT / "test4_climatedt_projection.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), int(sys.argv[2]) if len(sys.argv) > 2 else 256)
