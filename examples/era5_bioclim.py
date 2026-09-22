"""Test 3: are the covariates from the right period?

98.5% of the study's records postdate CHELSA's 1981-2010 climatology. This
builds the same six bioclim variables from ERA5 monthly means, for two periods:

- 1981-2010, the same period as CHELSA: isolates the dataset and support
  difference (ERA5 ~28 km reanalysis vs CHELSA ~1 km statistical downscaling);
- 2011-2025, the period the records actually come from: the time effect.

ERA5 is used rather than Climate DT because Climate DT is a free-running
projection: its 2021 is a plausible 2021 under SSP3-7.0, not the observed one.

Definitions follow CHELSA's file specification: bio1 annual mean temperature,
bio4 standard deviation of monthly means, bio5/bio6 warmest/coldest month mean
daily max/min, bio12 annual precipitation, bio15 coefficient of variation of
monthly precipitation (as a percentage).

Usage: python examples/era5_bioclim.py <repo> [nside]
"""

import calendar
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr

from healpix_connector.sources.gridded import sample_cells_from_grid

AREA = [45, -11, 34, 5]  # N, W, S, E: Iberia plus a margin
PERIODS = {"1981-2010": (1981, 2010), "2011-2025": (2011, 2025)}
# ERA5 monthly means offer no daily max/min temperature, so bio5 and bio6
# cannot be built faithfully from them (they would need the derived daily
# statistics, one year per request). Test 3 therefore uses the four variables
# that monthly means do support, and the CHELSA baseline is rerun with the same
# four, so only the dataset and the period differ.
VARS = {"2m_temperature": "t2m", "total_precipitation": "tp"}
VAR4 = ["bio1", "bio4", "bio12", "bio15"]
DATA = Path("data/era5")
OUT = Path(__file__).parent / "results"


def fetch(period: str) -> Path:
    y0, y1 = PERIODS[period]
    path = DATA / f"era5_monthly_{period}.nc"
    if path.exists():
        return path
    import cdsapi

    DATA.mkdir(parents=True, exist_ok=True)
    cdsapi.Client().retrieve("reanalysis-era5-single-levels-monthly-means", {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": list(VARS),
        "year": [str(y) for y in range(y0, y1 + 1)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": ["00:00"], "area": AREA,
        "data_format": "netcdf", "download_format": "unarchived",
    }, path)
    return path


def _open(path: Path) -> xr.Dataset:
    """Open a CDS result: one NetCDF, or a zip holding one file per stream."""
    import zipfile

    if not zipfile.is_zipfile(path):
        return xr.open_dataset(path)
    out = path.with_suffix("")
    out.mkdir(exist_ok=True)
    with zipfile.ZipFile(path) as z:
        z.extractall(out)
    parts = [xr.open_dataset(f) for f in sorted(out.glob("*.nc"))]
    return xr.merge(parts, compat="override", join="override")


def bioclim(path: Path) -> tuple[dict, np.ndarray, np.ndarray]:
    """The six variables from monthly means, in CHELSA's units."""
    ds = _open(path)
    time = "valid_time" if "valid_time" in ds.dims else "time"
    monthly = ds.groupby(f"{time}.month").mean(time)  # climatological month means
    days = np.array([calendar.monthrange(2001, m)[1] for m in monthly.month.values])

    t = monthly["t2m"] - 273.15
    # tp is the mean daily accumulation in metres: to mm per month.
    p = monthly["tp"] * xr.DataArray(days, dims="month", coords={"month": monthly.month}) * 1000.0

    out = {
        "bio1": t.mean("month"),
        "bio4": t.std("month"),
        "bio12": p.sum("month"),
        "bio15": 100.0 * p.std("month") / p.mean("month"),
    }
    lat = ds["latitude"].values
    lon = ds["longitude"].values
    return {k: v.values for k, v in out.items()}, lon, lat


def main(repo: Path, nside: int = 256):
    depth = int(np.log2(nside))
    theirs = xr.open_dataset(repo / "data/clean/predictors_cells.nc", group=f"nside_{nside}")
    cells = theirs["cell"].values.astype(np.uint64)
    chelsa_offset = {"bio1": -273.15, "bio4": 0.0, "bio12": 0.0, "bio15": 0.0}

    # The matched baseline: the study's own CHELSA values, same four variables.
    base = repo / "data/clean/predictors_cells_chelsa4.nc"
    if base.exists():
        base.unlink()
    theirs[VAR4].to_netcdf(base, mode="w", group=f"nside_{nside}", engine="netcdf4")
    print(f"wrote {base.name} (CHELSA, {VAR4})")
    report = {"nside": nside, "depth": depth, "periods": {}}

    for period in PERIODS:
        fields, lon, lat = bioclim(fetch(period))
        data, diffs = {}, {}
        for v, arr in fields.items():
            s = sample_cells_from_grid(arr, lon, lat, cells, depth, 0.25)
            data[v] = s.value
            chelsa = theirs[v].values * 0.1 + chelsa_offset[v]
            d = np.abs(chelsa - s.value)
            diffs[v] = {"median_abs_diff_vs_chelsa": round(float(np.nanmedian(d)), 2),
                        "p95": round(float(np.nanpercentile(d, 95)), 2),
                        "mean_era5": round(float(np.nanmean(s.value)), 2),
                        "mean_chelsa": round(float(np.nanmean(chelsa)), 2)}
            support = s
        report["periods"][period] = {
            "support_m": round(support.support_m), "cell_width_m": round(support.cell_width_m),
            "method": support.method, "per_variable": diffs}
        print(f"{period}: {support.method}")
        for v, d in diffs.items():
            print(f"   {v}: ERA5 mean {d['mean_era5']}, CHELSA mean {d['mean_chelsa']}, "
                  f"median |diff| {d['median_abs_diff_vs_chelsa']}")

        out = repo / f"data/clean/predictors_cells_era5_{period}.nc"
        if out.exists():
            out.unlink()
        xr.Dataset({v: (("cell",), data[v]) for v in fields},
                   coords={"cell": ("cell", cells.astype(np.int64)),
                           "lon": ("cell", theirs["lon"].values), "lat": ("cell", theirs["lat"].values)},
                   attrs={**theirs.attrs, "source": f"ERA5 monthly means {period} via healpix-connector",
                          "support_m": float(support.support_m), "sampling": support.method,
                          "period": period}).to_netcdf(
            out, mode="w", group=f"nside_{nside}", engine="netcdf4",
            encoding={v: {"zlib": True, "complevel": 4} for v in fields})
        print(f"   wrote {out.name}")

    # How much do the two ERA5 periods differ from each other? That is the time effect.
    a = xr.open_dataset(repo / "data/clean/predictors_cells_era5_1981-2010.nc", group=f"nside_{nside}")
    b = xr.open_dataset(repo / "data/clean/predictors_cells_era5_2011-2025.nc", group=f"nside_{nside}")
    report["time_effect_2011_2025_minus_1981_2010"] = {
        v: {"mean_change": round(float(np.nanmean(b[v].values - a[v].values)), 2),
            "max_abs_change": round(float(np.nanmax(np.abs(b[v].values - a[v].values))), 2)}
        for v in fields}
    print("\ntime effect (2011-2025 minus 1981-2010):")
    for v, d in report["time_effect_2011_2025_minus_1981_2010"].items():
        print(f"   {v}: mean change {d['mean_change']}, max |change| {d['max_abs_change']}")
    OUT.mkdir(exist_ok=True)
    (OUT / "test3_era5_periods.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), int(sys.argv[2]) if len(sys.argv) > 2 else 256)
