# Upstream issues (drafts, not filed)

Everything here is measured from healpix-connector's own use. Nothing is filed
until Anne has read it. Checked against `main` on 2026-09-22 (healpix-resample
`2026.7.0` installed, `main` newer; healpix-convert `2026.9.0`).

**Landscape, checked 2026-09-22.** healpix-resample: 5 open PRs (4 dependabot, #57
notebooks env), 5 open issues, none covering these; #27 is about resampler parameter
naming and #5 (user-defined output cell ids) is closed. healpix-convert: 1 open issue
(#9, licence links) and 7 open PRs, 6 of them ours from 2026-09-07, still unreviewed —
so file issues there before adding more code.

---

## 1. healpix-resample — `GroupByResampler`: `std`/`var` reductions, and cell counts

`reduce` accepts `sum`, `prod`, `mean`, `amax`, `amin` (`groupby.py:28`, passed to
`torch.scatter_reduce`). Binning a fine raster onto coarser cells needs the spread
*inside* the cell as well as its mean: that spread is the value's uncertainty, which
is the whole point of a declared support.

healpix-connector reimplements this in `binning.py` (mean, std, min, max, pixel count,
coverage per cell) purely because the reduction list stops at `mean`. With `std`/`var`
and the per-cell sample count exposed, that module disappears.

Measured use: CHELSA v2.1 (~1 km) onto depth 8 over Iberia, 1,849 cells; the within-cell
standard deviation of annual precipitation reaches 413 kg m-2, i.e. the same order as
the differences the science argues about.

Related, and worth keeping apart: `ConservativeResampler` (PR #45) takes per-sample
`area` and forwards `ellipsoid`, so area weighting is covered there — but it rejects
`out_cell_ids`, so it cannot answer "these cells, please".

## 1b. healpix-resample — `OverlapConservativeResampler` is spherical only

`overlap_conservative.py` (PR #62, merged 2026-09-11, not yet in the released
`2026.7.0`) hard-codes `ellipsoid="sphere"` (line 398) and says so: "the current
implementation is **spherical** ... an authalic-ellipsoid variant would only change
the lat -> z mapping".

That makes its cell ids a different grid from every other resampler's, which default to
WGS84. GRID4EARTH's own convention is WGS84, so overlap-conservative output cannot be
matched to it without a silent shift — ~20 km at depth 7, about 40 % of a cell.

Ask: an `ellipsoid` option, as the docstring already anticipates.

## 2. healpix-resample — `PSFResampler` has no non-negativity constraint

`resample()` solves a linear system, so it rings: on a non-negative field it returns
negative values with no warning.

Measured: Climate DT total precipitation rate, January 2030 monthly mean, depth 7,
480 Iberian cells as `out_cell_ids`, `threshold=0.5`, `ellipsoid="WGS84"`:

| call | cells < 0 |
|---|---|
| `lam=0.0` | 114 / 480 |
| `lam=0.01` | 107 / 480 |
| `lam=0.1` | 71 / 480 |
| `NearestResampler` | 0 / 480 |

Downstream that produced monthly precipitation totals of -193 mm and a coefficient of
variation of -2222 % where the mean crossed zero.

Ask: say so in the docstring, and ideally offer a constrained solve (or a documented
clip) for fields that cannot be negative. Regularisation alone does not fix it.

## 3. healpix-convert — `ClimateDTConverter`: accept a date range

`ClimateDTConverter(date=..., time=...)` takes one date and builds its cache name from
it (`_grib_path`). Polytope accepts a MARS date range, and the request the converter
already builds passes it straight through — but `"20300101/to/20300131"` makes an
unusable file name, so a range cannot be cached.

It matters at climate length: a decade of two variables twice a day is **120 requests**
as month-long ranges against **~90,000** as single timesteps. healpix-connector
subclasses the converter (`sources/climatedt.MonthConverter`) to do nothing but
sanitise that name.

Ask: accept a date range and build a filename-safe cache tag. A public `download()` /
`dataset()` pair would also help: we currently call `_download`/`_load_dataset`, since
`prepare()` writes a zarr skeleton we do not want.

## 4. healpix-convert — offer: a CHELSA (and WorldClim) converter

healpix-connector carries `sources/chelsa.py`: CHELSA v2.1 bioclimatic GeoTIFFs
(~1 km, CC0, DOI 10.16904/envidat.228) to zarr-conventions/dggs v1 with the CF grid
mapping and a STAC item — the same output shape as the ERA5, CAMS and Climate DT
converters, for the reference climatology biodiversity modelling actually uses.

One correction it carries that belongs upstream too: CHELSA's own specification gives
`bio4` in degC, but the published values are 100x the standard deviation (495.75
against 4.96 computed from ERA5 for the same cells and period).

Offer: move it in beside the others (it needs only rasterio; healpix-convert already
depends on rasterix), and add WorldClim on the same path.
