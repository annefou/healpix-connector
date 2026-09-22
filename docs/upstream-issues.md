# Upstream issues (drafts, not filed)

Five things healpix-connector had to work around. Each one says what we were doing,
what happened, and what we would like. Everything is measured from real use; nothing is
filed until Anne has read it.

Checked against `main` on 2026-09-22 (healpix-resample `2026.7.0` installed, `main`
newer; healpix-convert `2026.9.0`).

**Landscape.** healpix-resample: 5 open PRs (4 dependabot, #57 notebooks environment)
and 5 open issues, none covering these. healpix-convert: 1 open issue (#9, licence
links) and 7 open PRs, 6 of them ours from 2026-09-07, still unreviewed.

---

## 1. healpix-resample — a cell's mean without its spread

**What we were doing.** Putting CHELSA climate data (~1 km pixels) onto HEALPix cells
at depth 8 (~25 km) over Iberia: about 1,800 cells, each covering roughly 650 pixels.
A cell's value is the mean of its pixels, and we also need how much those pixels vary,
because that spread is the uncertainty of the value we hand back.

**What happened.** `GroupByResampler` can return the mean of the pixels in a cell, but
not their standard deviation: `reduce` accepts `sum`, `prod`, `mean`, `amax`, `amin`
(`groupby.py:28`, passed to `torch.scatter_reduce`). Nor does it say how many pixels
went into each cell, which is what tells you whether a cell is well covered or has two
pixels in a corner.

So healpix-connector computes all of it itself (`binning.py`: mean, standard deviation,
min, max, pixel count, coverage), duplicating the binning healpix-resample already does.

**Why it matters.** The spread is not a diagnostic, it is the answer: for annual
precipitation over Iberia the within-cell standard deviation reaches 413 kg m-2, the
same size as the differences these studies argue about. Anyone who reports a cell mean
without it is reporting a number that looks more precise than it is.

**What we would like.** `std` and `var` as reductions, and the per-cell sample count
exposed. Then `binning.py` disappears.

**Related, but a different class.** `ConservativeResampler` (PR #45) already takes
per-sample areas and honours `ellipsoid`, so area-weighted means are covered there.
What neither it nor `GroupByResampler` accepts is `out_cell_ids` — the parameter that
lets a caller say "return values for exactly these cells". They return whichever cells
happened to contain an input point, so a caller filling a region someone asked for gets
back a different list of cells and has to match it up; cells with no input simply do not
appear. That is a separate, smaller wish.

## 2. healpix-resample — `OverlapConservativeResampler` works on a sphere, everything else on WGS84

**What we were doing.** Nothing yet — we looked at it for precipitation, which needs a
conservative method, and stopped.

**What happened.** `overlap_conservative.py` (PR #62, merged 2026-09-11, not in the
released `2026.7.0`) hard-codes `ellipsoid="sphere"` (line 398). Its own docstring says
so: "the current implementation is **spherical** ... an authalic-ellipsoid variant would
only change the lat -> z mapping".

**Why it matters.** Every other resampler defaults to WGS84, and WGS84 is GRID4EARTH's
convention. The same cell id therefore means a different patch of ground depending on
which resampler produced it — about 20 km apart at depth 7, roughly 40 % of a cell.
Nothing warns you; the ids line up perfectly and the data is simply in the wrong place.

**What we would like.** An `ellipsoid` option, as the docstring already anticipates.

## 3. healpix-resample — `PSFResampler` returns negative rainfall

**What we were doing.** Climate DT arrives on HEALPix depth 7 but on a sphere, so we
used the PSF resampler — healpix-convert's own optional correction — to move 2 m
temperature and precipitation onto WGS84 cells.

**What happened.** Temperature was fine. Precipitation came back negative in a quarter
of the cells. PSF solves a linear system, so it overshoots at sharp edges, and nothing
stops the solution going below zero.

January 2030 monthly mean, depth 7, 480 Iberian cells as `out_cell_ids`,
`threshold=0.5`, `ellipsoid="WGS84"`:

| call | cells below zero |
|---|---|
| `lam=0.0` | 114 of 480 |
| `lam=0.01` | 107 of 480 |
| `lam=0.1` | 71 of 480 |
| `NearestResampler` | 0 of 480 |

Carried through to bioclimatic variables that gave monthly precipitation totals of
-193 mm, and a coefficient of variation of -2222 % in cells where the mean crossed zero.

**Why it matters.** It is silent. The output is a plausible-looking field, and the
failure only shows up later as an impossible number — or not at all, if nobody checks.

**What we would like.** A sentence in the docstring saying PSF is unsuitable for
quantities that cannot be negative, and, if it is cheap, a constrained solve or a
documented clip. Turning up `lam` is not a fix: it only halves the count.

## 4. healpix-convert — `ClimateDTConverter` can only ask for one timestep at a time

**What we were doing.** Twenty years of Climate DT (SSP3-7.0), two variables, twice a
day, to compare two decades of the same simulation.

**What happened.** `ClimateDTConverter(date=..., time=...)` takes a single date, and
builds its cache file name from it (`_grib_path`). Polytope itself accepts a MARS date
range, and the request the converter already builds passes one through untouched — but
`"20300101/to/20300131"` cannot be a file name, so a range cannot be cached.

Asked a month at a time, twenty years is **240 requests**. Asked a timestep at a time,
it is about **175,000**.

healpix-connector therefore subclasses the converter (`sources/climatedt.MonthConverter`)
for no reason other than to sanitise that file name.

**What we would like.** Accept a date range and build a file-name-safe cache tag. A
public `download()` and `dataset()` would help too: we currently call `_download` and
`_load_dataset`, because `prepare()` also writes a zarr skeleton we do not want.

## 5. healpix-convert — an offer: a CHELSA converter (and WorldClim after it)

**What we were doing.** Biodiversity models are mostly built on CHELSA or WorldClim
bioclimatic layers, not on reanalysis. We needed CHELSA on HEALPix cells.

**What we have.** `sources/chelsa.py` in healpix-connector converts CHELSA v2.1
GeoTIFFs (~1 km, CC0, DOI 10.16904/envidat.228) to zarr-conventions/dggs v1 with the CF
grid mapping and a STAC item — the same output shape as the ERA5, CAMS and Climate DT
converters. It needs rasterio; healpix-convert already depends on rasterix.

**One correction it carries.** CHELSA's own specification gives `bio4` in degC, but the
published values are 100x the standard deviation: 495.75 against 4.96 computed from ERA5
for the same cells and period. Anyone converting CHELSA will hit this.

**What we would like.** To move it in beside the other converters, and add WorldClim on
the same path — if you want it there. Say the word and we open the PR.
