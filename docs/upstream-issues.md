# Upstream issues

Five requests, written to stand on their own: a maintainer should not need to know what
healpix-connector is to act on them. **All five were filed on 2026-09-22:**
healpix-resample [#63](https://github.com/GRID4EARTH/healpix-resample/issues/63) (1),
[#64](https://github.com/GRID4EARTH/healpix-resample/issues/64) (2),
[#65](https://github.com/GRID4EARTH/healpix-resample/issues/65) (3); healpix-convert
[#12](https://github.com/GRID4EARTH/healpix-convert/issues/12) (4),
[#11](https://github.com/GRID4EARTH/healpix-convert/issues/11) (5). Checked against `main` on 2026-09-22
(healpix-resample `2026.7.0` released, `main` newer; healpix-convert `2026.9.0`).

**Status, checked 2026-11-04.**

| # | Issue | State |
|---|---|---|
| 1 | resample #63 | open, no reply |
| 2 | resample #64 | open, no reply |
| 3 | resample #65 | **closed 2026-10-15** by PR #72: `PSFResampler.resample(..., clip=(lo, hi), max_outer=10)`, a bounded (active-set) least-squares fit; `ResampleResults.n_clipped_low/high`. In releases >= 2026.10.2. Also closes healpix-convert #35 (use the bounded fit for bounded fields). Used here since PR #4 (`healpix-resample >= 2026.10.2`; `clip=(0.0, None)` for `tprate`). |
| 4 | convert #12 | open, no reply |
| 5 | convert #11 | **answered 2026-10-05** (benbovy): see below |
| 6 | convert [#104](https://github.com/GRID4EARTH/healpix-convert/issues/104) | **filed 2026-11-04** (reworded after #11; paramIds and base request added) |

**The answer to #11.** healpix-convert aims at a shared core (one internal data model in
`ConvertStagingCache`/`ConvertSettings`, generic group converters in `healpix_converters`,
STAC + Zarr-conventions output configured through settings) plus built-in or third-party
*settings*. The per-dataset modules in `healpix_convert/converters` (ERA5, CAMS, Climate DT)
are temporary and will be refactored onto the core. Built-in support for a product is
considered when the source is public, stable and mature; its conversion fits a small
notebook runnable in CI with little customisation; and it needs no dependency beyond the
core's. A mechanism to coerce custom input metadata into the data model is planned; a
registry for third-party settings is #18. Consequence for us: CHELSA stays downstream until
third-party settings exist; requests about Climate DT should describe the service, not
patch `converters/climatedt.py`.

**Landscape (2026-09-22).** healpix-resample: 5 open PRs (4 dependabot, #57 notebooks environment)
and 5 open issues, none covering these. healpix-convert: 1 open issue (#9, licence
links) and 7 open PRs, 6 of them ours from 2026-09-07, still unreviewed.

---

## 1. healpix-resample — add `std` and `var` reductions to `GroupByResampler`, and expose the per-cell sample count

**Request.** Let `reduce` take `std` and `var`, and make the number of samples that fell
into each cell available on the result.

**Today.** `reduce` accepts `sum`, `prod`, `mean`, `amax`, `amin` (`groupby.py:28`,
passed to `torch.scatter_reduce`), and the per-cell count is computed internally but not
returned.

**Why.** When a fine field is binned onto coarser cells, the spread of the samples
inside a cell is the uncertainty of the cell's value, and the count says whether the
cell is well covered or holds two samples in a corner. Both come free from the grouping
that has already been done; without them every caller recomputes the binning to get
them.

**Example.** CHELSA v2.1 (~1 km) onto depth 8 over Iberia, ~1,800 cells of ~650 pixels
each: the within-cell standard deviation of annual precipitation reaches 413 kg m-2 —
the same size as the differences such studies are arguing about, so a mean reported
without it looks more precise than it is.

**Smaller, related wish.** Neither `GroupByResampler` nor `ConservativeResampler`
accepts `out_cell_ids`, so a caller cannot ask for a fixed set of cells: the result
holds whichever cells happened to receive a sample, and cells with none are absent
rather than empty.

## 2. healpix-resample — give `OverlapConservativeResampler` an `ellipsoid` option

**Request.** Let `OverlapConservativeResampler` work on WGS84, not only on a sphere.

**Today.** `overlap_conservative.py` (PR #62, merged 2026-09-11, not in released
`2026.7.0`) hard-codes `ellipsoid="sphere"` (line 398). Its docstring already notes that
"an authalic-ellipsoid variant would only change the lat -> z mapping".

**Why.** The other resamplers default to WGS84, so the same cell id means a different
patch of ground depending on which resampler produced it — about 20 km at depth 7,
roughly 40 % of a cell. The ids line up, nothing warns, and the data is in the wrong
place.

## 3. healpix-resample — document that `PSFResampler` can return negative values, and ideally guard it

**Request.** State in the docstring that PSF is unsuitable for quantities that cannot be
negative, and, if cheap, offer a constrained solve or a documented clip.

**Today.** `resample()` solves a linear system with no non-negativity constraint, so it
overshoots at sharp edges and silently returns values below zero.

**Measured.** Precipitation rate, monthly mean, depth 7, 480 cells over Iberia passed as
`out_cell_ids`, `threshold=0.5`, `ellipsoid="WGS84"`:

| call | cells below zero |
|---|---|
| `lam=0.0` | 114 of 480 |
| `lam=0.01` | 107 of 480 |
| `lam=0.1` | 71 of 480 |
| `NearestResampler` | 0 of 480 |

Downstream this produced monthly precipitation totals of -193 mm and a coefficient of
variation of -2222 % where the mean crossed zero. Raising `lam` is not a fix: it only
halves the count.

**Why.** The failure is silent — the field looks plausible and only becomes visibly
wrong after further arithmetic, or never.

## 4. healpix-convert — let `ClimateDTConverter` take a date range

**Request.** Accept a MARS date range (`20300101/to/20300131`) and build a
file-name-safe cache tag from it. A public `download()` and `dataset()` alongside
`prepare()` would help too.

**Today.** The converter takes a single date and builds its cache file name from it
(`_grib_path`), so a range — which Polytope accepts, and which the request already
passes through untouched — cannot be cached. `prepare()` is the only public way in, and
it also writes a zarr skeleton, so code that only wants the GRIB calls `_download` and
`_load_dataset`.

**Why.** Climate work asks for years. Twenty years of two surface variables twice a day
is 240 requests a month at a time, against 14,610 a timestep at a time (7,305 days,
twice a day; a request already carries both variables).

## 5. healpix-convert — where should converters for new datasets live?

**Question first.** Is healpix-convert meant to grow one converter per dataset, so that
a new source is contributed here beside ERA5, CAMS and Climate DT? Or is it a core plus
a pattern, with each project carrying its own converters and only the conventions,
STAC and grid handling shared?

The answer decides where anyone bringing a new dataset should put their work, and it is
not written down anywhere we could find. If the answer is "here", it would help to say
what a converter has to provide to be accepted — output conventions, STAC item, tests,
which dependencies are acceptable.

**Concrete case.** If the answer is "here", we would like to contribute a converter for
CHELSA v2.1 bioclimatic layers (~1 km GeoTIFFs, CC0, DOI 10.16904/envidat.228), and
WorldClim on the same path afterwards. Species-distribution and biodiversity work is
built on those far more than on reanalysis, so they are the layers people arrive with.

Working code exists: it writes zarr-conventions/dggs v1 with the CF grid mapping and a
STAC item, the same output shape as the existing converters, and needs rasterio
(healpix-convert already depends on rasterix). Happy to open it as a PR, or to keep it
downstream and follow whatever pattern you prefer.


---

## 6. healpix-convert — Climate DT generation 2: what the service offers, for the settings (filed 2026-11-04 as healpix-convert #104)

**Title.** Climate DT generation 2 (incl. the monthly `clmn` stream): service facts for the settings

**Request.** When Climate DT moves onto the core, please cover generation 2, which is what
DestinE serves now. The settings today describe generation 1 only (`CDT_GENERATION = "1"`,
hourly `clte`, nside 128). We are not proposing a patch to `converters/climatedt.py`, since
#11 says it will be refactored; this is the information the new settings will need, measured
against the service on 2026-10-23.

**Generation 2, as served.**

| | |
|---|---|
| models | IFS-FESOM (LUMI), IFS-NEMO (MareNostrum 5) |
| streams | `clte` (hourly) and **`clmn` (monthly means)** |
| resolutions | `high` = nside 1024 (depth 10), `standard` = nside 128 (depth 7); native HEALPix NESTED on a sphere |
| activities | `baseline` / `hist`, 1990-2014; `projections` / `SSP3-7.0`, 2015-2049 |
| monthly fields | `avg_2t`, `avg_tprate`, `avg_2d`; **no monthly max/min temperature** |

**Two traps.**

- **Short names can resolve to a different field.** Requesting `mx2t24` on `clmn` returned
  `avg_snswrfcs` (paramId 235051), with no error. Request by numeric paramId, and check the
  `paramId` of every returned GRIB message against what was asked.
- **MareNostrum 5 drops large transfers** (~500 MB). One parameter and half a year per request
  goes through; a year of two fields at depth 10 is four requests of ~145 MB.

**Why.** Generation 1 has no monthly stream and covers 2020-2039 only, so decadal work means
hundreds of hourly requests (#12). Generation 2's `clmn` covers 1990-2049 in a handful.

**Related.** For the sphere -> WGS84 step on same-size grids, a conservative resampler
(bilinear in conservative mode, normalised) is within 0.07 degC of an exact polygon overlap
at depth 10, against 0.80 degC for taking the same cell id, and cannot go negative. Worth
offering as a resampler choice for this step beside `nearest` and PSF.
For PSF on precipitation, the bounded fit from healpix-resample #72 (`clip=(0.0, None)`) is
what #35 asks to use.
