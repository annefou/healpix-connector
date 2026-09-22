# CLAUDE.md — healpix-connector

Package that connects existing data repositories (GBIF, OBIS, Living Atlases, STAC
sources) to HEALPix. Design note: "GRID4EARTH HEALPix data service, with a
biodiversity extension" (shared Claude Doc; ask Anne for the link).

## Conventions (break these and results are silently wrong)

- NESTED ordering only; parent = cell >> 2.
- Always pass `ellipsoid="WGS84"` to healpix-geo: its default is `"sphere"`.
  Use `healpix_connector.conventions.ELLIPSOID`.
- Say *depth* in code and prose. healpix-resample's `level` and the zarr metadata's
  `refinement_level` are the same number; map at the boundary.
- Occurrences are read on request, never mirrored. Large requests use a GBIF
  download (record-level, with DOI).
- GBIF: always send `checklistKey` and record the taxonomy + version used. GBIF.org
  defaults to Catalogue of Life XR; the occurrence API still defaults to the legacy
  backbone. Use `gbifID` for record identity.
- Time matching keys on `eventDate`, never on publication date; report both.
- Library first: CLI and MCP only wrap library functions. No LLM in any pipeline.

## Commands

- Tests: `pixi run -e test test`

## Git

- Commit as `Anne Fouilloux <anne.fouilloux@lifewatch.eu>`. No Co-Authored-By trailer
  (destined for the GRID4EARTH organisation).

## Temporary code (to move upstream once test 1 settles the method)

- `binning.py` → healpix-resample: per-sample weights + `std`/`var` on `GroupByResampler`
  (it already supports WGS84; its reductions are unweighted mean/sum/min/max/prod).
- `sources/chelsa.py` (and a future WorldClim converter) → healpix-convert, beside ERA5,
  Climate DT and CAMS. Its converters handle non-Zarr inputs (ERA5 downloads GRIB) and it
  already depends on rasterix.
- `sources/climatedt.MonthConverter` → healpix-convert: `ClimateDTConverter` takes a single
  date and builds its cache name from it, so a MARS date range (which Polytope accepts, and
  which turns ~90,000 requests for a decade into 120) cannot be cached. Ask for a date range
  plus a filename-safe cache tag; everything else there is used unchanged.
- Also to report: `PSFResampler` has no non-negativity constraint, so on precipitation it
  rings into negative values (114 of 480 Iberian cells, a January mean; `lam` up to 0.1 only
  reduces it to 71). Worth a warning in its docstring, or a clipped/constrained variant.
- Drafted in `docs/upstream-issues.md` (2026-09-22), not filed until Anne has read them. After the upstream
  PRs land, delete the local converter; the connector reads converted data through STAC.

## Measured limits (do not rediscover)

- GBIF search: deep paging collapses, with or without a geometry filter. Measured
  2026-09-22: 0.5 s per page at offset 0, 0.4 s at 6,000, **361 s at 12,000**
  (394 s for the same query without geometry). Hence
  `gbif.PRACTICAL_PAGING_LIMIT = 5,000`; beyond that use a download.
- Climate DT (DestinE, SSP3-7.0 via Polytope): **no monthly stream** (`clmn`, `mnth`
  return 400); only hourly `clte`. Served for **2020-2039** only (2018 and 2040 are
  refused). `tp` is an hourly accumulation, too large for decadal work; use param
  **260048 `tprate`** (instantaneous rate, kg m-2 s-1). Date ranges and several
  params in one request work (a month, 2 params, 2 times ~ 46 MB in seconds).
  The GRIB is native `gridType=healpix, Nside=128, orderingConvention=nested`,
  i.e. depth 7 **on a sphere** - look study cells up with `ellipsoid="sphere"`
  rather than assuming ids match WGS84 ones.
  Measured 2026-09-22.
- The polytope client calls `sys.exit` on a refused request: catch `SystemExit`.
- Climate DT is **not** ours to re-implement: the request, endpoint, conventions and the
  sphere→WGS84 resampling come from healpix-convert (public, PyPI, byte-identical to the
  private GRID4EARTH `legacy-converters`). `healpix_resample.NearestResampler(ellipsoid=
  "WGS84")` - healpix-convert's declared resampler for these fields - reproduces a manual
  spherical-cell lookup exactly (max diff 0 K), so use it rather than hand-rolled indexing.
- healpix-convert pulls `torch`; the CUDA wheels do not fit on a normal working disk, so the
  `climatedt` feature pins the CPU wheel via the PyTorch CPU index.
- GBIF answers 429 to bursts; `_get` retries on 429 and 5xx, honouring `Retry-After`.
