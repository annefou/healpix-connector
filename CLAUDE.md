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
