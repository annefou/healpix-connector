# healpix-connector

Connect existing data repositories to the [HEALPix](https://healpix.sourceforge.io/) discrete global grid, the way the [GRID4EARTH](https://grid4earth.eu) project uses it: NESTED ordering on the WGS84 ellipsoid.

Every value it returns states its **support** (HEALPix depth), its **valid time** and its **uncertainty**, together with the identifiers needed to cite and reproduce it.

> **Status: early development.** The CHELSA converter (`sources.chelsa`), the Climate DT source (`sources.climatedt`), the GBIF connector (`connectors.gbif`) and sampling with support, valid time and uncertainty (`sample`) work; see `examples/`.
>
> **First scientific tests, done.** [sdm-hotspot-spatial-effort](https://github.com/annefou/sdm-hotspot-spatial-effort) was rerun four times with only its covariates changed: as declared-support cell means instead of centre pixels (test 1), from ERA5 for the records' own period instead of a 1981-2010 climatology (test 3), and from a Climate DT projection, one decade against another (test 4). Its *Contradicted* verdict survives all of them; misidentification moves between 94.3 % and 96.7 %, against a 47.8-68.6 % reference range, and every effect measured is the size of the 0.57 pp gap between our baseline and the published number. See `examples/results/test1_verdict.json`, `test3_verdict.json` and `test4_verdict.json`.
>
> **One finding worth its own line.** CHELSA's `bio4` is 100x the standard deviation, which neither its file specification nor the GeoTIFF's own GDAL tags say, and its `bio15` is a percentage documented as kg m-2. Both were checked against CHELSA's own monthly layers: `examples/verify_chelsa_scaling.py`.

## What it will do

- **Occurrence connectors.** Read records from GBIF, OBIS and the Living Atlases for a region (bounding box, polygon or HEALPix cells), and return each record with its cell, the cells its positional uncertainty covers, the download DOI and the taxonomy used. Records are read from the repository on request, never mirrored.
  - GBIF works today: name matching, region search and existing downloads. It always states the taxonomy (`checklistKey`), because GBIF.org now defaults to the Catalogue of Life while its API still defaults to the legacy backbone, and it retries when GBIF rate-limits or errors.
- **Environmental sources.** Any dataset published as Zarr following [zarr-conventions/dggs](https://github.com/zarr-conventions/dggs) and described by a STAC item can be used as a source. First converter: CHELSA climatologies (~1 km). Climate DT (DestinE) is read through GRID4EARTH's own [healpix-convert](https://github.com/GRID4EARTH/healpix-convert) converter, which defines the request, the conventions and the sphere-to-WGS84 resampling; we add only a month per request, monthly aggregation and the restriction to the cells asked for.
- **Matching.** Attach environmental values to records at a declared depth, with the record's date checked against the dataset's period, and uncertainty kept as separate components (the source's own spread inside the cell, and the spread across the cells the record's positional uncertainty covers). Depth differences are reconciled explicitly: a coarser source is *inherited*, a finer one *aggregated*, and the result says which.

## Design principles

- **Library first.** Everything is a Python function; a CLI and an MCP server are thin layers on top. Nothing requires an LLM.
- **Built on GRID4EARTH tools**: [healpix-geo](https://github.com/GRID4EARTH/healpix-geo), [healpix-resample](https://github.com/GRID4EARTH/healpix-resample), [healpix-convert](https://github.com/GRID4EARTH/healpix-convert). Generic improvements go upstream to them; what is held here temporarily, and the issues asking for it, are listed in [`docs/upstream-issues.md`](docs/upstream-issues.md).

## Development

```bash
pixi run -e test test
```

## Credits

Developed by [LifeWatch ERIC](https://www.lifewatch.eu) within GRID4EARTH, an ESA Digital Twin Earth project.

## Licence

[MIT](LICENSE).
