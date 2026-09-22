# healpix-connector

Connect existing data repositories to the [HEALPix](https://healpix.sourceforge.io/) discrete global grid, the way the [GRID4EARTH](https://grid4earth.eu) project uses it: NESTED ordering on the WGS84 ellipsoid.

Every value it returns states its **support** (HEALPix depth), its **valid time** and its **uncertainty**, together with the identifiers needed to cite and reproduce it.

> **Status: early development.** The CHELSA converter (`sources.chelsa`), the GBIF connector (`connectors.gbif`) and sampling with support, valid time and uncertainty (`sample`) work; see `examples/`.
>
> **First scientific test, done.** [sdm-hotspot-spatial-effort](https://github.com/annefou/sdm-hotspot-spatial-effort) was rerun with only its covariates replaced by declared-support cell means. Its *Contradicted* verdict survives: 96.67 % of hotspots misidentified with the study's cell-centre pixels, 96.09 % with cell means, against a 47.8-68.6 % reference range (`examples/results/test1_verdict.json`).

## What it will do

- **Occurrence connectors.** Read records from GBIF, OBIS and the Living Atlases for a region (bounding box, polygon or HEALPix cells), and return each record with its cell, the cells its positional uncertainty covers, the download DOI and the taxonomy used. Records are read from the repository on request, never mirrored.
  - GBIF works today: name matching, region search and existing downloads. It always states the taxonomy (`checklistKey`), because GBIF.org now defaults to the Catalogue of Life while its API still defaults to the legacy backbone, and it retries when GBIF rate-limits or errors.
- **Environmental sources.** Any dataset published as Zarr following [zarr-conventions/dggs](https://github.com/zarr-conventions/dggs) and described by a STAC item can be used as a source. First converter: CHELSA climatologies (~1 km).
- **Matching.** Attach environmental values to records at a declared depth, with the record's date checked against the dataset's period, and uncertainty kept as separate components (the source's own spread inside the cell, and the spread across the cells the record's positional uncertainty covers). Depth differences are reconciled explicitly: a coarser source is *inherited*, a finer one *aggregated*, and the result says which.

## Design principles

- **Library first.** Everything is a Python function; a CLI and an MCP server are thin layers on top. Nothing requires an LLM.
- **Built on GRID4EARTH tools**: [healpix-geo](https://github.com/GRID4EARTH/healpix-geo), [healpix-resample](https://github.com/GRID4EARTH/healpix-resample), [healpix-convert](https://github.com/GRID4EARTH/healpix-convert). Generic improvements go upstream to them.

## Development

```bash
pixi run -e test test
```

## Credits

Developed by [LifeWatch ERIC](https://www.lifewatch.eu) within GRID4EARTH, an ESA Digital Twin Earth project.

## Licence

[MIT](LICENSE).
