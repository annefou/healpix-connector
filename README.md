# healpix-connector

Connect existing data repositories to the [HEALPix](https://healpix.sourceforge.io/) discrete global grid, the way the [GRID4EARTH](https://grid4earth.eu) project uses it: NESTED ordering on the WGS84 ellipsoid.

Every value it returns states its **support** (HEALPix depth), its **valid time** and its **uncertainty**, together with the identifiers needed to cite and reproduce it.

> **Status: early development.** The CHELSA converter works for bio1 (`healpix_connector.sources.chelsa`); see `examples/` for the first comparison. Occurrence connectors are next.

## What it will do

- **Occurrence connectors.** Read records from GBIF, OBIS and the Living Atlases for a region (bounding box, polygon or HEALPix cells), and return each record with its cell, the cells its positional uncertainty covers, the download DOI and the taxonomy used. Records are read from the repository on request, never mirrored.
- **Environmental sources.** Any dataset published as Zarr following [zarr-conventions/dggs](https://github.com/zarr-conventions/dggs) and described by a STAC item can be used as a source. First converter: CHELSA climatologies (~1 km).
- **Matching.** Attach environmental values to records at declared depth, matched to each record's `eventDate`, with uncertainty reported as separate components.

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
