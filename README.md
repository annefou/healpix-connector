# healpix-connector

Connect existing data repositories to the [HEALPix](https://healpix.sourceforge.io/) discrete global grid, the way the [GRID4EARTH](https://grid4earth.eu) project uses it: NESTED ordering on the WGS84 ellipsoid.

Every value it returns states its **support** (HEALPix depth), its **valid time** and its **uncertainty**, together with the identifiers needed to cite and reproduce it.

> **Status: early development.** The CHELSA converter (`sources.chelsa`), the Climate DT sources (`sources.climatedt`, generation 1 hourly; `sources.climatedt_monthly`, generation 2 monthly means at depth 10, added 2026-10-23), the GBIF and OBIS connectors (`connectors.gbif`, `connectors.obis`) and sampling with support, valid time and uncertainty (`sample`) work; see `examples/`.
>
> **First scientific tests, done.** [sdm-hotspot-spatial-effort](https://github.com/annefou/sdm-hotspot-spatial-effort) was rerun four times with only its covariates changed: as declared-support cell means instead of centre pixels (test 1), from ERA5 for the records' own period instead of a 1981-2010 climatology (test 3), and from a Climate DT projection, one decade against another (test 4). Its *Contradicted* verdict survives all of them; misidentification moves between 94.3 % and 96.7 %, against a 47.8-68.6 % reference range, and every effect measured is the size of the 0.57 pp gap between our baseline and the published number. See `examples/results/test1_verdict.json`, `test3_verdict.json` and `test4_verdict.json`.
>
> **One finding worth its own line.** CHELSA's `bio4` is 100x the standard deviation, which neither its file specification nor the GeoTIFF's own GDAL tags say, and its `bio15` is a percentage documented as kg m-2. Both were checked against CHELSA's own monthly layers: `examples/verify_chelsa_scaling.py`.

## What it will do

- **Occurrence connectors.** Read records from GBIF, OBIS and the Living Atlases for a region (bounding box, polygon or HEALPix cells), and return each record with its cell, the cells its positional uncertainty covers, the download DOI and the taxonomy used. Records are read from the repository on request, never mirrored.
  - GBIF and OBIS work today. GBIF: name matching, region search and existing downloads. It always states the taxonomy (`checklistKey`), because GBIF.org now defaults to the Catalogue of Life while its API still defaults to the legacy backbone, and it retries when GBIF rate-limits or errors. OBIS: WoRMS AphiaIDs, cursor paging, region counts before paging, and the repository's own quality flags passed through; since OBIS mints no DOI for a query, a result is pinned by its dataset ids and retrieval date.
- **Environmental sources.** Any dataset published as Zarr following [zarr-conventions/dggs](https://github.com/zarr-conventions/dggs) and described by a STAC item can be used as a source. First converter: CHELSA climatologies (~1 km). Climate DT (DestinE) is read through GRID4EARTH's own [healpix-convert](https://github.com/GRID4EARTH/healpix-convert) converter, which defines the request, the conventions and the sphere-to-WGS84 resampling; we add only a month per request, monthly aggregation and the restriction to the cells asked for.
- **Matching.** Attach environmental values to records at a declared depth, with the record's date checked against the dataset's period, and uncertainty kept as separate components (the source's own spread inside the cell, and the spread across the cells the record's positional uncertainty covers). Depth differences are reconciled explicitly: a coarser source is *inherited*, a finer one *aggregated*, and the result says which.

## Design principles

- **Library first.** Everything is a Python function; a CLI and an MCP server are thin layers on top. Nothing requires an LLM.
- **Built on GRID4EARTH tools**: [healpix-geo](https://github.com/GRID4EARTH/healpix-geo), [healpix-resample](https://github.com/GRID4EARTH/healpix-resample), [healpix-convert](https://github.com/GRID4EARTH/healpix-convert). Generic improvements go upstream to them; what is held here temporarily, and the issues asking for it, are listed in [`docs/upstream-issues.md`](docs/upstream-issues.md).

## Development

```bash
pixi run -e test test
```

## How to cite

### This software

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22904850.svg)](https://doi.org/10.5281/zenodo.22904850)

Cite the **version you used**. Every result carries it as `producer` (e.g. `"healpix-connector 0.2.0"`) in its provenance.

> Fouilloux, A., & LifeWatch ERIC (2026). *healpix-connector: connect existing data repositories to the HEALPix grid, with declared support, valid time and uncertainty* (Version 0.2.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23603498

```bibtex
@software{healpix_connector_0_2_0,
  author    = {Fouilloux, Anne and {LifeWatch ERIC}},
  title     = {healpix-connector: connect existing data repositories to the HEALPix grid,
               with declared support, valid time and uncertainty},
  version   = {0.2.0},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23603498},
  url       = {https://doi.org/10.5281/zenodo.23603498}
}
```

| Release | DOI |
|---|---|
| all versions (concept DOI, resolves to the latest) | [10.5281/zenodo.22904850](https://doi.org/10.5281/zenodo.22904850) |
| v0.2.0 | [10.5281/zenodo.23603498](https://doi.org/10.5281/zenodo.23603498) |
| v0.1.0 | [10.5281/zenodo.22904851](https://doi.org/10.5281/zenodo.22904851) |

The metadata is in [`CITATION.cff`](CITATION.cff); GitHub's *Cite this repository* button reads it and exports APA and BibTeX.

### The data behind a result

Citing this library does not cite the data. Each repository has its own rules, and every connector and source returns what they ask for in its provenance, so a result can be cited without looking anything up again.

| Data | What to cite | Where healpix-connector puts it |
|---|---|---|
| **GBIF**, published results | The **download DOI**, as GBIF requires; also state the taxonomy (`checklistKey`). Make a download for anything you publish: search results have no DOI. | `gbif.download_records(...).provenance`: `doi`, `checklist_key`, `download_created`, `predicate` |
| GBIF, exploration | Not citable as it stands (no DOI). | `gbif.search(...).provenance`: `params`, `checklist_key`, `retrieved_at` |
| **OBIS** | OBIS mints no DOI for a query: cite **each dataset**, plus the retrieval date. | `obis.search(...).provenance`: `datasets`, `retrieved_at`; then `obis.dataset_metadata(id)` for each dataset's `citation` and `doi` (`None` when it has none) |
| **Darwin Core Archive** (IPT, OBIS export) | The archive's DOI where its EML states one, otherwise its title and package id. | `dwca.read(...).provenance`: `archive_doi`, `archive_title`, `retrieved_at` |
| **CHELSA** v2.1 | Karger et al., CHELSA v2.1, EnviDat, [10.16904/envidat.228](https://doi.org/10.16904/envidat.228) (CC0-1.0), and the CHELSA paper: Karger et al. (2017), *Scientific Data* 4, 170122, [10.1038/sdata.2017.122](https://doi.org/10.1038/sdata.2017.122). | `sample.provenance(src, depth)`: `dataset_doi`, `license`; the STAC item has a `cite-as` link |
| **Climate DT** (Destination Earth) | No DOI. Name the model, experiment, generation, resolution and period, and follow the terms of use of your DestinE access. | `climatedt.SOURCE`, `climatedt_gen2.json`; the dataset's `source`, `source_access` and `resampling_method` attributes |
| Any other DGGS dataset (e.g. ERA5) | Its own DOI and licence. Set them when you write the dataset, or the result cannot carry them. | the dataset's `source_doi` and `license` attributes (e.g. `CoarseSample.to_dataset(..., source_doi=, license=)`) → `sample.provenance`: `dataset_doi`, `license` |

When values were resampled or aggregated onto the grid, `sample.provenance` also states `source_depth`, `sampled_at_depth`, `resampling_method` and `period`. Report these in the methods section: the same cell id at a different depth, or from a different resampler, is a different number.

### The GRID4EARTH tools it builds on

If the grid handling matters to your result, cite these too: [healpix-resample](https://github.com/GRID4EARTH/healpix-resample) ([10.5281/zenodo.21723671](https://doi.org/10.5281/zenodo.21723671)), [healpix-geo](https://github.com/GRID4EARTH/healpix-geo) and [healpix-convert](https://github.com/GRID4EARTH/healpix-convert).

## Credits

Developed by [LifeWatch ERIC](https://www.lifewatch.eu) within GRID4EARTH, an ESA Digital Twin Earth project.

## Licence

[MIT](LICENSE).
