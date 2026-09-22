"""Records from OBIS for a region, with their cells, uncertainty and provenance.

Common dolphin (Delphinus delphis) off the Galician coast, 2015-2024. Shows the
three things the connector adds to a repository query: the region is translated
and then filtered exactly, each record carries the cells its uncertainty covers,
and the result says which taxonomy and which datasets it came from.

Usage: python examples/obis_region_records.py
"""

import json
from pathlib import Path

from healpix_connector.connectors import obis
from healpix_connector.region import Region

REGION = Region.from_bbox(-9.5, 42.0, -8.5, 43.0)
DEPTH = 8


def main():
    name = obis.match_name("Delphinus delphis")
    print(f"{name['query']} -> AphiaID {name['aphia_id']} ({name['backbone']}, rank {name['rank']})")

    counts = obis.statistics(REGION)
    print(f"region holds {counts['records']:,} records, {counts['species']} species, "
          f"{counts['datasets']} datasets, years {counts['year_range']}")

    res = obis.search(REGION, aphia_id=name["aphia_id"],
                      filters={"startdate": "2015-01-01", "enddate": "2024-12-31"})
    rows = obis.to_cells(res.records, DEPTH)
    known = [r for r in rows if r["footprint_known"]]
    flagged = [r for r in rows if r["flags"]]
    print(f"{res.provenance['matched_in_query_polygon']} in the query polygon, "
          f"{len(rows)} kept in the region, on {len({r['cell'] for r in rows})} cells at depth {DEPTH}")
    print(f"{len(known)} records state their positional uncertainty; {len(flagged)} carry OBIS flags")

    report = {"taxon": name, "region_statistics": counts, "provenance": res.provenance,
              "records": len(rows), "cells": sorted({r["cell"] for r in rows}),
              "with_known_uncertainty": len(known),
              "flags_seen": sorted({f for r in rows for f in r["flags"]}),
              "example": rows[0] if rows else None}
    out = Path(__file__).parent / "results/obis_region_records.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
