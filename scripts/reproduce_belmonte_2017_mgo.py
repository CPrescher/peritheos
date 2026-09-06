#!/usr/bin/env python3
"""Reproduce Belmonte (2017) Table 4 static MgO pressures."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "mgo_belmonte_2017_b3lyp_static_bm3"


def reproduce() -> dict[str, object]:
    document = get_material_document("mgo")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == IDENTIFIER
    )
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    path = ROOT / "peritheos/data/datasets/mgo-belmonte-2017-table4-300k.csv"
    rows = []
    with path.open(newline="", encoding="utf-8") as stream:
        for source in csv.DictReader(stream):
            ratio = float(source["volume_ratio"])
            observed = float(source["static_pressure_gpa"])
            calculated = float(eos.pressure(ratio * record["eos"]["parameters"]["V0"]))
            residual = calculated - observed
            assert abs(residual) < 0.005
            rows.append(
                {
                    "volume_ratio": ratio,
                    "published_static_pressure_gpa": observed,
                    "calculated_static_pressure_gpa": calculated,
                    "residual_gpa": residual,
                }
            )
    return {
        "record": IDENTIFIER,
        "observations": len(rows),
        "max_abs_residual_gpa": max(abs(row["residual_gpa"]) for row in rows),
        "rows": rows,
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
