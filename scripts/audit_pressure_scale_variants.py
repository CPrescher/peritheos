"""Inventory curated pressure-scale alternatives without inferring new groups."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/pressure-scale-variants.json"


def inventory() -> dict:
    groups = []
    materials = set()
    for path in sorted((ROOT / "peritheos/data/materials").glob("*.eosmat")):
        document = json.loads(path.read_text())
        records = {r["identifier"]: r for r in document["eos_records"]}
        variants = defaultdict(list)
        for record in records.values():
            if "pressure_scale_variant_of" in record:
                variants[record["pressure_scale_variant_of"]].append(
                    record["identifier"]
                )
        for source, references in sorted(variants.items()):
            materials.add(document["identifier"])
            groups.append(
                {
                    "material": document["identifier"],
                    "source_eos_record": source,
                    "reference_variants": sorted(references),
                    "primary_reference": records[source]["reference"],
                }
            )
    return {
        "policy": "One source per study/model/data selection; published pressure-scale alternatives remain explicit reference variants. No coefficients or calibration routes are inferred by grouping.",
        "source_groups": len(groups),
        "reference_variants": sum(len(g["reference_variants"]) for g in groups),
        "materials": len(materials),
        "groups": groups,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = json.dumps(inventory(), indent=2, ensure_ascii=False) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != payload:
            raise SystemExit("Pressure-scale variant inventory is stale")
    else:
        OUTPUT.write_text(payload)


if __name__ == "__main__":
    main()
