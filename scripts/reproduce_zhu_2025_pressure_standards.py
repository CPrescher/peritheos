#!/usr/bin/env python3
"""Execute the supported 300 K Vinet branches from Zhu et al. (2025)."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

RECORDS = {
    "platinum": ("platinum_zhu_2025_vinet_300k", 60.38, 277.3, 5.230),
    "gold": ("gold_zhu_2025_vinet_300k", 67.85, 167.0, 5.897),
    "mgo": ("mgo_zhu_2025_vinet_300k", 74.71, 160.3, 4.182),
}


def reproduce() -> dict[str, object]:
    results = {}
    for material_identifier, (identifier, v0, k0, k0_prime) in RECORDS.items():
        document = get_material_document(material_identifier)
        source = next(
            row for row in document["eos_records"] if row["identifier"] == identifier
        )
        assert source["eos"]["parameters"] == {
            "V0": v0,
            "K0": k0,
            "K0_prime": k0_prime,
        }
        model = Material.from_eosmat(
            document, record_identifiers=[identifier]
        ).get_eos_record(identifier)
        checkpoints = {}
        for fraction in (0.95, 0.80, 0.65):
            target = fraction * v0
            pressure = float(model.pressure(target))
            recovered = float(model.volume(pressure))
            assert isclose(recovered, target, rel_tol=2e-10)
            checkpoints[str(fraction)] = pressure
        results[identifier] = {"pressure_gpa": checkpoints}
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
