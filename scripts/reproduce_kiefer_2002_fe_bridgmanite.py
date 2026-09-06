#!/usr/bin/env python3
"""Reproduce Kiefer et al. (2002) Mg/Fe bridgmanite EOS checkpoints."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

RECORDS = {
    "bridgmanite": ("bridgmanite_kiefer_2002_gga_bm3", 159.6, 264.3, 3.94),
    "mg075fe025sio3_bridgmanite": (
        "mg075fe025sio3_bridgmanite_kiefer_2002_gga_bm3",
        160.5,
        270.0,
        3.96,
    ),
}


def reproduce() -> dict[str, object]:
    results = {}
    for material_identifier, (identifier, v0, k0, k0_prime) in RECORDS.items():
        document = get_material_document(material_identifier)
        record = next(
            row for row in document["eos_records"] if row["identifier"] == identifier
        )
        assert record["eos"]["parameters"] == {
            "V0": v0,
            "K0": k0,
            "K0_prime": k0_prime,
        }
        eos = Material.from_eosmat(
            document, record_identifiers=[identifier]
        ).get_eos_record(identifier)
        pressure = float(eos.pressure(0.8 * v0))
        assert isclose(float(eos.volume(pressure)), 0.8 * v0, rel_tol=2e-10)
        results[identifier] = {
            "pressure_at_0.8_v0_gpa": pressure,
            "roundtrip_volume_a3_cell": float(eos.volume(pressure)),
        }
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
