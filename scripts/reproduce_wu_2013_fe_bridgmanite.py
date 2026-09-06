#!/usr/bin/env python3
"""Execute Wu et al. (2013) Mg and Fe25 bridgmanite BM3 branches."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

RECORDS = {
    "bridgmanite": ("bridgmanite_wu_2013_gga_bm3", 162.378, 244.0, 4.0),
    "mg075fe025sio3_bridgmanite": (
        "mg075fe025sio3_bridgmanite_wu_2013_gga_bm3",
        164.977,
        248.0,
        4.03,
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
        recovered = float(eos.volume(pressure))
        assert isclose(recovered, 0.8 * v0, rel_tol=2e-10)
        results[identifier] = {"pressure_at_0.8_v0_gpa": pressure}
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
