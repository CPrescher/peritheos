#!/usr/bin/env python3
"""Reproduce analytical checkpoints for Li and Zeng (2009) MgSiO3 fits."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

RECORDS = {
    "bridgmanite": {
        "bridgmanite_li_zeng_2009_gga_bm3": (168.12, 234.8, 4.199),
        "bridgmanite_li_zeng_2009_gga_vinet": (168.08, 235.7, 4.294),
        "bridgmanite_li_zeng_2009_gga_natural_strain3": (168.04, 237.0, 4.384),
    },
    "mgsio3_post_perovskite": {
        "mgsio3_post_perovskite_li_zeng_2009_gga_bm3": (
            168.29,
            223.73,
            4.152,
        ),
        "mgsio3_post_perovskite_li_zeng_2009_gga_vinet": (
            168.05,
            224.2,
            4.406,
        ),
        "mgsio3_post_perovskite_li_zeng_2009_gga_natural_strain3": (
            167.75,
            226.7,
            4.684,
        ),
    },
}


def reproduce() -> dict[str, object]:
    results: dict[str, object] = {}
    for material_identifier, expected_records in RECORDS.items():
        document = get_material_document(material_identifier)
        for identifier, (v0, k0, k0_prime) in expected_records.items():
            record = next(
                row
                for row in document["eos_records"]
                if row["identifier"] == identifier
            )
            assert record["eos"]["parameters"] == {
                "V0": v0,
                "K0": k0,
                "K0_prime": k0_prime,
            }
            eos = Material.from_eosmat(
                document, record_identifiers=[identifier]
            ).get_eos_record(identifier)
            assert isclose(float(eos.pressure(v0)), 0.0, abs_tol=1e-10)
            pressure = float(eos.pressure(0.8 * v0))
            assert pressure > 0.0
            results[identifier] = {
                "V0_a3_cell": v0,
                "K0_gpa": k0,
                "K0_prime": k0_prime,
                "pressure_at_0.8_v0_gpa": pressure,
            }
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
