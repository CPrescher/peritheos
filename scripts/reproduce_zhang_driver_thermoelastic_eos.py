#!/usr/bin/env python3
"""Reproduce executable checkpoints for Zhang (2022) and Driver (2010)."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

RECORDS = {
    "bridgmanite": {
        "bridgmanite_zhang_wentzcovitch_2022_phq_lda_300k_bm3": (160.60, 100.0),
        "bridgmanite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3": (169.76, 100.0),
    },
    "mgsio3_post_perovskite": {
        "mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_lda_300k_bm3": (
            160.08,
            100.0,
        ),
        "mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3": (
            170.24,
            100.0,
        ),
    },
    "alpha_quartz": {
        "alpha_quartz_driver_2010_qmc_300k_vinet": (112.91675014, 20.0),
    },
    "sio2_stv_andr": {
        "sio2_stv_andr_driver_2010_qmc_300k_vinet": (47.12273825, 40.0),
    },
    "seifertite": {
        "seifertite_driver_2010_qmc_300k_vinet": (91.75597334, 150.0),
    },
}


def reproduce() -> dict[str, object]:
    results: dict[str, object] = {}
    for material_identifier, records in RECORDS.items():
        document = get_material_document(material_identifier)
        for record_identifier, (expected_v0, checkpoint_pressure) in records.items():
            record = next(
                row
                for row in document["eos_records"]
                if row["identifier"] == record_identifier
            )
            assert isclose(record["eos"]["parameters"]["V0"], expected_v0)
            eos = Material.from_eosmat(
                document, record_identifiers=[record_identifier]
            ).get_eos_record(record_identifier)
            volume = float(eos.volume(checkpoint_pressure))
            roundtrip = float(eos.pressure(volume))
            assert isclose(float(eos.pressure(expected_v0)), 0.0, abs_tol=1e-10)
            assert isclose(roundtrip, checkpoint_pressure, rel_tol=2e-10)
            results[record_identifier] = {
                "V0_a3_cell": expected_v0,
                "checkpoint_pressure_gpa": checkpoint_pressure,
                "checkpoint_volume_a3_cell": volume,
                "roundtrip_pressure_gpa": roundtrip,
            }
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
