#!/usr/bin/env python3
"""Verify the complete Holland et al. (2013) modified-Tait endmembers."""

from __future__ import annotations

import json
from math import isclose
from pathlib import Path

from peritheos.eos.rt import ModifiedTait

ANGSTROM3_PER_CM3_MOL = 1.0e24 / 6.02214076e23
ROOT = Path(__file__).parents[1]

ROWS = {
    "mpv": (
        "bridgmanite",
        "bridgmanite_holland_2013_mpv_modified_tait",
        2.445,
        4,
        2510.0,
        4.14,
        -0.00160,
    ),
    "fpv": (
        "fesio3_bridgmanite",
        "fesio3_bridgmanite_holland_2013_fpv_modified_tait",
        2.548,
        4,
        2810.0,
        4.14,
        -0.00160,
    ),
    "apv": (
        "al2o3_perovskite",
        "al2o3_perovskite_holland_2013_apv_modified_tait",
        2.540,
        4,
        2030.0,
        4.00,
        -0.00200,
    ),
    "cpv": (
        "ca_perovskite",
        "ca_perovskite_holland_2013_cpv_modified_tait",
        2.745,
        1,
        2360.0,
        3.90,
        -0.00160,
    ),
    "per": (
        "mgo",
        "mgo_holland_2013_per_modified_tait",
        1.125,
        4,
        1616.0,
        3.95,
        -0.00240,
    ),
    "fper": (
        "feo",
        "feo_holland_2013_fper_modified_tait",
        1.206,
        4,
        1520.0,
        4.90,
        -0.00320,
    ),
    "stv": (
        "sio2_stv_andr",
        "sio2_stv_andr_holland_2013_stv_modified_tait",
        1.401,
        2,
        3090.0,
        4.60,
        -0.00150,
    ),
}


def reproduce() -> dict[str, object]:
    results: dict[str, object] = {}
    for abbreviation, row in ROWS.items():
        (
            material_identifier,
            record_identifier,
            v0_j_bar_mol,
            z,
            k0_kbar,
            kp,
            kpp_kbar_inv,
        ) = row
        expected = {
            "V0": v0_j_bar_mol * 10.0 * ANGSTROM3_PER_CM3_MOL * z,
            "K0": k0_kbar / 10.0,
            "K0_prime": kp,
            "K0_double_prime": kpp_kbar_inv * 10.0,
        }
        document = json.loads(
            (
                ROOT / "peritheos/data/materials" / f"{material_identifier}.eosmat"
            ).read_text(encoding="utf-8")
        )
        stored = next(
            item
            for item in document["eos_records"]
            if item["identifier"] == record_identifier
        )
        for name, value in expected.items():
            assert isclose(stored["eos"]["parameters"][name], value, rel_tol=2e-10)
        eos = ModifiedTait(**stored["eos"]["parameters"])
        checkpoint_volume = float(eos.volume(30.0))
        roundtrip = float(eos.pressure(checkpoint_volume))
        assert isclose(roundtrip, 30.0, rel_tol=2e-10)
        results[abbreviation] = {
            "record": record_identifier,
            "source_parameters": {
                "V0_j_bar_mol": v0_j_bar_mol,
                "K0_kbar": k0_kbar,
                "K0_prime": kp,
                "K0_double_prime_kbar_inverse": kpp_kbar_inv,
            },
            "stored_parameters": expected,
            "volume_at_30_gpa_a3_cell": checkpoint_volume,
            "roundtrip_pressure_gpa": roundtrip,
        }
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
