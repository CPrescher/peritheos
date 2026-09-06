#!/usr/bin/env python3
"""Reproduce the two Karki et al. (1997) static MgO Birch fits."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

V0 = 4.2506**3
RECORDS = {
    "mgo_karki_1997_lda_static_bm4": {
        "type": "BM4",
        "parameters": {
            "V0": V0,
            "K0": 159.7,
            "K0_prime": 4.26,
            "K0_double_prime": -0.026,
        },
    },
    "mgo_karki_1997_lda_static_bm3": {
        "type": "BM3",
        "parameters": {"V0": V0, "K0": 159.4, "K0_prime": 4.28},
    },
}


def _bm3_pressure(volume: float, *, v0: float, k0: float, kp: float) -> float:
    eta = v0 / volume
    return (
        1.5
        * k0
        * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
        * (1.0 + 0.75 * (kp - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
    )


def _bm4_pressure(
    volume: float, *, v0: float, k0: float, kp: float, kpp: float
) -> float:
    # Karki et al. Equation (11), written with Eulerian strain f.
    f = 0.5 * ((v0 / volume) ** (2.0 / 3.0) - 1.0)
    a1 = 1.5 * (kp - 4.0)
    a2 = 1.5 * (k0 * kpp + (kp - 4.0) * (kp - 3.0) + 35.0 / 9.0)
    return 3.0 * k0 * f * (1.0 + 2.0 * f) ** 2.5 * (1.0 + a1 * f + a2 * f**2)


def reproduce() -> dict[str, object]:
    document = get_material_document("mgo")
    results: dict[str, object] = {}
    for identifier, expected in RECORDS.items():
        record = next(
            row for row in document["eos_records"] if row["identifier"] == identifier
        )
        assert record["eos"]["type"] == expected["type"]
        assert record["eos"]["parameters"] == expected["parameters"]
        eos = Material.from_eosmat(
            document, record_identifiers=[identifier]
        ).get_eos_record(identifier)

        residuals = []
        anchors = []
        for volume_ratio in (1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7):
            volume = V0 * volume_ratio
            if expected["type"] == "BM4":
                pressure_reference = _bm4_pressure(
                    volume,
                    v0=V0,
                    k0=159.7,
                    kp=4.26,
                    kpp=-0.026,
                )
            else:
                pressure_reference = _bm3_pressure(volume, v0=V0, k0=159.4, kp=4.28)
            pressure_peritheos = float(eos.pressure(volume))
            residuals.append(pressure_peritheos - pressure_reference)
            assert isclose(pressure_peritheos, pressure_reference, abs_tol=2e-11)
            assert isclose(float(eos.volume(pressure_peritheos)), volume, rel_tol=2e-10)
            anchors.append(
                {
                    "volume_ratio": volume_ratio,
                    "pressure_gpa": pressure_peritheos,
                }
            )
        results[identifier] = {
            "anchors": anchors,
            "max_abs_equation_residual_gpa": max(abs(value) for value in residuals),
        }
    return {"reference_volume_a3_cell": V0, "records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
