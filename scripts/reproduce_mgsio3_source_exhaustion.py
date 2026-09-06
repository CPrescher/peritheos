"""Reproduce exact equation mappings from the MgSiO3 source-exhaustion audit."""

from __future__ import annotations

import json
import math
from pathlib import Path

from peritheos import Material

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"


def _record(material_file: str, identifier: str):
    document = json.loads((MATERIALS / material_file).read_text(encoding="utf-8"))
    material = Material.from_eosmat(document, record_identifiers=[identifier])
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )
    return material.get_eos_record(identifier), source


def _bm3(volume: float, v0: float, k0: float, k0_prime: float) -> float:
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def _dorogokupets_reference(
    volume: float, v0: float, k0: float, k0_prime: float, k: float = 5.0
) -> float:
    x = (volume / v0) ** (1.0 / 3.0)
    eta = 1.5 * k0_prime - k + 0.5
    return 3.0 * k0 * x ** (-k) * (1.0 - x) * math.exp(eta * (1.0 - x))


def reproduce() -> dict[str, object]:
    """Return direct-formula residuals and density-volume conversions."""
    specifications = [
        (
            "bridgmanite.eosmat",
            "bridgmanite_dorogokupets_2015_298k_rydberg_stacey",
            "rydberg",
        ),
        (
            "akimotoite.eosmat",
            "akimotoite_dorogokupets_2015_298k_rydberg_stacey",
            "rydberg",
        ),
        (
            "mgsio3_post_perovskite.eosmat",
            "mgsio3_post_perovskite_dorogokupets_2015_298k_rydberg_stacey",
            "rydberg",
        ),
        (
            "bridgmanite.eosmat",
            "bridgmanite_hamahata_2000_md_300k_bm3",
            "bm3",
        ),
        (
            "mgsio3_liquid.eosmat",
            "mgsio3_liquid_akins_2004_adiabatic_bm3",
            "bm3",
        ),
    ]
    residuals: dict[str, float] = {}
    for material_file, identifier, equation in specifications:
        record, source = _record(material_file, identifier)
        parameters = source["eos"]["parameters"]
        v0 = float(parameters["V0"])
        maximum = 0.0
        for ratio in (0.98, 0.9, 0.8, 0.7):
            volume = ratio * v0
            actual = float(record.pressure(volume))
            if equation == "rydberg":
                expected = _dorogokupets_reference(
                    volume,
                    v0,
                    float(parameters["K0"]),
                    float(parameters["K0_prime"]),
                )
            else:
                expected = _bm3(
                    volume,
                    v0,
                    float(parameters["K0"]),
                    float(parameters["K0_prime"]),
                )
            maximum = max(maximum, abs(actual - expected))
        residuals[identifier] = maximum

    avogadro = 6.02214076e23
    return {
        "maximum_absolute_pressure_residual_gpa": residuals,
        "hamahata_v0_a3_z4": 4.0 * 100.387 / (4.030 * avogadro) * 1.0e24,
        "akins_v0_a3_per_formula_unit": 100.387 / (3.68 * avogadro) * 1.0e24,
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
