"""Reproduce the executable EOS checks from the final Ca/other source batch."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from peritheos import Material

ROOT = Path(__file__).parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
DATA = ROOT / "peritheos/data/datasets"


def _sun_2019_thermal_pressure(volume_cm3_mol: float, temperature_k: float) -> float:
    """Evaluate Sun et al. (2019) Equations (4) and (10)-(13)."""
    vx = 40.72
    cv = 204.347 + (-830.0 / (20.291 * math.sqrt(math.pi / 2.0))) * math.exp(
        -2.0 * ((volume_cm3_mol - 29.110) / 20.291) ** 2
    )
    reduced_volume = volume_cm3_mol / vx - 1.0
    gamma = 0.311 - 0.933 * reduced_volume + 1.144 * reduced_volume**2
    return gamma * cv * (temperature_k - 2500.0) / volume_cm3_mol / 1000.0


def _record(material_file: str, identifier: str):
    document = json.loads((MATERIALS / material_file).read_text(encoding="utf-8"))
    material = Material.from_eosmat(document, record_identifiers=[identifier])
    return material.get_eos_record(identifier)


def reproduce() -> dict[str, object]:
    """Return source checkpoints and exact equation-composition diagnostics."""
    karki = _record("ca_perovskite.eosmat", "ca_perovskite_karki_crain_1998_static_bm3")
    shim = _record(
        "casio3_perovskite_tetragonal.eosmat",
        "casio3_perovskite_tetragonal_shim_2002_bm3_1",
    )
    ono = _record("ca_perovskite.eosmat", "ca_perovskite_ono_2013_bm3_log_thermal")
    sun = _record("fesio3_liquid.eosmat", "fesio3_liquid_sun_2019_2500k_bm4_1")

    volume = 40.72 * (1.0e24 / 6.02214076e23)
    sun_pressure = float(sun.pressure(volume))
    with (DATA / "fesio3-liquid-sun-2019-table1-pvt.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = list(csv.DictReader(stream))
    liquid_rows = [row for row in rows if row["used_in_published_eos_fit"] == "yes"]
    nonliquid_rows = [row for row in rows if row["phase_state"] == "nonliquid"]
    sun_residuals = []
    for row in liquid_rows:
        calculated = float(sun.pressure(float(row["volume_a3_per_formula_unit"])))
        calculated += _sun_2019_thermal_pressure(
            float(row["volume_cm3_mol"]), float(row["temperature_k"])
        )
        sun_residuals.append(calculated - float(row["pressure_gpa"]))

    ono_volume = 40.0
    ono_temperature = 2000.0
    ono_cold = float(ono.eos.rt_eos.pressure(ono_volume))
    ono_expected_thermal = (0.0083 - 0.0031 * math.log(45.58 / ono_volume)) * (
        ono_temperature - 300.0
    )

    return {
        "karki_140_gpa_roundtrip_volume_a3": float(karki.volume(140.0)),
        "shim_45_8_gpa_roundtrip_volume_a3": float(shim.volume(45.8)),
        "ono_2000k_checkpoint": {
            "volume_a3": ono_volume,
            "cold_pressure_gpa": ono_cold,
            "thermal_pressure_gpa": float(
                ono.eos.thermal_pressure(ono_volume, ono_temperature)
            ),
            "expected_thermal_pressure_gpa": ono_expected_thermal,
            "total_pressure_gpa": float(ono.pressure(ono_volume, ono_temperature)),
        },
        "sun_2500k_table1_checkpoint": {
            "volume_cm3_mol": 40.72,
            "published_pressure_gpa": 1.07,
            "published_pressure_standard_error_gpa": 0.12,
            "calculated_pressure_gpa": sun_pressure,
            "absolute_difference_gpa": abs(sun_pressure - 1.07),
        },
        "sun_table1_states": len(rows),
        "sun_state_classification": {
            "liquid_states_used_in_published_eos_fit": len(liquid_rows),
            "nonliquid_states_excluded": len(nonliquid_rows),
            "nonliquid_volume_ratio_temperature_pairs": [
                [float(row["volume_ratio_to_vx"]), int(row["temperature_k"])]
                for row in nonliquid_rows
            ],
        },
        "sun_published_rounded_model_liquid_state_rmse_gpa": math.sqrt(
            sum(residual**2 for residual in sun_residuals) / len(sun_residuals)
        ),
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
