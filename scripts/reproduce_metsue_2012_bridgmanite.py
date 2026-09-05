"""Audit Metsue and Tsuchiya's (2012) seven static BM3 parameterizations."""

from __future__ import annotations

import json

import numpy as np

AVOGADRO = 6.02214076e23
MOLAR_TO_Z4_A3 = 4.0e24 / AVOGADRO

PUBLISHED = {
    "HS_model_1": (24.3325, 259.3412),
    "HS_model_2": (24.3326, 259.2916),
    "HS_model_3": (24.3317, 259.4608),
    "LS_model_1": (24.2583, 260.7987),
    "LS_model_2": (24.2605, 260.7980),
    "LS_model_3": (24.2596, 260.7441),
    "pure_MgSiO3": (24.2717, 258.7857),
}


def bm3_pressure(volume, v0, k0, kp=3.94):
    """Evaluate the standard Eulerian third-order Birch-Murnaghan EOS."""
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def reproduce() -> dict[str, object]:
    """Return conversion, configuration-spread, and composition checks."""
    converted = {
        name: {
            "V0_cm3_mol": values[0],
            "V0_a3_z4": values[0] * MOLAR_TO_Z4_A3,
            "K0_gpa": values[1],
            "K0_prime_fixed": 3.94,
            "pressure_at_v0_gpa": float(bm3_pressure(values[0], values[0], values[1])),
        }
        for name, values in PUBLISHED.items()
    }
    groups = {}
    for spin in ("HS", "LS"):
        values = [PUBLISHED[f"{spin}_model_{index}"] for index in (1, 2, 3)]
        volumes = np.array([value[0] for value in values])
        moduli = np.array([value[1] for value in values])
        groups[spin] = {
            "maximum_relative_volume_spread_percent": float(
                100.0 * np.ptp(volumes) / np.mean(volumes)
            ),
            "maximum_relative_modulus_spread_percent": float(
                100.0 * np.ptp(moduli) / np.mean(moduli)
            ),
        }
    pure_v, pure_k = PUBLISHED["pure_MgSiO3"]
    trends = {}
    for spin in ("HS", "LS"):
        volume, modulus = PUBLISHED[f"{spin}_model_1"]
        trends[spin] = {
            "dlnV0_dXFe": (volume / pure_v - 1.0) / 0.0625,
            "dlnK0_dXFe": (modulus / pure_k - 1.0) / 0.0625,
        }
    return {
        "molar_to_z4_a3": MOLAR_TO_Z4_A3,
        "records": converted,
        "configuration_spread": groups,
        "model1_composition_trends": trends,
    }


def main() -> None:
    """Print the audit as stable JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
