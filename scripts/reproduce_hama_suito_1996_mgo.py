"""Reproduce Hama and Suito's static-lattice MgO Vinet curve."""

from __future__ import annotations

import json

import numpy as np

BOHR_ANGSTROM = 0.529177210903
V0_SOURCE_BOHR3_PER_FORMULA = 123.747
Z = 4.0
V0 = V0_SOURCE_BOHR3_PER_FORMULA * BOHR_ANGSTROM**3 * Z
K0 = 157.0
K0_PRIME = 4.37


def vinet_pressure(volume):
    """Return pressure in GPa for conventional-cell volume in angstrom^3."""
    x = (np.asarray(volume, dtype=float) / V0) ** (1.0 / 3.0)
    return 3.0 * K0 * (1.0 - x) * np.exp(1.5 * (K0_PRIME - 1.0) * (1.0 - x)) / x**2


def reproduce():
    ratios = np.array([1.0, 0.9, 0.75, 0.5, 0.35])
    return {
        "v0_angstrom3_conventional_cell": V0,
        "pressures_gpa": {
            f"{ratio:.2f}": float(vinet_pressure(ratio * V0)) for ratio in ratios
        },
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
