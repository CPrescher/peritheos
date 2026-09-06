"""Reproduce Chen et al.'s tetragonal CaSiO3 second-order BM curve."""

from __future__ import annotations

import json

import numpy as np

V0_PER_FORMULA = 46.3
Z = 4.0
V0 = V0_PER_FORMULA * Z
K0 = 223.0


def bm2_pressure(volume):
    eta = (V0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return 1.5 * K0 * (eta**7 - eta**5)


def reproduce():
    ratios = np.array([1.0, 0.95, 0.9, 0.85, 0.8])
    return {
        "v0_conventional_cell_a3": V0,
        "pressures_gpa": {
            f"{ratio:.2f}": float(bm2_pressure(ratio * V0)) for ratio in ratios
        },
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
