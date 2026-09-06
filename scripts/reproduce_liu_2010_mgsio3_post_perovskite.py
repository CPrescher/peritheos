"""Reproduce Liu et al.'s static MgSiO3 post-perovskite BM3 curve."""

from __future__ import annotations

import json

import numpy as np

V0 = 163.3
K0 = 219.3
K0_PRIME = 4.4


def bm3_pressure(volume):
    eta = (V0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * K0 * (eta**7 - eta**5) * (1.0 + 0.75 * (K0_PRIME - 4.0) * (eta**2 - 1.0))
    )


def reproduce():
    ratios = np.array([1.0, 0.95, 0.9, 0.8, 0.7])
    return {
        "pressures_gpa": {
            f"{ratio:.2f}": float(bm3_pressure(ratio * V0)) for ratio in ratios
        }
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
