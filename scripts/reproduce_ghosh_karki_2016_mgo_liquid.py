"""Reproduce the Ghosh and Karki (2016) 3000 K liquid-MgO BM3."""

from __future__ import annotations

import json

import numpy as np

PUBLISHED = {"V0": 30.0, "K0": 16.6, "K0_prime": 6.14}


def bm3_pressure(volume, v0, k0, k0_prime):
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def reproduce():
    volumes = np.asarray([30.0, 27.0, 25.0, 22.0, 20.0, 18.0])
    pressures = bm3_pressure(volumes, *PUBLISHED.values())
    return {
        "temperature_k": 3000.0,
        "volumes_a3_per_formula": volumes.tolist(),
        "pressures_gpa": pressures.tolist(),
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
