"""Reproduce Zhang and Bukowinski's three source oxide EOS curves."""

from __future__ import annotations

import json

import numpy as np

PUBLISHED = {
    "MgO_B1": (74.60, 180.0, 4.04),
    "MgO_B2": (18.07, 175.0, 4.20),
    "SiO2_stishovite": (45.54, 378.0, 6.4),
}


def bm3_pressure(volume, v0, k0, kp):
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1 / 3)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def reproduce():
    return {
        "pressure_at_0.9_v0_gpa": {
            name: float(bm3_pressure(0.9 * p[0], *p)) for name, p in PUBLISHED.items()
        }
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
