"""Reproduce the three composition-specific Mao et al. (1991) BM2 curves."""

from __future__ import annotations

import json

import numpy as np

PUBLISHED = {
    "MgSiO3": 162.49,
    "Mg0.9Fe0.1SiO3": 162.79,
    "Mg0.8Fe0.2SiO3": 163.53,
}


def bm3_pressure(volume, v0, k0=261.0, kp=4.0):
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def reproduce() -> dict[str, object]:
    fractions = np.array([1.0, 0.95, 0.90])
    curves = {
        formula: dict(
            zip(
                map(str, fractions),
                map(float, bm3_pressure(fractions * v0, v0)),
            )
        )
        for formula, v0 in PUBLISHED.items()
    }
    return {"shared_K0_gpa": 261.0, "fixed_K0_prime": 4.0, "curves_gpa": curves}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
