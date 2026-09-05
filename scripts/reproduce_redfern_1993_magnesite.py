"""Reproduce the two Redfern et al. (1993) magnesite parameterizations."""

from __future__ import annotations

import json

import numpy as np

PUBLISHED = {
    "fixed_kp": {"V0": 279.4, "K0": 142.0, "K0_prime": 4.0},
    "free_kp": {"V0": 279.4, "K0": 151.0, "K0_prime": 2.5},
}


def bm3_pressure(volume, v0, k0, kp):
    """Evaluate the standard Eulerian finite-strain BM3 pressure."""
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def reproduce() -> dict[str, object]:
    """Return stable pressure checkpoints without inventing source P-V rows."""
    fractions = np.array([1.0, 0.98, 0.95, 0.90])
    curves = {}
    for name, pars in PUBLISHED.items():
        pressure = bm3_pressure(
            fractions * pars["V0"], pars["V0"], pars["K0"], pars["K0_prime"]
        )
        curves[name] = dict(zip(map(str, fractions), map(float, pressure)))
    return {"volume_fractions": fractions.tolist(), "curves_gpa": curves}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
