"""Reproduce Cohen and Lin's (2014) three static FeSiO3 Vinet fits.

Table III gives both the zero-pressure Vinet coefficients and independently
tabulated volume, bulk modulus, and pressure derivative at 100 GPa.  The
energy-volume grid in Figure 3 is not tabulated, so those high-pressure values
are used as deterministic convention checks without digitizing false precision.
"""

from __future__ import annotations

import json

import numpy as np
from scipy.optimize import brentq

PUBLISHED = {
    "Pv": {
        "V0": 44.31,
        "K0": 225.0,
        "K0_prime": 4.42,
        "V100": 34.27,
        "K100": 597.0,
        "K100_prime": 3.34,
    },
    "PPv": {
        "V0": 44.90,
        "K0": 189.0,
        "K0_prime": 4.73,
        "V100": 33.98,
        "K100": 579.0,
        "K100_prime": 3.47,
    },
    "PPv-II": {
        "V0": 45.45,
        "K0": 195.0,
        "K0_prime": 4.67,
        "V100": 34.49,
        "K100": 580.0,
        "K100_prime": 3.44,
    },
}


def vinet_pressure(volume: float | np.ndarray, v0: float, k0: float, kp: float):
    """Evaluate the standard Vinet pressure independently of Peritheos."""
    x = (np.asarray(volume, dtype=float) / v0) ** (1.0 / 3.0)
    eta = 1.5 * (kp - 1.0)
    return 3.0 * k0 * (1.0 - x) / x**2 * np.exp(eta * (1.0 - x))


def volume_at_pressure(pressure: float, v0: float, k0: float, kp: float) -> float:
    """Invert the Vinet pressure over the source's static compression range."""
    return float(
        brentq(
            lambda volume: vinet_pressure(volume, v0, k0, kp) - pressure,
            0.4 * v0,
            1.3 * v0,
        )
    )


def modulus_and_derivative(
    volume: float, v0: float, k0: float, kp: float
) -> tuple[float, float]:
    """Return K and dK/dP by symmetric numerical differentiation."""
    h = volume * 1.0e-5

    def bulk_modulus(at_volume: float) -> float:
        dh = at_volume * 1.0e-5
        derivative = (
            vinet_pressure(at_volume + dh, v0, k0, kp)
            - vinet_pressure(at_volume - dh, v0, k0, kp)
        ) / (2.0 * dh)
        return float(-at_volume * derivative)

    modulus = bulk_modulus(volume)
    derivative = (bulk_modulus(volume + h) - bulk_modulus(volume - h)) / (
        vinet_pressure(volume + h, v0, k0, kp) - vinet_pressure(volume - h, v0, k0, kp)
    )
    return modulus, float(derivative)


def reproduce() -> dict[str, object]:
    """Return stable 100 GPa checkpoint diagnostics for all three phases."""
    phases = {}
    for name, values in PUBLISHED.items():
        volume = volume_at_pressure(
            100.0, values["V0"], values["K0"], values["K0_prime"]
        )
        modulus, derivative = modulus_and_derivative(
            volume, values["V0"], values["K0"], values["K0_prime"]
        )
        phases[name] = {
            "calculated": {
                "V100_a3_per_formula": volume,
                "K100_gpa": modulus,
                "K100_prime": derivative,
            },
            "published": {
                "V100_a3_per_formula": values["V100"],
                "K100_gpa": values["K100"],
                "K100_prime": values["K100_prime"],
            },
            "absolute_differences": {
                "V100_a3_per_formula": abs(volume - values["V100"]),
                "K100_gpa": abs(modulus - values["K100"]),
                "K100_prime": abs(derivative - values["K100_prime"]),
            },
        }
    return {"pressure_gpa": 100.0, "phases": phases}


def main() -> None:
    """Print the reproduction report as JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
