"""Audit Cohen and Lin's (2014) three static FeSiO3 Vinet records.

The final article's Table III gives both the zero-pressure Vinet coefficients
and independently tabulated volume, bulk modulus, and pressure derivative at
100 GPa.  No inspected primary artifact tabulates the underlying E(V) rows.
This module therefore keeps exact table replay separate from a deliberately
approximate, volume-axis-only reading of the raster Figure 3 markers.
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

# The original arXiv v1 author PDF prints one extra decimal for K0 and K100.
# These superseded values are useful independent rounding diagnostics, but the
# production records correctly retain the final version-of-record Table III.
ARXIV_V1 = {
    "Pv": {
        "V0": 44.31,
        "K0": 224.8,
        "K0_prime": 4.42,
        "V100": 34.27,
        "K100": 596.5,
        "K100_prime": 3.34,
    },
    "PPv": {
        "V0": 44.90,
        "K0": 189.7,
        "K0_prime": 4.73,
        "V100": 33.98,
        "K100": 579.0,
        "K100_prime": 3.47,
    },
    "PPv-II": {
        "V0": 45.45,
        "K0": 194.7,
        "K0_prime": 4.67,
        "V100": 34.49,
        "K100": 580.1,
        "K100_prime": 3.44,
    },
}

PRESSURE_GRID_GPA = (-10.0, 0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0)

# Approximate marker centers read only from the x axis of the 1606 x 978 raster
# embedded in the arXiv v2 Figures_03.eps.  The plot frame maps x=231..1589 to
# V=30..50 A^3/FeSiO3.  These are graphical diagnostics, not primary rows and
# must never be passed to the primary-refit campaign as observations.
FIGURE_3_APPROXIMATE_VOLUMES = {
    "Pv": (46.3878, 44.1924, 40.3710, 37.8056, 35.8126, 34.2267, 32.9153, 31.7904),
    "PPv": (47.7393, 44.7695, 40.4073, 37.6030, 35.5508, 33.9255, 32.5843, 31.4585),
    "PPv-II": (47.9910, 45.2944, 40.9946, 38.1788, 36.1063, 34.4562, 33.0890, 31.9350),
}

TABLE_I_100_GPA_LATTICES = {
    "PPv": (2.508, 8.614, 6.283),
    "PPv-II": (10.082, 5.478, 2.495),
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
    """Return exact table checks plus explicitly approximate plot diagnostics."""
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
    arxiv_v1 = {}
    for name, values in ARXIV_V1.items():
        volume = volume_at_pressure(
            100.0, values["V0"], values["K0"], values["K0_prime"]
        )
        modulus, derivative = modulus_and_derivative(
            volume, values["V0"], values["K0"], values["K0_prime"]
        )
        arxiv_v1[name] = {
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
        }

    figure_diagnostic = {}
    for name, approximate_volumes in FIGURE_3_APPROXIMATE_VOLUMES.items():
        values = PUBLISHED[name]
        model_volumes = [
            volume_at_pressure(pressure, values["V0"], values["K0"], values["K0_prime"])
            for pressure in PRESSURE_GRID_GPA
        ]
        differences = np.asarray(approximate_volumes) - np.asarray(model_volumes)
        index_100 = PRESSURE_GRID_GPA.index(100.0)
        figure_diagnostic[name] = {
            "approximate_volumes_a3_per_formula": list(approximate_volumes),
            "published_curve_volumes_a3_per_formula": model_volumes,
            "maximum_absolute_volume_difference_a3_per_formula": float(
                np.max(np.abs(differences))
            ),
            "at_100_gpa": {
                "figure_3_approximate_volume_a3_per_formula": approximate_volumes[
                    index_100
                ],
                "table_iii_volume_a3_per_formula": values["V100"],
                "absolute_difference_a3_per_formula": abs(
                    approximate_volumes[index_100] - values["V100"]
                ),
            },
        }

    structure_check = {
        name: {
            "table_i_lattice_a3_per_formula": float(np.prod(lattice) / 4.0),
            "table_iii_volume_a3_per_formula": PUBLISHED[name]["V100"],
            "absolute_difference_a3_per_formula": abs(
                float(np.prod(lattice) / 4.0) - PUBLISHED[name]["V100"]
            ),
        }
        for name, lattice in TABLE_I_100_GPA_LATTICES.items()
    }

    return {
        "pressure_gpa": 100.0,
        "phases": phases,
        "arxiv_v1_rounding_diagnostic": arxiv_v1,
        "figure_3_approximate_diagnostic": {
            "classification": "approximate_plot_diagnostic_not_primary_observations",
            "pressure_grid_gpa": list(PRESSURE_GRID_GPA),
            "coordinate_uncertainty_a3_per_formula": 0.03,
            "phases": figure_diagnostic,
        },
        "table_i_structure_diagnostic": structure_check,
    }


def main() -> None:
    """Print the reproduction report as JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
