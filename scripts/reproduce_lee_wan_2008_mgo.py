#!/usr/bin/env python3
"""Reproduce the Lee and Wan (2008) MgO BM3 conversion and curve checks."""

from __future__ import annotations

import json

import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903
SUPER_CELL_FORMULA_UNITS = 32
CONVENTIONAL_CELL_FORMULA_UNITS = 4

SOURCE = {
    "LDA": {"V0_bohr3": 3782.30, "K0_gpa": 177.486, "K0_prime": 4.026},
    "GGA": {"V0_bohr3": 4046.43, "K0_gpa": 149.320, "K0_prime": 4.080},
}

# The five abscissae visible in Figures 1, 3, and 4. The source does not
# publish its raw simulation grid, so these values are used only for an
# internal curve/polynomial comparison, not as digitized fit observations.
FIGURE_CHECK_VOLUMES_BOHR3 = np.array([2740.0, 2980.0, 3240.0, 3500.0, 3800.0])


def bm3_pressure(volume, v0: float, k0: float, k0_prime: float):
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def conventional_z4_volume(volume_bohr3: float) -> float:
    return (
        volume_bohr3
        * BOHR_TO_ANGSTROM**3
        * CONVENTIONAL_CELL_FORMULA_UNITS
        / SUPER_CELL_FORMULA_UNITS
    )


def reproduce() -> dict[str, object]:
    lda = SOURCE["LDA"]
    gga = SOURCE["GGA"]
    lda_pressure = bm3_pressure(
        FIGURE_CHECK_VOLUMES_BOHR3,
        lda["V0_bohr3"],
        lda["K0_gpa"],
        lda["K0_prime"],
    )
    gga_pressure = bm3_pressure(
        FIGURE_CHECK_VOLUMES_BOHR3,
        gga["V0_bohr3"],
        gga["K0_gpa"],
        gga["K0_prime"],
    )
    pressure_difference = gga_pressure - lda_pressure

    # Figure 4 reports this visualization fit explicitly:
    # delta P = a1/V + a2/V^2, a1=1.53e4 GPa bohr^3,
    # a2=1.06e8 GPa bohr^6.
    polynomial = (
        1.53e4 / FIGURE_CHECK_VOLUMES_BOHR3 + 1.06e8 / FIGURE_CHECK_VOLUMES_BOHR3**2
    )
    return {
        "converted_parameters_z4": {
            method: {
                "V0_a3": conventional_z4_volume(values["V0_bohr3"]),
                "K0_gpa": values["K0_gpa"],
                "K0_prime": values["K0_prime"],
            }
            for method, values in SOURCE.items()
        },
        "figure_curve_check": {
            "volumes_bohr3": FIGURE_CHECK_VOLUMES_BOHR3.tolist(),
            "lda_pressure_gpa": lda_pressure.tolist(),
            "gga_pressure_gpa": gga_pressure.tolist(),
            "gga_minus_lda_gpa": pressure_difference.tolist(),
            "reported_polynomial_gpa": polynomial.tolist(),
            "maximum_absolute_difference_from_polynomial_gpa": float(
                np.max(np.abs(pressure_difference - polynomial))
            ),
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
