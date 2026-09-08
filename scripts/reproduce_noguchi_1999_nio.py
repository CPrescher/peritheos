#!/usr/bin/env python3
"""Reproduce the Noguchi et al. (1999) NiO shock-to-300 K reduction.

The Table 1 rows are shock Hugoniot states, not isothermal P-V observations.
This script first calculates the shock temperatures along the source-reported
linear ``Us-up`` Hugoniot and only then removes the Debye Mie-Gruneisen thermal
pressure before fitting BM3 to the reduced 300 K pressures.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import FitResult, fit_rt_eos
from peritheos.hugoniot import LinearUsUpHugoniot

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "nickel-oxide-noguchi-1999-table1-shock.csv"
)

REFERENCE_TEMPERATURE_K = 300.0
DEBYE_TEMPERATURE_K = 390.0
GRUNEISEN_0 = 1.38
GRUNEISEN_VOLUME_EXPONENT = 1.0
ATOMS_PER_FORMULA_UNIT = 2.0
MOLAR_MASS_G_MOL = 74.6928
BULK_DENSITY_G_CM3 = 6.781
SOURCE_C0_KM_S = 5.36
SOURCE_S = 1.19
SOURCE_CONFERENCE_K0_GPA = 184.0
SOURCE_CONFERENCE_K0_ERROR_GPA = 5.0
SOURCE_JOURNAL_K0_GPA = 191.0
SOURCE_JOURNAL_K0_PRIME = 3.9


@dataclass(frozen=True)
class ReducedStates:
    """Shock states and their calculated 300 K counterparts."""

    volume_ratio: NDArray[np.float64]
    hugoniot_pressure_gpa: NDArray[np.float64]
    shock_temperature_k: NDArray[np.float64]
    thermal_pressure_gpa: NDArray[np.float64]
    isothermal_pressure_gpa: NDArray[np.float64]


def _load_table1() -> dict[str, NDArray[np.float64]]:
    with DATASET.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {
        "volume_ratio": np.asarray(
            [float(row["final_normalized_volume"]) for row in rows]
        ),
        "pressure_gpa": np.asarray([float(row["final_pressure_gpa"]) for row in rows]),
        "particle_velocity_km_s": np.asarray(
            [float(row["final_particle_velocity_km_s"]) for row in rows]
        ),
        "shock_velocity_km_s": np.asarray(
            [float(row["final_shock_velocity_km_s"]) for row in rows]
        ),
    }


def _molar_reference_volume() -> float:
    """Return the initial molar volume in J bar^-1 mol^-1."""
    volume_cm3_mol = MOLAR_MASS_G_MOL / BULK_DENSITY_G_CM3
    return volume_cm3_mol / 10.0


def _models() -> tuple[LinearUsUpHugoniot, MieGruneisenDebye]:
    molar_volume = _molar_reference_volume()
    hugoniot = LinearUsUpHugoniot(
        V0=molar_volume,
        rho0=BULK_DENSITY_G_CM3,
        c0=SOURCE_C0_KM_S,
        s=SOURCE_S,
    )
    thermal = MieGruneisenDebye(
        BM3(
            molar_volume,
            SOURCE_JOURNAL_K0_GPA,
            SOURCE_JOURNAL_K0_PRIME,
        ),
        Tr=REFERENCE_TEMPERATURE_K,
        theta0=DEBYE_TEMPERATURE_K,
        gamma0=GRUNEISEN_0,
        q=GRUNEISEN_VOLUME_EXPONENT,
        n=ATOMS_PER_FORMULA_UNIT,
    )
    return hugoniot, thermal


def reduce_table1(*, rows: int = 8) -> ReducedStates:
    """Reduce the first *rows* final Hugoniot states to the 300 K isotherm.

    The temperature integration is the first-law relation along a Hugoniot,

    ``dT/dV = ((dE_H/dV + P_H) / C_V) - gamma*T/V``,

    with ``E_H-E_0 = (P_H+P_0)(V0-V)/2``.  Pressure-work units are converted
    by ``1 GPa * 1 J bar^-1 = 10^4 J``.  ``q=1`` is exactly the paper's
    ``gamma/V = constant`` convention; the integrated Gruneisen law therefore
    gives ``theta(V)=theta0*exp(gamma0*(1-V/V0))``.
    """
    if rows < 1 or rows > 8:
        raise ValueError("rows must be between one and eight")
    table = _load_table1()
    ratio = table["volume_ratio"][:rows]
    observed_pressure = table["pressure_gpa"][:rows]
    hugoniot, thermal = _models()
    volume = ratio * hugoniot.V0

    def temperature_derivative(
        current_volume: float, temperature: NDArray[np.float64]
    ) -> list[float]:
        pressure = float(hugoniot.pressure(current_volume))
        pressure_derivative = (
            -float(hugoniot.tangent_modulus(current_volume)) / current_volume
        )
        heat_capacity = float(
            thermal.molar_heat_capacity_v(current_volume, temperature[0])
        )
        gamma = float(thermal.gruneisen_parameter(current_volume))
        mechanical_heating = (
            0.5
            * ((hugoniot.V0 - current_volume) * pressure_derivative + pressure)
            * 1.0e4
            / heat_capacity
        )
        return [mechanical_heating - gamma * temperature[0] / current_volume]

    solution = solve_ivp(
        temperature_derivative,
        (hugoniot.V0, float(np.min(volume))),
        [REFERENCE_TEMPERATURE_K],
        dense_output=True,
        rtol=1.0e-10,
        atol=1.0e-10,
        max_step=hugoniot.V0 / 1000.0,
    )
    if not solution.success:
        raise RuntimeError(f"shock-temperature integration failed: {solution.message}")
    shock_temperature = np.asarray(solution.sol(volume)[0], dtype=float)
    thermal_pressure = np.asarray(
        thermal.thermal_pressure_increment(volume, shock_temperature), dtype=float
    )
    return ReducedStates(
        volume_ratio=ratio,
        hugoniot_pressure_gpa=observed_pressure,
        shock_temperature_k=shock_temperature,
        thermal_pressure_gpa=thermal_pressure,
        isothermal_pressure_gpa=observed_pressure - thermal_pressure,
    )


def fit_journal_isotherm() -> tuple[ReducedStates, FitResult]:
    """Fit BM3 to all eight thermally reduced journal Table 1 states."""
    states = reduce_table1(rows=8)
    molar_volume = _molar_reference_volume()
    fit = fit_rt_eos(
        BM3,
        states.volume_ratio * molar_volume,
        states.isothermal_pressure_gpa,
        initial={
            "K0": SOURCE_JOURNAL_K0_GPA,
            "K0_prime": SOURCE_JOURNAL_K0_PRIME,
        },
        fixed={"V0": molar_volume},
    )
    return states, fit


def fit_conference_crosscheck() -> tuple[ReducedStates, FitResult]:
    """Reproduce the seven-state conference result with K0' fixed to four."""
    states = reduce_table1(rows=7)
    molar_volume = _molar_reference_volume()
    fit = fit_rt_eos(
        BM3,
        states.volume_ratio * molar_volume,
        states.isothermal_pressure_gpa,
        initial={"K0": SOURCE_CONFERENCE_K0_GPA},
        fixed={"V0": molar_volume, "K0_prime": 4.0},
    )
    return states, fit


def main() -> None:
    states, journal_fit = fit_journal_isotherm()
    conference_states, conference_fit = fit_conference_crosscheck()

    print("Noguchi et al. NiO shock-to-300 K Debye Mie-Gruneisen reduction")
    print(
        "thermal parameters: Tr=300 K, theta0=390 K, gamma0=1.38, "
        "gamma/V=constant (q=1), n=2"
    )
    print("source Hugoniot: Us = 5.36 + 1.19 up (km/s)")
    print("eta       PH/GPa       TH/K    DeltaPth/GPa    P300/GPa")
    for index in range(states.volume_ratio.size):
        print(
            f"{states.volume_ratio[index]:.3f}  "
            f"{states.hugoniot_pressure_gpa[index]:10.3f}  "
            f"{states.shock_temperature_k[index]:9.3f}  "
            f"{states.thermal_pressure_gpa[index]:13.6f}  "
            f"{states.isothermal_pressure_gpa[index]:10.6f}"
        )
    print(
        "conference cross-check (first seven rows, K0'=4): "
        f"K0={conference_fit.parameters['K0']:.8f} GPa; "
        f"T(132.9 GPa)={conference_states.shock_temperature_k[-1]:.3f} K "
        "(source: 184+/-5 GPa and about 1600 K)"
    )
    print(
        "journal Table 1 reduced-state BM3 fit: "
        f"K0={journal_fit.parameters['K0']:.8f} GPa, "
        f"K0'={journal_fit.parameters['K0_prime']:.8f}; "
        f"RMSE={np.sqrt(np.mean(journal_fit.residuals**2)):.8f} GPa "
        "(source: 191 GPa, 3.9)"
    )


if __name__ == "__main__":
    main()
