#!/usr/bin/env python3
"""Audit Mosenfelder et al. (2009) MgSiO3 post-perovskite.

The paper fits the Guignot et al. (2007) static P-V-T observations together
with six Table 2 shock states.  The shock temperatures are not required for
the Hugoniot pressure-density reduction: Rankine-Hugoniot energy closes the
model, and temperature is an additional observable only for the three Luo
et al. pyrometry rows.

Sources
-------
Mosenfelder et al. (2009), doi:10.1029/2008JB005900, equations 2-12,
Tables 2-4; Mosenfelder et al. (2007), doi:10.1029/2006JB004364,
equations 5-11; Guignot et al. (2007), doi:10.1016/j.epsl.2007.01.025,
Table 1.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro
from scipy.integrate import quad
from scipy.optimize import least_squares, root_scalar

ROOT = Path(__file__).resolve().parents[1]
STATIC_DATA = ROOT / (
    "peritheos/data/datasets/mgsio3-post-perovskite-guignot-2007-table1-pvt.csv"
)
SHOCK_DATA = ROOT / (
    "peritheos/data/datasets/mgsio3-post-perovskite-mosenfelder-2009-table2-shock.csv"
)
DERIVED_DATA = ROOT / (
    "peritheos/data/datasets/"
    "mgsio3-post-perovskite-mosenfelder-2009-derived-reference-isentrope.csv"
)

MOLAR_MASS_G_MOL = 100.3875
FORMULA_UNITS_PER_CELL = 4.0
TR_K = 300.0
V0_A3 = 162.2
PUBLISHED = np.array([225.0, 4.21, 2.61, 2.1, 1.035, 990.0])
PARAMETER_NAMES = ("K0", "K0_prime", "gamma0", "q", "Cvm", "theta0")


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], name: str) -> float:
    value = row[name]
    return float(value) if value else np.nan


def cell_volume_from_density(density_g_cm3: np.ndarray | float) -> np.ndarray:
    """Convert MgSiO3 density to the conventional Cmcm Z=4 cell volume."""
    return (
        FORMULA_UNITS_PER_CELL
        * MOLAR_MASS_G_MOL
        / (Avogadro * np.asarray(density_g_cm3, dtype=float))
        * 1.0e24
    )


def density_from_cell_volume(volume_a3: np.ndarray | float) -> np.ndarray:
    return (
        FORMULA_UNITS_PER_CELL
        * MOLAR_MASS_G_MOL
        / (Avogadro * np.asarray(volume_a3, dtype=float))
        * 1.0e24
    )


def debye3(argument: np.ndarray | float) -> np.ndarray | float:
    """Evaluate the third Debye function used by the source."""
    values = np.asarray(argument, dtype=float)

    def scalar(value: float) -> float:
        if value < 1.0e-3:
            return 1.0 - 3.0 * value / 8.0 + value**2 / 20.0
        integral = quad(
            lambda x: x**3 * np.exp(-x) / (1.0 - np.exp(-x)),
            0.0,
            value,
            epsabs=1.0e-11,
            epsrel=1.0e-11,
        )[0]
        return 3.0 * integral / value**3

    result = np.array([scalar(value) for value in values.flat]).reshape(values.shape)
    return float(result) if result.ndim == 0 else result


def gamma(volume_a3: np.ndarray | float, gamma0: float, q: float):
    return gamma0 * (np.asarray(volume_a3, dtype=float) / V0_A3) ** q


def theta(volume_a3: np.ndarray | float, gamma0: float, q: float, theta0: float):
    current_gamma = gamma(volume_a3, gamma0, q)
    if q == 0.0:
        return theta0 * (np.asarray(volume_a3, dtype=float) / V0_A3) ** (-gamma0)
    return theta0 * np.exp((gamma0 - current_gamma) / q)


def isentropic_temperature(
    volume_a3: np.ndarray | float, gamma0: float, q: float, theta0: float
):
    """Equation 10: T_S/T0 equals theta(V)/theta0 for gamma=gamma(V)."""
    return TR_K * theta(volume_a3, gamma0, q, theta0) / theta0


def isentropic_pressure(volume_a3: np.ndarray | float, K0: float, K0_prime: float):
    strain = 0.5 * ((V0_A3 / np.asarray(volume_a3, dtype=float)) ** (2.0 / 3.0) - 1.0)
    return (
        3.0
        * K0
        * strain
        * (1.0 + 2.0 * strain) ** 2.5
        * (1.0 + 1.5 * (K0_prime - 4.0) * strain)
    )


def isentropic_energy_j_kg(volume_a3: np.ndarray | float, K0: float, K0_prime: float):
    density0 = float(density_from_cell_volume(V0_A3))
    strain = 0.5 * ((V0_A3 / np.asarray(volume_a3, dtype=float)) ** (2.0 / 3.0) - 1.0)
    return 4.5 * K0 / density0 * 1.0e6 * (strain**2 + (K0_prime - 4.0) * strain**3)


def specific_debye_energy_j_g(
    volume_a3: np.ndarray | float,
    temperature_k: np.ndarray | float,
    gamma0: float,
    q: float,
    Cvm: float,
    theta0: float,
):
    characteristic = theta(volume_a3, gamma0, q, theta0)
    return (
        Cvm
        * np.asarray(temperature_k, dtype=float)
        * debye3(characteristic / np.asarray(temperature_k, dtype=float))
    )


def static_pressure(
    volume_a3: np.ndarray | float,
    temperature_k: np.ndarray | float,
    parameters: np.ndarray,
):
    K0, K0_prime, gamma0, q, Cvm, theta0 = parameters
    volumes = np.asarray(volume_a3, dtype=float)
    density = density_from_cell_volume(volumes)
    reference_temperature = isentropic_temperature(volumes, gamma0, q, theta0)
    energy_difference = specific_debye_energy_j_g(
        volumes, temperature_k, gamma0, q, Cvm, theta0
    ) - specific_debye_energy_j_g(
        volumes, reference_temperature, gamma0, q, Cvm, theta0
    )
    # (g cm^-3)(J g^-1) = 10^-3 GPa.
    thermal = gamma(volumes, gamma0, q) * density * energy_difference / 1000.0
    return isentropic_pressure(volumes, K0, K0_prime) + thermal


def transition_energy_j_kg(starting_material: str) -> float:
    if starting_material == "Enstatite":
        return 1_189_000.0
    if starting_material == "Pv80Mj20":
        # Table 3 gives 189000 - 682000 X_Mj.  The only PPv aggregate row
        # reports 20 vol% majorite; the paper does not publish a separate
        # mole-fraction conversion used by the inversion.
        return 189_000.0 - 682_000.0 * 0.20
    raise ValueError(f"unsupported starting material {starting_material!r}")


def shock_pressure(
    density_g_cm3: np.ndarray | float,
    initial_density_g_cm3: np.ndarray | float,
    transition_energy: np.ndarray | float,
    parameters: np.ndarray,
):
    """Solve equations 5-9 of Mosenfelder et al. (2007) in closed form."""
    K0, K0_prime, gamma0, q, _Cvm, _theta0 = parameters
    density = np.asarray(density_g_cm3, dtype=float)
    initial_density = np.asarray(initial_density_g_cm3, dtype=float)
    volume = cell_volume_from_density(density)
    pressure_s = isentropic_pressure(volume, K0, K0_prime)
    energy_s = isentropic_energy_j_kg(volume, K0, K0_prime)
    current_gamma = gamma(volume, gamma0, q)
    compression_work_coefficient = (
        0.5 * density * current_gamma * (1.0 / initial_density - 1.0 / density)
    )
    numerator = pressure_s - (
        density
        * current_gamma
        * (np.asarray(transition_energy, dtype=float) + energy_s)
        / 1.0e6
    )
    return numerator / (1.0 - compression_work_coefficient)


def shock_temperature(
    density_g_cm3: float,
    initial_density_g_cm3: float,
    transition_energy: float,
    pressure_gpa: float,
    parameters: np.ndarray,
) -> float:
    """Infer T_H from equations 7 and 10-11 for one model shock state."""
    K0, K0_prime, gamma0, q, Cvm, theta0 = parameters
    volume = float(cell_volume_from_density(density_g_cm3))
    energy_s = float(isentropic_energy_j_kg(volume, K0, K0_prime))
    energy_h = (
        0.5 * pressure_gpa * (1.0 / initial_density_g_cm3 - 1.0 / density_g_cm3) * 1.0e6
    )
    energy_v_j_g = (energy_h - transition_energy - energy_s) / 1000.0
    temperature_s = float(isentropic_temperature(volume, gamma0, q, theta0))
    energy_s_debye = float(
        specific_debye_energy_j_g(volume, temperature_s, gamma0, q, Cvm, theta0)
    )

    def residual(temperature: float) -> float:
        return float(
            specific_debye_energy_j_g(volume, temperature, gamma0, q, Cvm, theta0)
            - energy_s_debye
            - energy_v_j_g
        )

    return float(root_scalar(residual, bracket=(1.0, 20_000.0)).root)


def arrays():
    static_rows = load_csv(STATIC_DATA)
    shock_rows = load_csv(SHOCK_DATA)
    static = {
        name: np.array([number(row, name) for row in static_rows])
        for name in static_rows[0]
    }
    shock = {
        name: np.array([number(row, name) for row in shock_rows])
        for name in shock_rows[0]
        if name
        not in {
            "starting_material",
            "reference",
            "shot",
            "flyer_driver_material",
            "phase_state",
        }
    }
    shock["transition_energy_j_kg"] = np.array(
        [transition_energy_j_kg(row["starting_material"]) for row in shock_rows]
    )
    shock["starting_material"] = np.array(
        [row["starting_material"] for row in shock_rows]
    )
    shock["shot"] = np.array([row["shot"] for row in shock_rows])
    return static, shock


def direct_residuals(parameters: np.ndarray, static: dict, shock: dict) -> np.ndarray:
    static_prediction = static_pressure(
        static["volume_a3"], static["temperature_k"], parameters
    )
    static_residual = (static_prediction - static["pressure_gpa"]) / static[
        "pressure_uncertainty_gpa"
    ]

    shock_prediction = shock_pressure(
        shock["shock_density_mg_m3"],
        shock["initial_density_mg_m3"],
        shock["transition_energy_j_kg"],
        parameters,
    )
    shock_residual = (shock_prediction - shock["pressure_gpa"]) / shock[
        "pressure_uncertainty_gpa"
    ]

    # Section 3 treats temperature as a fitted coordinate in principle, but
    # section 3 paragraph 19 explicitly says the shock-temperature observations
    # are fitted for MgSiO3 liquid.  They are not input constraints for this
    # solid PPv inversion and are printed below only as independent diagnostics.
    return np.concatenate((static_residual, shock_residual))


def orthogonal_residuals(values: np.ndarray, static: dict, shock: dict) -> np.ndarray:
    """Reproduce the source's errors-in-P,V,T projection objective."""
    parameters = values[:6]
    count_static = static["pressure_gpa"].size
    model_static_volume = values[6 : 6 + count_static]
    heated = np.isfinite(static["temperature_uncertainty_k"])
    count_heated = int(np.sum(heated))
    model_static_temperature = static["temperature_k"].copy()
    model_static_temperature[heated] = values[
        6 + count_static : 6 + count_static + count_heated
    ]
    model_shock_density = values[6 + count_static + count_heated :]

    static_pressure_residual = (
        static_pressure(model_static_volume, model_static_temperature, parameters)
        - static["pressure_gpa"]
    ) / static["pressure_uncertainty_gpa"]
    static_volume_residual = (model_static_volume - static["volume_a3"]) / static[
        "volume_uncertainty_a3"
    ]
    static_temperature_residual = (
        model_static_temperature[heated] - static["temperature_k"][heated]
    ) / static["temperature_uncertainty_k"][heated]

    shock_pressure_residual = (
        shock_pressure(
            model_shock_density,
            shock["initial_density_mg_m3"],
            shock["transition_energy_j_kg"],
            parameters,
        )
        - shock["pressure_gpa"]
    ) / shock["pressure_uncertainty_gpa"]
    shock_density_residual = (
        model_shock_density - shock["shock_density_mg_m3"]
    ) / shock["shock_density_uncertainty_mg_m3"]
    return np.concatenate(
        (
            static_pressure_residual,
            static_volume_residual,
            static_temperature_residual,
            shock_pressure_residual,
            shock_density_residual,
        )
    )


def diagnostics(parameters: np.ndarray, static: dict, shock: dict) -> None:
    static_delta = (
        static_pressure(static["volume_a3"], static["temperature_k"], parameters)
        - static["pressure_gpa"]
    )
    shock_model = shock_pressure(
        shock["shock_density_mg_m3"],
        shock["initial_density_mg_m3"],
        shock["transition_energy_j_kg"],
        parameters,
    )
    print(f"static pressure RMSE: {np.sqrt(np.mean(static_delta**2)):.4f} GPa")
    print(
        f"shock pressure RMSE: {np.sqrt(np.mean((shock_model - shock['pressure_gpa']) ** 2)):.4f} GPa"
    )
    for index in range(len(shock_model)):
        predicted_temperature = shock_temperature(
            shock["shock_density_mg_m3"][index],
            shock["initial_density_mg_m3"][index],
            shock["transition_energy_j_kg"][index],
            shock_model[index],
            parameters,
        )
        measured = shock["hugoniot_temperature_k"][index]
        measured_text = f"; measured {measured:.0f} K" if np.isfinite(measured) else ""
        print(
            f"shock {shock['shot'][index]}: Pmodel={shock_model[index]:.3f} GPa; "
            f"THmodel={predicted_temperature:.0f} K{measured_text}"
        )


def fit_orthogonal(static: dict, shock: dict):
    """Return the source-style errors-in-variables least-squares result."""
    heated = np.isfinite(static["temperature_uncertainty_k"])
    start = np.concatenate(
        (
            PUBLISHED,
            static["volume_a3"],
            static["temperature_k"][heated],
            shock["shock_density_mg_m3"],
        )
    )
    lower = np.concatenate(
        (
            [100.0, 2.0, 0.1, 0.01, 0.2, 100.0],
            np.full(static["volume_a3"].size, 100.0),
            np.full(int(np.sum(heated)), 300.0),
            np.full(shock["shock_density_mg_m3"].size, 4.0),
        )
    )
    upper = np.concatenate(
        (
            [400.0, 8.0, 8.0, 8.0, 2.0, 3000.0],
            np.full(static["volume_a3"].size, 140.0),
            np.full(int(np.sum(heated)), 5000.0),
            np.full(shock["shock_density_mg_m3"].size, 7.0),
        )
    )
    return least_squares(
        orthogonal_residuals,
        start,
        args=(static, shock),
        bounds=(lower, upper),
        x_scale="jac",
        max_nfev=20_000,
        ftol=1.0e-12,
        xtol=1.0e-12,
        gtol=1.0e-12,
    )


def write_derived_reference_states(shock: dict) -> None:
    """Write model-derived isentrope states separately from shock observations."""
    fields = [
        "shot",
        "reported_shock_density_mg_m3",
        "derived_cell_volume_a3",
        "derived_isentrope_pressure_gpa",
        "derived_isentrope_temperature_k",
        "derived_isentrope_energy_j_kg",
        "transition_energy_j_kg",
        "model_hugoniot_pressure_at_reported_density_gpa",
        "model_hugoniot_temperature_k",
    ]
    rows = []
    model_pressure = shock_pressure(
        shock["shock_density_mg_m3"],
        shock["initial_density_mg_m3"],
        shock["transition_energy_j_kg"],
        PUBLISHED,
    )
    for index, shot in enumerate(shock["shot"]):
        density = float(shock["shock_density_mg_m3"][index])
        volume = float(cell_volume_from_density(density))
        rows.append(
            {
                "shot": shot,
                "reported_shock_density_mg_m3": f"{density:.8g}",
                "derived_cell_volume_a3": f"{volume:.10g}",
                "derived_isentrope_pressure_gpa": f"{float(isentropic_pressure(volume, *PUBLISHED[:2])):.10g}",
                "derived_isentrope_temperature_k": f"{float(isentropic_temperature(volume, PUBLISHED[2], PUBLISHED[3], PUBLISHED[5])):.10g}",
                "derived_isentrope_energy_j_kg": f"{float(isentropic_energy_j_kg(volume, *PUBLISHED[:2])):.10g}",
                "transition_energy_j_kg": f"{float(shock['transition_energy_j_kg'][index]):.10g}",
                "model_hugoniot_pressure_at_reported_density_gpa": f"{float(model_pressure[index]):.10g}",
                "model_hugoniot_temperature_k": f"{shock_temperature(density, float(shock['initial_density_mg_m3'][index]), float(shock['transition_energy_j_kg'][index]), float(model_pressure[index]), PUBLISHED):.10g}",
            }
        )
    with DERIVED_DATA.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-derived",
        action="store_true",
        help="write the explicitly model-derived reference-isentrope table",
    )
    args = parser.parse_args()
    static, shock = arrays()
    if args.write_derived:
        write_derived_reference_states(shock)
    print("Published Table 4 coefficients")
    for name, value in zip(PARAMETER_NAMES, PUBLISHED):
        print(f"  {name}={value:g}")
    diagnostics(PUBLISHED, static, shock)

    fit = least_squares(
        direct_residuals,
        PUBLISHED,
        args=(static, shock),
        bounds=(
            [100.0, 2.0, 0.1, 0.01, 0.2, 100.0],
            [400.0, 8.0, 8.0, 8.0, 2.0, 3000.0],
        ),
        x_scale="jac",
        max_nfev=20_000,
    )
    dof = direct_residuals(fit.x, static, shock).size - fit.x.size
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    print("\nDirect-coordinate weighted reproduction")
    for name, value, error in zip(PARAMETER_NAMES, fit.x, errors):
        print(f"  {name}={value:.8g} +/- {error:.3g}")
    print(f"  reduced chi2={2.0 * fit.cost / dof:.6g}; dof={dof}")
    diagnostics(fit.x, static, shock)

    orthogonal = fit_orthogonal(static, shock)
    parameter_fit = orthogonal.x[:6]
    source_dof = static["pressure_gpa"].size + shock["pressure_gpa"].size - 6
    print("\nOrthogonal P-V-T / P-density reproduction")
    for name, value in zip(PARAMETER_NAMES, parameter_fit):
        print(f"  {name}={value:.8g}")
    print(
        f"  reduced chi2={2.0 * orthogonal.cost / source_dof:.6g}; "
        f"source-counting dof={source_dof}"
    )
    diagnostics(parameter_fit, static, shock)


if __name__ == "__main__":
    main()
