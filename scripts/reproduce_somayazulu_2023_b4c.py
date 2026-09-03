#!/usr/bin/env python3
"""Reproduce the Somayazulu et al. (2023) B4C thermal-EOS checks.

Data source: official supplement, Table 4, doi:10.6084/m9.figshare.c.6751752.
The publication DOI is 10.1098/rsta.2022.0331.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye, ThermalReferenceStateEOS
from peritheos.fitting import FitResult, fit_thermal_eos
from peritheos.units import cell_volume_to_molar_volume

DATASET_ID = "b4c_somayazulu_2023_table4_pvt"
CELL_TO_MOLAR_CM3 = float(cell_volume_to_molar_volume(1.0, 9.0))
FIXED_REFERENCE = BM3(
    V0=328.4 * CELL_TO_MOLAR_CM3,
    K0=221.0,
    K0_prime=3.3,
)
FIXED_THERMAL = {"Tr": 298.0, "theta0": 1425.0, "n": 5.0}


def load_data() -> np.ndarray:
    """Load the embedded 51-row primary P-V-T dataset from the material file."""
    document = get_material_document("b4c")
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    names = [column["name"] for column in dataset["columns"]]
    return np.array(
        [tuple(row) for row in dataset["rows"]],
        dtype=[(name, float) for name in names],
    )


def fit_q(data: np.ndarray, *, errors_in_variables: bool) -> FitResult:
    """Fit q with the article's gamma0=0.8 held fixed."""
    kwargs = {}
    if errors_in_variables:
        kwargs = {
            "volume_sigma": data["volume_sigma_a3"] * CELL_TO_MOLAR_CM3,
            "temperature_sigma": data["temperature_sigma_k"],
        }
    return fit_thermal_eos(
        MieGruneisenDebye,
        FIXED_REFERENCE,
        data["volume_a3"] * CELL_TO_MOLAR_CM3,
        data["temperature_k"],
        data["pressure_gpa"],
        initial={"q": 2.1},
        fixed={**FIXED_THERMAL, "gamma0": 0.8},
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={"q": (0.0, 8.0)},
        pressure_sigma=data["pressure_sigma_gpa"],
        absolute_sigma=True,
        **kwargs,
    )


def fit_gamma_q(data: np.ndarray) -> FitResult:
    """Fit gamma0 and q with all reported P, V, and T uncertainties."""
    return fit_thermal_eos(
        MieGruneisenDebye,
        FIXED_REFERENCE,
        data["volume_a3"] * CELL_TO_MOLAR_CM3,
        data["temperature_k"],
        data["pressure_gpa"],
        initial={"gamma0": 0.9, "q": 2.1},
        fixed=FIXED_THERMAL,
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={"gamma0": (0.1, 3.0), "q": (0.0, 8.0)},
        pressure_sigma=data["pressure_sigma_gpa"],
        volume_sigma=data["volume_sigma_a3"] * CELL_TO_MOLAR_CM3,
        temperature_sigma=data["temperature_sigma_k"],
        absolute_sigma=True,
    )


def fit_berman(data: np.ndarray, *, reference_volume_law: str = "berman") -> FitResult:
    """Refit a Berman or article-displayed thermal reference-state model."""
    return fit_thermal_eos(
        ThermalReferenceStateEOS,
        FIXED_REFERENCE,
        data["volume_a3"] * CELL_TO_MOLAR_CM3,
        data["temperature_k"],
        data["pressure_gpa"],
        initial={"alpha0": 1.94e-5, "dK_dT": -0.008},
        fixed={"Tr": 298.0, "alpha1": 0.0573e-8},
        configuration={
            "thermal_expansion_law": "linear_temperature",
            "reference_volume_law": reference_volume_law,
        },
        bounds={"alpha0": (-1.0e-4, 1.0e-4), "dK_dT": (-0.2, 0.2)},
        pressure_sigma=data["pressure_sigma_gpa"],
        volume_sigma=data["volume_sigma_a3"] * CELL_TO_MOLAR_CM3,
        temperature_sigma=data["temperature_sigma_k"],
        absolute_sigma=True,
    )


def berman_model(
    alpha0: float,
    dK_dT: float,
    *,
    reference_volume_law: str = "berman",
) -> ThermalReferenceStateEOS:
    """Build the fixed-reference model used for an EosFit-style comparison."""
    return ThermalReferenceStateEOS(
        FIXED_REFERENCE,
        Tr=298.0,
        alpha0=alpha0,
        dK_dT=dK_dT,
        alpha1=0.0573e-8,
        thermal_expansion_law="linear_temperature",
        reference_volume_law=reference_volume_law,
    )


def effective_pressure_sigma(
    model: ThermalReferenceStateEOS, data: np.ndarray
) -> np.ndarray:
    """Return EosFit7's effective pressure standard deviations.

    Angel et al. (2014) define sigma_eff^2 as sigma_P^2 plus the
    volume and temperature variances propagated through dP/dV and dP/dT.
    """
    volume = data["volume_a3"] * CELL_TO_MOLAR_CM3
    sigma_volume = data["volume_sigma_a3"] * CELL_TO_MOLAR_CM3
    temperature = data["temperature_k"]
    bulk_modulus = np.asarray(model.bulk_modulus(volume, temperature), dtype=float)
    dp_dv = -bulk_modulus / volume
    temperature_step = 0.1
    dp_dt = (
        np.asarray(model.pressure(volume, temperature + temperature_step), dtype=float)
        - np.asarray(
            model.pressure(volume, temperature - temperature_step), dtype=float
        )
    ) / (2.0 * temperature_step)
    return np.sqrt(
        data["pressure_sigma_gpa"] ** 2
        + (sigma_volume * dp_dv) ** 2
        + (data["temperature_sigma_k"] * dp_dt) ** 2
    )


def fit_berman_effective_variance(
    data: np.ndarray, *, reference_volume_law: str = "berman"
) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Emulate EosFit7's iteratively updated effective-variance weights."""
    volume = data["volume_a3"] * CELL_TO_MOLAR_CM3
    observed = data["pressure_gpa"]
    temperature = data["temperature_k"]
    values = np.array([1.94e-5, -0.008])
    optimization = None
    iterations = 0
    for iterations in range(1, 101):
        model = berman_model(*values, reference_volume_law=reference_volume_law)
        sigma = effective_pressure_sigma(model, data)

        def residual(parameters: np.ndarray) -> np.ndarray:
            candidate = berman_model(
                *parameters, reference_volume_law=reference_volume_law
            )
            return (
                np.asarray(candidate.pressure(volume, temperature), dtype=float)
                - observed
            ) / sigma

        optimization = least_squares(
            residual,
            values,
            bounds=([-1.0e-4, -0.2], [1.0e-4, 0.2]),
            x_scale="jac",
        )
        if np.allclose(optimization.x, values, rtol=1.0e-12, atol=1.0e-14):
            values = optimization.x
            break
        values = optimization.x
    assert optimization is not None
    final_model = berman_model(*values, reference_volume_law=reference_volume_law)
    final_sigma = effective_pressure_sigma(final_model, data)
    final_residual = (
        np.asarray(final_model.pressure(volume, temperature), dtype=float) - observed
    ) / final_sigma
    covariance = np.linalg.pinv(optimization.jac.T @ optimization.jac)
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    reduced_chi_square = float(np.sum(final_residual**2) / (len(data) - 2))
    return values, errors, reduced_chi_square, iterations


def print_effective_variance_fit(
    label: str, data: np.ndarray, *, reference_volume_law: str = "berman"
) -> None:
    """Print an EosFit-style effective-variance fit."""
    values, errors, reduced_chi_square, iterations = fit_berman_effective_variance(
        data, reference_volume_law=reference_volume_law
    )
    print(
        f"{label}: alpha0={values[0]:.8g} +/- {errors[0]:.3g}, "
        f"dK_dT={values[1]:.8g} +/- {errors[1]:.3g}; "
        f"reduced_chi2={reduced_chi_square:.6g}; IRLS iterations={iterations}"
    )


def published_berman_diagnostics(data: np.ndarray) -> tuple[float, float]:
    """Return RMSE and reduced effective chi-square for Table 2 parameters."""
    model = berman_model(1.94e-5, -0.008)
    volume = data["volume_a3"] * CELL_TO_MOLAR_CM3
    pressure_residual = (
        np.asarray(model.pressure(volume, data["temperature_k"]), dtype=float)
        - data["pressure_gpa"]
    )
    weighted = pressure_residual / effective_pressure_sigma(model, data)
    return (
        float(np.sqrt(np.mean(pressure_residual**2))),
        float(np.sum(weighted**2) / (len(data) - 2)),
    )


def published_parameter_diagnostics(data: np.ndarray) -> tuple[float, float]:
    """Return RMSE and chi-square per point at the printed V and T values."""
    model = MieGruneisenDebye(
        rt_eos=FIXED_REFERENCE,
        gamma0=0.8,
        q=2.1,
        debye_temperature_law="integrated_gruneisen",
        **FIXED_THERMAL,
    )
    residual = (
        model.pressure(
            data["volume_a3"] * CELL_TO_MOLAR_CM3,
            data["temperature_k"],
        )
        - data["pressure_gpa"]
    )
    rmse = float(np.sqrt(np.mean(residual**2)))
    chi_square_per_point = float(np.mean((residual / data["pressure_sigma_gpa"]) ** 2))
    return rmse, chi_square_per_point


def print_fit(label: str, result: FitResult) -> None:
    values = ", ".join(
        f"{name}={result.parameters[name]:.6g} +/- {result.standard_errors[name]:.3g}"
        for name in result.free_parameters
    )
    print(
        f"{label}: {values}; chi2={result.chi_square:.6g}; "
        f"reduced_chi2={result.reduced_chi_square:.6g}; "
        f"dof={result.degrees_of_freedom}"
    )


def main() -> None:
    data = load_data()
    heated = data[data["temperature_k"] > 300.0]
    rmse, chi_square_per_point = published_parameter_diagnostics(heated)
    print(
        f"Loaded {len(data)} rows ({len(heated)} heated, {len(data) - len(heated)} at 300 K)"
    )
    print(
        "Published gamma0=0.8, q=2.1 on heated rows at reported V,T: "
        f"RMSE={rmse:.6g} GPa; chi2/N={chi_square_per_point:.6g}"
    )
    print_fit(
        "Pressure-only q fit, heated rows", fit_q(heated, errors_in_variables=False)
    )
    print_fit(
        "Full P-V-T EIV q fit, heated rows", fit_q(heated, errors_in_variables=True)
    )
    print_fit("Full P-V-T EIV q fit, all rows", fit_q(data, errors_in_variables=True))
    gamma_q = fit_gamma_q(heated)
    print_fit("Full P-V-T EIV gamma0,q fit, heated rows", gamma_q)
    print(f"gamma0-q correlation={gamma_q.correlation[0, 1]:.6g}")
    berman_heated = fit_berman(heated)
    print_fit("Berman P-V-T EIV fit, heated rows", berman_heated)
    print(f"alpha0-dK_dT correlation={berman_heated.correlation[0, 1]:.6g}")
    print_fit("Berman P-V-T EIV fit, all rows", fit_berman(data))
    berman_rmse, berman_reduced_chi_square = published_berman_diagnostics(heated)
    print(
        "Published Berman parameters on heated rows: "
        f"RMSE={berman_rmse:.6g} GPa; "
        f"effective reduced_chi2={berman_reduced_chi_square:.6g}"
    )
    print_effective_variance_fit(
        "EosFit-style Berman effective-variance emulation, heated rows", heated
    )
    print_effective_variance_fit(
        "Article equation (Fei) effective-variance fit, heated rows",
        heated,
        reference_volume_law="integrated_expansivity",
    )


if __name__ == "__main__":
    main()
