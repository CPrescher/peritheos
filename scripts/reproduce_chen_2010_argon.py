"""Audit the unresolved Chen (2010) fcc argon EOS; no published-fit claim.

Run with python -m scripts.reproduce_chen_2010_argon. No source PDF or
optional PDF dependencies are required for the audit of bundled coordinates.
"""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from peritheos.eos.rt import BM3

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/data/chen-2010-argon-reproduction.json"
POWERS = np.array([3.0, 7 / 3, 5 / 3])
MASS_PER_CELL = 4 * 39.948 / 0.602214076


def heat_capacity_ratio(pressure):
    """Published expression following Eq. (4), pressure in GPa."""
    return 1 + 0.25 * np.exp(-0.197 * np.asarray(pressure))


def density_pressure_derivative(pressure, vp, vs):
    """Eq. (4): velocities km/s, output (g/cm^3)/GPa.

    This is a pointwise relation, not an invented velocity interpolation.
    Full original velocity samples and integration procedure are unavailable.
    """
    denominator = np.asarray(vp) ** 2 - 4 * np.asarray(vs) ** 2 / 3
    if np.any(denominator <= 0):
        raise ValueError("Positive acoustic bulk modulus is required")
    return heat_capacity_ratio(pressure) / denominator


def local_bm3_coefficients(p=2.0, bulk=15.1, derivative=5.4):
    """Unique BM3 density polynomial matching three local constraints.

    For x=rho/rho_ref, P=sum(c_i*x**n_i), K=D P and KK'=D^2 P,
    where D=x*d/dx and n=(3,7/3,5/3). This does NOT assume zero
    pressure at rho_ref, nor add a constant pressure to a zero-pressure BM3.
    """
    return np.linalg.solve(
        np.array([np.ones(3), POWERS, POWERS**2]), [p, bulk, bulk * derivative]
    )


def polynomial_pressure(density, coefficients, density_ref=2.18):
    x = np.asarray(density)[..., None] / density_ref
    return np.sum(coefficients * x**POWERS, axis=-1)


def read_points(suffix):
    path = ROOT / f"peritheos/data/datasets/argon-fcc-chen-2010-{suffix}.csv"
    with path.open() as stream:
        return list(csv.DictReader(stream))


def reproduce():
    c = local_bm3_coefficients()
    # The larger positive root is the mechanically stable zero-pressure branch.
    x0 = float(max(np.roots(c)) ** 1.5)
    rho0 = 2.18 * x0
    k0 = float(np.sum(POWERS * c * x0**POWERS))
    kp0 = float(np.sum(POWERS**2 * c * x0**POWERS) / k0)
    model = BM3(MASS_PER_CELL / rho0, k0, kp0)
    density_grid = np.linspace(2.02, 4, 80)
    native_delta = model.pressure(MASS_PER_CELL / density_grid) - polynomial_pressure(
        density_grid, c
    )
    curve = read_points("figure5-curve")
    pressure = np.array([float(r["pressure_gpa"]) for r in curve])
    density = np.array([float(r["density_g_cm3"]) for r in curve])
    near = (pressure > 1) & (pressure < 3)
    # A local diagnostic of the drawn line, not a refit of measured data.
    poly = np.polyfit(pressure[near], density[near], 3)
    curve_rho2 = float(np.polyval(poly, 2))
    slope = float(np.polyval(np.polyder(poly), 2))
    obs = read_points("figure5-brillouin")
    p_obs = np.array([float(r["pressure_gpa"]) for r in obs])
    rho_obs = np.array([float(r["density_g_cm3"]) for r in obs])
    residual = polynomial_pressure(rho_obs, c) - p_obs
    return {
        "study": "Chen et al. (2010), doi:10.1103/PhysRevB.81.144110",
        "reproduction_status": "not_reproduced",
        "outcome": "supporting_study_with_digitized_density_and_separate_curve",
        "reference_state": {
            "pressure_gpa": 2,
            "temperature_k": 290,
            "density_g_cm3": 2.18,
        },
        "standard_bm3_local_constraint_diagnostic": {
            "density_polynomial_coefficients_gpa": c.tolist(),
            "derived_zero_pressure_parameters_not_published": {
                "V0": MASS_PER_CELL / rho0,
                "K0": k0,
                "K0_prime": kp0,
            },
            "predicted_rho0_g_cm3": rho0,
            "printed_rho0_g_cm3": 1.52,
            "printed_rho0_uncertainty_g_cm3": 0.05,
            "rho0_difference_divided_by_printed_uncertainty": (rho0 - 1.52) / 0.05,
            "predicted_rho_at_1_3_gpa": 2.18
            * brentq(lambda x: np.sum(c * x**POWERS) - 1.3, x0, 1),
            "printed_integration_anchor_rho_g_cm3": 2.02,
            "predicted_K_Ksecond_at_2gpa": float(np.sum(POWERS**3 * c) / 15.1 - 5.4**2),
            "printed_K_Ksecond_at_2gpa": -7.3,
            "printed_K_Ksecond_uncertainty": 1.2,
            "native_pressure_max_difference_gpa": float(np.max(np.abs(native_delta))),
            "density_marker_pressure_rms_gpa": float(np.sqrt(np.mean(residual**2))),
            "density_marker_pressure_max_abs_gpa": float(np.max(np.abs(residual))),
        },
        "figure5_curve_diagnostic": {
            "rows": len(curve),
            "rho_at_2_gpa_g_cm3": curve_rho2,
            "local_cubic_1_to_3_gpa_bulk_modulus_at_2_gpa": curve_rho2 / slope,
            "note": "Derivative of digitized fitted line; not an independent measured elastic constant.",
        },
        "brillouin_density_markers": {
            "distinct_pdf_positions": len(obs),
            "pressure_range_gpa": [float(min(p_obs)), float(max(p_obs))],
            "rows_with_matched_pressure_errorbar": sum(
                bool(r["pressure_plot_halfwidth_gpa"]) for r in obs
            ),
            "rows_with_matched_density_errorbar": sum(
                bool(r["density_plot_halfwidth_g_cm3"]) for r in obs
            ),
            "note": "PDF render multiplicities are not sample counts; integrated densities, not XRD observations.",
        },
        "equation4_check": {
            "Cp_over_Cv_at_2gpa": float(heat_capacity_ratio(2)),
            "drho_dP_from_printed_figure4_velocities": float(
                density_pressure_derivative(2, 3.228, 1.153)
            ),
        },
        "limitation": "No unrounded EOS coefficients, original velocity/density table, covariance or integration weights supplied. Rounded-constraint inconsistency is not a significance test; errors are correlated and confidence convention is unspecified. No executable Chen EOS record is bundled.",
    }


if __name__ == "__main__":
    REPORT.write_text(json.dumps(reproduce(), indent=2) + "\n")
    print(REPORT)
