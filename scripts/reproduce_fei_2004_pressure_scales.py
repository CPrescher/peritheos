"""Independent Fei (2004) equations, table-output checks, and diagnostic refits."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos import get_eos_record

DATA = Path(__file__).resolve().parents[1] / "peritheos/data/datasets"
# V0, K0, K0', theta0, gamma0, q, mass, printed mass heat-capacity prefactor.
PARAMETERS = {
    "gold": (67.85, 167.0, 5.0, 170.0, 2.97, 0.7, 196.96657, 0.125),
    "platinum": (60.38, 273.0, 4.8, 230.0, 2.69, 0.5, 195.084, 0.12786),
}


def source_pressure(
    volume, temperature, material="gold", coefficients=None, reference_volume=None
):
    """BM3 + Debye quadrature in cell-volume units; no Peritheos calculations."""
    v0, k0, kp, theta0, gamma0, q, mass, c = PARAMETERS[material]
    if reference_volume is not None:
        v0 = reference_volume
    if coefficients is not None:
        k0, kp, gamma0, q = coefficients
    volume, temperature = np.broadcast_arrays(volume, temperature)
    ratio = volume / v0
    gamma = gamma0 * ratio**q
    theta = theta0 * np.exp(
        -gamma0 * np.log(ratio) if abs(q) < 1e-12 else (gamma0 - gamma) / q
    )

    def energy(t, th):
        y = th / t
        # Four atoms/cell. Printed C is J/(g K); 3 D3 integral / y^3.
        return (
            4
            * mass
            * c
            * t
            * 3
            / y**3
            * quad(lambda x: x**3 / np.expm1(x), 0, y, epsabs=1e-12, epsrel=1e-12)[0]
        )

    increment = np.array(
        [
            energy(t, th) - energy(300.0, th)
            for t, th in zip(temperature.flat, theta.flat)
        ]
    ).reshape(volume.shape)
    x = ratio ** (-1 / 3)
    return 1.5 * k0 * (x**7 - x**5) * (
        1 + 0.75 * (kp - 4) * (x * x - 1)
    ) + gamma * increment / (volume * 602.214076)


def load_rows(material):
    table = 1 if material == "gold" else 2
    with (DATA / f"{material}-fei-2004-table{table}.csv").open() as stream:
        return list(csv.DictReader(stream))


def diagnostic(material, room_only=False):
    rows = load_rows(material)
    if room_only:
        rows = [r for r in rows if float(r["temperature_k"]) == 300]
    v, t, p = [
        np.array([float(r[key]) for r in rows])
        for key in ("volume_a3", "temperature_k", "pressure_gpa")
    ]
    base = np.array(PARAMETERS[material][1:3] + PARAMETERS[material][4:6])
    free = [0, 1] if room_only else [3] if material == "gold" else [0, 1, 2, 3]
    if room_only:
        base[:2] = [290.0, 2.7]
    scales = np.array([100.0, 1.0, 1.0, 1.0])[free]

    def unpack(x):
        result = base.copy()
        result[free] = x * scales
        return result

    fitted = least_squares(
        lambda x: source_pressure(v, t, material, unpack(x)) - p,
        base[free] / scales,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=2000,
    )
    cov = (
        np.linalg.inv(fitted.jac.T @ fitted.jac)
        * np.sum(fitted.fun**2)
        / (len(rows) - len(free))
    )
    errors = np.sqrt(np.diag(cov)) * scales
    published = source_pressure(v, t, material, base)
    suffix = "bm3_300k" if room_only else "bm3_mgd"
    native = get_eos_record(f"{material}_fei_2004_{suffix}").pressure(v, temperature=t)
    output = np.array([float(r[f"{material}_fei_pressure_gpa"]) for r in rows])
    # Local finite-difference bound using independent +/- published errors.
    # Not a confidence interval: parameter correlations and confidence are absent.
    uncertainties = (
        [3.0, 0.2, 0.03, 0.3] if material == "gold" else [3.0, 0.3, 0.03, 0.5]
    )
    if room_only:
        uncertainties[:2] = [10.0, 0.9]
    bound = np.zeros(len(rows))
    for i, error in enumerate(uncertainties):
        plus, minus = base.copy(), base.copy()
        plus[i] += error
        minus[i] -= error
        bound += (
            np.abs(
                source_pressure(v, t, material, plus)
                - source_pressure(v, t, material, minus)
            )
            / 2
        )
    v0 = PARAMETERS[material][0]
    dv = 0.004 if material == "gold" else 0.01
    bound += (
        np.abs(
            source_pressure(v, t, material, base, v0 + dv)
            - source_pressure(v, t, material, base, v0 - dv)
        )
        / 2
    )
    return {
        "rows": len(rows),
        "run_ids": [r["run"] for r in rows],
        "native_equation_max_difference_gpa": float(np.max(np.abs(native - published))),
        "published_rmse_gpa": float(np.sqrt(np.mean((published - p) ** 2))),
        "published_mae_gpa": float(np.mean(np.abs(published - p))),
        "table_output_max_difference_gpa": None
        if room_only
        else float(np.max(np.abs(published - output))),
        "table_output_max_fraction_of_parameter_error_bound": None
        if room_only
        else float(np.max(np.abs(published - output) / bound)),
        "free_parameters": [["K0", "K0_prime", "gamma0", "q"][i] for i in free],
        "published_parameters": base[free].tolist(),
        "published_errors": [uncertainties[i] for i in free],
        "refit_parameters": (fitted.x * scales).tolist(),
        "diagnostic_standard_errors": errors.tolist(),
        "rmse_gpa": float(np.sqrt(np.mean(fitted.fun**2))),
        "solver_success": bool(fitted.success),
        "limitations": "Unweighted pressure residuals, printed volumes and pressures; residual-scaled linearized standard errors are diagnostic only. No source weights, covariance, confidence or calibration-error correlations. Pt V0 and theta0 fixed for this diagnostic; source does not enumerate its full fitted/fixed list. Au q fit excludes upstream shock observations and is not the source joint compromise.",
    }


def reproduce():
    return {
        "gold": diagnostic("gold"),
        "platinum": diagnostic("platinum"),
        "platinum_300k": diagnostic("platinum", True),
    }


def ledger_outcome(record):
    key = "gold" if record["identifier"].startswith("gold") else "platinum"
    if record["identifier"].endswith("300k"):
        key += "_300k"
    result = reproduce()[key]
    within = all(
        abs(a - b) <= e
        for a, b, e in zip(
            result["published_parameters"],
            result["refit_parameters"],
            result["published_errors"],
        )
    )
    # Au's published q deliberately compromises static and shock data.
    status = (
        "not_refittable"
        if key == "gold"
        else "similar"
        if within
        else "parity_not_achieved"
    )
    return {
        "status": status,
        "dataset_identifiers": record["fit_datasets"],
        "observations": result["rows"],
        "fit_kind": "diagnostic_unweighted_pressure",
        "published_rmse_gpa": result["published_rmse_gpa"],
        "rmse_gpa": result["rmse_gpa"],
        "solver_success": result["solver_success"],
        "reason": result["limitations"],
        "parameters": [
            {
                "parameter": n,
                "published": p,
                "refit": v,
                "published_error": e,
                "refit_error": se,
                "relative_difference": abs(v - p) / abs(p),
                "within_combined_2sigma": None,
                "similar": abs(v - p) <= e,
            }
            for n, p, v, e, se in zip(
                result["free_parameters"],
                result["published_parameters"],
                result["refit_parameters"],
                result["published_errors"],
                result["diagnostic_standard_errors"],
            )
        ],
        "reproduction": result,
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2))
