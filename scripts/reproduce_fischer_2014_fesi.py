#!/usr/bin/env python3
"""Independently evaluate Fischer (2014), Tables 1-2 and Equations (1)-(2).

The official S1/S2 observations are selected by their exact phase labels.
Unweighted pressure-residual fits are diagnostic because the authors do not
specify regression weights. No Peritheos evaluator is used. Volumes here are
cm^3/mol atoms; source lattice parameters are converted without rounding.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, Boltzmann
from scipy.integrate import quad
from scipy.optimize import least_squares

# V0 [cm^3/mol atoms], K0 [GPa], K0' [dimensionless].
STATIC = {
    "fesi_b20": (6.803, 192.2, 5.03),
    "fesi_b2": (6.414, 230.6, 4.17),
    "fe084si016_d03": (6.961, 183.0, 5.59),
    "fe084si016_hcp": (7.203, 129.1, 5.29),
}
# theta0 [K], gamma0, q. The reference temperature is 300 K.
THERMAL = {"fesi_b2": (417.0, 1.30, 1.7), "fe084si016_hcp": (417.0, 1.14, 1.0)}
ATOMS_PER_CELL = {
    "fe073si027_d03": 16,
    "fesi_b20": 8,
    "fesi_b2": 2,
    "fe084si016_d03": 16,
    "fe084si016_hcp": 2,
}

VINET_STATIC = {
    "fesi_b20": (6.803, 193.1, 5.03),
    "fe073si027_d03": (6.799, 193.9, 4.91),
    "fesi_b2": (6.414, 236.1, 4.17),
    "fe084si016_d03": (6.961, 182.0, 5.59),
    "fe084si016_hcp": (7.283, 111.4, 6.08),
}
VINET_THERMAL = {
    "fesi_b2": (417.0, 1.17, 1.41),
    "fe084si016_hcp": (417.0, 1.15, 1.0),
    "fe073si027_d03": (417.0, 1.90, 1.0),
}
ROOT = Path(__file__).resolve().parents[1]
PARAMETER_NAMES = ("V0", "K0", "K0_prime", "theta0", "gamma0", "q")
FREE = {
    "fesi_b20": (1, 2),
    "fesi_b2": (1, 4, 5),
    "fe084si016_d03": (1, 2),
    "fe084si016_hcp": (0, 1, 2, 4),
}


def parameters(branch, model="bm3"):
    """Published coefficients in the source's atomic molar normalization."""
    static, thermal = (
        (STATIC, THERMAL) if model == "bm3" else (VINET_STATIC, VINET_THERMAL)
    )
    return np.array(static[branch] + thermal.get(branch, (417.0, 0.0, 1.0)))


def load_data(branch):
    """Read the full source CSV, then apply only the stated single-phase selection."""
    table = 2 if branch.startswith("fesi_") else 1
    path = ROOT / f"peritheos/data/datasets/fischer-2014-table-s{table}-pvt.csv"
    if branch == "fe073si027_d03":
        path = ROOT / "peritheos/data/datasets/fischer-2012-table-s2-pvt.csv"
    phase = branch.rsplit("_", 1)[1]
    label = {"d03": "D03", "b20": "B20", "b2": "B2", "hcp": "hcp"}[phase]
    with path.open(newline="") as stream:
        return [row for row in csv.DictReader(stream) if row["phase"] == label]


def vinet(volume, v0, k0, kp):
    """Text S1, page 1: standard Vinet pressure, in GPa."""
    y = (np.asarray(volume, dtype=float) / v0) ** (1.0 / 3.0)
    return 3 * k0 * (1 - y) / y**2 * np.exp(1.5 * (kp - 1) * (1 - y))


def bm3(volume, v0, k0, kp):
    """Equation (1), expressed in finite strain, with pressure in GPa."""
    strain = ((v0 / np.asarray(volume, dtype=float)) ** (2.0 / 3.0) - 1) / 2
    return 3 * k0 * strain * (1 + 2 * strain) ** 2.5 * (1 + 1.5 * (kp - 4) * strain)


def pressure(branch, volume, temperature=300.0, model="bm3", coefficients=None):
    """Source pressure at molar atomic volume, using energy per physical atom."""
    v, t = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    coeff = parameters(branch, model) if coefficients is None else coefficients
    result = (bm3 if model == "bm3" else vinet)(v, *coeff[:3])
    if branch not in (THERMAL if model == "bm3" else VINET_THERMAL):
        if np.any(t != 300):
            raise ValueError("The source supplies only a 300 K fit for this branch")
        return result
    theta0, gamma0, q = coeff[3:]
    gamma = gamma0 * (v / coeff[0]) ** q
    theta = theta0 * np.exp((gamma0 - gamma) / q)

    def energy(temp, characteristic):
        integral = quad(lambda x: x**3 / np.expm1(x), 0, characteristic / temp)[0]
        return 9 * Boltzmann * temp * (temp / characteristic) ** 3 * integral

    delta = np.array(
        [energy(temp, th) - energy(300.0, th) for temp, th in zip(t.flat, theta.flat)]
    ).reshape(v.shape)
    # cm^3/mol atoms -> m^3/atom, then Pa -> GPa.
    return result + gamma * delta / (v * 1e-6 / Avogadro) / 1e9


def reproduce():
    """Independently refit all eight source branches with source-fixed choices."""
    output = {}
    for branch in STATIC:
        rows = load_data(branch)
        phase = branch.rsplit("_", 1)[1]
        v = (
            np.array([float(r[phase + "_volume_a3"]) for r in rows])
            * Avogadro
            / 1e24
            / ATOMS_PER_CELL[branch]
        )
        t = np.array([float(r["temperature_k"]) for r in rows])
        p = np.array([float(r["pressure_gpa"]) for r in rows])
        for model in ("bm3", "vinet"):
            published = parameters(branch, model)
            free = list(FREE[branch])

            def residual(values):
                coeff = published.copy()
                coeff[free] = values
                return pressure(branch, v, t, model, coeff) - p

            lower = np.array([3.0, 10.0, 0.0, 1.0, 0.01, -5.0])[free]
            upper = np.array([12.0, 500.0, 15.0, 3000.0, 5.0, 10.0])[free]
            fit = least_squares(
                residual,
                published[free],
                bounds=(lower, upper),
                x_scale="jac",
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
            observed_residual = residual(published[free])
            covariance = (
                np.linalg.pinv(fit.jac.T @ fit.jac)
                * np.sum(fit.fun**2)
                / (len(p) - len(free))
            )
            output[branch + "_" + model] = {
                "observations": len(rows),
                "selection": "phase=" + rows[0]["phase"],
                "objective": "unweighted_pressure_residual_diagnostic",
                "solver_success": bool(fit.success),
                "parameters": {
                    PARAMETER_NAMES[i]: float(x) for i, x in zip(free, fit.x)
                },
                "parameter_standard_errors": {
                    PARAMETER_NAMES[i]: float(x)
                    for i, x in zip(free, np.sqrt(np.diag(covariance)))
                },
                "published_rmse_gpa": float(np.sqrt(np.mean(observed_residual**2))),
                "published_residual_range_gpa": [
                    float(observed_residual.min()),
                    float(observed_residual.max()),
                ],
                "refit_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
                "pressure_range_gpa": [float(p.min()), float(p.max())],
                "temperature_range_k": [float(t.min()), float(t.max())],
                "source_checkpoints": [
                    {
                        "source_row": int(rows[i]["source_row"]),
                        "pressure_gpa": float(p[i]),
                        "calculated_pressure_gpa": float(p[i] + observed_residual[i]),
                    }
                    for i in (0, len(rows) // 2, len(rows) - 1)
                ],
            }
    # Source-inherited staged Fe–16Si protocol: static fit, then thermal fit
    # with the printed static coefficients, not our diagnostic static estimates.
    branch = "fe073si027_d03"
    rows = load_data(branch)
    v = np.array([float(r["d03_molar_atomic_volume_cm3_mol"]) for r in rows])
    t = np.array([float(r["temperature_k"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    published = parameters(branch, "vinet")
    mask = t == 300
    static_fit = least_squares(
        lambda x: vinet(v[mask], published[0], *x) - p[mask], published[1:3]
    )

    def thermal_residual(x):
        coeff = published.copy()
        coeff[4] = x[0]
        return pressure(branch, v[~mask], t[~mask], "vinet", coeff) - p[~mask]

    thermal_fit = least_squares(thermal_residual, [published[4]])
    residual = pressure(branch, v, t, "vinet") - p
    output[branch + "_vinet"] = {
        "observations": len(rows),
        "static_observations": int(mask.sum()),
        "thermal_observations": int((~mask).sum()),
        "static_parameters": dict(zip(["K0", "K0_prime"], static_fit.x.tolist())),
        "thermal_parameters": {"gamma0": float(thermal_fit.x[0])},
        "published_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_residual_range_gpa": [float(residual.min()), float(residual.max())],
        "static_refit_rmse_gpa": float(np.sqrt(np.mean(static_fit.fun**2))),
        "thermal_refit_rmse_gpa": float(np.sqrt(np.mean(thermal_fit.fun**2))),
        "source_checkpoints": [
            {
                "source_row": int(rows[i]["source_row"]),
                "pressure_gpa": float(p[i]),
                "calculated_pressure_gpa": float(p[i] + residual[i]),
            }
            for i in (0, len(rows) // 2, len(rows) - 1)
        ],
    }
    # This amplitude sensitivity diagnoses the B2 source inconsistency; it is
    # deliberately not installed as a chemically inconsistent n=4 FeSi model.
    rows = load_data("fesi_b2")
    v = np.array([float(r["b2_volume_a3"]) for r in rows]) * Avogadro / 1e24 / 2
    t = np.array([float(r["temperature_k"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    for model in ("bm3", "vinet"):
        cold = pressure("fesi_b2", v, 300.0, model)
        increment = pressure("fesi_b2", v, t, model) - cold
        residual = cold + 2 * increment - p
        output["fesi_b2_" + model][
            "unpublished_double_thermal_amplitude_diagnostic"
        ] = {
            "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "residual_range_gpa": [float(residual.min()), float(residual.max())],
            "production_accepted": False,
        }
    return output


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, allow_nan=False))
