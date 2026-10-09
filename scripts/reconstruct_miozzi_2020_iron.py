"""Conditional pressure reconstructions and staged Fe fits; no catalog mutations.

Run with the project's Python environment. Output is a diagnostic JSON and an
attributed row-level CSV alongside the original literature reproduction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/data/miozzi-2020-reconstruction"
R = 8.31446261815324
NA = 6.02214076e23
NODES, WEIGHTS = np.polynomial.legendre.leggauss(64)
PUBLISHED = np.array([22.81, 129.0, 6.24, 1.11, 0.3])
MODELS = {
    "speziale_variable_q_debye": {
        "source": "https://doi.org/10.1029/2000JB900318",
        "equations": "Speziale 2001 Eqs. 10-11; theta numerically integrated from gamma",
        "V0_a3": 74.71,
        "K0_gpa": 160.2,
        "K0_prime": 3.99,
        "gamma0": 1.524,
        "q0": 1.65,
        "q1": 11.8,
        "theta0_k": 773.0,
        "reference_temperature_k": 300.0,
        "oscillator": "Debye",
        "gas_constant": R,
    },
    "dorogokupets_2010_fit2_einstein": {
        "source": "https://doi.org/10.1007/s00269-010-0367-2",
        "equations": "Dorogokupets 2010 Eqs. 1 and 5, Table 1 Fit #2",
        "V0_a3": 11.248 * 4e24 / NA,
        "V0_cm3_mol": 11.248,
        "K0_gpa": 160.2,
        "K0_prime": 3.99,
        "gamma0": 1.524,
        "gamma_inf": 1.325,
        "beta": 11.8,
        "theta0_k": 599.0,
        "reference_temperature_k": 298.15,
        "oscillator": "Einstein",
        "gas_constant": 8.31451,
    },
    "speziale_constant_q_debye": {
        "source": "https://doi.org/10.1029/2000JB900318",
        "equations": "Speziale 2001 Eqs. 1-4; constant thermodynamic q sensitivity",
        "V0_a3": 74.71,
        "K0_gpa": 160.2,
        "K0_prime": 3.99,
        "gamma0": 1.524,
        "q": 1.65,
        "theta0_k": 773.0,
        "reference_temperature_k": 300.0,
        "oscillator": "Debye",
        "gas_constant": R,
    },
}


def bm3(v, v0, k0, kp):
    x = (v0 / np.asarray(v)) ** (1 / 3)
    return 1.5 * k0 * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def debye_energy(theta, temperature, atoms, gas_constant=R):
    theta, temperature = np.broadcast_arrays(theta, temperature)
    limit = theta / temperature
    z = limit[..., None] * (NODES + 1) / 2
    integral = limit / 2 * np.sum(WEIGHTS * z**3 / np.expm1(z), axis=-1)
    return 9 * atoms * gas_constant * temperature * integral / limit**3


def constant_q_gamma_theta(ratio, gamma0, q, theta0):
    """Integrated law with a stable q=0 limit, including negative q."""
    logarithm = np.log(ratio)
    gamma = gamma0 * np.exp(q * logarithm)
    integral = logarithm if abs(q) < 1e-8 else np.expm1(q * logarithm) / q
    return gamma, theta0 * np.exp(-gamma0 * integral)


def mgo_gamma_theta(v, model):
    m = MODELS[model]
    ratio = np.asarray(v) / m["V0_a3"]
    if "q" in m:
        return constant_q_gamma_theta(ratio, m["gamma0"], m["q"], m["theta0_k"])
    if "gamma_inf" in m:
        power = ratio ** m["beta"]
        gamma = m["gamma_inf"] + (m["gamma0"] - m["gamma_inf"]) * power
        theta = (
            m["theta0_k"]
            * ratio ** (-m["gamma_inf"])
            * np.exp((m["gamma0"] - m["gamma_inf"]) / m["beta"] * (1 - power))
        )
        return gamma, theta
    logarithm = np.log(ratio)
    gamma = m["gamma0"] * np.exp(m["q0"] / m["q1"] * (ratio ** m["q1"] - 1))
    u = logarithm[..., None] * (NODES + 1) / 2
    integrand = m["gamma0"] * np.exp(m["q0"] / m["q1"] * np.expm1(m["q1"] * u))
    integral = logarithm / 2 * np.sum(WEIGHTS * integrand, axis=-1)
    return gamma, m["theta0_k"] * np.exp(-integral)


def mgo_pressure(v, t, model):
    m = MODELS[model]
    v, t = np.broadcast_arrays(v, t)
    gamma, theta = mgo_gamma_theta(v, model)
    if m["oscillator"] == "Einstein":
        delta_energy = (
            6
            * m["gas_constant"]
            * theta
            * (
                1 / np.expm1(theta / t)
                - 1 / np.expm1(theta / m["reference_temperature_k"])
            )
        )
    else:
        delta_energy = debye_energy(theta, t, 2, m["gas_constant"]) - debye_energy(
            theta, m["reference_temperature_k"], 2, m["gas_constant"]
        )
    # MgO has four formula units per conventional cell, each with two atoms.
    return bm3(v, m["V0_a3"], m["K0_gpa"], m["K0_prime"]) + (
        gamma * delta_energy / (v * NA / 4e24) / 1000
    )


def fe_pressure(v, t, parameters):
    v0, k0, kp = parameters[:3]
    cold = bm3(v, v0, k0, kp)
    if len(parameters) == 3:
        return cold
    gamma, theta = constant_q_gamma_theta(
        np.asarray(v) / v0, parameters[3], parameters[4], 420.0
    )
    return (
        cold
        + gamma
        * (debye_energy(theta, t, 1) - debye_energy(theta, 300, 1))
        / (np.asarray(v) * NA / 2e24)
        / 1000
    )


def derivative(function, value, step):
    return (function(value + step) - function(value - step)) / (2 * step)


def read_data():
    rows, hashes = [], {}
    for medium in ("he", "mgo"):
        path = ROOT / f"peritheos/data/datasets/iron-miozzi-2020-{medium}-pvt.csv"
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        with path.open(encoding="utf-8") as stream:
            rows.extend(dict(row, medium=medium) for row in csv.DictReader(stream))
    names = (
        "pressure_gpa",
        "pressure_error_gpa",
        "temperature_k",
        "temperature_error_k",
        "volume_a3",
        "volume_error_a3",
        "mgo_volume_a3",
        "mgo_volume_error_a3",
    )
    data = {
        n: np.array([float(row[n]) if row.get(n) else np.nan for row in rows])
        for n in names
    }
    data["mgo"] = np.array([row["medium"] == "mgo" for row in rows])
    return rows, data, hashes


def pressure_column(data, model):
    p = data["pressure_gpa"].copy()
    pv, pt = np.zeros_like(p), np.zeros_like(p)
    mask = data["mgo"]
    if model != "printed":
        v, t = data["mgo_volume_a3"][mask], data["temperature_k"][mask]
        p[mask] = mgo_pressure(v, t, model)
        pv[mask] = derivative(lambda x: mgo_pressure(x, t, model), v, 1e-4)
        pt[mask] = derivative(lambda x: mgo_pressure(v, x, model), t, 1e-2)
    return p, pv, pt


def effective_sigma(data, model, parameters, thermal):
    """Residual covariance includes shared temperature of Fe and MgO.

    Assumes other coordinate errors independent; calibration coefficients fixed.
    Missing errors yield NaN and exclude that row from weighted fits.
    """
    _, pv, pt = pressure_column(data, model)
    v, t = data["volume_a3"], data["temperature_k"]
    fv = derivative(lambda x: fe_pressure(x, t, parameters), v, 1e-4)
    ft = derivative(lambda x: fe_pressure(v, x, parameters), t, 1e-2)
    variance = (fv * data["volume_error_a3"]) ** 2
    if model == "printed":
        variance += data["pressure_error_gpa"] ** 2
        if thermal:
            variance += (ft * data["temperature_error_k"]) ** 2
    else:
        mask = data["mgo"]
        variance[mask] += (pv[mask] * data["mgo_volume_error_a3"][mask]) ** 2
        variance[~mask] += data["pressure_error_gpa"][~mask] ** 2
        if thermal:
            variance += ((ft - pt) * data["temperature_error_k"]) ** 2
    return np.sqrt(variance)


def metrics(residual):
    return {
        "row_count": int(len(residual)),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "residual_range_gpa": [float(np.min(residual)), float(np.max(residual))],
    }


def conditional_covariance(result, scales, row_count, weighting):
    if not result.success or np.any(result.active_mask):
        raise ValueError("Fit failed or intersected bounds")
    _, singular, right = np.linalg.svd(result.jac, full_matrices=False)
    if singular[-1] <= singular[0] * max(result.jac.shape) * np.finfo(float).eps:
        raise ValueError("Rank-deficient conditional fit")
    dof = row_count - len(scales)
    if dof <= 0:
        raise ValueError("No covariance degrees of freedom")
    reduced = float(np.sum(result.fun**2) / dof)
    factor = max(1.0, reduced) if weighting == "effective_errors" else reduced
    covariance = (right.T / singular**2) @ right * factor
    covariance *= scales[:, None] * scales[None, :]
    return covariance, dof, reduced, factor, float(singular[0] / singular[-1])


def fit(data, model, selection, weighting, thermal):
    p, _, _ = pressure_column(data, model)
    v, t = data["volume_a3"], data["temperature_k"]
    names = ["K0", "K0_prime", "gamma0", "q"] if thermal else ["V0", "K0", "K0_prime"]
    scales = np.array([100, 5, 1, 1] if thermal else [23, 100, 5])
    lower = np.array([20, 1, 0.1, -5] if thermal else [20, 20, 1]) / scales
    upper = np.array([400, 10, 5, 10] if thermal else [25, 400, 10]) / scales
    initial = (
        np.array([129.0, 6.2, 1.11, 0.3] if thermal else [22.80, 129, 6.2]) / scales
    )

    def full(z):
        return np.r_[22.81, z * scales] if thermal else z * scales

    def sigma(z):
        return (
            effective_sigma(data, model, full(z), thermal)
            if weighting == "effective_errors"
            else np.ones_like(p)
        )

    errors = sigma(initial)
    usable = selection & np.isfinite(errors) & (errors > 0)
    if np.sum(usable) <= len(scales):
        raise ValueError("Insufficient observations for conditional covariance")
    stages = []
    if thermal:
        # Section 3.3 first varies gamma and q with the cold coefficients fixed.
        first = least_squares(
            lambda z: (
                (fe_pressure(v, t, full(np.r_[initial[:2], z]))[usable] - p[usable])
                / errors[usable]
            ),
            initial[2:],
            bounds=(lower[2:], upper[2:]),
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
            max_nfev=2000,
        )
        covariance, dof, reduced, factor, condition = conditional_covariance(
            first, scales[2:], int(np.sum(usable)), weighting
        )
        initial[2:] = first.x
        stages.append(
            {
                "fixed_cold_parameters": [22.81, 129.0, 6.2],
                "gamma0_q": first.x.tolist(),
                "free_parameter_names": ["gamma0", "q"],
                "standard_errors": np.sqrt(np.diag(covariance)).tolist(),
                "covariance": covariance.tolist(),
                "degrees_of_freedom": dof,
                "reduced_objective": reduced,
                "covariance_scale": factor,
                "jacobian_condition_number": condition,
                "weighting": "Initial-stage weights frozen at starting coefficients",
            }
        )
    converged = False
    for cycle in range(100):
        errors = sigma(initial)
        denominator = errors[usable]
        result = least_squares(
            lambda z: (fe_pressure(v, t, full(z))[usable] - p[usable]) / denominator,
            initial,
            bounds=(lower, upper),
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
            max_nfev=2000,
        )
        change = float(np.max(np.abs(result.x - initial)))
        initial = result.x
        if weighting == "unweighted" or change < 1e-7:
            converged = True
            break
    if not converged or not result.success or np.any(result.active_mask):
        raise ValueError(
            f"Fit failed/intersected bounds: {model}, {weighting}, {thermal}"
        )
    covariance, dof, reduced, factor, condition = conditional_covariance(
        result, scales, int(np.sum(usable)), weighting
    )
    fitted = full(initial)
    residual = fe_pressure(v, t, fitted) - p
    return {
        "parameters": fitted.tolist(),
        "free_parameter_names": names,
        "standard_errors": np.sqrt(np.diag(covariance)).tolist(),
        "covariance": covariance.tolist(),
        "degrees_of_freedom": dof,
        "reduced_objective": reduced,
        "covariance_scale": factor,
        "weighted_row_count": int(np.sum(usable)),
        "excluded_combined_row_numbers": (
            np.flatnonzero(selection & ~usable) + 1
        ).tolist(),
        "selection_metrics": metrics(residual[selection]),
        "fit_row_metrics": metrics(residual[usable]),
        "solver_success": bool(result.success),
        "weight_cycles": cycle + 1,
        "scaled_parameter_change": change,
        "jacobian_condition_number": condition,
        "initial_thermal_stage": stages,
    }


def reconstruct():
    rows, data, hashes = read_data()
    cold = data["temperature_k"] <= 300
    heated = ~cold
    outputs, derived = {}, []
    for model in ["printed", *MODELS]:
        p, pv, pt = pressure_column(data, model)
        fit_results = {}
        for selection_name, selection, thermal in [
            ("cold", cold, False),
            ("thermal_all", np.ones_like(cold), True),
            ("thermal_heated", heated, True),
        ]:
            fit_results[selection_name] = {
                mode: fit(data, model, selection, mode, thermal)
                for mode in ("unweighted", "effective_errors")
            }
        published_residual = (
            fe_pressure(data["volume_a3"], data["temperature_k"], PUBLISHED) - p
        )
        outputs[model] = {
            "mgo_pressure_shift": metrics((p - data["pressure_gpa"])[data["mgo"]]),
            "published_fe": {
                "all": metrics(published_residual),
                "heated": metrics(published_residual[heated]),
                "cold_bm3": metrics(
                    (
                        fe_pressure(
                            data["volume_a3"], data["temperature_k"], [22.80, 129, 6.2]
                        )
                        - p
                    )[cold]
                ),
            },
            "fits": fit_results,
            "doubled_fe_energy_hypothesis": {
                "qualification": "Deliberate normalization sensitivity only: doubles vibrational energy at unchanged gamma/theta and molar Fe volume. This is not the physical n=1 normalization and is not an adopted correction or evidence of authors' settings.",
                "published_coefficients": metrics(
                    2 * fe_pressure(data["volume_a3"], data["temperature_k"], PUBLISHED)
                    - fe_pressure(
                        data["volume_a3"], data["temperature_k"], PUBLISHED[:3]
                    )
                    - p
                ),
            },
        }
        if model == "printed":
            continue
        for i, row in enumerate(rows):
            if not data["mgo"][i]:
                continue
            sigma_p = np.sqrt(
                (pv[i] * data["mgo_volume_error_a3"][i]) ** 2
                + (pt[i] * data["temperature_error_k"][i]) ** 2
            )
            derived.append(
                dict(
                    model=model,
                    source_row=row["source_row"],
                    source_page=row["source_page"],
                    printed_pressure_gpa=row["pressure_gpa"],
                    printed_pressure_error_gpa=row["pressure_error_gpa"],
                    temperature_k=row["temperature_k"],
                    temperature_error_k=row["temperature_error_k"],
                    fe_volume_a3=row["volume_a3"],
                    fe_volume_error_a3=row["volume_error_a3"],
                    mgo_volume_a3=row["mgo_volume_a3"],
                    mgo_volume_error_a3=row["mgo_volume_error_a3"],
                    reconstructed_pressure_gpa=f"{p[i]:.10f}",
                    pressure_shift_gpa=f"{p[i] - data['pressure_gpa'][i]:.10f}",
                    conditional_pressure_error_gpa=f"{sigma_p:.10f}"
                    if np.isfinite(sigma_p)
                    else "",
                    pressure_temperature_covariance_gpa_k=f"{pt[i] * data['temperature_error_k'][i] ** 2:.10f}"
                    if np.isfinite(data["temperature_error_k"][i])
                    else "",
                    published_fe_residual_gpa=f"{published_residual[i]:.10f}",
                )
            )
    return {
        "status": "conditional_reconstruction_not_original_fit_recovered",
        "attribution": "Derived from Miozzi et al. (2020), DOI 10.3390/min10020100, supplementary observations, CC BY 4.0",
        "original_csv_sha256": hashes,
        "calibration_models": MODELS,
        "calibration_normalization": {
            "atoms_per_formula_unit": 2,
            "formula_units_per_cell": 4,
        },
        "fixed_fe_thermal": {
            "V0_a3": 22.81,
            "theta0_k": 420,
            "reference_temperature_k": 300,
            "atoms_per_formula_unit": 1,
            "formula_units_per_cell": 2,
        },
        "qualification": [
            "The exact Speziale thermal implementation, original fit input selection and error semantics are unknown. These three documented models are sensitivity scenarios, not identified original pressures.",
            "He pressures unchanged. Missing coordinate errors exclude rows only from weighted fits; all rows remain in unweighted diagnostics. Printed zero errors remain zero.",
            "After reconstruction, printed pressure errors are retained as source metadata but not reused. Residual variance uses MgO/Fe volume errors and (dP_Fe/dT-dP_MgO/dT)^2*sigma_T^2, accounting for shared temperature. Other coordinate errors are assumed independent.",
            "Printed-input control assumes independent P,V,T errors because its calibration covariance is unavailable. Objective and complete-error row selections therefore differ from reconstructed controls.",
            "Printed error widths provisionally treated as standard deviations; calibration parameter uncertainty, inter-row correlations, fixed-parameter uncertainty and source confidence definitions are not propagated.",
            "Initial thermal-only fit holds cold coefficients fixed; final fit varies K0,K0_prime,gamma0,q, updating effective weights between solves. This is not an exact replay of the EosFit console or unspecified authors' iteration sequence.",
            "Covariance is local at final weights, SVD-based; weighted covariance scales by max(1,reduced chi-square), unweighted by RSS/dof. Conditional errors are not reported publication uncertainties.",
            "All-row and heated-only thermal selections are retained to expose original-selection ambiguity. No inferred EOS record is added and published coefficients are unchanged.",
        ],
        "results": outputs,
    }, derived


def plot_results(summary):
    """Static scientific comparison; symbols represent supplied observations."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, data, _ = read_data()
    cold = data["temperature_k"] <= 300
    heated = ~cold
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    labels = {
        "printed": "Supplied pressure column",
        "speziale_variable_q_debye": "Variable-q Debye",
        "dorogokupets_2010_fit2_einstein": "Dorogokupets Fit #2",
        "speziale_constant_q_debye": "Constant-q Debye",
    }
    colors = ["#777777", "#0072B2", "#009E73", "#D55E00"]
    curve_v = np.linspace(16, 21.6, 200)
    axes[0].plot(
        curve_v, bm3(curve_v, 22.80, 129, 6.2), color="black", label="Published Fe BM3"
    )
    for (model, label), color in zip(labels.items(), colors):
        p, _, _ = pressure_column(data, model)
        axes[0].scatter(
            data["volume_a3"][cold], p[cold], s=15, color=color, alpha=0.65, label=label
        )
        residual = fe_pressure(data["volume_a3"], data["temperature_k"], PUBLISHED) - p
        axes[1].scatter(
            data["temperature_k"][heated],
            residual[heated],
            s=12,
            color=color,
            alpha=0.6,
        )
        pars = summary["results"][model]["fits"]["thermal_all"]["unweighted"][
            "parameters"
        ]
        fitted = fe_pressure(data["volume_a3"], data["temperature_k"], pars) - p
        axes[2].scatter(
            data["temperature_k"][heated], fitted[heated], s=12, color=color, alpha=0.6
        )
    axes[0].set(
        xlabel="Fe cell volume (Å³)",
        ylabel="Pressure (GPa)",
        title="Cold observations (36 rows)",
    )
    axes[0].legend(fontsize=7, loc="upper right")
    for ax in axes[1:]:
        ax.axhline(0, color="black", lw=0.7)
        ax.axhspan(-3, 3, color="gray", alpha=0.08)
        ax.set(
            xlabel="Temperature (K)", ylabel="Fe calculated − observed pressure (GPa)"
        )
    axes[1].set_title("Published Fe thermal coefficients")
    axes[2].set_title("Separate equal-weight thermal fits")
    fig.suptitle(
        "Miozzi (2020): conditional MgO reconstruction; original calibration unresolved",
        fontsize=12,
    )
    fig.savefig(OUT / "comparison.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true", help="Save comparison.png")
    args = parser.parse_args()
    summary, rows = reconstruct()
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    documents = {
        "summary.json": json.dumps(summary, indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        "pressures.csv": stream.getvalue(),
    }
    if args.check:
        # Numerical fit results can differ slightly between BLAS/SciPy builds.
        stored = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
        if (
            stored["original_csv_sha256"] != summary["original_csv_sha256"]
            or stored["calibration_models"] != summary["calibration_models"]
        ):
            raise SystemExit("Reconstruction inputs/model metadata are stale")
        for model, result in summary["results"].items():
            for selection, fits in result["fits"].items():
                for mode, fit_result in fits.items():
                    saved = stored["results"][model]["fits"][selection][mode]
                    if saved["weighted_row_count"] != fit_result[
                        "weighted_row_count"
                    ] or not np.allclose(
                        saved["parameters"],
                        fit_result["parameters"],
                        rtol=2e-6,
                        atol=2e-6,
                    ):
                        raise SystemExit(
                            f"Reconstruction fit is stale: {model}/{selection}/{mode}"
                        )
        from scripts.check_numerical_archive import check_csv

        check_csv(
            (OUT / "pressures.csv").read_text(encoding="utf-8"),
            documents["pressures.csv"],
            {
                "reconstructed_pressure_gpa",
                "pressure_shift_gpa",
                "conditional_pressure_error_gpa",
                "pressure_temperature_covariance_gpa_k",
                "published_fe_residual_gpa",
            },
        )
    else:
        OUT.mkdir(parents=True, exist_ok=True)
        for filename, document in documents.items():
            (OUT / filename).write_text(document, encoding="utf-8")
    if args.plot:
        plot_results(summary)


if __name__ == "__main__":
    main()
