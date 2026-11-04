"""Independent fits to the original published MgSiO3 observations and scales.

No Dorogokupets pressure rescaling, point weighting, or fitted synthetic data.
The source papers do not specify a fully reproducible residual/weighting rule;
unweighted pressure residuals are an explicit reconstruction choice, not a claim
about their original optimizer. Published parameters remain the catalog values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import least_squares
from scipy.stats import t as student_t

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/wang-zhou-komabayashi-mgsio3-refit.json"
WANG = "akimotoite-wang-2004-tables1-2-pvt.csv"
ZHOU = "akimotoite-zhou-2014-table1-elasticity.csv"
KOMABAYASHI = "mgsio3-post-perovskite-komabayashi-2008-table1-pvt.csv"
IDS = (
    "akimotoite_wang_2004_t0133_bm3",
    "akimotoite_wang_2004_t0150_bm3",
    "akimotoite_zhou_2014_300k_bm3",
    "mgsio3_post_perovskite_komabayashi_2008_thermal_bm3",
)
WANG_T0150 = IDS[1]


def static_pressure_uncertainty(volume, parameters, residuals):
    """Local OLS covariance, conditional on exact volumes and fixed K0=210.

    Estimate the common pressure-residual variance with n-2 degrees of freedom.
    The analytic BM3 Jacobian avoids optimizer finite-difference step dependence.
    This is a reconstruction error model, not the source's unpublished protocol.
    """
    volume = np.asarray(volume, dtype=float)
    residuals = np.asarray(residuals, dtype=float)
    v0, kp = parameters
    dof = len(volume) - 2
    if dof <= 0 or volume.shape != residuals.shape:
        raise ValueError("Covariance requires matching arrays and n > 2")
    x = (v0 / volume) ** (1 / 3)
    a = x**7 - x**5
    b = 1 + 0.75 * (kp - 4) * (x**2 - 1)
    jacobian = np.column_stack(
        (
            1.5
            * 210
            / (3 * v0)
            * ((7 * x**7 - 5 * x**5) * b + 1.5 * (kp - 4) * a * x**2),
            1.125 * 210 * a * (x**2 - 1),
        )
    )
    _, singular, vt = np.linalg.svd(jacobian, full_matrices=False)
    if singular[-1] <= np.finfo(float).eps * max(jacobian.shape) * singular[0]:
        raise ValueError("Cannot estimate covariance for a rank-deficient fit")
    variance = float(residuals @ residuals / dof)
    covariance = (vt.T / singular**2) @ vt * variance
    errors = np.sqrt(np.diag(covariance))
    return dict(
        parameters=["V0", "K0_prime"],
        parameter_units=["angstrom^3 per conventional cell (Z=6)", "1"],
        matrix=covariance.tolist(),
        standard_errors=dict(zip(("V0", "K0_prime"), map(float, errors))),
        correlation_matrix=(covariance / np.outer(errors, errors)).tolist(),
        degrees_of_freedom=dof,
        residual_variance_gpa2=variance,
        residual_standard_error_gpa=float(np.sqrt(variance)),
        local_student_t_95_intervals={
            name: [
                float(value - student_t.ppf(0.975, dof) * error),
                float(value + student_t.ppf(0.975, dof) * error),
            ]
            for name, value, error in zip(("V0", "K0_prime"), parameters, errors)
        },
        convention="s^2 (J^T J)^-1; s^2=sum(pressure residual^2)/(n-2); analytic BM3 Jacobian; errors are local 1-standard-error estimates, not published errors",
        assumptions=[
            "Independent, equal-variance pressure residuals; variance estimated from residual scatter, not supplied pressure errors.",
            "Observed volumes are conditioned on; printed volume errors are not weights in this pressure-residual reconstruction.",
            "K0=210 GPa is fixed; its uncertainty, shared pressure-calibration systematics and correlations between observations are not propagated.",
            "Source confidence convention and original weighting/covariance are unknown; no formal combined-two-sigma comparison is asserted.",
            f"Student-t intervals use {dof} residual degrees of freedom and a local linear approximation; they are conditional, marginal intervals.",
        ],
    )


def wang_t0150_volume_weighting_diagnostic(initial):
    """Check the printed V widths without inventing missing pressure errors."""
    selected = [
        r
        for r in rows(WANG)
        if r["run"] == "T0150" and float(r["temperature_k"]) == 298
    ]
    v, p, widths = np.array(
        [
            (
                float(r["unit_cell_volume_a3"]),
                float(r["pressure_original_gpa"]),
                float(r["unit_cell_volume_uncertainty_a3"]),
            )
            for r in selected
        ]
    ).T
    fit = solve(lambda x: (BM3(x[0], 210, x[1]).volume(p) - v) / widths, initial)
    return dict(
        objective="sum(((V_model(P)-V_observed)/printed_V_error)^2)",
        parameters=dict(zip(("V0", "K0_prime"), map(float, fit.x))),
        weighted_sum_squares=float(fit.fun @ fit.fun),
        degrees_of_freedom=len(v) - 2,
        weighted_sum_squares_per_dof=float(fit.fun @ fit.fun / (len(v) - 2)),
        solver_success=bool(fit.success),
        qualification="Sensitivity diagnostic only; source V error confidence is unspecified and P errors are unavailable. The large standardized scatter shows that the printed V widths alone cannot explain the data under exact P. This is not a recovered author objective or a replacement catalog EOS.",
    )


def rows(filename):
    with (DATA / filename).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def observations():
    """Explicit source masks; retain every selected row, including ambient rows."""
    wang = rows(WANG)
    result = {}
    for identifier, run in zip(IDS[:2], ("T0133", "T0150")):
        selected = [
            r for r in wang if r["run"] == run and float(r["temperature_k"]) == 298
        ]
        result[identifier] = [
            (
                int(r["source_order"]),
                float(r["unit_cell_volume_a3"]),
                float(r["pressure_original_gpa"]),
                298.0,
            )
            for r in selected
        ]
    result[IDS[2]] = [
        (
            int(r["source_order"]),
            6 * 100.387 * 1e24 / Avogadro / float(r["density_g_cm3"]),
            float(r["pressure_gpa"]),
            300.0,
        )
        for r in rows(ZHOU)
        if float(r["temperature_k"]) == 300
    ]
    result[IDS[3]] = [
        (
            int(r["source_order"]),
            float(r["ppv_volume_a3"]),
            float(r["pressure_mgo_gpa"]),
            float(r["temperature_k"]),
        )
        for r in rows(KOMABAYASHI)
        if r["ppv_fit_included"] == "true"
    ]
    return result


def solve(residual, initial):
    fit = least_squares(
        residual, initial, xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=10000
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    return fit


def thermal_model(k0, alpha1, dkdt):
    # Eq. 4: alpha(T)=alpha0+alpha1*T; alpha(300) is constrained, not alpha0.
    return ThermalReferenceStateEOS(
        BM3(163.813, k0, 4),
        300,
        1.7e-5 - 300 * alpha1,
        dkdt,
        alpha1,
        thermal_expansion_law="linear_temperature",
    )


def comparison(published, errors, fitted):
    return [
        dict(
            parameter=name,
            published=value,
            published_error=errors[name],
            refit=fitted[name],
            refit_error=None,
            difference=fitted[name] - value,
            relative_difference=abs(fitted[name] - value) / abs(value),
            within_reported_error=bool(abs(fitted[name] - value) <= errors[name]),
            within_combined_2sigma=None,
        )
        for name, value in published.items()
    ]


def reproduce():
    results = {}
    for identifier, data in observations().items():
        array = np.array(data)
        v, p, t = array[:, 1:].T
        stages = []
        sensitivity = None
        if identifier in IDS[:2]:
            published = (
                dict(V0=264.2, K0_prime=4.8)
                if identifier == IDS[0]
                else dict(V0=263.9, K0_prime=5.6)
            )
            errors = dict(V0=0.2, K0_prime=0.5 if identifier == IDS[0] else 0.8)
            fit = solve(lambda x: BM3(x[0], 210, x[1]).pressure(v) - p, [264, 5])
            fitted = dict(zip(published, map(float, fit.x)))
            published_model = BM3(published["V0"], 210, published["K0_prime"])
            predicted = published_model.pressure(v)
            fixed = dict(K0=210)
            selection = "All 298 K rows in " + (
                "Table 1, run T0133 (Au)"
                if identifier == IDS[0]
                else "Table 2, run T0150 (NaCl)"
            )
        elif identifier == IDS[2]:
            published, errors = dict(V0=262.45, K0=207), dict(V0=0.26, K0=3)
            fit = solve(lambda x: BM3(x[0], x[1], 4.6).pressure(v) - p, [262, 210])
            fitted = dict(zip(published, map(float, fit.x)))
            predicted = BM3(262.45, 207, 4.6).pressure(v)
            fixed = dict(K0_prime=4.6)
            selection = "All 17 akimotoite 300 K density rows, including three ambient density-only rows; V=6*M/(N_A*rho), M=100.387 g/mol"
        else:
            rt = t == 300
            room = solve(lambda x: BM3(163.813, x[0], 4).pressure(v[rt]) - p[rt], [220])
            k0 = float(room.x[0])

            def thermal_fit(k):
                return solve(
                    lambda x: (
                        thermal_model(k, x[0] * 1e-8, x[1] * 0.01).pressure(v, t) - p
                    ),
                    [0.1, -0.8],
                )

            fit = thermal_fit(k0)
            alpha1, dkdt = float(fit.x[0] * 1e-8), float(fit.x[1] * 0.01)
            published = dict(K0=223.2, alpha1=1.13e-9, dK_dT=-0.0085)
            errors = dict(K0=0.2, alpha1=1.14e-8, dK_dT=0.0011)
            fitted = dict(K0=k0, alpha1=alpha1, dK_dT=dkdt)
            fixed = dict(V0=163.813, K0_prime=4, Tr=300, alpha_at_300k=1.7e-5)
            stages = [
                dict(
                    stage="300 K BM3",
                    observations=int(rt.sum()),
                    free_parameters=["K0"],
                    fixed_parameters=["V0", "K0_prime"],
                    rmse_gpa=float(np.sqrt(np.mean(room.fun**2))),
                ),
                dict(
                    stage="P-V-T",
                    observations=len(p),
                    free_parameters=["alpha1", "dK_dT"],
                    fixed_parameters=["V0", "K0", "K0_prime", "Tr", "alpha_at_300k"],
                    dependent_parameters={"alpha0": "1.7e-5 - 300*alpha1"},
                ),
            ]
            rounded = thermal_fit(223.2)
            sensitivity = dict(
                description="Second stage with the published rounded K0=223.2 GPa fixed",
                alpha1=float(rounded.x[0] * 1e-8),
                dK_dT=float(rounded.x[1] * 0.01),
            )
            # Evaluate the actual printed coefficients, including rounded alpha0.
            predicted = ThermalReferenceStateEOS(
                BM3(163.813, 223.2, 4),
                300,
                1.667e-5,
                -0.0085,
                1.13e-9,
                thermal_expansion_law="linear_temperature",
            ).pressure(v, t)
            selection = "All 21 Table 1 rows with PPv volume, using original Speziale-MgO pressures; the Pv-only row is excluded"
        results[identifier] = dict(
            status="similar",
            observations=len(p),
            selection=selection,
            source_orders=[int(r[0]) for r in data],
            observed_pressure_range_gpa=[float(p.min()), float(p.max())],
            observed_temperature_range_k=[float(t.min()), float(t.max())],
            objective="unweighted pressure residuals",
            fixed_values=fixed,
            free_parameters=list(published),
            parameters=comparison(published, errors, fitted),
            rmse_gpa=float(np.sqrt(np.mean(fit.fun**2))),
            published_rmse_gpa=float(np.sqrt(np.mean((predicted - p) ** 2))),
            solver_success=bool(fit.success),
            stages=stages,
            qualification="All free coefficients agree within the printed source errors. Exact optimizer, weighting, covariance and confidence convention are not recovered; common-ledger classification is similar, not strict statistical parity.",
        )
        if sensitivity is not None:
            results[identifier].update(
                rounded_K0_sensitivity=sensitivity,
                derived_alpha0=1.7e-5 - 300 * fitted["alpha1"],
            )
        if identifier == WANG_T0150:
            uncertainty = static_pressure_uncertainty(v, fit.x, fit.fun)
            outcome = results[identifier]
            for parameter in outcome["parameters"]:
                parameter["refit_error"] = uncertainty["standard_errors"][
                    parameter["parameter"]
                ]
            outcome.update(
                status="parity"
                if all(p["within_reported_error"] for p in outcome["parameters"])
                else "parity_not_achieved",
                parity_basis="within_reported_parameter_errors",
                weighting="Equal weight per pressure residual; absolute_sigma=False; residual variance estimated with n-2 degrees of freedom.",
                absolute_sigma=False,
                degrees_of_freedom=uncertainty["degrees_of_freedom"],
                parameter_covariance=uncertainty["matrix"],
                covariance_convention=uncertainty["convention"],
                covariance={
                    **uncertainty,
                    "volume_weighting_sensitivity": wang_t0150_volume_weighting_diagnostic(
                        fit.x
                    ),
                },
                qualification="Conditional numerical parameter parity: both free coefficients agree within the printed source error widths, with reconstructed local standard errors comparable to those widths. Covariance uses residual-scaled unweighted pressure least squares with K0=210 GPa fixed. The source confidence convention, original weighting and covariance remain unknown; no formal combined-two-sigma or author-workflow reproduction is claimed. See docs/literature-reproductions/wang-zhou-komabayashi-mgsio3.md#t0150-uncertainty-reconstruction for the uncertainty budget.",
            )
    return dict(
        format="peritheos.wang-zhou-komabayashi-mgsio3-refit",
        format_version=1,
        generated_with="scripts/reproduce_wang_zhou_komabayashi_mgsio3.py",
        source_sha256={
            f: hashlib.sha256((DATA / f).read_bytes()).hexdigest()
            for f in (WANG, ZHOU, KOMABAYASHI)
        },
        density_conversion=dict(
            molar_mass_g_mol=100.387, formula_units=6, avogadro_mol_inverse=Avogadro
        ),
        fits=results,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--record",
        choices=IDS,
        help="Update/check only this fit, retaining the other saved results",
    )
    args = parser.parse_args()
    report = reproduce()
    if args.record:
        saved = json.loads(REPORT.read_text())
        if saved["source_sha256"] != report["source_sha256"]:
            raise SystemExit(
                "Source datasets changed; run the complete reconstruction before a selective update"
            )
        report["fits"] = {**saved["fits"], args.record: report["fits"][args.record]}
    if not all(
        p["within_reported_error"]
        for f in report["fits"].values()
        for p in f["parameters"]
    ):
        raise AssertionError("A recovered coefficient is outside its source error")
    # Retain full precision in the machine-readable diagnostics.
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.check:
        if json.loads(REPORT.read_text()) != report:
            raise SystemExit(f"Stale report: {REPORT}")
    else:
        REPORT.write_text(text, encoding="utf-8")
    print(json.dumps({k: v["status"] for k, v in report["fits"].items()}))


if __name__ == "__main__":
    main()
