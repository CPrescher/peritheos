"""Pd source-family audit retaining printed errors and a separate ambient datum.

No ambient pressure uncertainty is invented. BM3 treats the ambient volume as
a V0 observation. BM2 keeps the stored V0 fixed diagnostically, so the ambient
row contributes neither a residual nor a degree of freedom. The exact author
constraints and EoSFit7-GUI workflow are not established by these assumptions.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3

DATA = (
    Path(__file__).resolve().parents[1]
    / "peritheos/data/datasets/palladium-fedotenko-2020-table1-compression.csv"
)
RECORDS = {"palladium_fedotenko_2020_bm3_1", "palladium_fedotenko_2020_bm2_2"}
AUDIT = "docs/literature-reproductions/palladium-fedotenko-2020-audit.md"


def observations():
    rows = np.genfromtxt(DATA, delimiter=",", names=True)
    if (
        rows.size != 15
        or rows["pressure_gpa"][0] != 0
        or not np.isnan(rows["pressure_error_gpa"][0])
    ):
        raise ValueError("Unexpected Fedotenko ambient datum or row count")
    for name in ("pressure_error_gpa", "volume_error_a3_conventional_cell"):
        if not np.all(np.isfinite(rows[name][1:]) & (rows[name][1:] > 0)):
            raise ValueError(f"Fedotenko positive-pressure rows require {name}")
    return rows


def fit(record, rows):
    fixed_v0 = record["eos"]["type"] == "BM2"
    published = record["eos"]["parameters"]
    names = ["K0"] if fixed_v0 else ["V0", "K0", "K0_prime"]
    p, v = rows["pressure_gpa"], rows["volume_a3_conventional_cell"]
    sp, sv = rows["pressure_error_gpa"][1:], rows["volume_error_a3_conventional_cell"]
    n = len(names)

    def parameters(x):
        return {
            **published,
            "K0_prime": 4.0 if fixed_v0 else published["K0_prime"],
            **dict(zip(names, x[:n])),
        }

    def residual(x):
        values = parameters(x)
        return np.r_[
            (BM3(**values).pressure(x[n:]) - p[1:]) / sp,
            (x[n:] - v[1:]) / sv[1:],
            [] if fixed_v0 else [(values["V0"] - v[0]) / sv[0]],
        ]

    result = least_squares(
        residual,
        np.r_[[published[name] for name in names], v[1:]],
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-10,
        max_nfev=2000,
    )
    if not result.success:
        raise ValueError(f"Fedotenko audit fit failed: {result.message}")
    covariance = np.linalg.inv(result.jac.T @ result.jac)[:n, :n]
    values = parameters(result.x)
    dof = 14 + int(not fixed_v0) - n
    return SimpleNamespace(
        parameters=values,
        standard_errors=dict(zip(names, np.sqrt(np.diag(covariance)))),
        free_parameters=names,
        covariance=covariance,
        residuals=BM3(**values).pressure(v) - p,
        weighted_sum_squares=float(np.sum(result.fun**2)),
        degrees_of_freedom=dof,
        reduced_chi_square=float(np.sum(result.fun**2) / dof),
        success=bool(result.success),
        message=str(result.message),
    )


def ledger_outcome(record):
    from scripts.validate_primary_eos_refits import _compare

    rows = observations()
    result = fit(record, rows)
    status, comparisons = _compare(record, result, False, 1.0)
    fixed_v0 = record["eos"]["type"] == "BM2"
    source = {**record["eos"]["parameters"]}
    if fixed_v0:
        source["K0_prime"] = 4.0
    residuals = (
        BM3(**source).pressure(rows["volume_a3_conventional_cell"])
        - rows["pressure_gpa"]
    )
    outcome = {
        "status": status,
        "coefficient_comparison_status": status,
        "author_workflow_reproduced": False,
        "source_equation_status": "reproduced",
        "data_agreement_status": "poor"
        if fixed_v0
        else "consistent_with_printed_widths",
        "dataset_identifiers": ["palladium_fedotenko_2020_table1_compression"],
        "observations": 15,
        "selection": "all 15 Table 1 rows; ambient datum handled separately",
        "observed_pressure_range_gpa": [0.0, 31.2],
        "observed_volume_range": [52.411, 58.868],
        "observed_temperature_range_k": None,
        "columns": {
            "pressure": "pressure_gpa",
            "volume": "volume_a3_conventional_cell",
            "temperature": None,
        },
        "fit_kind": "isothermal_pv",
        "objective": "errors_in_variables_with_separate_ambient_datum",
        "weighting": {
            "positive_pressure_rows": 14,
            "pressure_sigma_gpa": 0.2,
            "volume_sigma_a3": 0.003,
            "ambient": "V0 fixed diagnostically at 58.868 A^3; ambient row contributes no residual or degree of freedom"
            if fixed_v0
            else "ambient volume observation 58.868 +/-0.002 A^3 at reference P=0; V0 free",
            "missing_ambient_pressure_sigma": "not imputed",
        },
        "absolute_sigma": True,
        "free_parameters": result.free_parameters,
        "fixed_parameters": ["V0", "K0_prime"] if fixed_v0 else [],
        "fixed_values": {"V0": source["V0"], "K0_prime": 4.0} if fixed_v0 else {},
        "parameters": comparisons,
        "parameter_covariance": result.covariance.tolist(),
        "covariance_convention": "inverse J^T J; absolute independent coordinate widths; no residual scaling",
        "degrees_of_freedom": result.degrees_of_freedom,
        "weighted_sum_squares": result.weighted_sum_squares,
        "reduced_chi_square": result.reduced_chi_square,
        "published_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "rmse_gpa": float(np.sqrt(np.mean(result.residuals**2))),
        "published_max_abs_pressure_residual_gpa": float(np.max(np.abs(residuals))),
        "pressure_rmse_coordinates": "measured volumes for both curves; all 15 rows",
        "solver_success": result.success,
        "solver_message": result.message,
        "reproduction_report": AUDIT,
        "qualification": "All 14 positive-pressure rows retain their printed 0.2 GPa pressure and 0.003 A^3 volume widths. No uncertainty is invented for ambient pressure. Printed widths are interpreted conditionally as independent absolute standard deviations; source confidence and covariance conventions remain unknown. The exact author constraints and EoSFit7-GUI workflow are not recovered.",
    }
    if fixed_v0:
        outcome.update(
            reproduction_status=status,
            usage_recommendation="not_recommended_for_quantitative_use",
            reason="Poor data agreement despite coefficient parity: the published fixed-K'=4 curve has 0.934246 GPa pressure RMS and 1.438559 GPa maximum residual on Table 1. The correctly weighted diagnostic refit remains poor (reduced chi-square about 25 under independent absolute-error assumptions). Published coefficients are preserved; agreement of K0 does not validate the curve.",
        )
        outcome["qualification"] += (
            " V0 fixed at 58.868 A^3 is a diagnostic constraint, not a documented author choice. Using Table 2's isolated 58.863 A^3 does not resolve the misfit."
        )
    return outcome
