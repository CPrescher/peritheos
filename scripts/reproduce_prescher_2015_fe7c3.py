"""Reconstruct the archived Fe7C3 fit on its author-supplied 28-row input.

Repository-only audit. The effective-variance objective and residual-scaled
covariance reproduce the saved export, but the original solver implementation
and source confidence convention are not established. Table S4 is archival only.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3

ROOT = Path(__file__).resolve().parents[1]
DATASET_ID = "iron_carbide_fe7c3_prescher_2015_author_fit_input"
RECORD_ID = "iron_carbide_fe7c3_prescher_2015_bm3_1"
DATA = ROOT / "peritheos/data/datasets/iron-carbide-fe7c3-prescher-2015-author-fit.csv"
EXPORT = (
    ROOT
    / "peritheos/data/datasets/iron-carbide-fe7c3-prescher-2015-author-fit-export.txt"
)
REPORT = "docs/literature-reproductions/iron-carbide-fe7c3-prescher-2015.md"
QUALIFICATION = (
    "Numerical parameter parity on the recovered 28-row author input, with all "
    "six archived coefficient/error values reproduced at their printed "
    "precision. Effective pressure-variance weighting and residual-scaled "
    "covariance are reconstructed interpretations; the original solver and "
    "confidence convention remain unrecovered. Pressure widths include "
    "before/after half-differences, not established standard deviations. "
    "Table S4 contains different cell values and is retained for provenance "
    "only: it should not be used for scientific calculations, EOS fitting or "
    "pressure calibration. The affected processing logs contain revised cells "
    "alongside inconsistent reduction/export cells. This does not establish "
    "the physically preferable reduction or exact author-workflow identity."
)


def observations():
    rows = np.genfromtxt(DATA, delimiter=",", names=True, encoding="utf-8")
    if rows.size != 28:
        raise ValueError("Prescher author input must contain all 28 observations")
    for name in ("pressure_gpa", "unit_cell_volume_a3"):
        if not np.all(np.isfinite(rows[name]) & (rows[name] > 0)):
            raise ValueError(f"Invalid author-input coordinate {name}")
    if not np.all(rows["pressure_uncertainty_gpa"] >= 0) or not np.all(
        np.isfinite(rows["pressure_uncertainty_gpa"])
    ):
        raise ValueError("Pressure widths must be finite and nonnegative")
    if not np.all(
        np.isfinite(rows["unit_cell_volume_uncertainty_a3"])
        & (rows["unit_cell_volume_uncertainty_a3"] > 0)
    ):
        raise ValueError("Volume widths must be finite and positive")
    return rows


def fit(rows):
    p, v = rows["pressure_gpa"], rows["unit_cell_volume_a3"]
    sp, sv = (
        rows["pressure_uncertainty_gpa"],
        rows["unit_cell_volume_uncertainty_a3"],
    )
    names = ["V0", "K0", "K0_prime"]

    def model(values):
        return BM3(**dict(zip(names, values)))

    def objective(values):
        eos = model(values)
        # dP/dV = -K/V. A printed zero pressure width is preserved, while the
        # positive volume width keeps the propagated denominator finite.
        sigma = np.sqrt(sp**2 + (eos.bulk_modulus(v) / v * sv) ** 2)
        return (eos.pressure(v) - p) / sigma

    optimization = least_squares(
        objective,
        [748.0, 168.0, 6.1],
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=2000,
    )
    if not optimization.success:
        raise ValueError(f"Prescher author-input replay failed: {optimization.message}")
    dof = len(p) - len(names)
    chi_square = float(np.sum(optimization.fun**2))
    covariance = np.linalg.inv(optimization.jac.T @ optimization.jac) * (
        chi_square / dof
    )
    return SimpleNamespace(
        parameters=dict(zip(names, optimization.x)),
        standard_errors=dict(zip(names, np.sqrt(np.diag(covariance)))),
        free_parameters=names,
        covariance=covariance,
        residuals=np.asarray(model(optimization.x).pressure(v)) - p,
        weighted_sum_squares=chi_square,
        degrees_of_freedom=dof,
        reduced_chi_square=chi_square / dof,
        success=bool(optimization.success),
        message=str(optimization.message),
    )


def ledger_outcome(record):
    from scripts.validate_primary_eos_refits import _compare

    rows = observations()
    result = fit(rows)
    status, comparisons = _compare(record, result, False, 1.0)
    source = BM3(**record["eos"]["parameters"])
    p, v = rows["pressure_gpa"], rows["unit_cell_volume_a3"]
    archived = [float(x) for x in EXPORT.read_text().splitlines()[1].split("\t")[3:9]]
    reconstructed = [
        value
        for name in result.free_parameters
        for value in (result.parameters[name], result.standard_errors[name])
    ]
    export_match = all(
        abs(actual - target) < 0.005 for actual, target in zip(reconstructed, archived)
    )
    if not export_match:
        raise ValueError("Author fit export is not reproduced at printed precision")
    return {
        "status": status,
        "coefficient_comparison_status": status,
        "parity_basis": "within_reported_parameter_errors",
        "author_workflow_reproduced": False,
        "author_fit_export_reproduced_at_printed_precision": export_match,
        "dataset_identifiers": [DATASET_ID],
        "observations": len(rows),
        "selection": "all 28 archived author EOS inputs; series a new and series d new, excluding series b",
        "observed_pressure_range_gpa": [float(np.min(p)), float(np.max(p))],
        "observed_volume_range": [float(np.min(v)), float(np.max(v))],
        "observed_temperature_range_k": None,
        "columns": {
            "pressure": "pressure_gpa",
            "volume": "unit_cell_volume_a3",
            "temperature": None,
        },
        "fit_kind": "isothermal_pv",
        "objective": "effective_variance_pressure_residuals",
        "weighting": {
            "denominator": "sqrt(sP^2 + (K(V)/V * sV)^2); updated with parameters",
            "zero_pressure_widths": "preserved as printed; positive volume widths keep the denominator finite",
            "source_width_semantics": "pressure widths include before/after half-differences; confidence/cross-row covariance unknown",
        },
        "absolute_sigma": False,
        "free_parameters": result.free_parameters,
        "fixed_parameters": [],
        "parameters": comparisons,
        "parameter_covariance": result.covariance.tolist(),
        "covariance_convention": "inverse(J.T J) times weighted residual sum of squares / 25; reconstructed convention",
        "degrees_of_freedom": result.degrees_of_freedom,
        "weighted_sum_squares": result.weighted_sum_squares,
        "reduced_chi_square": result.reduced_chi_square,
        "published_rmse_gpa": float(np.sqrt(np.mean((source.pressure(v) - p) ** 2))),
        "rmse_gpa": float(np.sqrt(np.mean(result.residuals**2))),
        "pressure_rmse_coordinates": "recorded author-input volumes for both curves; all 28 rows",
        "solver_success": result.success,
        "solver_message": result.message,
        "reproduction_report": REPORT,
        "qualification": QUALIFICATION,
    }
