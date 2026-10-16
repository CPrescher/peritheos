"""Expose the scoped rubidium-halide diagnostic refits to the common ledger."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _comparison(name, stored, fitted, error):
    relative = abs(fitted - stored) / abs(stored)
    limit = 0.05 if name == "V0" else 0.15 if name == "K0" else 0.30
    similar = relative <= limit
    if name == "K0_prime":
        similar = abs(fitted - stored) <= max(1.0, abs(stored) * 0.20)
    if name in {"gamma0", "q"}:
        similar = relative <= 0.25
    return {
        "parameter": name,
        "published": stored,
        "published_error": error,
        "refit": fitted,
        "refit_error": None,
        "difference": fitted - stored,
        "relative_difference": relative,
        "within_reported_error": abs(fitted - stored) <= error,
        "within_combined_2sigma": None,
        "similar": similar,
    }


def rubidium_outcome(record):
    from scripts.investigate_rubidium_halides_farla_2025 import weighted_refit

    compound = record["identifier"].split("_")[0]
    fit = weighted_refit(compound)
    names = ["K0", "K0_prime", "gamma0", "q"]
    values = {**record["eos"]["parameters"], **record["thermal"]["parameters"]}
    errors = {**record["parameter_errors"], **record["thermal"]["parameter_errors"]}
    comparisons = [
        _comparison(name, values[name], value, errors[name])
        for name, value in zip(names, fit["fit_parameters"])
    ]
    if not fit["converged"] or not all(p["similar"] for p in comparisons):
        raise AssertionError("Rubidium conditional refit no longer supports similarity")
    report = json.loads(
        (ROOT / "docs/data/rubidium-halides-farla-2025-reproduction.json").read_text()
    )
    published_rmse = report["compounds"][compound]["independent_fits"][0][
        "source_pressure_rmse_gpa"
    ]
    qualification = (
        "Conditional numerical similarity, not source-exact coefficient or workflow "
        "parity. Iterated frozen effective-variance weights use diagonal printed "
        "P/V/T errors without a pressure floor; all DAC and heated LVP observations "
        "are included, cold LVP rows excluded as a documented analyst choice. "
        "Fixed V0/theta0 and source-reported CsCl pressures are retained. No full "
        "calibration-error covariance or original author row flags are recovered. "
        "The independent fit does not round to every published coefficient; source "
        "benchmarks and CsCl replay discrepancies remain separate qualifications."
    )
    return {
        "status": "similar",
        "parameters": comparisons,
        "dataset_identifiers": record["fit_datasets"],
        "observations": fit["rows"],
        "fit_kind": "iterated_effective_variance_bm3_mgd",
        "objective": "iterated frozen diagonal effective P/V/T variance",
        "published_rmse_gpa": published_rmse,
        "rmse_gpa": fit["pressure_rmse_gpa"],
        "solver_success": fit["converged"],
        "original_author_workflow_parity": False,
        "original_author_coefficient_rounding_parity": False,
        "qualification": qualification,
        "reason": qualification,
        "reproduction_report": "docs/data/rubidium-halides-farla-2025-follow-up.json",
    }
