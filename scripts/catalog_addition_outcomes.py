"""Scoped reproduction outcomes for the October 2026 catalog integration.

The dedicated scientific routines remain authoritative for their objectives and
source corrections. These adapters expose their results to the common ledger.
"""

from __future__ import annotations

import json
from functools import lru_cache
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


@lru_cache(maxsize=1)
def _diopside_reports():
    from scripts.diopside_check_plots import inspect_plots
    from scripts.diopside_reproduce import reproduce

    return reproduce(), inspect_plots()


def diopside_outcome(record):
    report, plots = _diopside_reports()
    thermal = record["identifier"] == "diopside_zhao_1998_thermal_bm3"
    if thermal:
        fit = plots["graphical_interpretation_diagnostic"]
        checks = plots["parity_assessment"]["parameter_checks"]
        comparisons = [
            _comparison(
                name,
                p["published_value"],
                p["independent_refit_value"],
                p["reported_error"],
            )
            for name, p in checks.items()
        ]
        dataset = record["recommended_refit_dataset"]
        rows = 32
        published_rmse = fit["published_curve_pressure_rmse_gpa"]
        rmse = fit["refit_pressure_rmse_gpa"]
        qualification = (
            "Conditional parameter parity within reported error widths using the "
            "separate corrected analysis dataset: row 15 pressure is the rounded "
            "Figure 2 reading 4.58 GPa, not the printed 6.95 GPa. Its +/-0.05 GPa "
            "graphical bound is not a statistical error. All 32 rows, including "
            "the repeated ambient row, are retained; V0 and Tr are fixed. Original "
            "source tables and published coefficients are preserved. Exact original "
            "fit inputs, selection, weighting and covariance remain unrecovered."
        )
    else:
        source = report["thompson_2008"]
        fit = source["diagnostic_refit"]
        comparisons = [
            _comparison(
                name, value, fit["coefficients"][name], record["parameter_errors"][name]
            )
            for name, value in record["eos"]["parameters"].items()
        ]
        dataset = record["fit_datasets"][0]
        rows = source["observations"]
        published_rmse = source["published_curve_all_rows"]["pressure_rmse_gpa"]
        rmse = fit["pressure_rmse_gpa"]
        qualification = (
            "Conditional parameter parity within reported error widths on all 14 "
            "printed rows using effective P/V variance. Source confidence and exact "
            "author weights/covariance are not recovered; no combined-two-sigma "
            "or original-author workflow parity is asserted."
        )
    if not all(p["within_reported_error"] and p["similar"] for p in comparisons):
        raise AssertionError("Diopside dedicated comparison no longer supports parity")
    return {
        "status": "parity",
        "parity_basis": "within_reported_parameter_errors",
        "parameters": comparisons,
        "dataset_identifiers": [dataset],
        "original_source_dataset_identifiers": record["fit_datasets"],
        "observations": rows,
        "fit_kind": "graph_corroborated_corrected_pvt"
        if thermal
        else "effective_variance_pv",
        "published_rmse_gpa": published_rmse,
        "rmse_gpa": rmse,
        "solver_success": True,
        "original_author_workflow_parity": False,
        "qualification": qualification,
        "reason": qualification,
        "reproduction_report": "docs/data/diopside-plot-check.json"
        if thermal
        else "docs/data/diopside-reproduction.json",
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
