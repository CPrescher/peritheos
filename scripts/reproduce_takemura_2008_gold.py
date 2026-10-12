"""Conditional Takemura Au parameter parity from all 86 observed Table III rows.

Run with ``--check`` to verify the dedicated report. This script never writes
material records or shared ledgers. Equal weights are a diagnostic assumption.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/takemura-2008-gold-reproduction.json"
RECORD = "gold_takemura_2008_vinet_6"
DATASET = "gold_takemura_2008_table3_compression"
V0 = 4 * 16.962  # Published atomic reference volume, normalized to the fcc cell.
K0 = 167.0
PRINTED_WIDTH = 0.1  # No confidence convention is assigned to the source width.
QUALIFICATION = (
    "Conditional numerical parameter parity: all 86 observed Table III rows, "
    "with published V0/K0 fixed, recover K0-prime within the printed 0.1 width "
    "on both ruby coordinates under equal-pressure and equal-volume residual "
    "objectives. Original residual direction, weights, unrounded reductions, "
    "pointwise uncertainties, confidence convention and covariance are not "
    "recovered; no combined-two-sigma or exact author-fit claim. See the "
    "[dedicated audit](literature-reproductions/takemura-2008-gold.md)."
)


def pressure(volume, derivative):
    """Independent published Vinet equation in GPa and conventional cell Å³."""
    x = (np.asarray(volume) / V0) ** (1 / 3)
    return 3 * K0 * (1 - x) / x**2 * np.exp(1.5 * (derivative - 1) * (1 - x))


def volume(pressures, derivative):
    """Independent monotone compression-branch inversion, including P=0."""
    return np.array(
        [
            V0
            if p == 0
            else V0
            * brentq(
                lambda x: pressure(V0 * x**3, derivative) - p,
                0.4,
                1.0,
                xtol=1e-14,
            )
            ** 3
            for p in pressures
        ]
    )


def _metrics(v, p, derivative):
    residual = pressure(v, derivative) - p
    return {
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "pressure_max_abs_residual_gpa": float(np.max(abs(residual))),
        "cell_volume_rmse_angstrom3": float(
            np.sqrt(np.mean((volume(p, derivative) - v) ** 2))
        ),
    }


def _fit(v, p, published, objective):
    if objective == "equal_pressure_residuals":

        def residual(beta):
            return pressure(v, beta[0]) - p

        x = (v / V0) ** (1 / 3)

        def jacobian(beta):
            return (pressure(v, beta[0]) * 1.5 * (1 - x))[:, None]
    else:

        def residual(beta):
            return volume(p, beta[0]) - v

        jacobian = "3-point"
    result = least_squares(
        residual,
        [published],
        jac=jacobian,
        bounds=(3.0, 9.0),
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
    )
    if not result.success:
        raise RuntimeError(result.message)
    derivative = float(result.x[0])
    return {
        "objective": objective,
        "K0_prime": derivative,
        "difference_from_published": derivative - published,
        "within_reported_error": abs(derivative - published) <= PRINTED_WIDTH,
        "refit_error": None,
        "reduced_chi_square": None,
        "solver_success": bool(result.success),
        **_metrics(v, p, derivative),
    }


def reproduce():
    document = get_material_document("gold")
    material = Material.from_eosmat(document, record_identifiers=[RECORD])
    data = material.get_dataset(DATASET)
    if not data.checksum_verified or len(data) != 86:
        raise ValueError("Expected all 86 checksum-verified source observations")
    v = data["lattice_111_angstrom"] ** 3
    branches = {}
    for name, column, published in (
        ("dorogokupets_oganov", "ruby_pressure_dorogokupets_gpa", 5.9),
        ("zha", "ruby_pressure_zha_gpa", 5.5),
    ):
        p = data[column]
        branches[name] = {
            "pressure_column": column,
            "observed_pressure_range_gpa": [float(p.min()), float(p.max())],
            "published_K0_prime": published,
            "published_error_width": PRINTED_WIDTH,
            "published_curve": _metrics(v, p, published),
            "fits": {
                objective: _fit(v, p, published, objective)
                for objective in ("equal_pressure_residuals", "equal_volume_residuals")
            },
        }
    parity = all(
        fit["within_reported_error"]
        for branch in branches.values()
        for fit in branch["fits"].values()
    )
    return {
        "record_identifier": RECORD,
        "source_doi": "10.1103/PhysRevB.78.104119",
        "dataset_identifier": DATASET,
        "source_table_sha256": data.resource.sha256,
        "observations": len(data),
        "input_kind": "observed_pressure_and_111_lattice_parameter",
        "cell_volume_convention": "am(111)^3; four atoms per conventional fcc cell",
        "fixed_parameters": {"V0": V0, "K0": K0},
        "free_parameters": ["K0_prime"],
        "observed_cell_volume_range_angstrom3": [float(v.min()), float(v.max())],
        "branches": branches,
        "status": "parity" if parity else "parity_not_achieved",
        "parity_basis": "within_reported_parameter_errors",
        "qualification": QUALIFICATION,
        "author_fit_reproduced": False,
        "published_coefficients_changed": False,
    }


def ledger_outcome(record):
    """Record-specific conditional result for the common diagnostic ledger."""
    assert record["eos"]["parameters"] == {"V0": V0, "K0": K0, "K0_prime": 5.9}
    report = reproduce()
    branch = report["branches"]["dorogokupets_oganov"]
    fit = branch["fits"]["equal_pressure_residuals"]
    derivative = fit["K0_prime"]
    return {
        "status": report["status"],
        "parity_basis": report["parity_basis"],
        "qualification": report["qualification"],
        "author_fit_reproduced": False,
        "dataset_identifiers": [DATASET],
        "source_table_sha256": report["source_table_sha256"],
        "observations": report["observations"],
        "selection": "all 86 observed rows; no exclusions",
        "observed_pressure_range_gpa": branch["observed_pressure_range_gpa"],
        "experimental_pressure_range_gpa": record["experimental_pressure_range_gpa"],
        "observed_volume_range": report["observed_cell_volume_range_angstrom3"],
        "observed_temperature_range_k": None,
        "columns": {
            "pressure": branch["pressure_column"],
            "volume": "lattice_111_angstrom",
            "temperature": None,
        },
        "fit_kind": "isothermal_pv",
        "objective": fit["objective"],
        "absolute_sigma": False,
        "free_parameters": ["K0_prime"],
        "parameters": [
            {
                "parameter": "K0_prime",
                "published": 5.9,
                "published_error": PRINTED_WIDTH,
                "refit": derivative,
                "refit_error": None,
                "difference": derivative - 5.9,
                "relative_difference": abs(derivative - 5.9) / 5.9,
                "within_combined_2sigma": None,
                "within_reported_error": fit["within_reported_error"],
                "similar": abs(derivative - 5.9) <= 1.0,
            }
        ],
        "published_rmse_gpa": branch["published_curve"]["pressure_rmse_gpa"],
        "rmse_gpa": fit["pressure_rmse_gpa"],
        "reduced_chi_square": None,
        "degrees_of_freedom": report["observations"] - 1,
        "solver_success": fit["solver_success"],
        "solver_message": "dedicated equal-pressure refit; volume/scale sensitivities archived",
        "reproduction_report": "docs/data/takemura-2008-gold-reproduction.json",
    }


def check_saved(saved, current):
    if isinstance(current, dict):
        assert saved.keys() == current.keys()
        for key in current:
            check_saved(saved[key], current[key])
    elif isinstance(current, list):
        assert len(saved) == len(current)
        for old, new in zip(saved, current):
            check_saved(old, new)
    elif isinstance(current, float):
        np.testing.assert_allclose(saved, current, rtol=1e-7, atol=1e-9)
    else:
        assert saved == current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    if args.check:
        check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), report)
    else:
        OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(report["status"], report["parity_basis"])


if __name__ == "__main__":
    main()
