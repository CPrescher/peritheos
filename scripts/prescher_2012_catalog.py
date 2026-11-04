"""Repository-only catalog mapping of the separately preserved Fe3C recovery.

Published K coefficients/errors are never replaced by reconstructed values.
V0 is an explicit adoption from the drawn full-range curve, not a source value.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from scripts.recover_prescher_2012_fe3c import SOURCE, VECTORS, bm3, fit

ROOT = Path(__file__).resolve().parents[1]
V0 = 154.5678863704812
DOI = "10.1103/PhysRevB.85.140402"
REPORT = "docs/literature-reproductions/prescher-2012-recovery.md"
PREFIX = "cementite_prescher_2012_"
IDS = {"full": PREFIX + "full_range_bm3"}
RETIRED_IDS = {PREFIX + "below_22_gpa_bm3", PREFIX + "22_47_gpa_bm3"}
DATASETS = {"full": PREFIX + "figure1_digitized"}
REFERENCE = {
    "authors": [
        "Prescher",
        "Dubrovinsky",
        "McCammon",
        "Glazyrin",
        "Nakajima",
        "Kantor",
        "Merlini",
        "Hanfland",
    ],
    "year": 2012,
    "title": "Structurally hidden magnetic transitions in Fe3C at high pressures",
    "source": "Physical Review B",
    "volume": "85",
    "locator": "140402(R)",
    "doi": DOI,
}
COMMON = (
    "Conditional numerical parity reproduces K0, K0_prime and both error values "
    "at published printed precision, with V0 held fixed. V0=154.5678863704812 A^3 "
    "per conventional Fe3C cell (Z=4) is adopted from the Figure 1 full-range "
    "calculated curve; it is not a printed or author-confirmed reference volume. "
    "Published K coefficients and errors are preserved. The source confidence "
    "convention, original fit inputs, weighting and fixed/free V0 choice remain "
    "unconfirmed. Matching conditional residual-based errors does not establish "
    "the original experimental uncertainty calculation, author-workflow "
    "reproduction or a physical-validation/recommendation claim."
)
QUALIFICATIONS = {
    "full": (
        "The reconstruction uses all 22 digitized Figure 1 markers, including the "
        "ambient marker, with equal pressure weights and adopted V0 fixed. K0 and "
        "K0_prime are free. Covariance is inverse(J^T J) times SSE/(22-2); its errors "
        "condition on exact V0 and omit V0 uncertainty. The result is K0=160.757426 "
        "+/-2.353365 GPa and K0_prime=5.942125 +/-0.219843, which round to 161+/-2 "
        "and 5.9+/-0.2. Digitization/calibration errors are correlated and original "
        "experimental errors are unknown. The 19 final Table S1 rows remain a "
        "separate complete numeric dataset; three graph markers are absent from it. "
        "Published low/high split fits are retained only in the research audit, "
        "not offered as public EOS records."
    )
}


def selected_rows(side="full"):
    """Return all graph markers without modifying their recovered coordinates."""
    if side != "full":
        raise ValueError(
            "Split fits are audit-only; the catalog exposes the full range"
        )
    rows = json.loads(VECTORS.read_text())["markers"]
    return (
        np.array([r["pressure_gpa"] for r in rows]),
        np.array([r["volume_a3"] for r in rows]),
    )


def ledger_outcome(record):
    """Replay the full-range fixed-V0 coefficient and printed-error comparison."""
    if record["identifier"] != IDS["full"]:
        raise ValueError("Only the full-range Prescher EOS is in the catalog")
    p, v = selected_rows()
    result = fit(p, v, fixed_v0=V0)
    published, errors = SOURCE["full_range"]
    comparisons = [
        {
            "parameter": "V0",
            "published": V0,
            "published_error": None,
            "refit": V0,
            "refit_error": None,
            "difference": 0.0,
            "relative_difference": 0.0,
            "compared_for_parity": False,
            "role": "Adopted curve-derived reference held fixed; not a published measurement",
        }
    ]
    for i, key in enumerate(("K0", "K0_prime")):
        value = result["parameters"][key]
        error = result["residual_scaled_errors"][i]
        tolerance = (0.5, 0.05)[i]
        matched = abs(value - published[i]) < tolerance
        error_matches = abs(error - errors[i]) < tolerance
        if not (matched and error_matches):
            raise AssertionError(f"Prescher printed agreement diverged: {key}")
        comparisons.append(
            {
                "parameter": key,
                "published": published[i],
                "published_error": errors[i],
                "refit": value,
                "refit_error": error,
                "difference": value - published[i],
                "relative_difference": abs(value - published[i]) / abs(published[i]),
                "within_published_printing_precision": matched,
                "similar": matched,
                "within_combined_2sigma": None,
                "uncertainty_matches_printed_precision": error_matches,
                "uncertainty_finding": "matched_at_printed_precision",
                "display_decimal_places": i,
            }
        )
    return {
        "status": "parity",
        "parity_basis": "conditional_coefficients_and_errors_at_published_precision_fixed_V0",
        "comparison_scope": "parameter_agreement",
        "conditional": True,
        "compared_parameters": ["K0", "K0_prime"],
        "author_workflow_reproduced": False,
        "dataset_identifiers": [DATASETS["full"]],
        "fit_kind": "isothermal_pv",
        "objective": "pressure_residuals",
        "selection": "All 22 Figure 1 digitized markers, including ambient",
        "columns": {
            "pressure": "pressure_gpa",
            "volume": "volume_a3",
            "temperature": None,
        },
        "weighting": "Equal pressure weights; no experimental sigmas recovered",
        "absolute_sigma": False,
        "free_parameters": ["K0", "K0_prime"],
        "fixed_parameters": ["V0"],
        "fixed_values": {"V0": V0},
        "adopted_catalog_V0": V0,
        "parameters": comparisons,
        "observations": len(p),
        "observed_pressure_range_gpa": [float(min(p)), float(max(p))],
        "observed_volume_range": [float(min(v)), float(max(v))],
        "observed_temperature_range_k": [293, 293],
        "degrees_of_freedom": result["degrees_of_freedom"],
        "solver_success": result["success"],
        "rmse_gpa": result["pressure_rmse_gpa"],
        "published_rmse_gpa": float(
            np.sqrt(np.mean((bm3(v, V0, *published) - p) ** 2))
        ),
        "parameter_covariance": result["residual_scaled_covariance"],
        "covariance_convention": "Conditional residual-scaled one-sigma covariance: inverse(J^T J)*SSE/(22-2), V0 exact. Not recovered experimental uncertainty or confidence convention.",
        "reduced_chi_square": None,
        "qualification": COMMON + " " + QUALIFICATIONS["full"],
        "reproduction_report": REPORT,
    }


def catalog_records():
    """Define the sole public full-range source record with an explicitly adopted V0."""
    records = []
    for side, identifier in IDS.items():
        key = {"full": "full_range", "low": "below_22_gpa", "high": "above_22_gpa"}[
            side
        ]
        values, errors = SOURCE[key]
        location = (
            "page 140402-2, section III and Figure 1"
            if side == "full"
            else "page 140402-3, section III"
        )
        provenance = {
            "determination_method": "Single-crystal X-ray compression measurements at 293 K; adopted V0 reconstructed from the calculated curve of the same experimental publication.",
            "V0": {
                "value_kind": "adopted",
                "source": "Figure 1 full-range calculated curve; independent vector-graphics BM3 recovery",
                "report": "docs/data/prescher-2012-fe3c-recovery.json",
                "report_key": "source_drawn_curve_reconstruction.parameters.V0",
                "source_value_published": False,
                "source_constraint_recovered": False,
                "uncertainty": "not recovered; drawing-vertex covariance is not experimental V0 uncertainty",
                "applies_to": "Full-range record only; original fixed/free V0 constraint was not recovered",
            },
            "K0": location
            + "; published central coefficient and printed error unchanged",
            "K0_prime": location
            + "; published central coefficient and printed error unchanged",
            "reference_state": "P_ref=0 GPa, T_ref=293 K; conventional orthorhombic Fe3C cell, Z=4. Reference lattice of the material card is separate from this EOS V0.",
        }
        records.append(
            {
                "identifier": identifier,
                "label": "Prescher et al. (2012), "
                + {"full": "full-range", "low": "below 22 GPa", "high": "22-47 GPa"}[
                    side
                ]
                + " BM3 (adopted V0)",
                "model": "BM3",
                "temperature_ref": 293,
                "parameters": {
                    "V0": {"value": V0, "fixed": True},
                    "K0": {"value": values[0], "error": errors[0]},
                    "K0_prime": {"value": values[1], "error": errors[1]},
                },
                "parameter_error_confidence": None,
                "experimental_pressure_range_gpa": {
                    "full": [0, 47.4],
                    "low": [0, 22],
                    "high": [22, 47.4],
                }[side],
                "experimental_temperature_range_k": [293, 293],
                "pressure_range_status": "reported_qualitatively",
                "determination_method": "experimental",
                "reference": REFERENCE,
                "notes": COMMON
                + " "
                + QUALIFICATIONS[side]
                + " The main text rounds the upper limit to 47 GPa; Table S1 reaches 47.4 GPa. The nominal split ranges do not establish separate structural phases or a phase-stability envelope.",
                "fit_datasets": [PREFIX + "table_s1", DATASETS[side]]
                if side != "low"
                else [DATASETS[side]],
                "pressure_calibration": {
                    "status": "resolved",
                    "methods": [
                        {
                            "kind": "ruby_fluorescence",
                            "reference_calibration_record": "ruby_mao_1986",
                            "reference": {
                                "authors": ["Mao", "Xu", "Bell"],
                                "year": 1986,
                                "source": "Journal of Geophysical Research",
                                "doi": "10.1029/JB091iB05p04673",
                            },
                            "source_location": "Main article experimental section, reference 13",
                        }
                    ],
                    "recalculation": {
                        "status": "missing_calibrant_observations",
                        "notes": "Row-wise ruby wavelengths before/after measurement are unavailable; preserve published pressures and widths. No pressure recalibration applied.",
                    },
                    "audit_date": "2026-10-30",
                    "notes": "Neon medium. Pressure widths derive from before/after ruby-reading differences, not established independent standard deviations.",
                },
                "evidence": {
                    "parameter_provenance": provenance,
                    "scientific_validation": {
                        "status": "primary_source_validated",
                        "audit_date": "2026-10-30",
                        "note": "Published BM3 K coefficients/errors, units and reference conventions checked against the primary paper; V0 is an explicitly documented graph-derived adoption, not a published coefficient. Conditional agreement includes both printed error values; original confidence and author workflow remain unconfirmed.",
                        "verified_fields": [
                            "equation",
                            "published_K0_and_K0_prime",
                            "published_K_errors",
                            "units",
                            "volume_basis",
                            "temperature_ref",
                            "adopted_V0_derivation",
                            "conditional_parameter_parity",
                        ],
                        "primary_source_check": {
                            "doi": DOI,
                            "access_url": "https://doi.org/" + DOI,
                            "locations": [
                                location,
                                "Figure 1, drawn full-range curve (adopted V0)",
                                "Final supplementary Table S1, pages 1-2",
                            ],
                            "finding": "Published K coefficients/errors preserved. Unprinted V0 supplied by explicit Figure 1 full-curve adoption; source constraint and original workflow unresolved. See "
                            + REPORT,
                        },
                        "primary_data_check": {
                            "status": "primary_table",
                            "dataset_identifiers": [PREFIX + "table_s1", DATASETS[side]]
                            if side != "low"
                            else [DATASETS[side]],
                            "finding": "All 19 final Table S1 states retained separately from all 22 digitized graph markers. Split diagnostics remain outside the public catalog. The named comparison uses its own explicit row selection; these datasets are not asserted to be the original author fit inputs.",
                        },
                        "unresolved": [COMMON, QUALIFICATIONS[side]],
                    },
                },
            }
        )
    return records


def supporting_datasets():
    """Small inline evidence tables, separated by source kind."""
    columns = [
        {
            "name": "source_row",
            "quantity": "index",
            "unit": "dimensionless",
            "role": "value",
        },
        {
            "name": "pressure_gpa",
            "quantity": "pressure",
            "unit": "GPa",
            "role": "value",
        },
        {
            "name": "volume_a3",
            "quantity": "conventional_unit_cell_volume",
            "unit": "angstrom^3",
            "role": "value",
        },
    ]
    markers = json.loads(VECTORS.read_text())
    datasets = []
    for side in ("full",):
        rows = [
            [m["marker"], m["pressure_gpa"], m["volume_a3"]] for m in markers["markers"]
        ]
        datasets.append(
            {
                "identifier": DATASETS[side],
                "kind": "pressure_volume_digitized"
                if side == "full"
                else "archived_working_volume_with_published_pressure",
                "description": QUALIFICATIONS[side],
                "reference": REFERENCE,
                "source_location": "Figure 1 main panel, all 22 markers"
                if side == "full"
                else "Private Fe3C archive, own/SHELX p1-p19 CIFs; paired Table S1 pressures by run index (conditional)",
                "columns": columns,
                "rows": rows,
                "used_by_eos_records": [IDS[side]],
                "uncertainty": {
                    "pressure": "Experimental errors not recovered; no sigma inferred",
                    "volume": "Experimental errors not recovered; no sigma inferred",
                    "coordinate_precision": markers["precision_limit"]
                    if side == "full"
                    else "Original CIF tokens and hashes preserved in the linked evidence file",
                },
                "provenance": {
                    "type": markers["kind"]
                    if side == "full"
                    else "working_reduction_not_final_table",
                    "evidence_file": "research/cementite/prescher-2012-figure1-vectors.json"
                    if side == "full"
                    else "research/cementite/prescher-2012-own-shelx-cells.json",
                    "normalization": "Conventional Fe3C cell in A^3; source values retained",
                    "pressure_pairing": "Digitized coordinates"
                    if side == "full"
                    else "Table S1 pressure paired by run index; not confirmed original CIF pressure metadata",
                },
                "source_rights": "Unspecified; source figures and author archive are not relicensed.",
                "contributor_rights": "CC0 applies only to contributor-held transcription, normalization, metadata and arrangement rights.",
            }
        )
    return datasets
