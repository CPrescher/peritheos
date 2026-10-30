"""Repository-only catalog mapping of the separately preserved Fe3C recovery.

Published K coefficients/errors are never replaced by reconstructed values.
V0 is an explicit adoption from the drawn full-range curve, not a source value.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from scripts.recover_prescher_2012_fe3c import SOURCE, TABLE, VECTORS, bm3, fit

ROOT = Path(__file__).resolve().parents[1]
V0 = 154.5678863704812
DOI = "10.1103/PhysRevB.85.140402"
REPORT = "docs/literature-reproductions/prescher-2012-recovery.md"
PREFIX = "cementite_prescher_2012_"
IDS = {"full": PREFIX + "full_range_bm3"}
RETIRED_IDS = {PREFIX + "below_22_gpa_bm3", PREFIX + "22_47_gpa_bm3"}
DATASETS = {"full": PREFIX + "table_s1"}
GRAPH_DATASET = PREFIX + "figure1_digitized"
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
COMMON = "Conditional coefficient parity compares K0 and K0_prime at the published printed precision, with adopted V0=154.5678863704812 A^3 per conventional Fe3C cell (Z=4) held fixed. Published K0=161 +/-2 GPa and K0_prime=5.9 +/-0.2 are retained unchanged. V0 is adopted from the published calculated full-range curve, not a printed source coefficient or the ambient measurement. Its uncertainty and the original fixed/free V0 convention are unrecovered. Numerical coefficient agreement does not establish uncertainty agreement, author-workflow reproduction or physical validation."
QUALIFICATIONS = {
    "full": "The primary experimental dataset is all 19 numbered p1-p19 points from final Table S1, with published ruby pressures and reported widths, spanning 1.8-47.4 GPa at nominal 293 K. The recorded comparison uses equal pressure weights with V0 fixed and K0/K0_prime free. It gives K0=160.994011 +/-2.532193 GPa and K0_prime=5.919478 +/-0.240867, displayed as 161 +/-3 and 5.9 +/-0.2 at the paper precision. The K0 uncertainty is not reproduced: 2.532193 GPa rounds to 3, not the published 2 GPa. Covariance is inverse(J^T J)*SSE/(19-2), conditional on exact V0 and exact input volumes. Source coordinate widths are preserved but are not the weights in this comparison; source confidence and calibration covariance are unknown. The 22 digitized Figure 1 markers remain separate supplementary evidence, not the primary dataset. No ambient or p20-p22 point is added. For p16, the final source volume is 133.25 A^3, not the alternate SHELX cell 133.58 A^3. Published split fits remain audit-only."
}


def selected_rows(side="full"):
    """Return all final numeric p1-p19 rows without changing source coordinates."""
    if side != "full":
        raise ValueError(
            "Split fits are audit-only; the catalog exposes the full range"
        )
    rows = np.genfromtxt(TABLE, delimiter=",", names=True)
    return rows["pressure_gpa"], rows["volume_a3"]


def ledger_outcome(record):
    """Replay the numbered-data fixed-V0 coefficient comparison and error mismatch."""
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
        if not matched:
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
                "uncertainty_finding": "matched_at_printed_precision"
                if error_matches
                else "not_reproduced",
                "display_decimal_places": i,
            }
        )
    return {
        "status": "parity",
        "parity_basis": "conditional_K0_and_K0_prime_at_published_printed_precision",
        "comparison_scope": "parameter_agreement",
        "assessment_basis": "Conditional agreement of K0 and K0_prime at published printed precision for the 19 numbered p1-p19 observations, with adopted curve-derived V0 held fixed. K0 uncertainty is not reproduced: 2.532193 GPa rounds to 3, while the published error is 2 GPa.",
        "assessment_qualifications": [
            COMMON,
            QUALIFICATIONS["full"],
            "This is a conditional coefficient comparison, not recovery of source confidence, calibration covariance or the author fitting workflow.",
        ],
        "conditional": True,
        "compared_parameters": ["K0", "K0_prime"],
        "author_workflow_reproduced": False,
        "dataset_identifiers": [DATASETS["full"]],
        "fit_kind": "isothermal_pv",
        "objective": "pressure_residuals",
        "selection": "All 19 numbered p1-p19 final Table S1 rows; no ambient point",
        "columns": {
            "pressure": "pressure_gpa",
            "volume": "volume_a3",
            "temperature": None,
        },
        "weighting": "Equal pressure weights; reported coordinate widths retained but not used as weights",
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
        "covariance_convention": "Conditional residual-scaled one-sigma covariance: inverse(J^T J)*SSE/(19-2), V0 exact. Not recovered experimental uncertainty or confidence convention.",
        "reduced_chi_square": None,
        "qualification": COMMON + " " + QUALIFICATIONS["full"],
        "reproduction_report": "docs/literature-reproductions/prescher-2012-numbered-full-range-fit.md",
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
                "fit_datasets": [DATASETS[side]],
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
                        "note": "Published BM3 K coefficients/errors, units and reference conventions checked against the primary paper; V0 is an explicitly documented graph-derived adoption, not a published coefficient. Conditional agreement covers central K coefficients only; the K0 uncertainty mismatch is explicit and original confidence/workflow remain unconfirmed.",
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
                            "dataset_identifiers": [DATASETS[side]],
                            "digitized_dataset_identifiers": [GRAPH_DATASET],
                            "finding": "The primary dataset is all 19 numbered p1-p19 final Table S1 states with published ruby pressures. All 22 digitized graph markers remain separate supplementary evidence. Original author input selection and constraints remain unrecovered.",
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
                "identifier": GRAPH_DATASET,
                "kind": "pressure_volume_digitized"
                if side == "full"
                else "archived_working_volume_with_published_pressure",
                "description": "Supplementary evidence: all 22 digitized Figure 1 markers. Primary numerical dataset is final Table S1, numbered p1-p19; source/graph observations are not merged.",
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
