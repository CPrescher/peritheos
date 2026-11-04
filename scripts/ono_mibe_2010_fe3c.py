"""Repository-only audit and append-only integration of the original experiment.

Run ``uv run python -m scripts.ono_mibe_2010_fe3c`` to refresh the numerical
report; add ``--integrate`` once to append the source record and observations.
The 0 K calculations and Litasov's later recalibration are not installed here.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from decimal import Decimal
from pathlib import Path

import numpy as np
from jsonschema import Draft202012Validator
from scipy.optimize import brentq, least_squares

from peritheos import eosmat_schema, validate_eosmat_document
from peritheos._eosmat_format import format_eosmat
from scripts.standardize_reconstruction_metadata import assessment_for_record

ROOT = Path(__file__).resolve().parents[1]
IDENTIFIER = "cementite_ono_mibe_2010_bm3"
DATASET = "cementite_ono_mibe_2010_table1"
REPORT = "docs/literature-reproductions/ono-mibe-2010-fe3c.md"
NUMERICAL = "docs/data/ono-mibe-2010-fe3c.json"
TOKENS = ROOT / "peritheos/data/datasets/cementite-ono-mibe-2010-table1-tokens.csv"
CARD = ROOT / "peritheos/data/materials/cementite.eosmat"
REFERENCE = {
    "authors": ["Ono, S.", "Mibe, K."],
    "year": 2010,
    "title": "Magnetic transition of iron carbide at high pressures",
    "source": "Physics of the Earth and Planetary Interiors",
    "volume": "180",
    "locator": "1-6",
    "doi": "10.1016/j.pepi.2010.03.008",
}
QUALIFICATION = (
    "Original experimental BM3 at 300 K, restricted to the below-35-GPa fit "
    "(last measured fit point 34.3 GPa). V0=155.4(2) A^3 is adopted from the "
    "measured Table 1 ambient cell and checked against Figure 4; the paper does "
    "not explicitly state its V0 fitting constraint. K0=167 GPa and K0_prime=6.7 "
    "are the original Table 3 values, without reported coefficient errors. "
    "Source-equation agreement is checked at graphical resolution; the separate "
    "fixed-anchor equal-pressure refit is diagnostic, not coefficient parity "
    "or recovery of author weights. Original pressures and all 14 measurements "
    "are preserved. No Litasov recalibration, 0 K DFT branch, high-pressure "
    "magnetic transition model or thermal extrapolation is substituted."
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_token(token):
    value, sep, tail = token.partition("(")
    width = (
        float(
            Decimal(tail.rstrip(")"))
            * Decimal(10) ** Decimal(value).as_tuple().exponent
        )
        if sep
        else None
    )
    return float(value), width


def observations():
    with TOKENS.open(newline="", encoding="utf-8") as stream:
        raw = list(csv.DictReader(stream))
    rows = []
    for i, row in enumerate(raw):
        values = {
            "source_row": i + 1,
            "selection_code": 0 if i == 0 else 1 if i < 6 else 2,
        }
        for field in ("pressure_gpa", "a_a", "b_a", "c_a", "volume_a3"):
            values[field], values[field + "_error_width"] = parse_token(row[field])
        rows.append(values)
    return raw, rows


def pressure(volume, parameters=(167.0, 6.7)):
    k, kp = parameters
    x = (155.4 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def volume(pressures, parameters=(167.0, 6.7)):
    return np.array(
        [
            brentq(lambda v: pressure(v, parameters) - p, 120, 155.4)
            for p in np.atleast_1d(pressures)
        ]
    )


def reproduce():
    selected = [r for r in observations()[1] if r["selection_code"] == 1]
    p = np.array([r["pressure_gpa"] for r in selected])
    v = np.array([r["volume_a3"] for r in selected])
    result = least_squares(
        lambda x: pressure(v, x) - p, [167, 6.7], xtol=1e-13, ftol=1e-13, gtol=1e-13
    )
    # Manual readings of the dashed experimental curve, not experimental points.
    # Original embedded Figure 4 JPEG: 632 x 517 px. Axes x=90..615 map
    # 0..100 GPa; y=11..403 map 160..120 A^3. Reading envelope is +/-4 px.
    pixels = np.array([[142.5, 131], [195, 186], [247.5, 229]])
    gp = (pixels[:, 0] - 90) * 100 / 525
    gv = 160 - (pixels[:, 1] - 11) * 40 / 392
    predicted = volume(gp)
    graph = {
        "source": "Figure 4 dashed experimental BM3, p. 3; manual curve readings, not observations",
        "source_image_sha256": "416d2b579379daff74c33ad37ecca9d294e4968454a014c6180ab41546706cbd",
        "image_dimensions_px": [632, 517],
        "axis_calibration": {
            "x_px": [90, 615],
            "pressure_gpa": [0, 100],
            "y_px": [11, 403],
            "volume_a3": [160, 120],
        },
        "pixel_checkpoints": pixels.tolist(),
        "pressure_gpa": gp.tolist(),
        "read_volume_a3": gv.tolist(),
        "predicted_volume_a3": predicted.tolist(),
        "difference_a3": (predicted - gv).tolist(),
        "tolerance_a3": 0.41,
        "tolerance_basis": "Four source-image pixels in volume, rounded up; includes dashed-line thickness, axis/readout and pressure-location ambiguity. Not a measurement uncertainty.",
        "passed": bool(np.all(np.abs(predicted - gv) <= 0.41)),
    }
    rp, rv = pressure(v) - p, volume(p) - v
    fit_rp = pressure(v, result.x) - p
    return {
        "record_identifier": IDENTIFIER,
        "source_equation_check": graph,
        "source_observation_comparison": {
            "source_rows": [r["source_row"] for r in selected],
            "pressure_residuals_gpa": rp.tolist(),
            "volume_residuals_a3": rv.tolist(),
            "rms_pressure_gpa": float(np.sqrt(np.mean(rp**2))),
            "max_abs_pressure_gpa": float(max(abs(rp))),
            "max_abs_volume_a3": float(max(abs(rv))),
            "notes": "No pressure errors reported; volume-error confidence unspecified. No chi-square or formal significance inferred.",
        },
        "independent_diagnostic_refit": {
            "selection": "Five nonambient Table 1 rows below 35 GPa; measured ambient row supplies fixed V0 only",
            "fixed_V0_a3": 155.4,
            "objective": "Equal pressure residuals, no supplied sigmas; author objective and weights unknown",
            "parameters_K0_K0_prime": result.x.tolist(),
            "difference_from_published": (result.x - [167, 6.7]).tolist(),
            "rms_pressure_gpa": float(np.sqrt(np.mean(fit_rp**2))),
            "solver_success": bool(result.success),
            "coefficient_parity_claimed": False,
            "author_workflow_reproduced": False,
            "notes": "Published coefficients unchanged. Missing source coefficient errors and unknown constraints preclude an uncertainty-based parameter comparison.",
        },
        "qualification": QUALIFICATION,
    }


def ledger_outcome(record):
    result = reproduce()
    if not result["source_equation_check"]["passed"]:
        raise AssertionError("Original Ono-Mibe experimental curve check failed")
    return {
        "status": "source_reconstruction",
        "parity_basis": "source_curve_reading_and_print_precision",
        "comparison_scope": "source_equation",
        "fit_kind": "source_curve_validation_only",
        "author_workflow_reproduced": False,
        "dataset_identifiers": [DATASET],
        "observations": 5,
        "selection": "Original experimental fit <35 GPa; adopted measured ambient V0; five nonambient rows through 34.3 GPa",
        "published_rmse_gpa": result["source_observation_comparison"][
            "rms_pressure_gpa"
        ],
        "independent_diagnostic_refit": result["independent_diagnostic_refit"],
        "source_equation_check": result["source_equation_check"],
        "qualification": QUALIFICATION,
        "reproduction_report": REPORT,
    }


def catalog_dataset():
    raw, rows = observations()
    columns = [
        {"name": n, "quantity": "index", "unit": "1", "role": "flag"}
        for n in ("source_row", "selection_code")
    ]
    for field in ("pressure_gpa", "a_a", "b_a", "c_a", "volume_a3"):
        quantity = (
            "pressure"
            if field == "pressure_gpa"
            else "conventional_unit_cell_volume"
            if field == "volume_a3"
            else "lattice_parameter"
        )
        unit = (
            "GPa"
            if field == "pressure_gpa"
            else "angstrom^3"
            if field == "volume_a3"
            else "angstrom"
        )
        columns.extend(
            [
                {"name": field, "quantity": quantity, "unit": unit, "role": "value"},
                {
                    "name": field + "_error_width",
                    "quantity": quantity,
                    "unit": unit,
                    "role": "uncertainty",
                    "of": field,
                },
            ]
        )
    return {
        "identifier": DATASET,
        "kind": "pressure_volume_lattice",
        "description": "Complete Table 1: 14 measured states at 300 K, including all axes/volume errors and the excluded >35 GPa states. Codes: 0 measured ambient anchor, 1 five nonambient fit-domain states, 2 outside fit domain.",
        "reference": copy.deepcopy(REFERENCE),
        "source_location": "Table 1, p. 3; Figure 4 caption establishes 300 K",
        "columns": columns,
        "rows": [[r[c["name"]] for c in columns] for r in rows],
        "source_tokens": raw,
        "row_count": 14,
        "used_by_eos_records": [IDENTIFIER],
        "uncertainty_convention": "Parenthetical last-digit lattice/volume errors; confidence and covariance unspecified. Pressure errors not reported and retained as null, not zero.",
        "provenance": {
            "source_manifest": "datasets/cementite-ono-mibe-2010-source.json",
            "source_pdf_sha256": "3d254d8de9bcf45e50eec1665f3b2edf7f424c057af780a574937bec5dea7f5a",
            "source_tokens_resource": "datasets/" + TOKENS.name,
            "source_tokens_sha256": digest(TOKENS),
            "normalization": "Pnma conventional-cell volumes used literally, Z=4; products a*b*c are not substituted for measured printed V.",
            "calibrant_checkpoint": {
                "source": "Figure 2 caption, p. 3; outside represented low-pressure fit",
                "pressure_gpa": 58.5,
                "NaCl_phase": "B2",
                "NaCl_a_a": 2.872,
                "NaCl_a_error": None,
            },
        },
        "source_rights": "Elsevier 2010 article copyright; source numerical license unspecified. Article text, figures and PDF are not relicensed.",
        "contributor_rights": "CC0-1.0 only for contributor-held transcription, normalization, metadata and arrangement rights.",
        "observation_views": [
            {
                "identifier": IDENTIFIER + "_observations",
                "eos_records": [IDENTIFIER],
                "evidence_kind": "experimental",
                "notes": "All 14 original measured states retained on their published NaCl pressure scale; only the below-35-GPa subset belongs to the represented experimental branch.",
                "coordinates": {
                    "pressure": {"kind": "column", "column": "pressure_gpa"},
                    "volume": {
                        "kind": "column",
                        "column": "volume_a3",
                        "unit": "angstrom^3/conventional_unit_cell",
                        "reference": "Table 1 conventional Pnma cell, Z=4; printed volumes retained directly",
                    },
                },
                "uncertainties": {"volume": "volume_a3_error_width"},
                "fit_selection": {
                    "row_selection": [{"column": "selection_code", "values": [1]}],
                    "notes": "Conditional fixed-anchor diagnostic selection; ambient measurement supplies V0. High-pressure rows preserved but excluded from this low-pressure fit.",
                },
            }
        ],
    }


def catalog_record():
    record = {
        "identifier": IDENTIFIER,
        "label": "Ono and Mibe (2010), original experimental 300 K BM3, below 35 GPa",
        "model": "BM3",
        "record_kind": "published",
        "determination_method": "experimental",
        "temperature_ref": 300,
        "parameters": {
            "V0": {"value": 155.4, "error": 0.2},
            "K0": {"value": 167, "error": None},
            "K0_prime": {"value": 6.7, "error": None},
        },
        "parameter_error_confidence": None,
        "experimental_pressure_range_gpa": [0, 34.3],
        "experimental_temperature_range_k": [300, 300],
        "pressure_range_status": "reported_exactly",
        "validity": {
            "pressure_gpa": [0, 34.3],
            "temperature_k": [300, 300],
            "notes": [
                "Source Table 3 says <=35 GPa; actual last fit-domain measurement is 34.3 GPa. All 14 measurements extend to 66.5 GPa but do not validate this branch there.",
                "Source-assigned low-pressure ferromagnetic branch; no magnetic stability boundary is implemented. Annealing at 1500-2000 K precedes room-temperature XRD and is not thermal EOS coverage.",
            ],
        },
        "reference": copy.deepcopy(REFERENCE),
        "range_provenance": {
            "pressure": "Figure 4 experimental fit <35 GPa; Table 3 <=35 GPa; Table 1 actual included span 0-34.3 GPa",
            "temperature": "Figure 4 experimental circles at 300 K; Methods/Results post-anneal XRD",
        },
        "fit_datasets": [DATASET],
        "notes": QUALIFICATION,
        "pressure_calibration": {
            "status": "partially_resolved",
            "audit_date": "2026-11-01",
            "methods": [
                {
                    "kind": "equation_of_state",
                    "reference": {
                        "authors": ["Brown, J. M."],
                        "year": 1999,
                        "title": "The NaCl pressure standard",
                        "source": "Journal of Applied Physics",
                        "volume": "86",
                        "locator": "5801-5808",
                    },
                    "source_location": "Methods p. 2 and Table 1 footnote: NaCl EOS of Brown (1999) and Ono et al. (2006)",
                    "notes": "NaCl is the primary calibrant and pressure medium. Exact per-row calibrant observations and phase/scale assignment are not tabulated. No bundled scale substituted.",
                }
            ],
            "recalculation": {
                "status": "missing_calibrant_observations",
                "notes": "No NaCl cell readings for the five fit-domain states. Figure 2 gives one B2 NaCl a=2.872 A at 58.5 GPa outside this fit, retained as provenance.",
            },
            "notes": "Preliminary ruby calibration cites Dorogokupets and Oganov (2007); final pressures are NaCl based. Original pressures retained; later Litasov (2013) scale/fit is separate.",
        },
        "evidence": {
            "parameter_provenance": {
                "determination_method": "Table 3 footnote e and Figure 4: experimental post-anneal XRD at 300 K, BM3 fit below 35 GPa; no theoretical coefficients included.",
                "V0": "Table 1 measured ambient volume 155.4(2) A^3 adopted as reference anchor; independently consistent with Figure 4 dashed curve. Original fit constraint is not explicit, so no fixed=true source claim.",
                "K0": "Table 3, ferromagnetic experimental this-study row e: 167 GPa; no error reported.",
                "K0_prime": "Table 3 same row: 6.7; no error reported.",
                "equation": "Figure 4 caption: third-order Birch-Murnaghan; explicit pressure-volume expression p. 4 with x=(V/V0)^(1/3), identical to native BM3.",
                "reference_state": "Figure 4: experiment at 300 K. Table 1 ambient cell at P=0; conventional Pnma Fe12C4, Z=4, 16 atoms (p. 2). Existing material structure is separate from this measured EOS anchor.",
            },
            "scientific_validation": {
                "status": "primary_source_validated",
                "audit_date": "2026-11-01",
                "report_path": REPORT,
                "note": QUALIFICATION,
                "usage_recommendation": "qualified_use",
                "recommendation_reason": "Use original experimental curve within 300 K and 0-34.3 GPa; ambient anchor adoption and unavailable author fit details remain explicit.",
                "primary_source_check": {
                    "access_url": "https://doi.org/10.1016/j.pepi.2010.03.008",
                    "doi": REFERENCE["doi"],
                    "locations": [
                        "Table 1 and Figure 4, p. 3",
                        "BM3 pressure equation, p. 4",
                        "Table 3, experimental row e, p. 5",
                        "Methods/Results, p. 2",
                    ],
                    "finding": QUALIFICATION,
                },
                "verified_fields": [
                    "original_experimental_coefficients",
                    "adopted_measured_reference_volume",
                    "BM3_equation",
                    "source_curve_checkpoints",
                    "all_Table1_rows_and_errors",
                    "phase_cell_and_limits",
                    "pressure_scale_provenance",
                ],
                "primary_data_check": {
                    "status": "bundled",
                    "dataset_identifiers": [DATASET],
                    "finding": "All 14 original Table 1 states and errors retained. Five nonambient rows used for a separate conditional fixed-anchor diagnostic; published curve checked independently against Figure 4.",
                },
                "unresolved": [
                    "Original V0 fitting constraint, objective, weights, coefficient errors and covariance are not reported.",
                    "NaCl readings needed for a fit-domain pressure recalculation are not supplied.",
                ],
            },
            "reproduction": {
                "author_workflow_reproduced": False,
                "source_equation": reproduce()["source_equation_check"],
                "independent_diagnostic_refit": reproduce()[
                    "independent_diagnostic_refit"
                ],
            },
        },
    }
    outcome = {
        "material_identifier": "cementite",
        "record_identifier": IDENTIFIER,
        "fixed_parameters": [],
        "usage_recommendation": record["evidence"]["scientific_validation"][
            "usage_recommendation"
        ],
        "recommendation_reason": record["evidence"]["scientific_validation"][
            "recommendation_reason"
        ],
        **ledger_outcome(record),
    }
    outcome["qualification"] += (
        " " + record["evidence"]["scientific_validation"]["recommendation_reason"]
    )
    record["evidence"]["reproduction"]["assessment"] = assessment_for_record(
        record, outcome
    )
    return record


def integrate():
    original = json.loads(CARD.read_text())
    document = copy.deepcopy(original)
    if any(r["identifier"] == IDENTIFIER for r in document["eos_records"]):
        raise ValueError(
            "Original Ono-Mibe record already installed; no replacement authorized"
        )
    document["eos_records"].append(catalog_record())
    document.setdefault("datasets", []).append(catalog_dataset())
    validate_eosmat_document(document)
    Draft202012Validator(eosmat_schema(4)).validate(document)
    rendered = format_eosmat(document)
    restored = json.loads(rendered)
    assert restored["eos_records"][:-1] == original["eos_records"]
    assert restored["datasets"][:-1] == original.get("datasets", [])
    assert all(
        restored[k] == v
        for k, v in original.items()
        if k not in {"eos_records", "datasets"}
    )
    CARD.write_text(rendered)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrate", action="store_true")
    args = parser.parse_args()
    (ROOT / NUMERICAL).write_text(
        json.dumps(reproduce(), indent=2, allow_nan=False) + "\n"
    )
    if args.integrate:
        integrate()
