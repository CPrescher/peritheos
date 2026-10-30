"""Copy existing reconstruction findings into native EOSMAT's standard assessment.

Run ``python -m scripts.standardize_reconstruction_metadata [paths ...]``;
``--check`` is read-only. No fitting, reclassification or coefficient replacement
is performed. The common ledger and original record evidence remain intact.
"""

from __future__ import annotations

import argparse
import copy
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from peritheos import validate_eosmat_document
from peritheos._eosmat_format import format_eosmat
from scripts.reconstruction_links import report_links

ROOT = Path(__file__).resolve().parents[1]
REPORT = "docs/data/primary-eos-refits.json"
STATUS = {
    "parity": "parity",
    "similar": "similar",
    "parity_not_achieved": "not_reproduced",
    "not_refittable": "not_refittable",
    "source_reconstruction": "source_reconstruction",
}
BASIS = {
    "conditional_coefficients_and_errors_at_published_precision_fixed_V0": "Conditional agreement of K0, K0_prime and both error values at published printed precision with curve-derived V0 held fixed. Errors are residual-based and conditional on exact V0; the source confidence convention and original author workflow remain unconfirmed.",
    "conditional_K0_and_K0_prime_at_published_printed_precision": "Conditional agreement of K0 and K0_prime at published printed precision with an explicitly adopted, fixed reference volume. Parameter-error reproduction is assessed separately; coefficient agreement does not establish source confidence, calibration covariance or author-workflow reproduction.",
    "conditional_reported_widths_and_print_precision": "Conditional numerical agreement within reported source error widths where available, and explicitly verified half-last-digit printing tolerances otherwise. Printing precision is not a standard error or confidence interval; fixed constraints and source-equation versus independent-data scope remain explicit.",
    "conditional_combined_error_parameter_agreement": "Conditional parameter agreement under both the combined-error comparison and numerical similarity screen; source confidence conventions remain explicit, and unknown source confidence or shared-data covariance does not become a formal confidence test or author-workflow reproduction.",
    "within_reported_parameter_errors": "Agreement within the source's reported parameter error widths; this is not an inferred confidence convention.",
    "conditional_fixed_V0_diagnostic_only": "Conditional fixed-V0 comparison under the documented diagnostic protocol; original source weighting and workflow are not established.",
    "registered_peritheos_refit_reproduced": "Reproduction of the registered Peritheos refit's stored coefficients under its declared protocol.",
    "reproduction_of_stored_peritheos_refit": "Reproduction of the stored independent Peritheos refit under its declared protocol.",
    "conditional_central_value_at_published_precision": "Conditional agreement of the fitted central value at published precision; this does not establish uncertainty or calibration-network parity.",
    "published_decimal_precision": "Agreement at the source's printed decimal precision.",
    "source_curve_reading_and_print_precision": "Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.",
}
FIT_FIELDS = (
    "fit_kind",
    "objective",
    "selection",
    "columns",
    "weighting",
    "absolute_sigma",
    "free_parameters",
    "fixed_parameters",
    "validation_fixed_parameters",
    "fixed_values",
    "tied_parameters",
    "observations",
    "observed_pressure_range_gpa",
    "observed_temperature_range_k",
    "observed_volume_range",
    "degrees_of_freedom",
    "solver_success",
    "solver_message",
    "rmse_gpa",
    "published_rmse_gpa",
    "published_pressure_rmse_gpa",
    "refitted_pressure_rmse_gpa",
    "rmse_volume_a3",
    "rmse_shock_velocity_km_s",
    "reduced_chi_square",
    "parameter_covariance",
    "covariance",
    "covariance_convention",
    "comparison_scope",
    "reconstruction_kind",
    "pressure_reconstruction",
    "dataset_identifiers",
    "reproduction_report",
    "published_ruby_parameters",
    "sensitivity_fit_ruby_parameters",
    "stages",
)


def _basis(row: Mapping[str, Any]) -> str:
    explicit = row.get("assessment_basis")
    if explicit is not None:
        if not isinstance(explicit, str) or not explicit.strip():
            raise ValueError("An explicit assessment basis must be nonempty prose")
        return explicit
    code = row.get("parity_basis")
    if isinstance(code, str) and code:
        return BASIS.get(code, code)
    comparisons = row.get("parameters", [])
    if (
        row["status"] == "parity"
        and comparisons
        and all(
            p.get("within_combined_2sigma") is True and p.get("similar") is True
            for p in comparisons
        )
    ):
        return "Every compared free coefficient meets both the recorded combined-two-sigma comparison and parameter-specific numerical similarity screen. Error conventions and fit qualifications remain applicable."
    if (
        row["status"] == "similar"
        and comparisons
        and all(
            p.get("within_combined_2sigma") is True or p.get("similar") is True
            for p in comparisons
        )
    ):
        return "Every compared free coefficient meets either the recorded combined-two-sigma comparison or parameter-specific numerical similarity screen; the audit does not classify the complete comparison as parity."
    for key in ("reason", "qualification"):
        if isinstance(row.get(key), str) and row[key].strip():
            return row[key]
    return {
        "parity": "Numerical agreement under the recorded audit protocol; consult its parameter comparisons and source qualifications for the criterion.",
        "similar": "Qualified numerical similarity under the recorded audit protocol; no stronger result is inferred.",
        "parity_not_achieved": "The recorded parameter reconstruction did not meet the audit's agreement criterion; the parameter comparisons retain the actual reconstructed values.",
        "not_refittable": "The existing audit does not establish a complete independent refit from the available observations and source-reduction inputs.",
        "source_reconstruction": "Reconstruction of a source calculation or calibration relationship without an independent full EOS coefficient fit.",
    }[row["status"]]


def assessment_for_record(
    record: Mapping[str, Any], row: Mapping[str, Any]
) -> dict[str, Any]:
    """Organize a recorded result; never compute or promote a scientific outcome."""
    if record["identifier"] != row["record_identifier"]:
        raise ValueError("Reconstruction ledger identifier does not match EOS record")
    status = STATUS[row["status"]]
    scope = (
        "source_equation"
        if status == "source_reconstruction"
        or row.get("fit_kind") == "source_curve_validation_only"
        else "independent_refit"
        if record.get("record_kind") == "refit"
        else "parameter_agreement"
    )
    qualifications: list[str] = []

    def retain(value: Any, label: str = "") -> None:
        if isinstance(value, list):
            for item in value:
                retain(item, label)
        elif isinstance(value, str) and value.strip():
            item = f"{label}: {value}" if label else value
            if item not in qualifications:
                qualifications.append(item)

    basis = _basis(row)
    for key in (
        "qualification",
        "reason",
        "recommendation_reason",
        "usage_recommendation",
    ):
        if row.get(key) != basis:
            retain(
                row.get(key),
                "Usage recommendation" if key == "usage_recommendation" else "",
            )
    validation = record.get("evidence", {}).get("scientific_validation", {})
    for key in ("unresolved", "recommendation_reason", "usage_recommendation"):
        retain(
            validation.get(key),
            "Usage recommendation" if key == "usage_recommendation" else "",
        )
    legacy = record.get("evidence", {}).get("reproduction", {})
    if isinstance(legacy, Mapping):
        retain(legacy.get("reason"))
        retain(
            legacy.get("original_publication_reproduction_status"),
            "Original publication reproduction status",
        )
    retain(
        row.get("original_publication_reproduction_status"),
        "Original publication reproduction status",
    )
    if scope == "independent_refit":
        qualifications.insert(
            0,
            "This comparison concerns the stored independent Peritheos refit, not recovery of the original publication's coefficients or author workflow.",
        )
    else:
        qualifications.append(
            "This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use."
        )
    explicit_qualifications = row.get("assessment_qualifications")
    if explicit_qualifications is not None:
        if (
            not isinstance(explicit_qualifications, list)
            or not explicit_qualifications
            or any(
                not isinstance(value, str) or not value.strip()
                for value in explicit_qualifications
            )
        ):
            raise ValueError(
                "Explicit assessment qualifications must be a nonempty list of prose"
            )
        # Preserve the audit's reviewed narrative verbatim. This changes no
        # classification, numerical comparison, constraint or confidence field.
        qualifications = copy.deepcopy(explicit_qualifications)
    result: dict[str, Any] = {
        "status": status,
        "scope": scope,
        "basis": basis,
        "qualifications": qualifications,
    }
    if row.get("parameters"):
        # Copy the original comparison rows, names, values, errors and flags
        # verbatim. In a refit assessment, 'published' denotes the stored refit.
        result["parameters"] = copy.deepcopy(row["parameters"])
    else:
        qualifications.append(
            "No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred."
        )
    result["fit"] = {key: copy.deepcopy(row[key]) for key in FIT_FIELDS if key in row}
    if "data_agreement_status" in row:
        # A data-agreement assessment retains its distinct coefficient result
        # and measured-coordinate diagnostics alongside the fitting objective.
        result["fit"].update(
            {
                key: copy.deepcopy(row[key])
                for key in (
                    "coefficient_comparison_status",
                    "data_agreement_status",
                    "weighted_sum_squares",
                    "published_max_abs_pressure_residual_gpa",
                    "pressure_rmse_coordinates",
                )
                if key in row
            }
        )
    result["report"] = REPORT
    if row.get("material_identifier"):
        result.update(report_links(row["material_identifier"], record["identifier"]))
    result["source_status"] = row["status"]
    if row.get("parity_basis"):
        result["basis_code"] = copy.deepcopy(row["parity_basis"])
    return result


def standardize_document(
    document: Mapping[str, Any], rows: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    """Add one assessment per record, retaining every pre-existing field/value."""
    result = copy.deepcopy(dict(document))
    if result.get("format_version") != 4:
        raise ValueError("Reconstruction standardization requires native EOSMAT v4")
    for record in result["eos_records"]:
        identifier = record["identifier"]
        if identifier not in rows:
            raise ValueError(
                f"No existing reconstruction ledger entry for {identifier}"
            )
        if (
            rows[identifier].get("material_identifier", result["identifier"])
            != result["identifier"]
        ):
            raise ValueError(
                f"Reconstruction ledger material does not match {identifier}"
            )
        assessment = assessment_for_record(record, rows[identifier])
        reproduction = record.setdefault("evidence", {}).setdefault("reproduction", {})
        if not isinstance(reproduction, dict):
            raise ValueError(
                f"Historical reproduction metadata is not an object: {identifier}"
            )
        existing = reproduction.get("assessment")
        # Add navigation to an otherwise identical existing assessment. Never
        # overwrite scientific findings or an explicitly different report link.
        previous = {
            k: v
            for k, v in assessment.items()
            if k not in {"report_url", "report_path"}
        }
        if existing is not None and existing not in (assessment, previous):
            raise ValueError(
                f"Existing assessment differs for {identifier}; review it explicitly before replacing"
            )
        reproduction["assessment"] = assessment
    validate_eosmat_document(result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--ledger", type=Path, default=ROOT / REPORT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    ledger = json.loads(args.ledger.read_text(encoding="utf-8"))
    rows = {row["record_identifier"]: row for row in ledger["records"]}
    if len(rows) != len(ledger["records"]):
        raise ValueError("Duplicate reconstruction ledger identifiers")
    paths = args.paths or [ROOT / "peritheos/data/materials"]
    files = sorted(
        {
            p
            for path in paths
            for p in (path.glob("*.eosmat") if path.is_dir() else [path])
        }
    )
    changed = []
    for path in files:
        original = path.read_bytes()
        result = standardize_document(json.loads(original), rows)
        rendered = format_eosmat(result).encode("utf-8")
        if rendered == original:
            continue
        changed.append(path)
        if not args.check:
            if path.read_bytes() != original:
                raise RuntimeError(f"Material changed while standardizing: {path}")
            path.write_bytes(rendered)
    if args.check and changed:
        for path in changed:
            print(f"Needs reconstruction metadata: {path}")
        return 1
    print(
        f"{'Checked' if args.check else 'Standardized'} {len(files)} files; {len(changed)} changed."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
