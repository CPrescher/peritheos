"""Human-oriented JSON layout, without changing any EOSMAT field values."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any, NamedTuple

_MATERIAL_KEYS = (
    "format",
    "format_version",
    "identifier",
    "name",
    "formula",
    "phase",
    "aliases",
    "family_id",
    "units",
    "formula_units_per_cell",
    "cell_contents",
    "symmetry",
    "space_group",
    "space_group_number",
    "lattice",
    "atom_sites",
    "peaks",
    "notes",
    "recommendation_hold",
    "eos_records",
    "datasets",
    "source",
)
_RECORD_KEYS = (
    "identifier",
    "model",
    "parameters",
    "equation_options",
    "label",
    "aliases",
    "record_kind",
    "determination_method",
    "equation_kind",
    "default_for",
    "default",
    "catalog_access",
    "pressure_scale_variant_of",
    "reference",
    "temperature_ref",
    "temperature_ref_provenance",
    "eos",
    "thermal",
    "parameter_errors",
    "parameter_error_confidence",
    "fixed_parameters",
    "parameter_covariance",
    "volume_basis",
    "volume",
    "experimental_pressure_range_gpa",
    "pressure_range_status",
    "experimental_temperature_range_k",
    "validity",
    "loading_path",
    "branch_kind",
    "initial_state",
    "branch_domain",
    "notes",
    "pressure_calibration",
    "scientific_validation",
    "parameter_provenance",
    "source_lineage",
    "derived_from_record",
    "fit_provenance",
    "derivation",
    "evidence",
)

# EOSMAT v4 records (and v5 records written before the v5 reading-order
# revision) put the information needed to read/use the EOS ahead of long
# citations and audit history. Unrecognized extensions stay ahead of that tail.
_V4_RECORD_KEYS = (
    "identifier",
    "label",
    "model",
    "temperature_ref",
    "reference_pressure_gpa",
    "parameters",
    "equation_options",
    "thermal",
    "parameter_error_confidence",
    "parameter_error_volume_basis",
    "experimental_pressure_range_gpa",
    "experimental_temperature_range_k",
    "pressure_range_status",
    "validity",
    "determination_method",
    "record_kind",
    "equation_kind",
    "default_for",
    "default",
    "catalog_access",
    "recommendation_hold",
    "phase_scope",
    "phase_selection",
    "loading_path",
    "branch_kind",
    "initial_state",
    "volume_basis",
    "volume",
    "branch_domain",
    "sample_composition",
    "experimental_configuration",
    "computational_configuration",
    "notes",
    "range_provenance",
    "temperature_ref_provenance",
    "pressure_scale_variant_of",
    "derived_from_record",
    "aliases",
)
_V4_RECORD_TAIL = (
    "parameter_covariance",
    "reference",
    "check",
    "pressure_calibration",
    "evidence",
    "superseded_record_text",
)

# Native v5 records read in this order: what the record is and how it was
# checked, the model and its reference state, coefficients, error convention
# and covariance, the stated envelope, how the values were determined, notes,
# linked datasets and the citation. The remaining keys keep their earlier
# relative order. Unrecognized extensions stay ahead of the audit tail: a
# source-specific scientific qualification must not disappear below it.
_V5_RECORD_KEYS = (
    "identifier",
    "label",
    "check",
    "model",
    "temperature_ref",
    "reference_pressure_gpa",
    "parameters",
    "equation_options",
    "thermal",
    "parameter_error_confidence",
    "parameter_error_volume_basis",
    "parameter_covariance",
    "validity",
    "experimental_pressure_range_gpa",
    "experimental_temperature_range_k",
    "pressure_range_status",
    "range_provenance",
    "determination_method",
    "record_kind",
    "notes",
    "fit_datasets",
    "supporting_datasets",
    "validation_datasets",
    "comparison_datasets",
    "diagnostic_datasets",
    "analysis_datasets",
    "context_datasets",
    "reference",
    "equation_kind",
    "default_for",
    "default",
    "catalog_access",
    "recommendation_hold",
    "phase_scope",
    "phase_selection",
    "loading_path",
    "branch_kind",
    "initial_state",
    "volume_basis",
    "volume",
    "branch_domain",
    "sample_composition",
    "experimental_configuration",
    "computational_configuration",
    "temperature_ref_provenance",
    "pressure_scale_variant_of",
    "derived_from_record",
    "aliases",
)
_V5_RECORD_TAIL = (
    "pressure_calibration",
    "evidence",
    "superseded_record_text",
)
_V4_MATERIAL_FIRST = (
    "format_version",
    "identifier",
    "name",
    "formula",
    "phase",
    "units",
    "formula_units_per_cell",
    "recommendation_hold",
    "eos_records",
)
_V5_MATERIAL_FIRST = (
    "format_version",
    "identifier",
    "name",
    "formula",
    "phase",
    "units",
    "formula_units_per_cell",
    "recommendation_hold",
    "default_records",
    "eos_records",
)
_CHECK_KEYS = (
    "reproduction",
    "compared",
    "why",
    "warning",
    "source_check",
    "source_checked_on",
    "source_note",
    "audit",
)
_THERMAL_KEYS = (
    "model",
    "parameters",
    "thermal_pressure_reference",
    "debye_temperature_law",
    "thermal_expansion_law",
    "reference_volume_law",
    "bulk_modulus_law",
)
_ASSESSMENT_KEYS = (
    "status",
    "scope",
    "source_status",
    "basis",
    "basis_code",
    "qualifications",
    "report_url",
    "report_path",
    "parameters",
    "fit",
    "report",
)
_COMPARISON_KEYS = (
    "parameter",
    "published",
    "published_error",
    "published_sigma",
    "refit",
    "refit_error",
    "refit_sigma",
    "difference",
    "delta",
    "relative_difference",
    "relative_delta",
)


class Layout(NamedTuple):
    """Key-order tables of one native layout (v4/v5 records, check blocks)."""

    material_first: tuple[str, ...]
    record_keys: tuple[str, ...]
    record_tail: tuple[str, ...]
    check_keys: tuple[str, ...]


#: The EOSMAT v5 layout written by ``save_eosmat`` and ``format``.
V5_LAYOUT = Layout(
    material_first=_V5_MATERIAL_FIRST,
    record_keys=_V5_RECORD_KEYS,
    record_tail=_V5_RECORD_TAIL,
    check_keys=_CHECK_KEYS,
)
#: The EOSMAT v4 layout. Bundled v5 files used it before the v5 reading-order
#: revision; repository receipts use it to recover those files byte for byte.
V4_LAYOUT = Layout(
    material_first=_V4_MATERIAL_FIRST,
    record_keys=_V4_RECORD_KEYS,
    record_tail=_V4_RECORD_TAIL,
    check_keys=_CHECK_KEYS,
)


def format_eosmat(document: Mapping[str, Any], *, layout: Layout | None = None) -> str:
    """Format a validated document; retain array order and extension metadata.

    Material/record fields, native thermal components and recognized reconstruction
    assessment/comparison fields are reordered.
    Parameter maps, covariance metadata, arrays and arbitrary nested extensions
    retain their order. Historical v2/v3 layout is unchanged. Short scalar groups
    fit on one line; long strings remain intact even above 120 columns.
    ``layout`` defaults to :data:`V5_LAYOUT` for v5 and :data:`V4_LAYOUT`
    otherwise.
    """
    if layout is None:
        layout = V5_LAYOUT if document.get("format_version") == 5 else V4_LAYOUT
    return _format(document, audit=False, layout=layout)


def format_audit(sidecar: Mapping[str, Any]) -> str:
    """Format an audit sidecar; each entry uses the native record ordering."""
    return _format(sidecar, audit=True, layout=V5_LAYOUT)


def _format(document: Mapping[str, Any], *, audit: bool, layout: Layout) -> str:
    """Render ``document``; ``layout`` selects the native key-order tables."""
    native = audit or document.get("format_version") in (4, 5)

    def scalar(value: Any) -> bool:
        return not isinstance(value, (Mapping, list, tuple))

    def encode(value: Any) -> str:
        return json.dumps(value, ensure_ascii=False, allow_nan=False)

    def render(value: Any, depth: int, context: str = "") -> str:
        if scalar(value):
            return encode(value)
        if isinstance(value, Mapping) and any(
            not isinstance(key, str) for key in value
        ):
            raise TypeError("EOSMAT JSON object keys must be strings")
        if not value:
            return "{}" if isinstance(value, Mapping) else "[]"
        if context not in {"material", "record", "thermal", "audit", "audit_records"}:
            items = value.values() if isinstance(value, Mapping) else value
            if all(scalar(item) for item in items):
                compact = encode(dict(value) if isinstance(value, Mapping) else value)
                if depth * 2 + len(compact) <= 120:
                    return compact
        indent = "  " * depth
        inner = indent + "  "
        if isinstance(value, Mapping):
            priority: tuple[str, ...] = (
                _MATERIAL_KEYS
                if context == "material"
                else _RECORD_KEYS
                if context == "record"
                else ()
            )
            if context == "material" and native:
                first: tuple[str, ...] = layout.material_first
                priority = (
                    first
                    + tuple(key for key in priority if key not in first)
                    + ("shared_metadata",)
                )
            if context == "record" and native:
                priority = layout.record_keys
            elif context == "record" and "model" in value:
                # Retain the historical formatter's treatment of extensions.
                first = (
                    "identifier",
                    "label",
                    "model",
                    "temperature_ref",
                    "parameters",
                    "thermal",
                )
                priority = first + tuple(key for key in priority if key not in first)
            if context == "thermal":
                priority = _THERMAL_KEYS
            if context == "check":
                priority = layout.check_keys
            if context == "evidence":
                priority = ("reproduction",)
            if context == "reproduction":
                priority = ("assessment",)
            if context == "assessment":
                priority = _ASSESSMENT_KEYS
            if context == "comparison":
                priority = _COMPARISON_KEYS
            keys = [key for key in priority if key in value]
            tail = layout.record_tail if context == "record" and native else ()
            keys.extend(key for key in value if key not in priority and key not in tail)
            keys.extend(key for key in tail if key in value)
            lines = []
            for key in keys:
                child = key if context == "material" and key == "eos_records" else ""
                if context == "audit" and key == "records":
                    child = "audit_records"
                if context == "audit_records":
                    child = "record"
                if context == "record" and native and key == "thermal":
                    child = "thermal"
                if context == "record" and native and key == "evidence":
                    child = "evidence"
                if context == "record" and native and key == "check":
                    child = "check"
                if context == "evidence" and key == "reproduction":
                    child = "reproduction"
                if context == "reproduction" and key == "assessment":
                    child = "assessment"
                if context == "assessment" and key == "parameters":
                    child = "comparisons"
                lines.append(
                    inner + encode(key) + ": " + render(value[key], depth + 1, child)
                )
            return "{\n" + ",\n".join(lines) + "\n" + indent + "}"
        child = (
            "record"
            if context == "eos_records"
            else "comparison"
            if context == "comparisons"
            else ""
        )
        return (
            "[\n"
            + ",\n".join(inner + render(item, depth + 1, child) for item in value)
            + "\n"
            + indent
            + "]"
        )

    return render(document, 0, "audit" if audit else "material") + "\n"
