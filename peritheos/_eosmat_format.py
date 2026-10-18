"""Human-oriented JSON layout, without changing any EOSMAT field values."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

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

# Native records put the information needed to read/use the EOS ahead of long
# citations and audit history. Unrecognized extensions stay ahead of that tail
# too: a source-specific scientific qualification must not disappear below it.
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
    "pressure_calibration",
    "evidence",
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


def format_eosmat(document: Mapping[str, Any]) -> str:
    """Format a validated document; retain array order and extension metadata.

    Material/record fields and native thermal-component fields are reordered.
    Parameter maps, covariance metadata, arrays and arbitrary nested extensions
    retain their order. Historical v2/v3 layout is unchanged. Short scalar groups
    fit on one line; long strings remain intact even above 120 columns.
    """

    native = document.get("format_version") == 4

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
        if context not in {"material", "record", "thermal"}:
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
                first: tuple[str, ...] = (
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
                priority = (
                    first
                    + tuple(key for key in priority if key not in first)
                    + ("shared_metadata",)
                )
            if context == "record" and native:
                priority = _V4_RECORD_KEYS
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
            keys = [key for key in priority if key in value]
            tail = _V4_RECORD_TAIL if context == "record" and native else ()
            keys.extend(key for key in value if key not in priority and key not in tail)
            keys.extend(key for key in tail if key in value)
            lines = []
            for key in keys:
                child = key if context == "material" and key == "eos_records" else ""
                if context == "record" and native and key == "thermal":
                    child = "thermal"
                lines.append(
                    inner + encode(key) + ": " + render(value[key], depth + 1, child)
                )
            return "{\n" + ",\n".join(lines) + "\n" + indent + "}"
        child = "record" if context == "eos_records" else ""
        return (
            "[\n"
            + ",\n".join(inner + render(item, depth + 1, child) for item in value)
            + "\n"
            + indent
            + "]"
        )

    return render(document, 0, "material") + "\n"
