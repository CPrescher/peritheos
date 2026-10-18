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


def format_eosmat(document: Mapping[str, Any]) -> str:
    """Format a validated document; retain array order and extension metadata.

    Only material and record fields are reordered. Other mappings retain their
    order, including parameter and covariance metadata. Short scalar groups fit
    on one line; long strings remain intact even when they exceed 120 columns.
    """

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
        if context not in {"material", "record"}:
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
            if context == "material" and value.get("format_version") == 4:
                first = (
                    "format_version",
                    "identifier",
                    "name",
                    "formula",
                    "phase",
                    "units",
                    "formula_units_per_cell",
                    "eos_records",
                )
                priority = (
                    first
                    + tuple(key for key in priority if key not in first)
                    + ("shared_metadata",)
                )
            if context == "record" and "model" in value:
                first = (
                    "identifier",
                    "label",
                    "model",
                    "temperature_ref",
                    "parameters",
                    "thermal",
                )
                priority = first + tuple(key for key in priority if key not in first)
            keys = [key for key in priority if key in value]
            keys.extend(key for key in value if key not in priority)
            lines = []
            for key in keys:
                child = key if context == "material" and key == "eos_records" else ""
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
