"""Native EOSMAT v4 structure and explicit migration to the historical v3 layout.

No files or secondary authoring representations are used during conversion.
The numerical engine can reuse its v3-shaped in-memory construction data.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping
from typing import Any

from peritheos.errors import EosmatError

UNITS = {
    "pressure": "GPa",
    "temperature": "K",
    "volume": "angstrom^3/conventional_unit_cell",
}
EVIDENCE_FIELDS = (
    "scientific_validation",
    "parameter_provenance",
    "source_lineage",
    "fit_provenance",
    "fit_protocol",
    "reproduction",
    "audit_corrections",
)

EVIDENCE_RESERVED = frozenset(
    {
        "derivation",
        "notes",
        "aliases",
        "identifier",
        "label",
        "record_kind",
        "determination_method",
        "reference",
        "default",
        "default_for",
        "equation_kind",
        "loading_path",
        "branch_kind",
        "initial_state",
        "volume_basis",
        "branch_domain",
        "derived_from_record",
        "pressure_scale_variant_of",
        "eos",
        "thermal",
        "parameter_errors",
        "fixed_parameters",
        "parameter_error_confidence",
        "parameter_covariance",
        "volume",
        "validity",
        "experimental_pressure_range_gpa",
        "pressure_range_status",
        "experimental_temperature_range_k",
        "temperature_ref",
        "temperature_ref_provenance",
        "pressure_calibration",
        "catalog_access",
        "model",
        "parameters",
        "equation_options",
        "evidence",
    }
)


def _object(value: Any, location: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise EosmatError(f"{location} must be an object")
    return copy.deepcopy(dict(value))


def _models(thermal: bool = False) -> Mapping[str, str]:
    from peritheos.eosmat import _RT_MODELS, _THERMAL_MODELS

    return _THERMAL_MODELS if thermal else _RT_MODELS


def _unpack_parameters(value: Any) -> tuple[dict, dict, list]:
    parameters = _object(value, "parameters")
    values, errors, fixed = {}, {}, []
    for name, parameter in parameters.items():
        if isinstance(parameter, Mapping):
            if set(parameter) - {"value", "error", "fixed"}:
                raise EosmatError(f"Unsupported parameter annotation for {name!r}")
            if "value" in parameter:
                values[name] = copy.deepcopy(parameter["value"])
            if "error" in parameter:
                errors[name] = copy.deepcopy(parameter["error"])
            if "fixed" in parameter:
                if not isinstance(parameter["fixed"], bool):
                    raise EosmatError(f"parameters.{name}.fixed must be boolean")
                if parameter["fixed"]:
                    fixed.append(name)
            if not parameter or set(parameter) == {"fixed"} and not parameter["fixed"]:
                raise EosmatError(
                    f"parameters.{name} has no value or parameter information"
                )
        else:
            values[name] = copy.deepcopy(parameter)
    return values, errors, fixed


def _pack_parameters(values: Any, errors: Any, fixed: Any) -> dict:
    values = _object(values, "parameters")
    errors = _object(errors, "parameter_errors")
    if not isinstance(fixed, list) or not all(isinstance(k, str) for k in fixed):
        raise EosmatError("fixed_parameters must be an array of names")
    if len(fixed) != len(set(fixed)):
        raise EosmatError("Duplicate fixed parameters cannot be migrated")
    names = list(dict.fromkeys([*values, *errors, *fixed]))
    # Object order has no scientific meaning. Preserve the original fixed-array
    # order through the order of its annotations, including metadata-only names.
    fixed_positions = [i for i, key in enumerate(names) if key in fixed]
    for index, name in zip(fixed_positions, fixed):
        names[index] = name
    result: dict[str, Any] = {}
    for name in names:
        annotation: dict[str, Any] = {}
        if name in values:
            annotation["value"] = values[name]
        if name in errors:
            annotation["error"] = errors[name]
        if name in fixed:
            annotation["fixed"] = True
        result[name] = (
            annotation["value"] if set(annotation) == {"value"} else annotation
        )
    return result


def to_v3(document: Mapping[str, Any]) -> dict[str, Any]:
    """Expand a native v4 document for historical consumers, without inference."""
    from peritheos._eosmat_metadata import expand_metadata

    result = expand_metadata(_object(document, "document"))
    if result.get("format_version") != 4:
        raise EosmatError("Expected native EOSMAT format_version 4")
    if "format" in result:
        raise EosmatError("EOSMAT v4 uses format_version alone; remove format")
    if not isinstance(result.get("identifier"), str) or not result["identifier"]:
        raise EosmatError("EOSMAT v4 requires identifier")
    result["format"] = "peritheos.material"
    result["format_version"] = 3
    result.setdefault("name", result["identifier"])
    result.setdefault("units", copy.deepcopy(UNITS))
    records = result.get("eos_records")
    if not isinstance(records, list):
        raise EosmatError("eos_records must be an array")
    expanded = []
    for raw in records:
        record = _object(raw, "EOS record")
        forbidden = {"eos", "parameter_errors", "fixed_parameters", *EVIDENCE_FIELDS}
        if forbidden & set(record):
            raise EosmatError("EOSMAT v4 record mixes native and historical v3 fields")
        model = record.pop("model", None)
        if not isinstance(model, str) or model not in _models():
            raise EosmatError(f"Unknown EOSMAT v4 model {model!r}")
        parameters, errors, fixed = _unpack_parameters(record.pop("parameters", None))
        options = _object(record.pop("equation_options", {}), "equation_options")
        if set(options) & {"type", "model", "parameters"}:
            raise EosmatError("equation_options cannot replace model or parameters")
        record["eos"] = {
            **options,
            "type": model,
            "model": _models()[model],
            "parameters": parameters,
        }
        record["parameter_errors"] = errors
        record["fixed_parameters"] = fixed
        evidence = _object(record.pop("evidence", {}), "evidence")
        if set(evidence) & (set(record) | EVIDENCE_RESERVED):
            raise EosmatError("Evidence fields collide with EOS record fields")
        record.update(evidence)
        record.setdefault("label", record.get("identifier", ""))
        record.setdefault("reference", "User-supplied parameters")
        record.setdefault(
            "scientific_validation", {"status": "pending_primary_source_check"}
        )
        if "thermal" in record:
            thermal = _object(record["thermal"], "thermal")
            if set(thermal) & {"type", "parameter_errors", "fixed_parameters"}:
                raise EosmatError("EOSMAT v4 thermal component mixes v3 fields")
            kind = thermal.get("model")
            if not isinstance(kind, str) or kind not in _models(True):
                raise EosmatError(f"Unknown EOSMAT v4 thermal model {kind!r}")
            p, e, f = _unpack_parameters(thermal.get("parameters"))
            thermal.update(type=kind, model=_models(True)[kind], parameters=p)
            # Empty thermal metadata is optional in v3; no scientific
            # information is lost when those empty containers are omitted.
            if e:
                thermal["parameter_errors"] = e
            if f:
                thermal["fixed_parameters"] = f
            record["thermal"] = thermal
        expanded.append(record)
    result["eos_records"] = expanded
    return result


def to_v4(document: Mapping[str, Any]) -> dict[str, Any]:
    """Migrate a v3 archive to the native v4 structure, retaining scientific data."""
    result = _object(document, "document")
    if (
        result.get("format_version") != 3
        or result.get("format") != "peritheos.material"
    ):
        raise EosmatError("Conversion requires a Peritheos EOSMAT v3 or v4 document")
    result.pop("format")
    result["format_version"] = 4
    if result.get("name") == result.get("identifier"):
        result.pop("name")
    if result.get("units") == UNITS:
        result.pop("units")
    for record in result["eos_records"]:
        if set(record) & {"model", "parameters", "equation_options", "evidence"}:
            raise EosmatError("Existing extension collides with a native v4 field")
        eos = _object(record.pop("eos"), "eos")
        model = eos.pop("type")
        if eos.pop("model") != _models().get(model):
            raise EosmatError("EOS model/type mismatch")
        record["model"] = model
        record["parameters"] = _pack_parameters(
            eos.pop("parameters"),
            record.pop("parameter_errors"),
            record.pop("fixed_parameters"),
        )
        if eos:
            record["equation_options"] = eos
        evidence = {key: record.pop(key) for key in EVIDENCE_FIELDS if key in record}
        if evidence:
            record["evidence"] = evidence
        if record.get("label") == record.get("identifier"):
            record.pop("label")
        if "thermal" in record:
            thermal = _object(record["thermal"], "thermal")
            kind = thermal.pop("type")
            if thermal.pop("model") != _models(True).get(kind):
                raise EosmatError("Thermal model/type mismatch")
            thermal["model"] = kind
            thermal["parameters"] = _pack_parameters(
                thermal["parameters"],
                thermal.pop("parameter_errors", {}),
                thermal.pop("fixed_parameters", []),
            )
            record["thermal"] = thermal
    return result
