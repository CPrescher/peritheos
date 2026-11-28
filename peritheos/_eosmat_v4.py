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
        "parameter_error_volume_basis",
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
        "check",
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


DEFAULT_CATEGORIES = ("equilibrium", "hugoniot")
_HUGONIOT_MODELS = frozenset({"LinearUsUpHugoniot"})


def default_category(record: Mapping[str, Any]) -> str:
    """``hugoniot`` for a Hugoniot model, otherwise ``equilibrium``."""
    eos = record.get("eos")
    model = eos.get("type") if isinstance(eos, Mapping) else record.get("model")
    return "hugoniot" if model in _HUGONIOT_MODELS else "equilibrium"


def _flags_default(record: Mapping[str, Any]) -> bool:
    return record.get("default") is True or record.get("default_for") is not None


def expand_default_records(document: dict[str, Any]) -> dict[str, Any]:
    """Write a material ``default_records`` map as per-record flags, in place.

    The named record gets ``default_for: <category>``, the form historical
    v3/v4 readers and Material exports use; the map is removed. Per-record
    flags that are already present must agree with the map: a named record
    flagged ``default: false`` and any other record flagged as a default are
    errors, never silently resolved.
    """
    if "default_records" not in document:
        return document
    mapping = document.pop("default_records")
    if not isinstance(mapping, Mapping) or not mapping:
        raise EosmatError("default_records must be a non-empty object")
    unknown = set(mapping) - set(DEFAULT_CATEGORIES)
    if unknown:
        raise EosmatError(
            f"default_records keys must be {' or '.join(DEFAULT_CATEGORIES)}; "
            f"got {sorted(unknown)}"
        )
    records = document.get("eos_records")
    if not isinstance(records, list):
        raise EosmatError("eos_records must be an array")
    by_identifier = {
        record.get("identifier"): record
        for record in records
        if isinstance(record, Mapping)
    }
    for category, identifier in mapping.items():
        if not isinstance(identifier, str) or not identifier:
            raise EosmatError(f"default_records.{category} must be a record identifier")
        record = by_identifier.get(identifier)
        if record is None:
            raise EosmatError(
                f"default_records.{category} names unknown record {identifier!r}"
            )
        if default_category(record) != category:
            raise EosmatError(
                f"default_records.{category} names {identifier!r}, which is not a "
                f"{category} record"
            )
    for record in records:
        if not isinstance(record, dict):
            continue
        category = default_category(record)
        if mapping.get(category) == record.get("identifier"):
            if record.get("default") is False:
                raise EosmatError(
                    f"default_records.{category} names {record['identifier']!r}, "
                    "which is flagged default: false"
                )
            record["default_for"] = category
        elif _flags_default(record):
            named = mapping.get(category)
            raise EosmatError(
                f"Record {record.get('identifier')!r} is flagged as the default "
                f"{category} record, but default_records names "
                + (repr(named) if named else f"no {category} record")
            )
    return document


def write_explicit_default_flags(document: dict[str, Any]) -> dict[str, Any]:
    """Give every record an explicit ``default`` boolean, in place.

    Used for v3/v4 exports of a v5 document, so historical readers see the
    selection on every record: the default of each category gets ``default:
    true`` and ``default_for: <category>``, every other record ``default:
    false``. A ``default_records`` map is expanded first.
    """
    expand_default_records(document)
    for record in document.get("eos_records", []):
        if not isinstance(record, dict):
            continue
        category = default_category(record)
        selected = (
            record.get("default") is True or record.get("default_for") == category
        )
        record["default"] = selected
        if selected:
            record["default_for"] = category
    return document


def collect_default_records(document: dict[str, Any]) -> dict[str, Any]:
    """Move per-record default flags into the material ``default_records`` map.

    In place. Every ``default`` and ``default_for`` record key is removed;
    records flagged ``default: true`` or ``default_for`` are named in the map,
    which is created only when a record is flagged. An existing map must
    agree with the flags.
    """
    records = document.get("eos_records", [])
    mapping: dict[str, str] = {}
    if "default_records" in document:
        expand_default_records(document)
    for record in records:
        if not isinstance(record, dict):
            continue
        if _flags_default(record):
            category = default_category(record)
            if (
                mapping.setdefault(category, record["identifier"])
                != record["identifier"]
            ):
                raise EosmatError(
                    f"A material may have at most one default {category} EOS record"
                )
        record.pop("default", None)
        record.pop("default_for", None)
    if mapping:
        document["default_records"] = {
            category: mapping[category]
            for category in DEFAULT_CATEGORIES
            if category in mapping
        }
    return document


def to_v3(
    document: Mapping[str, Any], *, explicit_defaults: bool = False
) -> dict[str, Any]:
    """Expand a native v4 document for historical consumers, without inference.

    With ``explicit_defaults``, a v5 document's default selection is written
    as an explicit per-record flag (see :func:`write_explicit_default_flags`).
    """
    from peritheos._eosmat_audit import (
        validate_check,
        validate_check_consistency,
        validate_check_text,
        validate_superseded_material_text,
        validate_superseded_record_text,
    )
    from peritheos._eosmat_metadata import expand_metadata

    result = expand_metadata(_object(document, "document"))
    if "superseded_material_text" in result:
        validate_superseded_material_text(result["superseded_material_text"])
    native_version = result.get("format_version")
    if native_version not in (4, 5):
        raise EosmatError("Expected native EOSMAT format_version 4 or 5")
    if "default_records" in result and native_version != 5:
        raise EosmatError("default_records requires EOSMAT format_version 5")
    if explicit_defaults and native_version == 5:
        write_explicit_default_flags(result)
    else:
        expand_default_records(result)
    if native_version == 5:
        raw_records = result.get("eos_records")
        if not isinstance(raw_records, list):
            raise EosmatError("eos_records must be an array")
        for record in raw_records:
            if isinstance(record, Mapping) and (
                "parameter_error_volume_basis" in record
                or isinstance(record.get("parameter_covariance"), Mapping)
                and "volume_basis" in record["parameter_covariance"]
            ):
                raise EosmatError(
                    "EOSMAT v5 uncertainties always use parameter units; remove legacy unit flags"
                )
    if "format" in result:
        raise EosmatError(
            f"EOSMAT v{native_version} uses format_version alone; remove format"
        )
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
        location = f"eos_records[{record.get('identifier', '?')!r}]"
        if "check" in record:
            validate_check(record["check"], location)
            validate_check_consistency(record, location)
        if "superseded_record_text" in record:
            validate_superseded_record_text(record["superseded_record_text"], location)
        validate_check_text(record, document, location)
        evidence = _object(record.pop("evidence", {}), "evidence")
        if set(evidence) & (set(record) | EVIDENCE_RESERVED):
            raise EosmatError("Evidence fields collide with EOS record fields")
        record.update(evidence)
        record.setdefault("label", record.get("identifier", ""))
        record.setdefault("reference", "User-supplied parameters")
        check = record.get("check")
        if (
            "scientific_validation" not in record
            and isinstance(check, Mapping)
            and "source_check" in check
        ):
            # A record file whose evidence lives in an audit sidecar keeps the
            # stored source-check status in its headline (a fixed renaming of
            # the sidecar code); nothing is inferred.
            from peritheos._eosmat_audit import SOURCE_CHECK_CODES

            record["scientific_validation"] = {
                "status": SOURCE_CHECK_CODES[check["source_check"]]
            }
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
    if native_version == 5:
        from peritheos._eosmat_uncertainty import mark_parameter_units

        mark_parameter_units(result)
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
