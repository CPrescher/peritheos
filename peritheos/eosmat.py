"""Shared Peritheos/Dioptas ``.eosmat`` material documents.

The format is owned by Peritheos and deliberately includes optional
crystallographic fields that Peritheos preserves but does not interpret.
Dioptas 0.10.0 can read the format-2 document shape directly.
"""

from __future__ import annotations

import copy
import json
from collections.abc import Mapping
from importlib import resources
from pathlib import Path
from typing import Any

from peritheos.errors import EosmatError, MaterialLookupError

EOSMAT_FORMAT = "peritheos.material"
EOSMAT_FORMAT_VERSION = 3
_MATERIAL_PACKAGE = "peritheos.data.materials"
_RT_TYPES = {
    "BM2",
    "BM3",
    "BM4",
    "Holzapfel",
    "ModifiedTait",
    "Murnaghan",
    "NaturalStrain2",
    "NaturalStrain3",
    "NaturalStrain4",
    "Vinet",
}
_THERMAL_TYPES = {
    "AlphaKT",
    "AsymptoticPowerLawMieGruneisenDebye",
    "DoubleDebyeHelmholtz",
    "DoubleDebyeLogMomentHelmholtz",
    "LinearThermalPressure",
    "LogVolumeThermalPressure",
    "MieGruneisenDebye",
    "MieGruneisenEinstein",
    "MultiOscillatorGruneisen",
    "Sokolova2016",
    "ThermalModifiedTait",
}
_RT_MODELS = {
    "BM2": "birch_murnaghan_2",
    "BM3": "birch_murnaghan_3",
    "BM4": "birch_murnaghan_4",
    "Holzapfel": "holzapfel",
    "ModifiedTait": "modified_tait",
    "Murnaghan": "murnaghan",
    "NaturalStrain2": "natural_strain_2",
    "NaturalStrain3": "natural_strain_3",
    "NaturalStrain4": "natural_strain_4",
    "Vinet": "vinet",
}
_THERMAL_MODELS = {
    "AlphaKT": "thermal_reference_state",
    "AsymptoticPowerLawMieGruneisenDebye": ("asymptotic_power_law_mie_gruneisen_debye"),
    "DoubleDebyeHelmholtz": "double_debye_helmholtz",
    "DoubleDebyeLogMomentHelmholtz": "double_debye_log_moment_helmholtz",
    "LinearThermalPressure": "linear_thermal_pressure",
    "LogVolumeThermalPressure": "log_volume_thermal_pressure",
    "MieGruneisenDebye": "mie_gruneisen_debye",
    "MieGruneisenEinstein": "mie_gruneisen_einstein",
    "MultiOscillatorGruneisen": ("multi_oscillator_gruneisen_thermal_pressure"),
    "Sokolova2016": "multi_oscillator_gruneisen_thermal_pressure",
    "ThermalModifiedTait": "thermal_modified_tait",
}


def _require_mapping(value: Any, location: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise EosmatError(f"{location} must be a JSON object")
    return value


def _finite_number(value: Any, location: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EosmatError(f"{location} must be a number")
    number = float(value)
    if not (-float("inf") < number < float("inf")):
        raise EosmatError(f"{location} must be finite")
    if positive and number <= 0.0:
        raise EosmatError(f"{location} must be greater than zero")
    return number


def validate_eosmat_document(document: Mapping[str, Any]) -> None:
    """Validate the shared format-2 structure without executing EOS code.

    The absent ``format`` discriminator used by native Dioptas 0.10.0 files is
    accepted for backward compatibility. Scientific validation status is
    metadata and is not inferred from structural validity.
    """
    document = _require_mapping(document, "document")
    format_name = document.get("format")
    version = document.get("format_version")
    is_legacy_dioptas = format_name is None and version == 2
    is_canonical = format_name == EOSMAT_FORMAT and version == EOSMAT_FORMAT_VERSION
    if not (is_legacy_dioptas or is_canonical):
        if format_name not in (None, EOSMAT_FORMAT):
            raise EosmatError(f"Unsupported material format {format_name!r}")
        raise EosmatError(
            "Supported eosmat documents are Peritheos format 3 or legacy "
            "Dioptas format 2"
        )

    for key in ("name", "formula"):
        if not isinstance(document.get(key), str):
            raise EosmatError(f"{key} must be a string")
    for key in ("phase", "symmetry", "notes"):
        if key in document and not isinstance(document[key], str):
            raise EosmatError(f"{key} must be a string")
    if document.get("cell_contents") is not None and not isinstance(
        document["cell_contents"], str
    ):
        raise EosmatError("cell_contents must be a string")
    if "identifier" in document and not isinstance(document["identifier"], str):
        raise EosmatError("identifier must be a string")
    if is_canonical:
        if not document.get("identifier"):
            raise EosmatError("Canonical eosmat requires a material identifier")
        units = _require_mapping(document.get("units"), "units")
        if units != {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        }:
            raise EosmatError("Canonical eosmat units must be GPa, K, and cell Å³")

    if "lattice" in document:
        lattice = _require_mapping(document["lattice"], "lattice")
        for key in ("a", "alpha", "beta", "gamma"):
            _finite_number(lattice.get(key), f"lattice.{key}")
        for key in ("b", "c"):
            if lattice.get(key) is not None:
                _finite_number(lattice[key], f"lattice.{key}")

    formula_units = document.get("formula_units_per_cell")
    if formula_units is not None:
        _finite_number(formula_units, "formula_units_per_cell", positive=True)
    for key in ("aliases", "atom_sites", "peaks", "eos_records"):
        if not isinstance(document.get(key, []), list):
            raise EosmatError(f"{key} must be a JSON array")

    record_identifiers: set[str] = set()
    default_count = 0
    for index, raw_record in enumerate(document.get("eos_records", [])):
        location = f"eos_records[{index}]"
        record = _require_mapping(raw_record, location)
        if not isinstance(record.get("label"), str):
            raise EosmatError(f"{location}.label must be a string")
        identifier = record.get("identifier")
        if is_canonical and identifier is None:
            raise EosmatError(f"{location} requires an identifier")
        if identifier is not None:
            if not isinstance(identifier, str) or not identifier:
                raise EosmatError(f"{location}.identifier must be a non-empty string")
            if identifier in record_identifiers:
                raise EosmatError(f"Duplicate EOS record identifier {identifier!r}")
            record_identifiers.add(identifier)
        default_count += record.get("default") is True

        reference = record.get("reference")
        if not isinstance(reference, (str, Mapping)):
            raise EosmatError(f"{location}.reference must be a string or object")
        if isinstance(reference, Mapping):
            if not isinstance(reference.get("authors"), list):
                raise EosmatError(f"{location}.reference.authors must be an array")
            if isinstance(reference.get("year"), bool) or not isinstance(
                reference.get("year"), int
            ):
                raise EosmatError(f"{location}.reference.year must be an integer")

        eos = _require_mapping(record.get("eos"), f"{location}.eos")
        eos_type = eos.get("type")
        if eos_type not in _RT_TYPES:
            raise EosmatError(f"{location}.eos.type {eos_type!r} is unsupported")
        if is_canonical and eos.get("model") != _RT_MODELS[eos_type]:
            raise EosmatError(f"{location}.eos.model does not match eos.type")
        parameters = _require_mapping(
            eos.get("parameters"), f"{location}.eos.parameters"
        )
        for name, value in parameters.items():
            _finite_number(value, f"{location}.eos.parameters.{name}")
        if "V0" not in parameters:
            raise EosmatError(f"{location}.eos.parameters requires V0")

        errors = _require_mapping(
            record.get("parameter_errors"), f"{location}.parameter_errors"
        )
        for name, value in errors.items():
            if value is not None:
                _finite_number(value, f"{location}.parameter_errors.{name}")
        error_confidence = record.get("parameter_error_confidence")
        if error_confidence is not None:
            error_confidence = _finite_number(
                error_confidence, f"{location}.parameter_error_confidence"
            )
            if not 0.0 < error_confidence < 1.0:
                raise EosmatError(
                    f"{location}.parameter_error_confidence must lie between zero "
                    "and one"
                )
        if not isinstance(record.get("fixed_parameters"), list):
            raise EosmatError(f"{location}.fixed_parameters must be an array")

        volume = record.get("volume")
        if volume is not None:
            volume = _require_mapping(volume, f"{location}.volume")
            if volume.get("reference_value") is not None:
                _finite_number(
                    volume["reference_value"],
                    f"{location}.volume.reference_value",
                    positive=True,
                )
            if volume.get("public_to_model_scale") is not None:
                _finite_number(
                    volume["public_to_model_scale"],
                    f"{location}.volume.public_to_model_scale",
                    positive=True,
                )
            if volume.get("model_unit") is not None and not isinstance(
                volume["model_unit"], str
            ):
                raise EosmatError(f"{location}.volume.model_unit must be a string")

        covariance = record.get("parameter_covariance")
        if covariance is not None:
            covariance = _require_mapping(
                covariance, f"{location}.parameter_covariance"
            )
            matrix = covariance.get("matrix")
            order = covariance.get("parameter_order")
            if not isinstance(matrix, list) or not matrix:
                raise EosmatError(
                    f"{location}.parameter_covariance.matrix must be a non-empty array"
                )
            if not isinstance(order, list) or not all(
                isinstance(name, str) and name for name in order
            ):
                raise EosmatError(
                    f"{location}.parameter_covariance.parameter_order must be "
                    "an array of non-empty strings"
                )
            if len(matrix) != len(order) or any(
                not isinstance(row, list) or len(row) != len(order) for row in matrix
            ):
                raise EosmatError(
                    f"{location}.parameter_covariance matrix must be square and "
                    "match parameter_order"
                )
            for row_index, row in enumerate(matrix):
                for column_index, value in enumerate(row):
                    _finite_number(
                        value,
                        f"{location}.parameter_covariance.matrix"
                        f"[{row_index}][{column_index}]",
                    )

        thermal = record.get("thermal")
        if thermal is not None:
            thermal = _require_mapping(thermal, f"{location}.thermal")
            thermal_type = thermal.get("type")
            if thermal_type not in _THERMAL_TYPES:
                raise EosmatError(
                    f"{location}.thermal.type {thermal_type!r} is unsupported"
                )
            if is_canonical and thermal.get("model") != _THERMAL_MODELS[thermal_type]:
                raise EosmatError(
                    f"{location}.thermal.model does not match thermal.type"
                )
            debye_temperature_law = thermal.get("debye_temperature_law")
            if thermal_type == "MieGruneisenDebye":
                if debye_temperature_law is not None and (
                    not isinstance(debye_temperature_law, str)
                    or debye_temperature_law
                    not in {"integrated_gruneisen", "variable_exponent"}
                ):
                    raise EosmatError(
                        f"{location}.thermal.debye_temperature_law is invalid"
                    )
            elif debye_temperature_law is not None:
                raise EosmatError(
                    f"{location}.thermal.debye_temperature_law requires "
                    "MieGruneisenDebye"
                )
            thermal_expansion_law = thermal.get("thermal_expansion_law")
            reference_volume_law = thermal.get("reference_volume_law")
            if thermal_type == "AlphaKT":
                if thermal_expansion_law is not None and (
                    not isinstance(thermal_expansion_law, str)
                    or thermal_expansion_law not in {"constant", "linear_temperature"}
                ):
                    raise EosmatError(
                        f"{location}.thermal.thermal_expansion_law is invalid"
                    )
                if reference_volume_law is not None and (
                    not isinstance(reference_volume_law, str)
                    or reference_volume_law
                    not in {"integrated_expansivity", "linear_temperature"}
                ):
                    raise EosmatError(
                        f"{location}.thermal.reference_volume_law is invalid"
                    )
            elif thermal_expansion_law is not None:
                raise EosmatError(
                    f"{location}.thermal.thermal_expansion_law requires AlphaKT"
                )
            elif reference_volume_law is not None:
                raise EosmatError(
                    f"{location}.thermal.reference_volume_law requires AlphaKT"
                )
            thermal_parameters = _require_mapping(
                thermal.get("parameters"), f"{location}.thermal.parameters"
            )
            for name, value in thermal_parameters.items():
                _finite_number(value, f"{location}.thermal.parameters.{name}")
            if thermal_expansion_law == "linear_temperature" and (
                "alpha1" not in thermal_parameters
            ):
                raise EosmatError(
                    f"{location}.thermal.parameters requires alpha1 for "
                    "linear_temperature thermal expansion"
                )
            if (
                thermal_expansion_law in {None, "constant"}
                and thermal_parameters.get("alpha1", 0.0) != 0.0
            ):
                raise EosmatError(
                    f"{location}.thermal.parameters.alpha1 must be zero for "
                    "constant thermal expansion"
                )
            if reference_volume_law == "linear_temperature" and (
                thermal_expansion_law not in {None, "constant"}
                or thermal_parameters.get("alpha1", 0.0) != 0.0
            ):
                raise EosmatError(
                    f"{location}.thermal linear_temperature reference volume "
                    "requires constant thermal expansion and alpha1=0"
                )
            thermal_errors = _require_mapping(
                thermal.get("parameter_errors", {}),
                f"{location}.thermal.parameter_errors",
            )
            for name, value in thermal_errors.items():
                if value is not None:
                    _finite_number(value, f"{location}.thermal.parameter_errors.{name}")
            if not isinstance(thermal.get("fixed_parameters", []), list):
                raise EosmatError(
                    f"{location}.thermal.fixed_parameters must be an array"
                )

        for range_name in (
            "experimental_pressure_range_gpa",
            "experimental_temperature_range_k",
        ):
            values = record.get(range_name)
            if values is not None:
                if not isinstance(values, list) or len(values) != 2:
                    raise EosmatError(f"{location}.{range_name} must have two values")
                low = _finite_number(values[0], f"{location}.{range_name}[0]")
                high = _finite_number(values[1], f"{location}.{range_name}[1]")
                if low > high:
                    raise EosmatError(f"{location}.{range_name} must be ordered")

        validity = record.get("validity")
        if validity is not None:
            validity = _require_mapping(validity, f"{location}.validity")
            for range_name in ("pressure_gpa", "temperature_k", "volume_ratio"):
                values = validity.get(range_name)
                if values is None:
                    continue
                if not isinstance(values, list) or len(values) != 2:
                    raise EosmatError(
                        f"{location}.validity.{range_name} must have two values"
                    )
                low = _finite_number(values[0], f"{location}.validity.{range_name}[0]")
                high = _finite_number(values[1], f"{location}.validity.{range_name}[1]")
                if low > high:
                    raise EosmatError(
                        f"{location}.validity.{range_name} must be ordered"
                    )
            if not isinstance(validity.get("notes", []), list):
                raise EosmatError(f"{location}.validity.notes must be an array")

        validation = record.get("scientific_validation")
        if is_canonical:
            validation = _require_mapping(
                validation, f"{location}.scientific_validation"
            )
            if validation.get("status") not in {
                "primary_source_validated",
                "pending_primary_source_check",
                "deferred",
            }:
                raise EosmatError(f"{location}.scientific_validation.status is invalid")

    if default_count > 1:
        raise EosmatError("A material may have at most one default EOS record")


def load_eosmat(path: str | Path) -> dict[str, Any]:
    """Load and structurally validate a Peritheos or Dioptas `.eosmat` file."""
    source = Path(path)
    with source.open(encoding="utf-8") as stream:
        try:
            document = json.load(stream)
        except json.JSONDecodeError as error:
            raise EosmatError(
                f"Invalid eosmat JSON at line {error.lineno}, column {error.colno}: "
                f"{error.msg}",
                code="eosmat.json",
                operation="load",
                context={
                    "path": str(source),
                    "line": error.lineno,
                    "column": error.colno,
                },
            ) from error
    validate_eosmat_document(document)
    return document


def save_eosmat(path: str | Path, document: Mapping[str, Any]) -> None:
    """Validate and save a Dioptas-compatible `.eosmat` document."""
    validate_eosmat_document(document)
    with Path(path).open("w", encoding="utf-8") as stream:
        json.dump(document, stream, indent=1, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def list_material_documents() -> tuple[str, ...]:
    """Return all material identifiers in the bundled shared library."""
    root = resources.files(_MATERIAL_PACKAGE)
    return tuple(
        sorted(
            path.name.removesuffix(".eosmat")
            for path in root.iterdir()
            if path.name.endswith(".eosmat")
        )
    )


def get_material_document(identifier: str) -> dict[str, Any]:
    """Return a defensive copy of one bundled `.eosmat` document."""
    if identifier not in list_material_documents():
        raise MaterialLookupError(
            f"Unknown material document {identifier!r}; available: "
            f"{list_material_documents()}"
        )
    resource = resources.files(_MATERIAL_PACKAGE).joinpath(f"{identifier}.eosmat")
    document = json.loads(resource.read_text(encoding="utf-8"))
    validate_eosmat_document(document)
    return copy.deepcopy(document)


def eosmat_schema() -> dict[str, Any]:
    """Return the bundled normative JSON Schema document."""
    resource = resources.files("peritheos.data").joinpath("eosmat-v3.schema.json")
    return json.loads(resource.read_text(encoding="utf-8"))


__all__ = [
    "EOSMAT_FORMAT",
    "EOSMAT_FORMAT_VERSION",
    "eosmat_schema",
    "get_material_document",
    "list_material_documents",
    "load_eosmat",
    "save_eosmat",
    "validate_eosmat_document",
]
