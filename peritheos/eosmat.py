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
from peritheos.pressure_calibrations import list_pressure_calibrations

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
_PRESSURE_CALIBRATION_STATUSES = {
    "resolved",
    "partially_resolved",
    "not_applicable",
    "unresolved",
}
_PRESSURE_CALIBRATION_METHODS = {
    "equation_of_state",
    "ruby_fluorescence",
    "other_optical_gauge",
    "shock_wave",
    "ultrasonic",
    "ab_initio",
    "self_consistent",
    "ambient_pressure",
    "other",
}
_PRESSURE_RECALCULATION_STATUSES = {
    "ready",
    "missing_calibrant_observations",
    "reference_eos_not_bundled",
    "reference_model_not_supported",
    "not_applicable",
    "not_possible",
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
    for key in ("aliases", "atom_sites", "peaks", "datasets", "eos_records"):
        if not isinstance(document.get(key, []), list):
            raise EosmatError(f"{key} must be a JSON array")

    dataset_identifiers: set[str] = set()
    dataset_record_links: list[tuple[str, str]] = []
    for index, raw_dataset in enumerate(document.get("datasets", [])):
        location = f"datasets[{index}]"
        dataset = _require_mapping(raw_dataset, location)
        identifier = dataset.get("identifier")
        if not isinstance(identifier, str) or not identifier:
            raise EosmatError(f"{location}.identifier must be a non-empty string")
        if identifier in dataset_identifiers:
            raise EosmatError(f"Duplicate dataset identifier {identifier!r}")
        dataset_identifiers.add(identifier)
        for key in ("kind", "source_location"):
            if not isinstance(dataset.get(key), str) or not dataset[key]:
                raise EosmatError(f"{location}.{key} must be a non-empty string")
        if not isinstance(dataset.get("reference"), (str, Mapping)):
            raise EosmatError(f"{location}.reference must be a string or object")

        columns = dataset.get("columns")
        if not isinstance(columns, list) or not columns:
            raise EosmatError(f"{location}.columns must be a non-empty array")
        column_names: set[str] = set()
        for column_index, raw_column in enumerate(columns):
            column_location = f"{location}.columns[{column_index}]"
            column = _require_mapping(raw_column, column_location)
            name = column.get("name")
            if not isinstance(name, str) or not name:
                raise EosmatError(f"{column_location}.name must be a non-empty string")
            if name in column_names:
                raise EosmatError(f"Duplicate column name {name!r} in {location}")
            column_names.add(name)
            for key in ("quantity", "unit", "role"):
                if not isinstance(column.get(key), str) or not column[key]:
                    raise EosmatError(
                        f"{column_location}.{key} must be a non-empty string"
                    )
            if column["role"] not in {
                "value",
                "standard_deviation",
                "standard_error",
                "bound",
                "flag",
            }:
                raise EosmatError(f"{column_location}.role is invalid")
            if column.get("of") is not None and column["of"] not in column_names:
                raise EosmatError(
                    f"{column_location}.of must reference an earlier column"
                )

        has_rows = "rows" in dataset
        has_resource = "resource" in dataset
        if has_rows == has_resource:
            raise EosmatError(
                f"{location} must contain exactly one of rows or resource"
            )
        if has_rows:
            rows = dataset["rows"]
            if not isinstance(rows, list):
                raise EosmatError(f"{location}.rows must be an array")
            for row_index, row in enumerate(rows):
                if not isinstance(row, list) or len(row) != len(columns):
                    raise EosmatError(
                        f"{location}.rows[{row_index}] must match the column count"
                    )
                for column_index, value in enumerate(row):
                    if value is not None:
                        _finite_number(
                            value,
                            f"{location}.rows[{row_index}][{column_index}]",
                        )
        else:
            resource = _require_mapping(dataset["resource"], f"{location}.resource")
            for key in ("path", "sha256", "media_type"):
                if not isinstance(resource.get(key), str) or not resource[key]:
                    raise EosmatError(
                        f"{location}.resource.{key} must be a non-empty string"
                    )
            path = Path(resource["path"])
            if path.is_absolute() or ".." in path.parts:
                raise EosmatError(
                    f"{location}.resource.path must be relative and local"
                )
            sha256 = resource["sha256"]
            if len(sha256) != 64 or any(
                character not in "0123456789abcdef" for character in sha256
            ):
                raise EosmatError(
                    f"{location}.resource.sha256 must be 64 lowercase hex characters"
                )

        used_by = dataset.get("used_by_eos_records")
        if not isinstance(used_by, list) or not used_by:
            raise EosmatError(
                f"{location}.used_by_eos_records must be a non-empty array"
            )
        for record_identifier in used_by:
            if not isinstance(record_identifier, str) or not record_identifier:
                raise EosmatError(
                    f"{location}.used_by_eos_records must contain non-empty strings"
                )
            dataset_record_links.append((identifier, record_identifier))

    record_identifiers: set[str] = set()
    derived_record_links: list[tuple[str, str]] = []
    fit_dataset_links: list[tuple[str, str]] = []
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

        record_kind = record.get("record_kind", "published")
        if record_kind not in {"published", "refit", "diagnostic"}:
            raise EosmatError(f"{location}.record_kind is invalid")
        derived_from = record.get("derived_from_record")
        if derived_from is not None:
            if not isinstance(derived_from, str) or not derived_from:
                raise EosmatError(
                    f"{location}.derived_from_record must be a non-empty string"
                )
            if identifier is not None:
                derived_record_links.append((identifier, derived_from))
        fit_provenance = record.get("fit_provenance")
        if record_kind == "refit" and (
            derived_from is None or not isinstance(fit_provenance, Mapping)
        ):
            raise EosmatError(
                f"{location} refit records require derived_from_record and "
                "fit_provenance"
            )
        if fit_provenance is not None:
            fit_provenance = _require_mapping(
                fit_provenance, f"{location}.fit_provenance"
            )
            dataset_identifier = fit_provenance.get("dataset")
            if not isinstance(dataset_identifier, str) or not dataset_identifier:
                raise EosmatError(
                    f"{location}.fit_provenance.dataset must be a non-empty string"
                )
            if identifier is not None:
                fit_dataset_links.append((identifier, dataset_identifier))

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
                    not in {
                        "integrated_expansivity",
                        "linear_temperature",
                        "berman",
                    }
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
                parameter_location = f"{location}.thermal.parameters.{name}"
                if (
                    name == "Tr"
                    and value is None
                    and thermal_type
                    in {"DoubleDebyeHelmholtz", "DoubleDebyeLogMomentHelmholtz"}
                ):
                    continue
                _finite_number(
                    value,
                    parameter_location,
                    positive=name == "Tr",
                )
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
            if reference_volume_law == "berman" and (
                thermal_expansion_law != "linear_temperature"
            ):
                raise EosmatError(
                    f"{location}.thermal berman reference volume requires "
                    "linear_temperature thermal expansion"
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

        calibration = record.get("pressure_calibration")
        if calibration is not None:
            calibration = _require_mapping(
                calibration, f"{location}.pressure_calibration"
            )
            if calibration.get("status") not in _PRESSURE_CALIBRATION_STATUSES:
                raise EosmatError(f"{location}.pressure_calibration.status is invalid")
            methods = calibration.get("methods")
            if not isinstance(methods, list):
                raise EosmatError(
                    f"{location}.pressure_calibration.methods must be an array"
                )
            if (
                calibration["status"]
                in {
                    "resolved",
                    "partially_resolved",
                    "not_applicable",
                }
                and not methods
            ):
                raise EosmatError(
                    f"{location}.pressure_calibration.methods must not be empty "
                    f"for status {calibration['status']!r}"
                )
            for method_index, raw_method in enumerate(methods):
                method_location = (
                    f"{location}.pressure_calibration.methods[{method_index}]"
                )
                method = _require_mapping(raw_method, method_location)
                if method.get("kind") not in _PRESSURE_CALIBRATION_METHODS:
                    raise EosmatError(f"{method_location}.kind is invalid")
                if (
                    not isinstance(method.get("source_location"), str)
                    or not method["source_location"]
                ):
                    raise EosmatError(
                        f"{method_location}.source_location must be a non-empty string"
                    )
                method_reference = method.get("reference")
                if method_reference is not None and not isinstance(
                    method_reference, (str, Mapping)
                ):
                    raise EosmatError(
                        f"{method_location}.reference must be a string or object"
                    )
                reference_record = method.get("reference_eos_record")
                if reference_record is not None and (
                    not isinstance(reference_record, str) or not reference_record
                ):
                    raise EosmatError(
                        f"{method_location}.reference_eos_record must be a "
                        "non-empty string"
                    )
                if (
                    reference_record is not None
                    and method["kind"] != "equation_of_state"
                ):
                    raise EosmatError(
                        f"{method_location}.reference_eos_record requires an "
                        "equation_of_state method"
                    )
                reference_calibration = method.get(
                    "reference_calibration_record"
                )
                if reference_calibration is not None and (
                    not isinstance(reference_calibration, str)
                    or not reference_calibration
                ):
                    raise EosmatError(
                        f"{method_location}.reference_calibration_record must be a "
                        "non-empty string"
                    )
                if (
                    reference_calibration is not None
                    and method["kind"] != "ruby_fluorescence"
                ):
                    raise EosmatError(
                        f"{method_location}.reference_calibration_record requires a "
                        "ruby_fluorescence method"
                    )
                if method["kind"] == "equation_of_state" and method_reference is None:
                    raise EosmatError(
                        f"{method_location}.reference is required for an EOS method"
                    )
            recalculation = _require_mapping(
                calibration.get("recalculation"),
                f"{location}.pressure_calibration.recalculation",
            )
            if recalculation.get("status") not in _PRESSURE_RECALCULATION_STATUSES:
                raise EosmatError(
                    f"{location}.pressure_calibration.recalculation.status is invalid"
                )

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
    for dataset_identifier, record_identifier in dataset_record_links:
        if record_identifier not in record_identifiers:
            raise EosmatError(
                f"Dataset {dataset_identifier!r} references unknown EOS record "
                f"{record_identifier!r}"
            )
    for record_identifier, parent_identifier in derived_record_links:
        if parent_identifier == record_identifier:
            raise EosmatError(
                f"EOS record {record_identifier!r} cannot derive from itself"
            )
        if parent_identifier not in record_identifiers:
            raise EosmatError(
                f"EOS record {record_identifier!r} derives from unknown EOS record "
                f"{parent_identifier!r}"
            )
    for record_identifier, dataset_identifier in fit_dataset_links:
        if dataset_identifier not in dataset_identifiers:
            raise EosmatError(
                f"EOS record {record_identifier!r} fit provenance references unknown "
                f"dataset {dataset_identifier!r}"
            )


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


def list_eos_record_documents() -> tuple[str, ...]:
    """Return all globally unique EOS-record identifiers in the library."""
    return tuple(
        sorted(
            record["identifier"]
            for material_identifier in list_material_documents()
            for record in get_material_document(material_identifier)["eos_records"]
        )
    )


def get_eos_record_document(identifier: str) -> dict[str, Any]:
    """Return one executable EOS record from the shared material library."""
    matches = [
        record
        for material_identifier in list_material_documents()
        for record in get_material_document(material_identifier)["eos_records"]
        if record["identifier"] == identifier
    ]
    if len(matches) != 1:
        if not matches:
            raise MaterialLookupError(
                f"Unknown EOS record document {identifier!r}; available: "
                f"{list_eos_record_documents()}"
            )
        raise EosmatError(f"Duplicate EOS record identifier {identifier!r}")
    return copy.deepcopy(matches[0])


def validate_pressure_calibration_references() -> None:
    """Verify that every bundled EOS and calibration link resolves uniquely."""
    identifiers = list_eos_record_documents()
    if len(identifiers) != len(set(identifiers)):
        raise EosmatError("Bundled EOS record identifiers must be globally unique")
    available = set(identifiers)
    available_calibrations = set(list_pressure_calibrations())
    for material_identifier in list_material_documents():
        for record in get_material_document(material_identifier)["eos_records"]:
            for method in record["pressure_calibration"]["methods"]:
                reference_identifier = method.get("reference_eos_record")
                if (
                    reference_identifier is not None
                    and reference_identifier not in available
                ):
                    raise EosmatError(
                        f"EOS record {record['identifier']!r} references missing "
                        f"pressure EOS {reference_identifier!r}"
                    )
                calibration_identifier = method.get(
                    "reference_calibration_record"
                )
                if (
                    calibration_identifier is not None
                    and calibration_identifier not in available_calibrations
                ):
                    raise EosmatError(
                        f"EOS record {record['identifier']!r} references missing "
                        f"pressure calibration {calibration_identifier!r}"
                    )


def eosmat_schema() -> dict[str, Any]:
    """Return the bundled normative JSON Schema document."""
    resource = resources.files("peritheos.data").joinpath("eosmat-v3.schema.json")
    return json.loads(resource.read_text(encoding="utf-8"))


__all__ = [
    "EOSMAT_FORMAT",
    "EOSMAT_FORMAT_VERSION",
    "eosmat_schema",
    "get_eos_record_document",
    "get_material_document",
    "list_eos_record_documents",
    "list_material_documents",
    "load_eosmat",
    "save_eosmat",
    "validate_eosmat_document",
    "validate_pressure_calibration_references",
]
