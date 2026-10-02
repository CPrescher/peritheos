"""Explicit observation pressure coordinates; never infer a scale from a name."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

from peritheos.datasets import _PRESSURE_QUANTITIES, Dataset, _plain_mapping
from peritheos.errors import DatasetError, EosmatError
from peritheos.pressure_calibrations import (
    RubyFluorescenceCalibration,
    get_pressure_calibration,
    get_pressure_calibration_document,
    list_pressure_calibrations,
    recalculate_ruby_pressure,
)

CONVENTION = "same_corrected_ruby_r1_ratio"
_PRESSURE_UNITS = {"Pa", "kPa", "MPa", "GPa", "bar", "kbar", "Mbar"}


def validate_dataset_pressure_metadata(document: Mapping[str, Any]) -> None:
    """Validate optional pressure provenance and reductions in a material.

    Absence means unresolved. A dataset link alone never authorizes a reduction.
    Called by the Python eosmat validator; JSON Schema handles the same shape.
    """
    calibrations: set[str] | None = None
    records = {
        r["identifier"]: r for r in document.get("eos_records", []) if "identifier" in r
    }

    def provenance(value: Mapping[str, Any], location: str) -> None:
        if (
            not isinstance(value.get("reference"), (str, Mapping))
            or not value.get("reference")
            or (isinstance(value["reference"], str) and not value["reference"].strip())
        ):
            raise EosmatError(f"{location}.reference is required")
        if (
            not isinstance(value.get("source_location"), str)
            or not value["source_location"].strip()
        ):
            raise EosmatError(f"{location}.source_location is required")

    def calibration(identifier: Any) -> None:
        nonlocal calibrations
        if calibrations is None:
            calibrations = set(list_pressure_calibrations())
        if not isinstance(identifier, str) or identifier not in calibrations:
            raise EosmatError(f"Unknown pressure calibration {identifier!r}")

    for dataset in document.get("datasets", []):
        columns = {c["name"]: c for c in dataset["columns"]}
        for column in columns.values():
            if "pressure_scale" not in column:
                continue
            scale = column["pressure_scale"]
            if not isinstance(scale, Mapping) or set(scale) != {
                "calibration_record",
                "reference",
                "source_location",
            }:
                raise EosmatError("pressure_scale requires calibration and provenance")
            if (
                column["quantity"] not in _PRESSURE_QUANTITIES
                or column["role"] != "value"
                or column["unit"] not in _PRESSURE_UNITS
            ):
                raise EosmatError("pressure_scale requires a pressure value column")
            calibration(scale["calibration_record"])
            provenance(scale, "pressure_scale")
        reductions = dataset.get("pressure_reductions", [])
        if not isinstance(reductions, list):
            raise EosmatError("pressure_reductions must be an array")
        seen: set[str] = set()
        for reduction in reductions:
            if not isinstance(reduction, Mapping):
                raise EosmatError("pressure reduction must be an object")
            identifier = reduction.get("eos_record")
            if not isinstance(identifier, str) or identifier not in records:
                raise EosmatError(
                    f"Unknown pressure reduction EOS record {identifier!r}"
                )
            if identifier not in dataset["used_by_eos_records"]:
                raise EosmatError("Pressure reduction EOS must be linked to dataset")
            if identifier in seen:
                raise EosmatError(
                    "Duplicate/ambiguous pressure reduction for EOS record"
                )
            seen.add(identifier)
            status = reduction.get("status")
            if not isinstance(status, str):
                raise EosmatError("Invalid pressure reduction status")
            if "notes" in reduction and (
                not isinstance(reduction["notes"], str)
                or not reduction["notes"].strip()
            ):
                raise EosmatError("Pressure reduction notes must be a non-empty string")
            keys = {"eos_record", "status", "reference", "source_location", "notes"}
            if status in {"as_reported", "transformed"}:
                keys |= {"pressure_column", "target_calibration_record", "row_scope"}
                if reduction.get("row_scope") != "all_rows":
                    raise EosmatError("Pressure reduction row_scope must be all_rows")
            if status == "transformed":
                keys.add("convention")
            if set(reduction) - keys:
                raise EosmatError("Unexpected or conflicting pressure reduction fields")
            provenance(reduction, "pressure_reduction")
            if status == "unresolved":
                if (
                    not isinstance(reduction.get("notes"), str)
                    or not reduction["notes"].strip()
                ):
                    raise EosmatError("Unresolved pressure reduction requires notes")
                continue
            if status not in {"as_reported", "transformed"}:
                raise EosmatError("Invalid pressure reduction status")
            name = reduction.get("pressure_column")
            if not isinstance(name, str) or name not in columns:
                raise EosmatError(f"Unknown pressure column {name!r}")
            scale = columns[name].get("pressure_scale")
            if scale is None:
                raise EosmatError("Pressure reduction source scale is unresolved")
            source = scale["calibration_record"]
            target = reduction.get("target_calibration_record")
            calibration(target)
            eos_calibration = records[identifier].get("pressure_calibration") or {}
            methods = eos_calibration.get("methods", [])
            if (
                eos_calibration.get("status") != "resolved"
                or len(methods) != 1
                or methods[0].get("reference_calibration_record") != target
                or methods[0].get("kind")
                != get_pressure_calibration_document(target)["kind"]
            ):
                raise EosmatError("Reduction conflicts with EOS pressure_calibration")
            if status == "as_reported" and source != target:
                raise EosmatError(
                    "as_reported requires identical source and target scales"
                )
            if status == "transformed":
                if source == target or reduction.get("convention") != CONVENTION:
                    raise EosmatError(
                        "Invalid pressure transformation convention or scales"
                    )
                if not all(
                    isinstance(get_pressure_calibration(c), RubyFluorescenceCalibration)
                    for c in (source, target)
                ):
                    raise EosmatError("Ruby ratio convention requires two ruby scales")


@dataclass(frozen=True)
class DatasetPressure:
    """Resolved pressures with original observations and auditable qualifications.

    Row identity is (dataset.identifier, zero-based row_indices). No filtering,
    sorting, uncertainty propagation, or observation-temperature inference occurs.
    Calibration documents retain validity and reference-temperature qualifications;
    extrapolation masks concern pressure only, not full experimental validity.
    """

    dataset: Dataset
    eos_record: str
    reduction: Mapping[str, Any]
    pressure_gpa: NDArray[np.float64] = field(repr=False)
    reported_pressure: NDArray[np.float64] = field(repr=False)
    reported_pressure_unit: str
    row_indices: NDArray[np.int64] = field(repr=False)
    source_calibration: Mapping[str, Any]
    target_calibration: Mapping[str, Any]
    source_pressure_extrapolated: NDArray[np.bool_] | None = field(repr=False)
    target_pressure_extrapolated: NDArray[np.bool_] | None = field(repr=False)
    eos_reference_temperature_k: float | None
    uncertainty_treatment: str = "not_propagated"

    def __post_init__(self) -> None:
        for name in ("reduction", "source_calibration", "target_calibration"):
            object.__setattr__(self, name, _plain_mapping(getattr(self, name)))


def resolve_dataset_pressure(
    document: Mapping[str, Any],
    eos_record: str,
    dataset_identifier: str,
    *,
    resource_root: Any | None = None,
    check_validity: bool = False,
) -> DatasetPressure:
    """Resolve an explicitly declared dataset pressure coordinate in GPa.

    Raises DatasetError for missing/unresolved reductions, never guesses from
    column names or EOS suffixes. With check_validity, reject pressures outside
    either calibration's declared pressure range. Otherwise retain them and
    return per-row extrapolation masks. Temperature validity is not inferred.
    """
    from peritheos.eosmat import validate_eosmat_document

    validate_eosmat_document(document)
    metadata = next(
        (
            d
            for d in document.get("datasets", [])
            if d["identifier"] == dataset_identifier
        ),
        None,
    )
    if metadata is None:
        raise DatasetError(f"Unknown dataset {dataset_identifier!r}")
    reduction = next(
        (
            r
            for r in metadata.get("pressure_reductions", [])
            if r["eos_record"] == eos_record
        ),
        None,
    )
    if reduction is None or reduction["status"] == "unresolved":
        reason = "No explicit reduction" if reduction is None else reduction["notes"]
        raise DatasetError(f"Unresolved dataset pressure for {eos_record!r}: {reason}")
    dataset = Dataset.from_mapping(metadata, resource_root=resource_root)
    column = dataset.get_column(reduction["pressure_column"])
    source_id = column.metadata["pressure_scale"]["calibration_record"]
    target_id = reduction["target_calibration_record"]
    source = get_pressure_calibration_document(source_id)
    target = get_pressure_calibration_document(target_id)
    reported = dataset.values(column.name)
    try:
        source_pressure = np.asarray(
            dataset.values(column.name, unit="GPa"), dtype=float
        )
    except (TypeError, ValueError) as error:
        raise DatasetError("Pressure reduction requires numeric pressures") from error
    if not np.all(np.isfinite(source_pressure)) or np.any(source_pressure < 0):
        raise DatasetError("Pressure reduction requires finite non-negative pressures")
    pressure = (
        source_pressure
        if reduction["status"] == "as_reported"
        else np.asarray(
            recalculate_ruby_pressure(source_pressure, source_id, target_id),
            dtype=float,
        )
    )

    pressure.setflags(write=False)

    def extrapolated(calibration: Mapping[str, Any], values: NDArray[Any]):
        limits = calibration.get("validity", {}).get("pressure_gpa")
        if limits is None:
            return None
        mask = (values < limits[0]) | (values > limits[1])
        mask.setflags(write=False)
        if check_validity and np.any(mask):
            raise DatasetError(
                f"Pressure outside calibration {calibration['identifier']!r} validity"
            )
        return mask

    indices = np.arange(len(dataset), dtype=np.int64)
    indices.setflags(write=False)
    record = next(
        r for r in document["eos_records"] if r.get("identifier") == eos_record
    )
    return DatasetPressure(
        dataset=dataset,
        eos_record=eos_record,
        reduction=reduction,
        pressure_gpa=pressure,
        reported_pressure=reported,
        reported_pressure_unit=column.unit,
        row_indices=indices,
        source_calibration=source,
        target_calibration=target,
        source_pressure_extrapolated=extrapolated(source, source_pressure),
        target_pressure_extrapolated=extrapolated(target, pressure),
        eos_reference_temperature_k=record.get("temperature_ref"),
    )


__all__ = ["DatasetPressure", "resolve_dataset_pressure"]
