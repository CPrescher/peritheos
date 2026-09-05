"""Typed discovery over the canonical executable material catalog."""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping
from numbers import Real
from typing import Any, Literal, Union

from peritheos.eos import ThermalEOS
from peritheos.errors import MaterialError
from peritheos.materials import EOSRecord, Material

RangeQuery = Union[float, tuple[float, float]]
RangeSemantics = Literal["contains", "overlaps"]
ValidationStatus = Literal[
    "primary_source_validated",
    "pending_primary_source_check",
    "deferred",
]


def _normalized(value: str) -> str:
    return " ".join(re.split(r"[^\w]+", value.casefold())).strip()


def _flatten_text(value: Any) -> tuple[str, ...]:
    if isinstance(value, Mapping):
        return tuple(text for item in value.values() for text in _flatten_text(item))
    if isinstance(value, (tuple, list)):
        return tuple(text for item in value for text in _flatten_text(item))
    if value is None:
        return ()
    return (str(value),)


def _matches_text(query: str | None, values: Iterable[str]) -> bool:
    if query is None:
        return True
    terms = _normalized(query).split()
    if not terms:
        return True
    haystack = _normalized(" ".join(values))
    return all(term in haystack for term in terms)


def _record_reference_text(record: EOSRecord) -> tuple[str, ...]:
    return (
        record.reference.authors,
        str(record.reference.year),
        record.reference.title,
        record.reference.doi,
        *record.reference.locations,
        *_flatten_text(record.eosmat_metadata.get("reference")),
        *_flatten_text(record.eosmat_metadata.get("scientific_validation")),
    )


def _record_models(record: EOSRecord) -> tuple[str, ...]:
    values = [type(record.eos).__name__]
    if record.is_thermal:
        values.append(type(record.eos.rt_eos).__name__)
    for component_name in ("eos", "thermal"):
        component = record.eosmat_metadata.get(component_name)
        if isinstance(component, Mapping):
            values.extend(
                str(component[key])
                for key in ("model", "type")
                if component.get(key) is not None
            )
    return tuple(values)


def _has_caloric_model(record: EOSRecord) -> bool:
    if not isinstance(record.eos, ThermalEOS):
        return False
    return (
        type(record.eos).molar_heat_capacity_v is not ThermalEOS.molar_heat_capacity_v
    )


def _range_query(value: RangeQuery | None, field: str) -> tuple[float, float] | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise MaterialError(f"{field} must be a finite number or two-value range")
    if isinstance(value, tuple):
        try:
            if len(value) != 2:
                raise MaterialError(f"{field} range must contain exactly two values")
            bounds = (float(value[0]), float(value[1]))
        except (TypeError, ValueError) as error:
            raise MaterialError(
                f"{field} must be a finite number or two-value range"
            ) from error
    elif isinstance(value, Real):
        bounds = (float(value), float(value))
    else:
        raise MaterialError(f"{field} must be a finite number or two-value range")
    if not all(math.isfinite(item) for item in bounds) or bounds[0] > bounds[1]:
        raise MaterialError(f"{field} must be finite and ordered low-to-high")
    if field == "temperature_k" and bounds[0] <= 0.0:
        raise MaterialError("temperature_k must be greater than zero")
    return bounds


def _calibration_range(
    record: EOSRecord, field: Literal["pressure_gpa", "temperature_k"]
) -> tuple[float, float] | None:
    """Return a reported range, without turning missing metadata into coverage."""
    metadata = record.eosmat_metadata
    validity = metadata.get("validity")
    value: Any = validity.get(field) if isinstance(validity, Mapping) else None
    if value is None:
        legacy_name = {
            "pressure_gpa": "experimental_pressure_range_gpa",
            "temperature_k": "experimental_temperature_range_k",
        }[field]
        value = metadata.get(legacy_name)
    if value is not None:
        return (float(value[0]), float(value[1]))
    if field == "temperature_k" and not record.is_thermal:
        # An isothermal fit is calibrated on its declared reference isotherm,
        # even when a separate one-point temperature interval was not repeated.
        return (record.reference_temperature, record.reference_temperature)
    return None


def _range_matches(
    available: tuple[float, float] | None,
    requested: tuple[float, float] | None,
    semantics: RangeSemantics,
) -> bool:
    if requested is None:
        return True
    if available is None:
        return False
    if semantics == "contains":
        return available[0] <= requested[0] and available[1] >= requested[1]
    return available[0] <= requested[1] and available[1] >= requested[0]


def _validation_statuses(
    value: ValidationStatus | Iterable[ValidationStatus] | None,
) -> frozenset[str] | None:
    if value is None:
        return None
    statuses: frozenset[str] = (
        frozenset((value,)) if isinstance(value, str) else frozenset(value)
    )
    allowed = {
        "primary_source_validated",
        "pending_primary_source_check",
        "deferred",
    }
    invalid = statuses - allowed
    if invalid:
        raise MaterialError(f"Unknown validation status: {sorted(invalid)}")
    return statuses


def _normalize_doi(value: str) -> str:
    normalized = value.strip().casefold()
    for prefix in (
        "https://doi.org/",
        "http://doi.org/",
        "https://dx.doi.org/",
        "http://dx.doi.org/",
        "doi:",
    ):
        if normalized.startswith(prefix):
            return normalized.removeprefix(prefix).strip()
    return normalized


def _validate_optional_bool(value: bool | None, field: str) -> None:
    if value is not None and not isinstance(value, bool):
        raise MaterialError(f"{field} must be true, false, or None")


def _record_matches(
    record: EOSRecord,
    *,
    text: str | None,
    name: str | None,
    alias: str | None,
    formula: str | None,
    phase: str | None,
    model_family: str | None,
    doi: str | None,
    author: str | None,
    reference: str | None,
    thermal: bool | None,
    caloric: bool | None,
    pressure: tuple[float, float] | None,
    temperature: tuple[float, float] | None,
    range_semantics: RangeSemantics,
    uncertainty: bool | None,
    statuses: frozenset[str] | None,
    material: Material | None,
) -> bool:
    reference_text = _record_reference_text(record)
    material_text: tuple[str, ...] = ()
    if material is not None:
        material_text = (
            material.identifier,
            material.name or "",
            *material.aliases,
        )
    text_values = (
        record.identifier,
        record.name,
        record.material,
        record.phase,
        *record.aliases,
        *material_text,
        *reference_text,
    )
    if not _matches_text(text, text_values):
        return False
    if not _matches_text(name, (record.name, *(material_text[1:2]))):
        return False
    if not _matches_text(
        alias, (*record.aliases, *(material.aliases if material else ()))
    ):
        return False
    if formula is not None and record.material.casefold() != formula.strip().casefold():
        return False
    phase_values: tuple[str, ...] = (record.phase,)
    if material is not None:
        phase_values = (
            *phase_values,
            material.name or "",
            material.symmetry or "",
            *material.aliases,
        )
    if not _matches_text(phase, phase_values):
        return False
    if not _matches_text(model_family, _record_models(record)):
        return False
    if doi is not None and _normalize_doi(doi) != _normalize_doi(record.reference.doi):
        return False
    if not _matches_text(author, (record.reference.authors,)):
        return False
    if not _matches_text(reference, reference_text):
        return False
    if thermal is not None and record.is_thermal is not thermal:
        return False
    if caloric is not None and _has_caloric_model(record) is not caloric:
        return False
    if uncertainty is not None:
        available = record.parameter_covariance is not None or bool(
            record.parameter_errors
        )
        if available is not uncertainty:
            return False
    if statuses is not None and record.scientific_validation_status not in statuses:
        return False
    return _range_matches(
        _calibration_range(record, "pressure_gpa"), pressure, range_semantics
    ) and _range_matches(
        _calibration_range(record, "temperature_k"), temperature, range_semantics
    )


def _prepare_filters(
    *,
    pressure_gpa: RangeQuery | None,
    temperature_k: RangeQuery | None,
    range_semantics: RangeSemantics,
    thermal: bool | None,
    caloric: bool | None,
    uncertainty: bool | None,
    validation_status: ValidationStatus | Iterable[ValidationStatus] | None,
) -> tuple[
    tuple[float, float] | None,
    tuple[float, float] | None,
    frozenset[str] | None,
]:
    if range_semantics not in ("contains", "overlaps"):
        raise MaterialError("range_semantics must be 'contains' or 'overlaps'")
    for field, value in (
        ("thermal", thermal),
        ("caloric", caloric),
        ("uncertainty", uncertainty),
    ):
        _validate_optional_bool(value, field)
    return (
        _range_query(pressure_gpa, "pressure_gpa"),
        _range_query(temperature_k, "temperature_k"),
        _validation_statuses(validation_status),
    )


def search_eos_records(
    text: str | None = None,
    *,
    name: str | None = None,
    alias: str | None = None,
    formula: str | None = None,
    phase: str | None = None,
    model_family: str | None = None,
    doi: str | None = None,
    author: str | None = None,
    reference: str | None = None,
    thermal: bool | None = None,
    caloric: bool | None = None,
    pressure_gpa: RangeQuery | None = None,
    temperature_k: RangeQuery | None = None,
    range_semantics: RangeSemantics = "contains",
    uncertainty: bool | None = None,
    validation_status: ValidationStatus | Iterable[ValidationStatus] | None = None,
) -> tuple[EOSRecord, ...]:
    """Search the complete catalog and return executable EOS records.

    Scalar range filters select records whose closed, reported calibration
    interval contains that point. A two-value query uses ``range_semantics``:
    ``"contains"`` requires complete coverage and ``"overlaps"`` requires
    any closed-interval overlap. A missing published range never matches a
    range query.
    """
    from peritheos.catalog import _catalog_index, list_eos_records

    pressure, temperature, statuses = _prepare_filters(
        pressure_gpa=pressure_gpa,
        temperature_k=temperature_k,
        range_semantics=range_semantics,
        thermal=thermal,
        caloric=caloric,
        uncertainty=uncertainty,
        validation_status=validation_status,
    )
    material_by_record = _catalog_index().material_by_record
    return tuple(
        record
        for record in list_eos_records()
        if _record_matches(
            record,
            text=text,
            name=name,
            alias=alias,
            formula=formula,
            phase=phase,
            model_family=model_family,
            doi=doi,
            author=author,
            reference=reference,
            thermal=thermal,
            caloric=caloric,
            pressure=pressure,
            temperature=temperature,
            range_semantics=range_semantics,
            uncertainty=uncertainty,
            statuses=statuses,
            material=material_by_record.get(record.identifier),
        )
    )


def search_materials(
    text: str | None = None,
    *,
    name: str | None = None,
    alias: str | None = None,
    formula: str | None = None,
    phase: str | None = None,
    model_family: str | None = None,
    doi: str | None = None,
    author: str | None = None,
    reference: str | None = None,
    thermal: bool | None = None,
    caloric: bool | None = None,
    pressure_gpa: RangeQuery | None = None,
    temperature_k: RangeQuery | None = None,
    range_semantics: RangeSemantics = "contains",
    uncertainty: bool | None = None,
    validation_status: ValidationStatus | Iterable[ValidationStatus] | None = None,
) -> tuple[Material, ...]:
    """Search the complete catalog and return executable materials.

    Material identity filters apply to the material itself. All criteria that
    depend on records, including record-only free-text matches, must be
    satisfied by one record belonging to the returned material.
    """
    from peritheos.catalog import list_materials

    pressure, temperature, statuses = _prepare_filters(
        pressure_gpa=pressure_gpa,
        temperature_k=temperature_k,
        range_semantics=range_semantics,
        thermal=thermal,
        caloric=caloric,
        uncertainty=uncertainty,
        validation_status=validation_status,
    )
    results = []
    for material in list_materials():
        material_text = (
            material.identifier,
            material.name or "",
            material.formula,
            material.phase,
            *material.aliases,
        )
        material_text_matches = _matches_text(text, material_text)
        if not material_text_matches and not any(
            _matches_text(
                text,
                (
                    record.identifier,
                    record.name,
                    *record.aliases,
                    *_record_reference_text(record),
                ),
            )
            for record in material.eos_records
        ):
            continue
        if not _matches_text(name, (material.name or "",)):
            continue
        if not _matches_text(alias, material.aliases):
            continue
        if (
            formula is not None
            and material.formula.casefold() != formula.strip().casefold()
        ):
            continue
        if not _matches_text(
            phase,
            (
                material.phase,
                material.name or "",
                material.symmetry or "",
                *material.aliases,
            ),
        ):
            continue
        if any(
            _record_matches(
                record,
                text=None if material_text_matches else text,
                name=None,
                alias=None,
                formula=None,
                phase=None,
                model_family=model_family,
                doi=doi,
                author=author,
                reference=reference,
                thermal=thermal,
                caloric=caloric,
                pressure=pressure,
                temperature=temperature,
                range_semantics=range_semantics,
                uncertainty=uncertainty,
                statuses=statuses,
                material=material,
            )
            for record in material.eos_records
        ):
            results.append(material)
    return tuple(results)


__all__ = [
    "RangeQuery",
    "RangeSemantics",
    "ValidationStatus",
    "search_eos_records",
    "search_materials",
]
