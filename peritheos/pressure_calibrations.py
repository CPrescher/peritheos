"""Executable pressure calibrations bundled with Peritheos.

Ruby pressure-scale conversion is performed through the common measured
quantity, the corrected R1 wavelength ratio ``lambda / lambda0``.  This makes
it possible to convert published pressures between scales even when the paper
reports pressure rather than the original wavelength shift.
"""

from __future__ import annotations

import copy
import json
from collections.abc import Mapping
from dataclasses import dataclass
from importlib import resources
from types import MappingProxyType
from typing import Any, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from peritheos.errors import MaterialLookupError, ValidationError

NumericResult = Union[float, NDArray[np.float64]]
_CALIBRATION_DATA_FILE = "pressure-calibrations.json"


def _as_finite_array(value: ArrayLike, name: str) -> NDArray[np.float64]:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValidationError(f"{name} must contain only finite values", field=name)
    return array


def _result(value: NDArray[np.float64]) -> NumericResult:
    if value.ndim == 0:
        return float(value)
    return value


@dataclass(frozen=True)
class RubyFluorescenceCalibration:
    """An executable room-temperature ruby R1 pressure scale."""

    identifier: str
    label: str
    model: str
    parameters: Mapping[str, float]
    reference: Mapping[str, Any]
    validity: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "parameters", MappingProxyType(dict(self.parameters))
        )
        object.__setattr__(self, "reference", MappingProxyType(dict(self.reference)))
        object.__setattr__(self, "validity", MappingProxyType(dict(self.validity)))

    @property
    def reference_wavelength_nm(self) -> float:
        """Ambient-pressure R1 wavelength stored with the published scale."""
        return float(self.parameters["reference_wavelength_nm"])

    def pressure_from_ratio(self, wavelength_ratio: ArrayLike) -> NumericResult:
        """Return pressure in GPa from corrected ``lambda / lambda0``."""
        ratio = _as_finite_array(wavelength_ratio, "wavelength_ratio")
        if np.any(ratio < 1.0):
            raise ValidationError(
                "wavelength_ratio must be at least one for non-negative pressure",
                field="wavelength_ratio",
            )

        a = float(self.parameters["A_gpa"])
        if self.model == "power_law":
            b = float(self.parameters["B"])
            pressure = (a / b) * (np.power(ratio, b) - 1.0)
        elif self.model == "quadratic_shift":
            shift_ratio = ratio - 1.0
            m = float(self.parameters["m"])
            pressure = a * shift_ratio * (1.0 + m * shift_ratio)
        elif self.model == "holzapfel_freund_ingalls":
            b = float(self.parameters["B"])
            c = float(self.parameters["C"])
            exponent = ((b + c) / c) * (1.0 - np.power(ratio, -c))
            pressure = (a / (b + c)) * np.expm1(exponent)
        else:  # pragma: no cover - protected by validation of bundled data
            raise ValidationError(f"Unsupported ruby calibration model {self.model!r}")
        return _result(np.asarray(pressure, dtype=float))

    def wavelength_ratio(self, pressure_gpa: ArrayLike) -> NumericResult:
        """Invert pressure in GPa to the corrected R1 wavelength ratio."""
        pressure = _as_finite_array(pressure_gpa, "pressure_gpa")
        if np.any(pressure < 0.0):
            raise ValidationError(
                "pressure_gpa must be non-negative", field="pressure_gpa"
            )

        a = float(self.parameters["A_gpa"])
        if self.model == "power_law":
            b = float(self.parameters["B"])
            ratio = np.power(1.0 + b * pressure / a, 1.0 / b)
        elif self.model == "quadratic_shift":
            m = float(self.parameters["m"])
            shift_ratio = np.expm1(0.5 * np.log1p(4.0 * m * pressure / a)) / (
                2.0 * m
            )
            ratio = 1.0 + shift_ratio
        elif self.model == "holzapfel_freund_ingalls":
            b = float(self.parameters["B"])
            c = float(self.parameters["C"])
            base = 1.0 - (c / (b + c)) * np.log1p((b + c) * pressure / a)
            if np.any(base <= 0.0):
                raise ValidationError(
                    "pressure_gpa is outside the invertible Holzapfel scale domain",
                    field="pressure_gpa",
                )
            ratio = np.power(base, -1.0 / c)
        else:  # pragma: no cover - protected by validation of bundled data
            raise ValidationError(f"Unsupported ruby calibration model {self.model!r}")
        return _result(np.asarray(ratio, dtype=float))

    def pressure_from_wavelength(
        self,
        wavelength_nm: ArrayLike,
        *,
        reference_wavelength_nm: float | None = None,
    ) -> NumericResult:
        """Return pressure from a temperature-corrected R1 wavelength."""
        wavelength = _as_finite_array(wavelength_nm, "wavelength_nm")
        reference = (
            self.reference_wavelength_nm
            if reference_wavelength_nm is None
            else float(reference_wavelength_nm)
        )
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValidationError(
                "reference_wavelength_nm must be finite and positive",
                field="reference_wavelength_nm",
            )
        return self.pressure_from_ratio(wavelength / reference)

    def pressure_from_shift(
        self,
        wavelength_shift_nm: ArrayLike,
        *,
        reference_wavelength_nm: float | None = None,
    ) -> NumericResult:
        """Return pressure from a temperature-corrected R1 wavelength shift."""
        shift = _as_finite_array(wavelength_shift_nm, "wavelength_shift_nm")
        reference = (
            self.reference_wavelength_nm
            if reference_wavelength_nm is None
            else float(reference_wavelength_nm)
        )
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValidationError(
                "reference_wavelength_nm must be finite and positive",
                field="reference_wavelength_nm",
            )
        return self.pressure_from_ratio(1.0 + shift / reference)

    def wavelength_from_pressure(
        self,
        pressure_gpa: ArrayLike,
        *,
        reference_wavelength_nm: float | None = None,
    ) -> NumericResult:
        """Return the corrected R1 wavelength implied by pressure."""
        reference = (
            self.reference_wavelength_nm
            if reference_wavelength_nm is None
            else float(reference_wavelength_nm)
        )
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValidationError(
                "reference_wavelength_nm must be finite and positive",
                field="reference_wavelength_nm",
            )
        return _result(np.asarray(self.wavelength_ratio(pressure_gpa)) * reference)


@dataclass(frozen=True)
class XrdPressureRecalculation:
    """Paired source and recalculated pressures for an XRD standard."""

    source_pressure_gpa: NumericResult
    source_wavelength_ratio: NumericResult
    target_pressure_gpa: NumericResult
    pressure_difference_gpa: NumericResult
    source_calibration_record: str
    target_eos_record: str


def pressure_calibration_library() -> dict[str, Any]:
    """Return a copy of the bundled calibration-library document."""
    resource = resources.files("peritheos.data").joinpath(_CALIBRATION_DATA_FILE)
    return json.loads(resource.read_text(encoding="utf-8"))


def list_pressure_calibrations() -> tuple[str, ...]:
    """Return identifiers of all bundled executable pressure calibrations."""
    return tuple(
        entry["identifier"]
        for entry in pressure_calibration_library()["calibrations"]
    )


def list_xrd_pressure_standards() -> tuple[str, ...]:
    """Return bundled EOS records explicitly used as XRD pressure standards."""
    # Local imports avoid a cycle while eosmat validates calibration links.
    from peritheos.eosmat import get_material_document, list_material_documents

    return tuple(
        sorted(
            {
                method["reference_eos_record"]
                for material_identifier in list_material_documents()
                for record in get_material_document(material_identifier)["eos_records"]
                for method in record["pressure_calibration"]["methods"]
                if "reference_eos_record" in method
            }
        )
    )


def get_pressure_calibration_document(identifier: str) -> dict[str, Any]:
    """Return one pressure-calibration data record."""
    matches = [
        entry
        for entry in pressure_calibration_library()["calibrations"]
        if entry["identifier"] == identifier
    ]
    if len(matches) != 1:
        if not matches:
            raise MaterialLookupError(
                f"Unknown pressure calibration {identifier!r}; available: "
                f"{list_pressure_calibrations()}"
            )
        raise ValidationError(f"Duplicate pressure calibration {identifier!r}")
    return copy.deepcopy(matches[0])


def get_pressure_calibration(identifier: str) -> RubyFluorescenceCalibration:
    """Construct an executable bundled pressure calibration."""
    document = get_pressure_calibration_document(identifier)
    if document["kind"] != "ruby_fluorescence":
        raise ValidationError(
            f"Pressure calibration {identifier!r} is not a ruby scale"
        )
    return RubyFluorescenceCalibration(
        identifier=document["identifier"],
        label=document["label"],
        model=document["model"],
        parameters=document["parameters"],
        reference=document["reference"],
        validity=document["validity"],
    )


def recalculate_ruby_pressure(
    pressure_gpa: ArrayLike,
    source_calibration: str | RubyFluorescenceCalibration,
    target_calibration: str | RubyFluorescenceCalibration,
) -> NumericResult:
    """Convert ruby-derived pressure between scales via the R1 ratio.

    This conversion assumes that both scales are applied to the same corrected
    R1 measurement.  It cannot undo temperature, stress, composition, peak-fit,
    or spatial-gradient effects in the original experiment.
    """
    source = (
        get_pressure_calibration(source_calibration)
        if isinstance(source_calibration, str)
        else source_calibration
    )
    target = (
        get_pressure_calibration(target_calibration)
        if isinstance(target_calibration, str)
        else target_calibration
    )
    return target.pressure_from_ratio(source.wavelength_ratio(pressure_gpa))


def xrd_standard_pressure(
    reference_eos_record: str,
    volume: ArrayLike,
    temperature_k: ArrayLike | None = None,
    *,
    check_validity: bool = False,
) -> NumericResult:
    """Calculate pressure from a bundled XRD standard EOS.

    ``volume`` uses the conventional-unit-cell volume documented by the EOS
    record.  ``temperature_k`` defaults to that record's reference
    temperature.  Thermal standards require temperature for non-reference
    conditions; isothermal standards reject other temperatures.
    """
    # Local imports avoid a cycle while eosmat validates calibration links.
    from peritheos.eosmat import get_material_document, list_material_documents
    from peritheos.materials import Material

    matches = []
    for material_identifier in list_material_documents():
        document = get_material_document(material_identifier)
        if any(
            record["identifier"] == reference_eos_record
            for record in document["eos_records"]
        ):
            matches.append(document)
    if len(matches) != 1:
        raise MaterialLookupError(
            f"Expected exactly one bundled EOS record {reference_eos_record!r}; "
            f"found {len(matches)}"
        )
    material = Material.from_eosmat(
        matches[0], record_identifiers=[reference_eos_record]
    )
    standard = material.get_eos_record(reference_eos_record)
    result = standard.pressure(
        volume,
        temperature_k,
        check_validity=check_validity,
    )
    return _result(np.asarray(result, dtype=float))


def recalculate_ruby_to_xrd_pressure(
    source_pressure_gpa: ArrayLike,
    source_calibration_record: str,
    target_eos_record: str,
    target_volume: ArrayLike,
    target_temperature_k: ArrayLike | None = None,
    *,
    check_validity: bool = False,
) -> XrdPressureRecalculation:
    """Replace ruby-derived pressures with a paired XRD-standard pressure.

    Each target volume (and temperature for a thermal standard) must describe
    the XRD calibrant measured at the same experimental state as the source
    pressure.  Pressure alone contains no information from which an Au, Pt,
    MgO, or other diffraction-standard volume can be reconstructed.
    """
    source = get_pressure_calibration(source_calibration_record)
    source_pressure = _as_finite_array(source_pressure_gpa, "source_pressure_gpa")
    source_ratio = np.asarray(source.wavelength_ratio(source_pressure), dtype=float)
    target_pressure = np.asarray(
        xrd_standard_pressure(
            target_eos_record,
            target_volume,
            target_temperature_k,
            check_validity=check_validity,
        ),
        dtype=float,
    )
    try:
        source_pressure, source_ratio, target_pressure = np.broadcast_arrays(
            source_pressure, source_ratio, target_pressure
        )
    except ValueError as error:
        raise ValidationError(
            "source pressures and target XRD observations must be broadcast-compatible"
        ) from error
    difference = target_pressure - source_pressure
    return XrdPressureRecalculation(
        source_pressure_gpa=_result(source_pressure),
        source_wavelength_ratio=_result(source_ratio),
        target_pressure_gpa=_result(target_pressure),
        pressure_difference_gpa=_result(difference),
        source_calibration_record=source.identifier,
        target_eos_record=target_eos_record,
    )


__all__ = [
    "RubyFluorescenceCalibration",
    "XrdPressureRecalculation",
    "get_pressure_calibration",
    "get_pressure_calibration_document",
    "list_pressure_calibrations",
    "list_xrd_pressure_standards",
    "pressure_calibration_library",
    "recalculate_ruby_pressure",
    "recalculate_ruby_to_xrd_pressure",
    "xrd_standard_pressure",
]
