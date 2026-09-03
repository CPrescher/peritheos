"""Executable pressure calibrations and scale transformations.

Ruby pressure-scale conversion is performed through the common measured
quantity, the corrected R1 wavelength ratio ``lambda / lambda0``.  This makes
it possible to convert published pressures between scales even when the paper
reports pressure rather than the original wavelength shift.

XRD scale conversion uses the unit-cell volume of a common standard as an
internal coordinate: invert the source-standard EOS, then evaluate the target
EOS at the same virtual volume.  Material EOS records can follow either edge
from their audited pressure-calibration provenance.
"""

from __future__ import annotations

import copy
import json
from collections import deque
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
        object.__setattr__(self, "parameters", MappingProxyType(dict(self.parameters)))
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
            shift_ratio = np.expm1(0.5 * np.log1p(4.0 * m * pressure / a)) / (2.0 * m)
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
class DiamondRamanCalibration:
    """An executable diamond-anvil high-frequency Raman-edge scale."""

    identifier: str
    label: str
    model: str
    parameters: Mapping[str, float]
    reference: Mapping[str, Any]
    validity: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameters", MappingProxyType(dict(self.parameters)))
        object.__setattr__(self, "reference", MappingProxyType(dict(self.reference)))
        object.__setattr__(self, "validity", MappingProxyType(dict(self.validity)))

    @property
    def reference_wavenumber_cm1(self) -> float:
        """Ambient high-frequency-edge wavenumber stored with the scale."""
        return float(self.parameters["reference_wavenumber_cm1"])

    def pressure_from_ratio(self, wavenumber_ratio: ArrayLike) -> NumericResult:
        """Return pressure in GPa from ``omega / omega0``."""
        ratio = _as_finite_array(wavenumber_ratio, "wavenumber_ratio")
        if np.any(ratio < 1.0):
            raise ValidationError(
                "wavenumber_ratio must be at least one for non-negative pressure",
                field="wavenumber_ratio",
            )
        shift = ratio - 1.0
        if self.model == "normalized_quadratic":
            a = float(self.parameters["A_gpa"])
            b = float(self.parameters["B_gpa"])
            pressure = a * shift + b * shift**2
        elif self.model == "akahama_quadratic":
            k0 = float(self.parameters["K0_gpa"])
            k0_prime = float(self.parameters["K0_prime"])
            pressure = k0 * shift * (1.0 + 0.5 * (k0_prime - 1.0) * shift)
        else:  # pragma: no cover - protected by bundled-data tests
            raise ValidationError(
                f"Unsupported diamond Raman calibration model {self.model!r}"
            )
        return _result(np.asarray(pressure, dtype=float))

    def wavenumber_ratio(self, pressure_gpa: ArrayLike) -> NumericResult:
        """Invert pressure in GPa to ``omega / omega0``."""
        pressure = _as_finite_array(pressure_gpa, "pressure_gpa")
        if np.any(pressure < 0.0):
            raise ValidationError(
                "pressure_gpa must be non-negative", field="pressure_gpa"
            )
        if self.model == "normalized_quadratic":
            linear = float(self.parameters["A_gpa"])
            quadratic = float(self.parameters["B_gpa"])
        elif self.model == "akahama_quadratic":
            linear = float(self.parameters["K0_gpa"])
            quadratic = 0.5 * linear * (
                float(self.parameters["K0_prime"]) - 1.0
            )
        else:  # pragma: no cover - protected by bundled-data tests
            raise ValidationError(
                f"Unsupported diamond Raman calibration model {self.model!r}"
            )
        shift = (
            np.sqrt(linear**2 + 4.0 * quadratic * pressure) - linear
        ) / (2.0 * quadratic)
        return _result(np.asarray(1.0 + shift, dtype=float))

    def pressure_from_wavenumber(
        self,
        wavenumber_cm1: ArrayLike,
        *,
        reference_wavenumber_cm1: float | None = None,
    ) -> NumericResult:
        """Return pressure from the measured diamond-edge wavenumber."""
        wavenumber = _as_finite_array(wavenumber_cm1, "wavenumber_cm1")
        reference = (
            self.reference_wavenumber_cm1
            if reference_wavenumber_cm1 is None
            else float(reference_wavenumber_cm1)
        )
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValidationError(
                "reference_wavenumber_cm1 must be finite and positive",
                field="reference_wavenumber_cm1",
            )
        return self.pressure_from_ratio(wavenumber / reference)

    def wavenumber_from_pressure(
        self,
        pressure_gpa: ArrayLike,
        *,
        reference_wavenumber_cm1: float | None = None,
    ) -> NumericResult:
        """Return the diamond-edge wavenumber implied by pressure."""
        reference = (
            self.reference_wavenumber_cm1
            if reference_wavenumber_cm1 is None
            else float(reference_wavenumber_cm1)
        )
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValidationError(
                "reference_wavenumber_cm1 must be finite and positive",
                field="reference_wavenumber_cm1",
            )
        return _result(np.asarray(self.wavenumber_ratio(pressure_gpa)) * reference)


@dataclass(frozen=True)
class XrdPressureRecalculation:
    """Pressures re-reduced from paired, measured XRD-standard volumes."""

    source_pressure_gpa: NumericResult
    source_wavelength_ratio: NumericResult
    target_pressure_gpa: NumericResult
    pressure_difference_gpa: NumericResult
    source_calibration_record: str
    target_eos_record: str


@dataclass(frozen=True)
class PressureScaleRecalculation:
    """A pressure-scale transformation through a virtual standard state.

    ``implied_standard_volume`` is calculated by inverting the source-standard
    EOS; it is not a measured calibrant volume from the sample experiment.
    """

    source_pressure_gpa: NumericResult
    implied_standard_volume: NumericResult
    target_pressure_gpa: NumericResult
    pressure_difference_gpa: NumericResult
    source_standard_eos_record: str
    target_standard_eos_record: str


@dataclass(frozen=True)
class RubyPressureScaleRecalculation:
    """A ruby-to-XRD transformation through a ruby-linked standard EOS."""

    source_pressure_gpa: NumericResult
    source_wavelength_ratio: NumericResult
    bridge_pressure_gpa: NumericResult
    implied_standard_volume: NumericResult
    target_pressure_gpa: NumericResult
    pressure_difference_gpa: NumericResult
    source_calibration_record: str
    bridge_calibration_record: str
    bridge_eos_record: str
    target_standard_eos_record: str


@dataclass(frozen=True)
class CalibrationPathRecalculation:
    """A recursive pressure transformation through the calibration graph."""

    source_pressure_gpa: NumericResult
    target_pressure_gpa: NumericResult
    pressure_difference_gpa: NumericResult
    source_node: str
    target_node: str
    calibration_path: tuple[str, ...]
    edge_identifiers: tuple[str, ...]
    intermediate_states: tuple[Mapping[str, Any], ...]


def pressure_calibration_library() -> dict[str, Any]:
    """Return a copy of the bundled calibration-library document."""
    resource = resources.files("peritheos.data").joinpath(_CALIBRATION_DATA_FILE)
    return json.loads(resource.read_text(encoding="utf-8"))


def list_pressure_calibrations() -> tuple[str, ...]:
    """Return identifiers of all bundled executable pressure calibrations."""
    return tuple(
        entry["identifier"] for entry in pressure_calibration_library()["calibrations"]
    )


def list_ruby_pressure_calibrations() -> tuple[str, ...]:
    """Return identifiers of bundled ruby R1 fluorescence scales."""
    return tuple(
        entry["identifier"]
        for entry in pressure_calibration_library()["calibrations"]
        if entry["kind"] == "ruby_fluorescence"
    )


def list_diamond_raman_calibrations() -> tuple[str, ...]:
    """Return identifiers of bundled diamond-anvil Raman-edge scales."""
    return tuple(
        entry["identifier"]
        for entry in pressure_calibration_library()["calibrations"]
        if entry["kind"] == "diamond_raman"
    )


def list_xrd_pressure_standards() -> tuple[str, ...]:
    """Return bundled EOS records explicitly used as XRD pressure standards."""
    # Local imports avoid a cycle while eosmat validates calibration links.
    from peritheos.eosmat import get_material_document, list_material_documents

    standards = set(
        pressure_calibration_library().get("pressure_standard_records", [])
    )
    standards.update(
        method["reference_eos_record"]
        for material_identifier in list_material_documents()
        for record in get_material_document(material_identifier)["eos_records"]
        for method in record["pressure_calibration"]["methods"]
        if "reference_eos_record" in method
    )
    return tuple(sorted(standards))


def _resolve_bundled_eos_record(identifier: str) -> tuple[str, dict[str, Any], Any]:
    """Return material id, document, and executable record for one EOS id."""
    # Local imports avoid a cycle while eosmat validates calibration links.
    from peritheos.eosmat import get_material_document, list_material_documents
    from peritheos.materials import Material

    matches: list[tuple[str, dict[str, Any]]] = []
    for material_identifier in list_material_documents():
        document = get_material_document(material_identifier)
        if any(
            record["identifier"] == identifier for record in document["eos_records"]
        ):
            matches.append((material_identifier, document))
    if len(matches) != 1:
        raise MaterialLookupError(
            f"Expected exactly one bundled EOS record {identifier!r}; "
            f"found {len(matches)}"
        )
    material_identifier, document = matches[0]
    material = Material.from_eosmat(document, record_identifiers=[identifier])
    return material_identifier, document, material.get_eos_record(identifier)


def _ruby_calibration_for_eos(identifier: str) -> str:
    """Return the single executable ruby calibration linked by an EOS record."""
    from peritheos.eosmat import get_eos_record_document

    record = get_eos_record_document(identifier)
    calibrations = {
        method["reference_calibration_record"]
        for method in record["pressure_calibration"]["methods"]
        if method["kind"] == "ruby_fluorescence"
        and "reference_calibration_record" in method
    }
    if len(calibrations) != 1:
        raise ValidationError(
            f"EOS record {identifier!r} must link exactly one executable ruby "
            f"calibration; found {sorted(calibrations)}",
            field="bridge_eos_record",
        )
    return next(iter(calibrations))


def list_ruby_xrd_bridges(target_eos_record: str | None = None) -> tuple[str, ...]:
    """List ruby-linked EOS records usable as virtual XRD bridges.

    If ``target_eos_record`` is supplied, only records for the same material
    and conventional unit-cell definition are returned.
    """
    from peritheos.eosmat import get_material_document, list_material_documents

    target_material = None
    if target_eos_record is not None:
        target_material, _, _ = _resolve_bundled_eos_record(target_eos_record)

    bridges = []
    for material_identifier in list_material_documents():
        if target_material is not None and material_identifier != target_material:
            continue
        document = get_material_document(material_identifier)
        for record in document["eos_records"]:
            if any(
                method["kind"] == "ruby_fluorescence"
                and "reference_calibration_record" in method
                for method in record["pressure_calibration"]["methods"]
            ):
                bridges.append(record["identifier"])
    return tuple(sorted(bridges))


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


def get_pressure_calibration(
    identifier: str,
) -> RubyFluorescenceCalibration | DiamondRamanCalibration:
    """Construct an executable bundled pressure calibration."""
    document = get_pressure_calibration_document(identifier)
    if document["kind"] == "ruby_fluorescence":
        return RubyFluorescenceCalibration(
            identifier=document["identifier"],
            label=document["label"],
            model=document["model"],
            parameters=document["parameters"],
            reference=document["reference"],
            validity=document["validity"],
        )
    if document["kind"] == "diamond_raman":
        return DiamondRamanCalibration(
            identifier=document["identifier"],
            label=document["label"],
            model=document["model"],
            parameters=document["parameters"],
            reference=document["reference"],
            validity=document["validity"],
        )
    raise ValidationError(
        f"Pressure calibration {identifier!r} has unsupported kind "
        f"{document['kind']!r}"
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
    if not isinstance(source, RubyFluorescenceCalibration) or not isinstance(
        target, RubyFluorescenceCalibration
    ):
        raise ValidationError("source and target must both be ruby calibrations")
    return target.pressure_from_ratio(source.wavelength_ratio(pressure_gpa))


def recalculate_diamond_raman_pressure(
    pressure_gpa: ArrayLike,
    source_calibration: str | DiamondRamanCalibration,
    target_calibration: str | DiamondRamanCalibration,
) -> NumericResult:
    """Convert diamond-edge pressure scales through a common Raman ratio."""
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
    if not isinstance(source, DiamondRamanCalibration) or not isinstance(
        target, DiamondRamanCalibration
    ):
        raise ValidationError(
            "source and target must both be diamond Raman calibrations"
        )
    return target.pressure_from_ratio(source.wavenumber_ratio(pressure_gpa))


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
    _, _, standard = _resolve_bundled_eos_record(reference_eos_record)
    result = standard.pressure(
        volume,
        temperature_k,
        check_validity=check_validity,
    )
    return _result(np.asarray(result, dtype=float))


def recalculate_xrd_pressure_scale(
    source_pressure_gpa: ArrayLike,
    source_standard_eos_record: str,
    target_standard_eos_record: str,
    temperature_k: ArrayLike | None = None,
    *,
    check_validity: bool = False,
) -> PressureScaleRecalculation:
    """Transform pressure between two EOSs of the same physical standard.

    The source-standard EOS is inverted to a virtual standard volume and the
    target-standard EOS is evaluated at that volume.  The two records may be
    members of different pressure scales, but must use the same material and
    conventional unit cell.
    """
    source_material, source_document, source_standard = _resolve_bundled_eos_record(
        source_standard_eos_record
    )
    target_material, target_document, target_standard = _resolve_bundled_eos_record(
        target_standard_eos_record
    )
    if source_material != target_material or (
        source_document.get("formula_units_per_cell")
        != target_document.get("formula_units_per_cell")
    ):
        raise ValidationError(
            "source and target EOS records must describe the same XRD standard "
            "and conventional unit cell",
            field="target_standard_eos_record",
        )

    source_pressure = _as_finite_array(source_pressure_gpa, "source_pressure_gpa")
    source_temperature = (
        source_standard.reference_temperature
        if temperature_k is None or not source_standard.is_thermal
        else temperature_k
    )
    target_temperature = (
        target_standard.reference_temperature
        if not target_standard.is_thermal
        else (
            source_standard.reference_temperature
            if temperature_k is None
            else temperature_k
        )
    )
    implied_volume = np.asarray(
        source_standard.volume(
            source_pressure,
            source_temperature,
            check_validity=check_validity,
        ),
        dtype=float,
    )
    target_pressure = np.asarray(
        target_standard.pressure(
            implied_volume,
            target_temperature,
            check_validity=check_validity,
        ),
        dtype=float,
    )
    try:
        source_pressure, implied_volume, target_pressure = np.broadcast_arrays(
            source_pressure, implied_volume, target_pressure
        )
    except ValueError as error:
        raise ValidationError(
            "source pressures and standard temperature must be broadcast-compatible"
        ) from error
    return PressureScaleRecalculation(
        source_pressure_gpa=_result(source_pressure),
        implied_standard_volume=_result(implied_volume),
        target_pressure_gpa=_result(target_pressure),
        pressure_difference_gpa=_result(target_pressure - source_pressure),
        source_standard_eos_record=source_standard_eos_record,
        target_standard_eos_record=target_standard_eos_record,
    )


def recalculate_ruby_with_measured_xrd_standard(
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
    if not isinstance(source, RubyFluorescenceCalibration):
        raise ValidationError(
            "source_calibration_record must identify a ruby calibration",
            field="source_calibration_record",
        )
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


def recalculate_ruby_to_xrd_pressure(
    source_pressure_gpa: ArrayLike,
    source_calibration_record: str,
    bridge_eos_record: str,
    target_eos_record: str,
    temperature_k: ArrayLike | None = None,
    *,
    check_validity: bool = False,
) -> RubyPressureScaleRecalculation:
    """Transform ruby-based pressure through two EOSs of the same standard.

    The ruby-linked ``bridge_eos_record`` is inverted to obtain the virtual
    standard volume implied by each source pressure.  ``target_eos_record`` is
    evaluated at that same volume.  No measured XRD-standard volume is needed.

    If the source and bridge EOS use different ruby calibrations, the source
    pressure is first converted through their common corrected R1 ratio.
    """
    source = get_pressure_calibration(source_calibration_record)
    if not isinstance(source, RubyFluorescenceCalibration):
        raise ValidationError(
            "source_calibration_record must identify a ruby calibration",
            field="source_calibration_record",
        )
    bridge_calibration_record = _ruby_calibration_for_eos(bridge_eos_record)
    bridge_calibration = get_pressure_calibration(bridge_calibration_record)
    if not isinstance(bridge_calibration, RubyFluorescenceCalibration):
        raise ValidationError(
            "bridge_eos_record must link to a ruby calibration",
            field="bridge_eos_record",
        )
    bridge_material, bridge_document, bridge = _resolve_bundled_eos_record(
        bridge_eos_record
    )
    target_material, target_document, target = _resolve_bundled_eos_record(
        target_eos_record
    )
    if bridge_material != target_material or (
        bridge_document.get("formula_units_per_cell")
        != target_document.get("formula_units_per_cell")
    ):
        raise ValidationError(
            "bridge and target EOS records must describe the same XRD standard "
            "and conventional unit cell",
            field="target_eos_record",
        )

    source_pressure = _as_finite_array(source_pressure_gpa, "source_pressure_gpa")
    source_ratio = np.asarray(source.wavelength_ratio(source_pressure), dtype=float)
    bridge_pressure = np.asarray(
        bridge_calibration.pressure_from_ratio(source_ratio), dtype=float
    )
    bridge_temperature = (
        bridge.reference_temperature
        if temperature_k is None or not bridge.is_thermal
        else temperature_k
    )
    target_temperature = (
        target.reference_temperature
        if not target.is_thermal
        else (
            bridge.reference_temperature if temperature_k is None else temperature_k
        )
    )
    implied_volume = np.asarray(
        bridge.volume(
            bridge_pressure,
            bridge_temperature,
            check_validity=check_validity,
        ),
        dtype=float,
    )
    target_pressure = np.asarray(
        target.pressure(
            implied_volume,
            target_temperature,
            check_validity=check_validity,
        ),
        dtype=float,
    )
    try:
        (
            source_pressure,
            source_ratio,
            bridge_pressure,
            implied_volume,
            target_pressure,
        ) = np.broadcast_arrays(
            source_pressure,
            source_ratio,
            bridge_pressure,
            implied_volume,
            target_pressure,
        )
    except ValueError as error:
        raise ValidationError(
            "source pressures and reference temperature must be broadcast-compatible"
        ) from error
    difference = target_pressure - source_pressure
    return RubyPressureScaleRecalculation(
        source_pressure_gpa=_result(source_pressure),
        source_wavelength_ratio=_result(source_ratio),
        bridge_pressure_gpa=_result(bridge_pressure),
        implied_standard_volume=_result(implied_volume),
        target_pressure_gpa=_result(target_pressure),
        pressure_difference_gpa=_result(difference),
        source_calibration_record=source.identifier,
        bridge_calibration_record=bridge_calibration.identifier,
        bridge_eos_record=bridge_eos_record,
        target_standard_eos_record=target_eos_record,
    )


def list_cross_calibration_edges() -> tuple[str, ...]:
    """Return identifiers from the explicit cross-calibration edge registry."""
    return tuple(
        edge["identifier"]
        for edge in pressure_calibration_library().get(
            "cross_calibration_edges", []
        )
    )


def get_cross_calibration_edge(identifier: str) -> dict[str, Any]:
    """Return one explicit, source-documented cross-calibration edge."""
    matches = [
        edge
        for edge in pressure_calibration_library().get(
            "cross_calibration_edges", []
        )
        if edge["identifier"] == identifier
    ]
    if len(matches) != 1:
        if not matches:
            raise MaterialLookupError(
                f"Unknown cross-calibration edge {identifier!r}; available: "
                f"{list_cross_calibration_edges()}"
            )
        raise ValidationError(f"Duplicate cross-calibration edge {identifier!r}")
    return copy.deepcopy(matches[0])


def _calibration_graph_edges() -> tuple[dict[str, Any], ...]:
    """Build explicit and mechanically valid calibration-graph edges."""
    from peritheos.eosmat import (
        get_material_document,
        list_material_documents,
    )

    library = pressure_calibration_library()
    edges: list[dict[str, Any]] = [
        copy.deepcopy(edge)
        for edge in library.get("cross_calibration_edges", [])
        if edge.get("executable", True)
    ]

    # An EOS pressure is expressed on the pressure basis used to reduce its
    # observations. These ancestry links preserve pressure numerically.
    documents = {
        material_identifier: get_material_document(material_identifier)
        for material_identifier in list_material_documents()
    }
    for document in documents.values():
        for record in document["eos_records"]:
            record_identifier = record["identifier"]
            for method_index, method in enumerate(
                record.get("pressure_calibration", {}).get("methods", [])
            ):
                reference_node = method.get("reference_eos_record") or method.get(
                    "reference_calibration_record"
                )
                if reference_node is None:
                    continue
                edges.append(
                    {
                        "identifier": (
                            f"ancestry__{record_identifier}__{method_index}"
                        ),
                        "kind": "calibration_ancestry",
                        "source_node": record_identifier,
                        "target_node": reference_node,
                        "transformation": "pressure_identity",
                        "bidirectional": True,
                        "generated": True,
                    }
                )

    # Every pair of executable EOSs for the same registered physical marker
    # can be crossed through that marker's virtual conventional-cell volume.
    standards = set(library.get("pressure_standard_records", []))
    standards.update(list_xrd_pressure_standards())
    metadata: dict[str, tuple[str, Any]] = {}
    standard_keys: set[tuple[str, Any]] = set()
    for material_identifier, document in documents.items():
        for record in document["eos_records"]:
            if record["identifier"] in standards:
                standard_keys.add(
                    (
                        material_identifier,
                        document.get("formula_units_per_cell"),
                    )
                )
    for material_identifier, document in documents.items():
        key = (
            material_identifier,
            document.get("formula_units_per_cell"),
        )
        if key not in standard_keys:
            continue
        for record in document["eos_records"]:
            metadata[record["identifier"]] = (
                material_identifier,
                document.get("formula_units_per_cell"),
            )
    grouped: dict[tuple[str, Any], list[str]] = {}
    for record_identifier, key in metadata.items():
        grouped.setdefault(key, []).append(record_identifier)
    for marker_records in grouped.values():
        ordered = sorted(marker_records)
        for index, source_node in enumerate(ordered):
            for target_node in ordered[index + 1 :]:
                edges.append(
                    {
                        "identifier": f"same_marker__{source_node}__{target_node}",
                        "kind": "same_marker_eos",
                        "source_node": source_node,
                        "target_node": target_node,
                        "transformation": "same_marker_volume",
                        "bidirectional": True,
                        "generated": True,
                    }
                )

    # Optical scales of the same measured feature have an exact common
    # coordinate even when their pressure equations differ.
    calibrations_by_kind: dict[str, list[str]] = {}
    for calibration in library["calibrations"]:
        calibrations_by_kind.setdefault(calibration["kind"], []).append(
            calibration["identifier"]
        )
    optical_transformations = {
        "ruby_fluorescence": "ruby_wavelength_ratio",
        "diamond_raman": "diamond_wavenumber_ratio",
    }
    for kind, transformation in optical_transformations.items():
        ordered = sorted(calibrations_by_kind.get(kind, []))
        for index, source_node in enumerate(ordered):
            for target_node in ordered[index + 1 :]:
                edges.append(
                    {
                        "identifier": f"same_signal__{source_node}__{target_node}",
                        "kind": "same_optical_signal",
                        "source_node": source_node,
                        "target_node": target_node,
                        "transformation": transformation,
                        "bidirectional": True,
                        "generated": True,
                    }
                )
    return tuple(edges)


def find_pressure_calibration_path(
    source_node: str, target_node: str
) -> tuple[dict[str, Any], ...]:
    """Find the shortest executable path between two pressure-scale nodes.

    Nodes are EOS-record or optical-calibration identifiers. Explicit
    literature edges are combined with generated same-marker, same-signal,
    and audited EOS-ancestry edges. Returned edge documents are oriented in
    traversal order.
    """
    if source_node == target_node:
        return ()
    adjacency: dict[str, list[dict[str, Any]]] = {}
    for edge in _calibration_graph_edges():
        adjacency.setdefault(edge["source_node"], []).append(edge)
        if edge.get("bidirectional", False):
            reverse = copy.deepcopy(edge)
            reverse["source_node"], reverse["target_node"] = (
                edge["target_node"],
                edge["source_node"],
            )
            reverse["reverse"] = True
            adjacency.setdefault(reverse["source_node"], []).append(reverse)

    queue: deque[tuple[str, tuple[dict[str, Any], ...]]] = deque(
        [(source_node, ())]
    )
    visited = {source_node}
    while queue:
        node, path = queue.popleft()
        for edge in sorted(
            adjacency.get(node, []),
            key=lambda item: (item.get("generated", False), item["identifier"]),
        ):
            next_node = edge["target_node"]
            next_path = path + (edge,)
            if next_node == target_node:
                return tuple(copy.deepcopy(item) for item in next_path)
            if next_node not in visited:
                visited.add(next_node)
                queue.append((next_node, next_path))
    raise ValidationError(
        f"No executable pressure-calibration path from {source_node!r} to "
        f"{target_node!r}",
        field="target_node",
    )


def recalculate_pressure_calibration_path(
    source_pressure_gpa: ArrayLike,
    source_node: str,
    target_node: str,
    temperature_k: ArrayLike | None = None,
    *,
    check_validity: bool = False,
) -> CalibrationPathRecalculation:
    """Execute a recursively discovered pressure-calibration path."""
    original = _as_finite_array(source_pressure_gpa, "source_pressure_gpa")
    current = np.asarray(original, dtype=float)
    current_temperature = (
        None
        if temperature_k is None
        else _as_finite_array(temperature_k, "temperature_k")
    )
    path = find_pressure_calibration_path(source_node, target_node)
    nodes = [source_node]
    states: list[Mapping[str, Any]] = []
    for edge in path:
        transformation = edge["transformation"]
        source = edge["source_node"]
        target = edge["target_node"]
        state: dict[str, Any] = {
            "edge_identifier": edge["identifier"],
            "source_node": source,
            "target_node": target,
            "transformation": transformation,
            "source_pressure_gpa": _result(np.asarray(current, dtype=float)),
        }
        if current_temperature is not None:
            state["source_temperature_k"] = _result(current_temperature)
        if transformation == "pressure_identity":
            recalculated = current
        elif transformation == "same_marker_volume":
            result = recalculate_xrd_pressure_scale(
                current,
                source,
                target,
                current_temperature,
                check_validity=check_validity,
            )
            recalculated = np.asarray(result.target_pressure_gpa, dtype=float)
            state["implied_standard_volume"] = result.implied_standard_volume
        elif transformation == "ruby_wavelength_ratio":
            recalculated = np.asarray(
                recalculate_ruby_pressure(current, source, target), dtype=float
            )
        elif transformation == "diamond_wavenumber_ratio":
            recalculated = np.asarray(
                recalculate_diamond_raman_pressure(current, source, target),
                dtype=float,
            )
        else:  # pragma: no cover - protected by bundled-data tests
            raise ValidationError(
                f"Unsupported calibration-edge transformation {transformation!r}"
            )
        temperature_transformation = edge.get("temperature_transformation")
        if temperature_transformation is not None:
            if current_temperature is None:
                raise ValidationError(
                    f"Calibration edge {edge['identifier']!r} requires temperature_k",
                    field="temperature_k",
                )
            scale = float(temperature_transformation["scale"])
            offset = float(temperature_transformation["offset_k"])
            if edge.get("reverse", False):
                if scale == 0.0:  # pragma: no cover - bundled-data validation
                    raise ValidationError(
                        f"Calibration edge {edge['identifier']!r} has a zero "
                        "temperature scale"
                    )
                current_temperature = (current_temperature - offset) / scale
            else:
                current_temperature = scale * current_temperature + offset
            if not np.all(np.isfinite(current_temperature)) or np.any(
                current_temperature <= 0.0
            ):
                raise ValidationError(
                    f"Calibration edge {edge['identifier']!r} produced an invalid "
                    "temperature",
                    field="temperature_k",
                )
            state["target_temperature_k"] = _result(current_temperature)
        state["target_pressure_gpa"] = _result(
            np.asarray(recalculated, dtype=float)
        )
        states.append(MappingProxyType(state))
        current = np.asarray(recalculated, dtype=float)
        nodes.append(target)
    original, current = np.broadcast_arrays(original, current)
    return CalibrationPathRecalculation(
        source_pressure_gpa=_result(original),
        target_pressure_gpa=_result(current),
        pressure_difference_gpa=_result(current - original),
        source_node=source_node,
        target_node=target_node,
        calibration_path=tuple(nodes),
        edge_identifiers=tuple(edge["identifier"] for edge in path),
        intermediate_states=tuple(states),
    )


def recalculate_eos_pressure_scale(
    source_eos_record: str,
    sample_volume: ArrayLike,
    target_standard_eos_record: str,
    sample_temperature_k: ArrayLike | None = None,
    *,
    bridge_eos_record: str | None = None,
    standard_temperature_k: ArrayLike | None = None,
    check_validity: bool = False,
) -> (
    PressureScaleRecalculation
    | RubyPressureScaleRecalculation
    | CalibrationPathRecalculation
):
    """Evaluate a sample EOS and transform it to a target pressure scale.

    The recursive calibration graph follows the source record's audited
    ancestry and may traverse same-marker EOS pairs, optical scales sharing a
    measured signal, and explicit cross-calibration-family edges. Supplying
    ``bridge_eos_record`` retains the direct, user-selected ruby-to-XRD route.
    """
    from peritheos.eosmat import get_eos_record_document

    _, _, source_eos = _resolve_bundled_eos_record(source_eos_record)
    source_pressure = source_eos.pressure(
        sample_volume,
        sample_temperature_k,
        check_validity=check_validity,
    )
    record = get_eos_record_document(source_eos_record)
    ruby_sources = {
        method["reference_calibration_record"]
        for method in record["pressure_calibration"]["methods"]
        if method["kind"] == "ruby_fluorescence"
        and "reference_calibration_record" in method
    }
    if bridge_eos_record is not None and len(ruby_sources) == 1:
        return recalculate_ruby_to_xrd_pressure(
            source_pressure,
            next(iter(ruby_sources)),
            bridge_eos_record,
            target_standard_eos_record,
            standard_temperature_k,
            check_validity=check_validity,
        )
    graph_temperature = (
        sample_temperature_k
        if standard_temperature_k is None
        else standard_temperature_k
    )
    return recalculate_pressure_calibration_path(
        source_pressure,
        source_eos_record,
        target_standard_eos_record,
        graph_temperature,
        check_validity=check_validity,
    )


__all__ = [
    "CalibrationPathRecalculation",
    "DiamondRamanCalibration",
    "PressureScaleRecalculation",
    "RubyFluorescenceCalibration",
    "RubyPressureScaleRecalculation",
    "XrdPressureRecalculation",
    "get_pressure_calibration",
    "get_pressure_calibration_document",
    "get_cross_calibration_edge",
    "find_pressure_calibration_path",
    "list_cross_calibration_edges",
    "list_diamond_raman_calibrations",
    "list_pressure_calibrations",
    "list_ruby_xrd_bridges",
    "list_ruby_pressure_calibrations",
    "list_xrd_pressure_standards",
    "pressure_calibration_library",
    "recalculate_eos_pressure_scale",
    "recalculate_diamond_raman_pressure",
    "recalculate_pressure_calibration_path",
    "recalculate_ruby_pressure",
    "recalculate_ruby_to_xrd_pressure",
    "recalculate_ruby_with_measured_xrd_standard",
    "recalculate_xrd_pressure_scale",
    "xrd_standard_pressure",
]
