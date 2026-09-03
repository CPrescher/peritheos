"""Materials and EOS records with provenance and calibration coverage.

The public pressure unit is GPa. Catalog volumes are conventional unit-cell
volumes in angstrom^3/cell because that is the directly measured quantity in
the source diffraction studies. Thermal EOS volume conversions are private
implementation details and are recorded in serialized metadata.
"""

from __future__ import annotations

import copy
import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np
from scipy.constants import Avogadro, electron_volt
from scipy.special import ndtri

from peritheos.eos import EosBase, NumericType, ThermalEOS
from peritheos.eos.rt import (
    BM2,
    BM3,
    BM4,
    Holzapfel,
    ModifiedTait,
    Murnaghan,
    NaturalStrain2,
    NaturalStrain3,
    NaturalStrain4,
    Vinet,
)
from peritheos.eos.thermal import (
    DoubleDebyeHelmholtz,
    DoubleDebyeLogMomentHelmholtz,
    LinearThermalPressure,
    LogVolumeThermalPressure,
    MieGruneisenDebye,
    MieGruneisenEinstein,
    MultiOscillatorGruneisenThermalEOS,
    Tange2009Debye,
    ThermalModifiedTait,
    ThermalReferenceStateEOS,
)
from peritheos.errors import (
    ConfigurationError,
    MaterialError,
    MaterialLookupError,
)
from peritheos.uncertainty import EOSUncertainty, PredictionUncertainty


def _validated_aliases(
    aliases: Iterable[str], *, identifier: str, kind: str, strict: bool = True
) -> tuple[str, ...]:
    if isinstance(aliases, str):
        raise MaterialError(f"{kind} aliases must be an iterable of strings")
    try:
        normalized = tuple(aliases)
    except TypeError as error:
        raise MaterialError(f"{kind} aliases must be an iterable of strings") from error
    if any(not isinstance(alias, str) for alias in normalized):
        raise MaterialError(f"{kind} aliases must contain only strings")
    if strict:
        if any(not alias.strip() for alias in normalized):
            raise MaterialError(f"{kind} aliases must be non-empty strings")
        if any(alias != alias.strip() for alias in normalized):
            raise MaterialError(f"{kind} aliases must not have surrounding whitespace")
        if len(normalized) != len(set(normalized)):
            raise MaterialError(f"{kind} aliases must be unique")
        if identifier in normalized:
            raise MaterialError(f"{kind} aliases must not repeat its identifier")
    return normalized


@dataclass(frozen=True)
class LiteratureReference:
    """A primary publication and the exact locations used by an entry."""

    authors: str
    year: int
    title: str
    doi: str
    locations: tuple[str, ...]


@dataclass(frozen=True)
class ValidityRange:
    """Published calibration/data envelope, not an extrapolation limit.

    The historical public name is retained for compatibility.  These bounds
    describe the states used to fit or assess a parameterization; they do not
    delimit where its equations may be evaluated or model phase stability.
    """

    pressure_gpa: tuple[float, float]
    temperature_k: tuple[float, float]
    volume_ratio: tuple[float, float] | None = None
    notes: tuple[str, ...] = ()

    def contains(
        self,
        pressure: Any,
        temperature: Any,
        volume_ratio: Any,
    ) -> np.ndarray:
        """Return whether broadcast states lie inside the marginal envelope."""
        pressure, temperature, volume_ratio = np.broadcast_arrays(
            np.asarray(pressure, dtype=float),
            np.asarray(temperature, dtype=float),
            np.asarray(volume_ratio, dtype=float),
        )

        def closed_interval(values, bounds):
            lower, upper = bounds
            lower_tolerance = (
                16.0 * np.finfo(float).eps * max(1.0, abs(lower))
                if np.isfinite(lower)
                else 0.0
            )
            upper_tolerance = (
                16.0 * np.finfo(float).eps * max(1.0, abs(upper))
                if np.isfinite(upper)
                else 0.0
            )
            return (values >= lower - lower_tolerance) & (
                values <= upper + upper_tolerance
            )

        valid = closed_interval(pressure, self.pressure_gpa) & closed_interval(
            temperature, self.temperature_k
        )
        if self.volume_ratio is not None:
            valid &= closed_interval(volume_ratio, self.volume_ratio)
        return valid


@dataclass(frozen=True)
class DeferredEOSRecord:
    """A catalog candidate excluded from the validated implementation."""

    material: str
    phase: str
    references: tuple[str, ...]
    reason: str


@dataclass(frozen=True)
class EOSRecord:
    """A literature-specific EOS parameterization for one material phase.

    Parameters are stored in the EOS's working volume unit. ``volume_scale``
    converts the public unit-cell volume to that working unit. The public API
    always returns pressure in GPa and accepts temperature in K.
    """

    identifier: str
    name: str
    material: str
    phase: str
    cell_contents: str
    eos: EosBase
    reference_temperature: float
    reference: LiteratureReference
    validity: ValidityRange
    parameter_provenance: Mapping[str, str]
    parameter_errors: Mapping[str, float] | None = None
    parameter_error_confidence: float | None = None
    parameter_covariance: tuple[tuple[float, ...], ...] | None = None
    covariance_parameters: tuple[str, ...] | None = None
    notes: tuple[str, ...] = ()
    volume_unit: str = "angstrom^3/conventional_unit_cell"
    volume_scale: float = 1.0
    scientific_validation_status: str = "primary_source_validated"
    scientific_validation_note: str = ""
    eosmat_metadata: Mapping[str, Any] = field(
        default_factory=dict, compare=False, repr=False
    )
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        allowed = {
            "primary_source_validated",
            "pending_primary_source_check",
            "deferred",
        }
        if self.scientific_validation_status not in allowed:
            raise MaterialError(
                f"scientific_validation_status must be one of {sorted(allowed)}"
            )
        if self.parameter_error_confidence is not None and not (
            0.0 < self.parameter_error_confidence < 1.0
        ):
            raise MaterialError(
                "parameter_error_confidence must lie between zero and one"
            )
        object.__setattr__(
            self,
            "eosmat_metadata",
            MappingProxyType(copy.deepcopy(dict(self.eosmat_metadata))),
        )
        object.__setattr__(
            self,
            "aliases",
            _validated_aliases(
                self.aliases, identifier=self.identifier, kind="EOS record"
            ),
        )

    @property
    def reference_volume(self) -> float:
        """Reference unit-cell volume in ``volume_unit``."""
        reference_eos = (
            self.eos.rt_eos if isinstance(self.eos, ThermalEOS) else self.eos
        )
        return float(reference_eos.V0 / self.volume_scale)

    @property
    def is_thermal(self) -> bool:
        """Whether temperature contributes to the pressure equation."""
        return isinstance(self.eos, ThermalEOS)

    def _temperature(self, temperature: Any | None) -> np.ndarray:
        if temperature is None:
            temperature = self.reference_temperature
        values = np.asarray(temperature, dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise MaterialError("Temperature must be finite and greater than zero")
        if not self.is_thermal and not np.allclose(
            values, self.reference_temperature, rtol=0.0, atol=1.0e-8
        ):
            raise MaterialError(
                f"{self.identifier} is an isothermal {self.reference_temperature:g} K "
                "EOS record and has no thermal correction"
            )
        return values

    def _validate_range(
        self, pressure: Any, temperature: Any, public_volume: Any
    ) -> None:
        ratio = np.asarray(public_volume, dtype=float) / self.reference_volume
        if not np.all(self.validity.contains(pressure, temperature, ratio)):
            raise MaterialError(
                f"State is outside the published calibration/data envelope for "
                f"{self.identifier}; the EOS can still be evaluated by leaving "
                "check_validity=False"
            )

    def pressure(
        self,
        volume: NumericType,
        temperature: Any | None = None,
        *,
        check_validity: bool = False,
    ) -> NumericType:
        """Calculate pressure in GPa from unit-cell volume and temperature."""
        temperatures = self._temperature(temperature)
        internal_volume = np.asarray(volume, dtype=float) * self.volume_scale
        if isinstance(self.eos, ThermalEOS):
            result = self.eos.pressure(internal_volume, temperatures)
        else:
            result = self.eos.pressure(internal_volume)
        if check_validity:
            self._validate_range(result, temperatures, volume)
        return result

    def volume(
        self,
        pressure: NumericType,
        temperature: Any | None = None,
        *,
        check_validity: bool = False,
    ) -> NumericType:
        """Invert unit-cell volume from pressure in GPa and temperature."""
        temperatures = self._temperature(temperature)
        if isinstance(self.eos, ThermalEOS):
            internal_volume = self.eos.volume(pressure, temperatures)
        else:
            internal_volume = self.eos.volume(pressure)
        result = np.asarray(internal_volume, dtype=float) / self.volume_scale
        if check_validity:
            self._validate_range(pressure, temperatures, result)
        if result.ndim == 0:
            return float(result)
        return result

    def thermal_pressure_increment(
        self,
        volume: NumericType,
        temperature: Any,
        *,
        check_validity: bool = False,
    ) -> NumericType:
        """Thermal pressure above the reference isotherm in GPa.

        ``volume`` uses this record's public conventional-cell volume unit.
        Unlike an absolute free-energy model's raw thermal contribution, this
        increment is zero at the record's reference temperature.
        """
        if not isinstance(self.eos, ThermalEOS):
            raise MaterialError(
                f"{self.identifier} is isothermal and has no thermal pressure"
            )
        temperatures = self._temperature(temperature)
        internal_volume = np.asarray(volume, dtype=float) * self.volume_scale
        result = self.eos.thermal_pressure_increment(internal_volume, temperatures)
        if check_validity:
            hot_pressure = self.eos.pressure(internal_volume, temperatures)
            self._validate_range(hot_pressure, temperatures, volume)
        return result

    def dac_thermal_pressure(
        self,
        volume: NumericType,
        temperature: Any,
        *,
        f_dac: float,
        check_validity: bool = False,
    ) -> NumericType:
        """Retained DAC pressure increment in GPa at a heated state."""
        if not isinstance(self.eos, ThermalEOS):
            raise MaterialError(
                f"{self.identifier} is isothermal and has no thermal pressure"
            )
        temperatures = self._temperature(temperature)
        internal_volume = np.asarray(volume, dtype=float) * self.volume_scale
        result = self.eos.dac_thermal_pressure(internal_volume, temperatures, f_dac)
        if check_validity:
            hot_pressure = self.eos.pressure(internal_volume, temperatures)
            self._validate_range(hot_pressure, temperatures, volume)
        return result

    def volume_with_dac_confinement(
        self,
        cold_pressure: NumericType,
        temperature: Any,
        *,
        f_dac: float,
        check_validity: bool = False,
    ) -> NumericType:
        """Predict heated cell volume from cold pressure and confinement.

        ``cold_pressure`` is the GPa pressure established at the record's
        reference temperature. The returned volume uses the record's public
        conventional-cell volume unit.
        """
        if not isinstance(self.eos, ThermalEOS):
            raise MaterialError(
                f"{self.identifier} is isothermal and cannot apply DAC confinement"
            )
        temperatures = self._temperature(temperature)
        internal_volume = self.eos.volume_with_dac_confinement(
            cold_pressure,
            temperatures,
            f_dac=f_dac,
        )
        result = np.asarray(internal_volume, dtype=float) / self.volume_scale
        if check_validity:
            hot_pressure = self.eos.pressure(internal_volume, temperatures)
            self._validate_range(hot_pressure, temperatures, result)
        if result.ndim == 0:
            return float(result)
        return result

    def temperature_from_volumes(
        self,
        ambient_volume: NumericType,
        heated_volume: NumericType,
        *,
        f_dac: float,
        check_validity: bool = False,
    ) -> NumericType:
        """Infer DAC temperature from two public unit-cell volumes.

        The volumes use this record's public ``volume_unit``. ``f_dac`` is the
        fraction of the EOS thermal-pressure increment retained as an increase
        above the reference-temperature pressure. See the underlying thermal
        EOS :meth:`temperature_from_volumes` method for model-specific details.
        """
        if not isinstance(self.eos, ThermalEOS):
            raise MaterialError(
                f"{self.identifier} is isothermal and cannot invert temperature"
            )
        ambient_internal = np.asarray(ambient_volume, dtype=float) * self.volume_scale
        heated_internal = np.asarray(heated_volume, dtype=float) * self.volume_scale
        result = self.eos.temperature_from_volumes(
            ambient_internal,
            heated_internal,
            f_dac=f_dac,
        )
        if check_validity:
            hot_pressure = self.eos.pressure(heated_internal, result)
            self._validate_range(hot_pressure, result, heated_volume)
        return result

    def within_validity(
        self, volume: NumericType, temperature: Any | None = None
    ) -> bool | np.ndarray:
        """Compatibility alias for :meth:`within_calibration_range`."""
        return self.within_calibration_range(volume, temperature)

    def within_calibration_range(
        self, volume: NumericType, temperature: Any | None = None
    ) -> bool | np.ndarray:
        """Test whether states lie in the published calibration/data envelope."""
        temperatures = self._temperature(temperature)
        pressure = self.pressure(volume, temperatures, check_validity=False)
        ratio = np.asarray(volume, dtype=float) / self.reference_volume
        result = self.validity.contains(pressure, temperatures, ratio)
        if result.ndim == 0:
            return bool(result)
        return result

    def _uncertainty(self) -> EOSUncertainty:
        if self.parameter_covariance is not None:
            return EOSUncertainty(
                self.eos,
                covariance=self.parameter_covariance,
                parameter_names=self.covariance_parameters,
            )
        if self.parameter_errors is not None:
            errors = self.parameter_errors
            if self.parameter_error_confidence is not None:
                normal_half_width = ndtri(0.5 * (1.0 + self.parameter_error_confidence))
                errors = {
                    name: value / normal_half_width
                    for name, value in self.parameter_errors.items()
                }
            return EOSUncertainty(self.eos, parameter_errors=errors)
        return EOSUncertainty.state_only(self.eos)

    def pressure_with_uncertainty(
        self,
        volume: Any,
        temperature: Any | None = None,
        *,
        volume_sigma: Any | None = None,
        temperature_sigma: Any | None = None,
        check_validity: bool = False,
        **options: Any,
    ) -> PredictionUncertainty:
        """Propagate published parameter and measured V/T uncertainty."""
        temperatures = self._temperature(temperature)
        if check_validity:
            self.pressure(volume, temperatures, check_validity=True)
        internal_volume = np.asarray(volume, dtype=float) * self.volume_scale
        internal_sigma = (
            None
            if volume_sigma is None
            else np.asarray(volume_sigma, dtype=float) * self.volume_scale
        )
        uncertainty = self._uncertainty()
        if self.is_thermal:
            return uncertainty.pressure(
                internal_volume,
                temperatures,
                volume_sigma=internal_sigma,
                temperature_sigma=temperature_sigma,
                **options,
            )
        if temperature_sigma is not None:
            raise MaterialError(
                "Temperature uncertainty cannot be propagated through an "
                "isothermal EOS record without a published thermal model"
            )
        return uncertainty.pressure(
            internal_volume, volume_sigma=internal_sigma, **options
        )

    def volume_with_uncertainty(
        self,
        pressure: Any,
        temperature: Any | None = None,
        *,
        pressure_sigma: Any | None = None,
        temperature_sigma: Any | None = None,
        check_validity: bool = False,
        **options: Any,
    ) -> PredictionUncertainty:
        """Propagate uncertainty into an inverted unit-cell volume."""
        temperatures = self._temperature(temperature)
        if check_validity:
            self.volume(pressure, temperatures, check_validity=True)
        uncertainty = self._uncertainty()
        if self.is_thermal:
            prediction = uncertainty.volume(
                pressure,
                temperatures,
                pressure_sigma=pressure_sigma,
                temperature_sigma=temperature_sigma,
                **options,
            )
        else:
            if temperature_sigma is not None:
                raise MaterialError(
                    "Temperature uncertainty cannot be propagated through an "
                    "isothermal EOS record without a published thermal model"
                )
            prediction = uncertainty.volume(
                pressure, pressure_sigma=pressure_sigma, **options
            )
        scale = self.volume_scale
        covariance = (
            None if prediction.covariance is None else prediction.covariance / scale**2
        )

        def public_volume(value: Any) -> NumericType:
            result = np.asarray(value, dtype=float) / scale
            if result.ndim == 0:
                return float(result)
            return result

        return PredictionUncertainty(
            value=public_volume(prediction.value),
            standard_error=public_volume(prediction.standard_error),
            lower=public_volume(prediction.lower),
            upper=public_volume(prediction.upper),
            covariance=covariance,
            method=prediction.method,
            confidence=prediction.confidence,
            assumptions=prediction.assumptions,
            rejected_fraction=prediction.rejected_fraction,
        )


@dataclass(frozen=True)
class Material:
    """A material phase with one or more literature-specific EOS records."""

    identifier: str
    formula: str
    phase: str
    cell_contents: str
    eos_records: tuple[EOSRecord, ...]
    volume_unit: str = "angstrom^3/conventional_unit_cell"
    name: str | None = None
    symmetry: str | None = None
    lattice: Mapping[str, float | None] | None = None
    formula_units_per_cell: float | None = None
    space_group: str | None = None
    space_group_number: int | None = None
    atom_sites: tuple[Mapping[str, Any], ...] = ()
    peaks: tuple[tuple[float, ...], ...] = ()
    eosmat_metadata: Mapping[str, Any] = field(
        default_factory=dict, compare=False, repr=False
    )
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.identifier:
            raise MaterialError("Material identifier must not be empty")
        if not self.eos_records:
            raise MaterialError("A material requires at least one EOS record")
        identifiers = [record.identifier for record in self.eos_records]
        if len(identifiers) != len(set(identifiers)):
            raise MaterialError(
                "EOS record identifiers must be unique within a material"
            )
        record_names = [
            name
            for record in self.eos_records
            for name in (record.identifier, *record.aliases)
        ]
        if len(record_names) != len(set(record_names)):
            raise MaterialError(
                "EOS record identifiers and aliases must be unique within a material"
            )
        expected = (self.formula, self.phase, self.cell_contents, self.volume_unit)
        for record in self.eos_records:
            actual = (
                record.material,
                record.phase,
                record.cell_contents,
                record.volume_unit,
            )
            if actual != expected:
                raise MaterialError(
                    "Every EOS record must match its material identity and volume unit"
                )
        if self.formula_units_per_cell is not None:
            if (
                not np.isfinite(self.formula_units_per_cell)
                or self.formula_units_per_cell <= 0.0
            ):
                raise MaterialError(
                    "formula_units_per_cell must be positive and finite"
                )
        if self.space_group_number is not None and not (
            1 <= self.space_group_number <= 230
        ):
            raise MaterialError("space_group_number must be between 1 and 230")
        if self.lattice is not None:
            object.__setattr__(
                self,
                "lattice",
                MappingProxyType(copy.deepcopy(dict(self.lattice))),
            )
        object.__setattr__(
            self,
            "atom_sites",
            tuple(
                MappingProxyType(copy.deepcopy(dict(site))) for site in self.atom_sites
            ),
        )
        object.__setattr__(
            self,
            "eosmat_metadata",
            MappingProxyType(copy.deepcopy(dict(self.eosmat_metadata))),
        )
        object.__setattr__(
            self,
            "aliases",
            _validated_aliases(
                self.aliases,
                identifier=self.identifier,
                kind="Material",
                strict=False,
            ),
        )

    def get_eos_record(self, identifier: str) -> EOSRecord:
        """Return one EOS parameterization by its stable identifier."""
        for record in self.eos_records:
            if identifier == record.identifier or identifier in record.aliases:
                return record
        from difflib import get_close_matches

        choices = {
            choice
            for record in self.eos_records
            for choice in (record.identifier, *record.aliases)
        }
        suggestions = get_close_matches(identifier, sorted(choices), n=3, cutoff=0.5)
        hint = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
        raise MaterialLookupError(
            f"Unknown EOS record {identifier!r} for {self.identifier!r}; "
            f"available: {[record.identifier for record in self.eos_records]}."
            f"{hint}",
            operation="lookup_eos_record",
            field="identifier",
            context={"identifier": identifier, "material": self.identifier},
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-safe ``.eosmat`` format-3 document."""
        return self.to_eosmat()

    def to_eosmat(self) -> dict[str, Any]:
        """Return one executable and Dioptas-compatible material document."""
        return _material_to_eosmat(self)

    def to_snapshot_dict(self) -> dict[str, Any]:
        """Return the deprecated executable snapshot-v2 representation."""
        return _eos_records_to_document(self.eos_records, self.identifier)

    @classmethod
    def from_dict(cls, document: Mapping[str, Any]) -> "Material":
        """Load canonical format 3 or the deprecated snapshot format 2."""
        return material_from_dict(document)

    @classmethod
    def from_eosmat(
        cls,
        document: Mapping[str, Any],
        *,
        require_primary_validation: bool = True,
        record_identifiers: Iterable[str] | None = None,
    ) -> "Material":
        """Construct an executable material from canonical ``.eosmat`` data.

        By default every record must have been checked against its primary
        source. Pass ``require_primary_validation=False`` only when explicitly
        accepting migrated or otherwise unaudited parameters.
        """
        return _material_from_eosmat(
            document,
            require_primary_validation=require_primary_validation,
            record_identifiers=record_identifiers,
        )


_MODEL_IDENTIFIERS = MappingProxyType(
    {
        "BM2": "birch_murnaghan_2",
        "BM3": "birch_murnaghan_3",
        "BM4": "birch_murnaghan_4",
        "Vinet": "vinet",
        "Holzapfel": "holzapfel",
        "ModifiedTait": "modified_tait",
        "Murnaghan": "murnaghan",
        "NaturalStrain2": "natural_strain_2",
        "NaturalStrain3": "natural_strain_3",
        "NaturalStrain4": "natural_strain_4",
        "DoubleDebyeHelmholtz": "double_debye_helmholtz",
        "DoubleDebyeLogMomentHelmholtz": ("double_debye_log_moment_helmholtz"),
        "LinearThermalPressure": "linear_thermal_pressure",
        "LogVolumeThermalPressure": "log_volume_thermal_pressure",
        "ThermalReferenceStateEOS": "thermal_reference_state",
        "MieGruneisenDebye": "mie_gruneisen_debye",
        "MieGruneisenEinstein": "mie_gruneisen_einstein",
        "Tange2009Debye": "asymptotic_power_law_mie_gruneisen_debye",
        "MultiOscillatorGruneisenThermalEOS": (
            "multi_oscillator_gruneisen_thermal_pressure"
        ),
        "ThermalModifiedTait": "thermal_modified_tait",
    }
)

_MODEL_CLASSES = MappingProxyType(
    {
        _MODEL_IDENTIFIERS[model.__name__]: model
        for model in (
            BM2,
            BM3,
            BM4,
            Holzapfel,
            ModifiedTait,
            Murnaghan,
            NaturalStrain2,
            NaturalStrain3,
            NaturalStrain4,
            Vinet,
            DoubleDebyeHelmholtz,
            DoubleDebyeLogMomentHelmholtz,
            LinearThermalPressure,
            LogVolumeThermalPressure,
            ThermalReferenceStateEOS,
            MieGruneisenDebye,
            MieGruneisenEinstein,
            MultiOscillatorGruneisenThermalEOS,
            Tange2009Debye,
            ThermalModifiedTait,
        )
    }
)

_EOSMAT_TYPES = MappingProxyType(
    {
        "birch_murnaghan_2": "BM2",
        "birch_murnaghan_3": "BM3",
        "birch_murnaghan_4": "BM4",
        "vinet": "Vinet",
        "holzapfel": "Holzapfel",
        "modified_tait": "ModifiedTait",
        "murnaghan": "Murnaghan",
        "natural_strain_2": "NaturalStrain2",
        "natural_strain_3": "NaturalStrain3",
        "natural_strain_4": "NaturalStrain4",
        "double_debye_helmholtz": "DoubleDebyeHelmholtz",
        "double_debye_log_moment_helmholtz": ("DoubleDebyeLogMomentHelmholtz"),
        "linear_thermal_pressure": "LinearThermalPressure",
        "log_volume_thermal_pressure": "LogVolumeThermalPressure",
        "thermal_reference_state": "AlphaKT",
        "mie_gruneisen_debye": "MieGruneisenDebye",
        "mie_gruneisen_einstein": "MieGruneisenEinstein",
        "asymptotic_power_law_mie_gruneisen_debye": (
            "AsymptoticPowerLawMieGruneisenDebye"
        ),
        "multi_oscillator_gruneisen_thermal_pressure": ("MultiOscillatorGruneisen"),
        "thermal_modified_tait": "ThermalModifiedTait",
    }
)

_MOLAR_VOLUME_THERMAL_MODELS = frozenset(
    {
        "mie_gruneisen_debye",
        "double_debye_helmholtz",
        "double_debye_log_moment_helmholtz",
        "mie_gruneisen_einstein",
        "asymptotic_power_law_mie_gruneisen_debye",
        "multi_oscillator_gruneisen_thermal_pressure",
        "thermal_modified_tait",
    }
)


def _model_identifier(eos: EosBase) -> str:
    """Return a stable mechanism-oriented model identifier."""
    try:
        return _MODEL_IDENTIFIERS[type(eos).__name__]
    except KeyError as error:
        raise MaterialError(
            f"EOS class {type(eos).__name__!r} is not registered for catalog export"
        ) from error


def _plain_data(value: Any) -> Any:
    """Return nested mappings/tuples as independently mutable JSON-like data."""
    if isinstance(value, Mapping):
        return {str(key): _plain_data(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain_data(item) for item in value]
    return copy.deepcopy(value)


def _eosmat_component(eos: EosBase) -> dict[str, Any]:
    model = _model_identifier(eos)
    parameters: dict[str, Any] = eos.parameter_values(include_reference=False)
    if (
        model
        in {
            "double_debye_helmholtz",
            "double_debye_log_moment_helmholtz",
        }
        and getattr(eos, "_reference_isotherm_temperature", None) is None
    ):
        parameters["Tr"] = None
    component = {
        "type": _EOSMAT_TYPES[model],
        "model": model,
        "parameters": parameters,
    }
    configuration = eos.configuration_values()
    if configuration:
        component["configuration"] = configuration
        # Kept at the established location for Dioptas format-3 readers.
        for name in (
            "debye_temperature_law",
            "thermal_expansion_law",
            "reference_volume_law",
        ):
            if name in configuration:
                component[name] = configuration[name]
    return component


def _reference_to_eosmat(reference: LiteratureReference) -> dict[str, Any]:
    return {
        "authors": [reference.authors],
        "year": reference.year,
        "title": reference.title,
        "source": "primary publication",
        "doi": reference.doi,
        "details": "; ".join(reference.locations),
        "locations": list(reference.locations),
    }


def _merge_eosmat_component(
    metadata: Any, component: Mapping[str, Any]
) -> dict[str, Any]:
    """Update executable fields while retaining component extensions."""
    result = _plain_data(metadata) if isinstance(metadata, Mapping) else {}
    existing_parameters = result.get("parameters")
    parameters = (
        _plain_data(existing_parameters)
        if isinstance(existing_parameters, Mapping)
        else {}
    )
    parameters.update(_plain_data(component["parameters"]))
    for key, value in component.items():
        if key == "parameters":
            continue
        if key in {"type", "model"}:
            result.setdefault(key, _plain_data(value))
        else:
            result[key] = _plain_data(value)
    result["parameters"] = parameters
    return result


def _record_to_eosmat(record: EOSRecord) -> dict[str, Any]:
    reference_eos = (
        record.eos.rt_eos if isinstance(record.eos, ThermalEOS) else record.eos
    )
    thermal_eos = record.eos if record.is_thermal else None
    reference_errors, thermal_errors, _ = _component_parameter_mapping(
        record, record.parameter_errors
    )
    component_provenance = _component_parameter_mapping(
        record, record.parameter_provenance
    )
    result: dict[str, Any] = _plain_data(record.eosmat_metadata)
    reference_component = _merge_eosmat_component(
        result.get("eos"), _eosmat_component(reference_eos)
    )
    reference_component["parameters"]["V0"] = record.reference_volume
    stored_errors = result.get("parameter_errors")
    merged_errors = (
        _plain_data(stored_errors) if isinstance(stored_errors, Mapping) else {}
    )
    merged_errors.update(reference_errors)
    stored_volume = result.get("volume")
    volume = _plain_data(stored_volume) if isinstance(stored_volume, Mapping) else {}
    volume.update(
        {
            "reference_value": record.reference_volume,
            "public_to_model_scale": record.volume_scale,
            "model_unit": (
                record.volume_unit if record.volume_scale == 1.0 else "J bar^-1 mol^-1"
            ),
        }
    )
    stored_validation = result.get("scientific_validation")
    validation = (
        _plain_data(stored_validation) if isinstance(stored_validation, Mapping) else {}
    )
    validation.update(
        {
            "status": record.scientific_validation_status,
            "note": record.scientific_validation_note,
        }
    )
    result.update(
        {
            "identifier": record.identifier,
            "label": record.name,
            "eos": reference_component,
            "parameter_errors": merged_errors,
            "parameter_error_confidence": record.parameter_error_confidence,
            "temperature_ref": record.reference_temperature,
            "volume": volume,
            "scientific_validation": validation,
        }
    )
    if record.aliases:
        result["aliases"] = list(record.aliases)
    else:
        result.pop("aliases", None)
    result.setdefault("reference", _reference_to_eosmat(record.reference))
    result.setdefault("fixed_parameters", [])
    result.setdefault(
        "parameter_provenance",
        {
            "reference_isotherm": component_provenance[0],
            "thermal_correction": component_provenance[1],
            "additional": component_provenance[2],
        },
    )
    result.setdefault("notes", "\n".join(record.notes))
    if record.parameter_covariance is not None:
        result["parameter_covariance"] = {
            "matrix": [list(row) for row in record.parameter_covariance],
            "parameter_order": list(record.covariance_parameters or ()),
        }
    else:
        result.setdefault("parameter_covariance", None)
    stored_validity = result.get("validity")
    validity = (
        _plain_data(stored_validity) if isinstance(stored_validity, Mapping) else {}
    )
    validity["notes"] = list(record.validity.notes)
    if np.all(np.isfinite(record.validity.pressure_gpa)):
        pressure_range = list(record.validity.pressure_gpa)
        result.setdefault("experimental_pressure_range_gpa", pressure_range)
        validity["pressure_gpa"] = pressure_range
    if np.all(np.isfinite(record.validity.temperature_k)):
        temperature_range = list(record.validity.temperature_k)
        result.setdefault("experimental_temperature_range_k", temperature_range)
        validity["temperature_k"] = temperature_range
    if record.validity.volume_ratio is not None:
        validity["volume_ratio"] = list(record.validity.volume_ratio)
    if len(validity) > 1 or validity["notes"]:
        result["validity"] = validity
    if thermal_eos is not None:
        thermal = _merge_eosmat_component(
            result.get("thermal"), _eosmat_component(thermal_eos)
        )
        stored_thermal_errors = thermal.get("parameter_errors")
        merged_thermal_errors = (
            _plain_data(stored_thermal_errors)
            if isinstance(stored_thermal_errors, Mapping)
            else {}
        )
        merged_thermal_errors.update(thermal_errors)
        thermal["parameter_errors"] = merged_thermal_errors
        thermal.setdefault("fixed_parameters", [])
        result["thermal"] = thermal
    return result


def _material_to_eosmat(material: Material) -> dict[str, Any]:
    document = _plain_data(material.eosmat_metadata)
    document.update(
        {
            "format": "peritheos.material",
            "format_version": 3,
            "identifier": material.identifier,
            "name": material.name or material.formula,
            "formula": material.formula,
            "phase": material.phase,
            "cell_contents": material.cell_contents,
            "units": {
                "pressure": "GPa",
                "temperature": "K",
                "volume": material.volume_unit,
            },
            "eos_records": [
                _record_to_eosmat(record) for record in material.eos_records
            ],
        }
    )
    optional = {
        "symmetry": material.symmetry,
        "lattice": material.lattice,
        "formula_units_per_cell": material.formula_units_per_cell,
        "space_group": material.space_group,
        "space_group_number": material.space_group_number,
    }
    if material.aliases:
        document["aliases"] = list(material.aliases)
    else:
        document.pop("aliases", None)
    for key, value in optional.items():
        if value is not None:
            document[key] = _plain_data(value)
        else:
            document.pop(key, None)
    if material.atom_sites:
        document["atom_sites"] = _plain_data(material.atom_sites)
    else:
        document.pop("atom_sites", None)
    if material.peaks:
        document["peaks"] = _plain_data(material.peaks)
    else:
        document.pop("peaks", None)
    from peritheos.eosmat import validate_eosmat_document

    validate_eosmat_document(document)
    return document


def _equation_component(eos: EosBase) -> dict[str, Any]:
    """Serialize one independently selectable EOS component."""
    component = {
        "model": _model_identifier(eos),
        "implementation": f"{type(eos).__module__}.{type(eos).__name__}",
        "parameters": eos.parameter_values(include_reference=False),
    }
    configuration = eos.configuration_values()
    if configuration:
        component["configuration"] = configuration
    return component


def _component_parameter_mapping(
    record: EOSRecord, values: Mapping[str, Any] | None
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Split dotted EOS parameter metadata between equation components."""
    reference: dict[str, Any] = {}
    thermal: dict[str, Any] = {}
    additional: dict[str, Any] = {}
    if values is None:
        return reference, thermal, additional

    reference_eos = (
        record.eos.rt_eos if isinstance(record.eos, ThermalEOS) else record.eos
    )
    reference_names = set(reference_eos.parameter_values(include_reference=False))
    thermal_names = (
        set(record.eos.parameter_values(include_reference=False))
        if record.is_thermal
        else set()
    )
    for name, value in values.items():
        if record.is_thermal and name.startswith("rt_eos."):
            component_name = name.removeprefix("rt_eos.")
            if component_name in reference_names:
                reference[component_name] = value
            else:
                additional[name] = value
        elif name in thermal_names:
            thermal[name] = value
        elif name in reference_names:
            reference[name] = value
        else:
            additional[name] = value
    return reference, thermal, additional


def _qualified_parameter_name(record: EOSRecord, name: str) -> str:
    """Qualify covariance names with their version-2 component path."""
    if record.is_thermal and name.startswith("rt_eos."):
        return f"reference_isotherm.{name.removeprefix('rt_eos.')}"
    component = "thermal_correction" if record.is_thermal else "reference_isotherm"
    return f"{component}.{name}"


def _eos_record_dict(record: EOSRecord) -> dict[str, Any]:
    """Serialize one scale inside a material-oriented document."""
    reference_eos = (
        record.eos.rt_eos if isinstance(record.eos, ThermalEOS) else record.eos
    )
    thermal_eos = record.eos if record.is_thermal else None
    reference_provenance, thermal_provenance, additional_provenance = (
        _component_parameter_mapping(record, record.parameter_provenance)
    )
    reference_errors, thermal_errors, additional_errors = _component_parameter_mapping(
        record, record.parameter_errors
    )
    model_volume_unit = (
        record.volume_unit
        if np.isclose(record.volume_scale, 1.0, rtol=0.0, atol=0.0)
        else "J bar^-1 mol^-1"
    )
    return {
        "identifier": record.identifier,
        "name": record.name,
        "reference_temperature": {"value": record.reference_temperature, "unit": "K"},
        "volume": {
            "reference_value": record.reference_volume,
            "public_unit": record.volume_unit,
            "model_unit": model_volume_unit,
            "public_to_model_scale": record.volume_scale,
        },
        "equation": {
            "reference_isotherm": _equation_component(reference_eos),
            "thermal_correction": (
                None if thermal_eos is None else _equation_component(thermal_eos)
            ),
            "combination": {
                "status": "source_parameterization",
                "validated_as_composed": True,
                "note": (
                    "The catalog record pins the component combination reported by "
                    "the cited primary source. Substituting either component creates "
                    "a user-composed model outside this record's validation claim."
                ),
            },
        },
        "provenance": {
            "primary_reference": asdict(record.reference),
            "parameters": {
                "reference_isotherm": reference_provenance,
                "thermal_correction": thermal_provenance,
                "additional": additional_provenance,
            },
        },
        "validity": asdict(record.validity),
        "uncertainty": {
            "standard_errors": {
                "reference_isotherm": reference_errors,
                "thermal_correction": thermal_errors,
                "additional": additional_errors,
            },
            "covariance": record.parameter_covariance,
            "parameter_order": (
                None
                if record.covariance_parameters is None
                else tuple(
                    _qualified_parameter_name(record, name)
                    for name in record.covariance_parameters
                )
            ),
            "error_confidence": record.parameter_error_confidence,
        },
        "notes": record.notes,
    }


def _eos_records_to_document(
    records: Iterable[EOSRecord],
    material_identifier: str,
) -> dict[str, Any]:
    """Serialize compatible EOS records as one version-2 material document.

    Every record must describe the same material, phase, cell contents, and
    public volume convention. Equation components remain record-local because
    published scales may use different reference isotherms or thermal models.
    """
    records = tuple(records)
    if not records:
        raise MaterialError("At least one EOS record is required")
    if not all(isinstance(record, EOSRecord) for record in records):
        raise ConfigurationError("records must contain only EOSRecord objects")
    first = records[0]
    material_key = (
        first.material,
        first.phase,
        first.cell_contents,
        first.volume_unit,
    )
    for record in records[1:]:
        if (
            record.material,
            record.phase,
            record.cell_contents,
            record.volume_unit,
        ) != material_key:
            raise MaterialError(
                "All records in one material document must use the same material, "
                "phase, cell contents, and public volume unit"
            )
    return {
        "format": "peritheos.material-snapshot",
        "format_version": 2,
        "material": {
            "identifier": material_identifier,
            "formula": first.material,
            "phase": first.phase,
            "cell_contents": first.cell_contents,
        },
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": first.volume_unit,
        },
        "eos_records": [_eos_record_dict(record) for record in records],
    }


def _load_component(
    component: Mapping[str, Any], *, rt_eos: EosBase | None = None
) -> EosBase:
    """Construct a whitelisted EOS component from a version-2 record."""
    try:
        model_identifier = component["model"]
        parameters = dict(component["parameters"])
        configuration = dict(component.get("configuration", {}))
    except (KeyError, TypeError, ValueError) as error:
        raise MaterialError("Invalid equation component") from error
    try:
        model_class = _MODEL_CLASSES[model_identifier]
    except KeyError as error:
        raise MaterialError(f"Unknown equation model {model_identifier!r}") from error
    if rt_eos is not None:
        parameters["rt_eos"] = rt_eos
    parameters.update(configuration)
    try:
        return model_class(**parameters)
    except (TypeError, ValueError) as error:
        raise MaterialError(
            f"Invalid parameters for equation model {model_identifier!r}"
        ) from error


def _flatten_component_mapping(
    values: Mapping[str, Any], *, thermal: bool
) -> dict[str, Any]:
    """Convert component-separated metadata back to EOS dotted names."""
    try:
        reference = dict(values.get("reference_isotherm", {}))
        thermal_values = dict(values.get("thermal_correction", {}))
        additional = dict(values.get("additional", {}))
    except (TypeError, ValueError) as error:
        raise MaterialError("Invalid component parameter metadata") from error
    flattened = {
        (f"rt_eos.{name}" if thermal else name): value
        for name, value in reference.items()
    }
    flattened.update(thermal_values)
    flattened.update(additional)
    return flattened


def _unqualify_parameter_name(name: str, *, thermal: bool) -> str:
    """Convert a version-2 covariance path to an EOS uncertainty name."""
    reference_prefix = "reference_isotherm."
    thermal_prefix = "thermal_correction."
    if name.startswith(reference_prefix):
        parameter = name.removeprefix(reference_prefix)
        return f"rt_eos.{parameter}" if thermal else parameter
    if thermal and name.startswith(thermal_prefix):
        return name.removeprefix(thermal_prefix)
    raise MaterialError(f"Invalid covariance parameter path {name!r}")


def _eosmat_reference(value: Any) -> LiteratureReference:
    if isinstance(value, str):
        return LiteratureReference(
            authors=value,
            year=0,
            title=value,
            doi="",
            locations=("legacy eosmat reference string",),
        )
    if not isinstance(value, Mapping):
        raise MaterialError("reference must be a string or object")
    authors_value = value.get("authors", ())
    if not isinstance(authors_value, list):
        raise MaterialError("reference.authors must be an array")
    locations_value = value.get("locations")
    if locations_value is None:
        locations_value = tuple(
            str(value[key]) for key in ("locator", "details") if value.get(key)
        )
    return LiteratureReference(
        authors=", ".join(str(author) for author in authors_value),
        year=int(value["year"]),
        title=str(value.get("title") or value.get("source") or "Untitled reference"),
        doi=str(value.get("doi", "")),
        locations=tuple(str(location) for location in locations_value),
    )


def _eosmat_validity(record: Mapping[str, Any]) -> ValidityRange:
    validity = record.get("validity")
    if isinstance(validity, Mapping):
        pressure = validity.get("pressure_gpa")
        temperature = validity.get("temperature_k")
        volume_ratio = validity.get("volume_ratio")
        notes = validity.get("notes", ())
    else:
        pressure = record.get("experimental_pressure_range_gpa")
        temperature = record.get("experimental_temperature_range_k")
        volume_ratio = None
        notes = ()
    # Legacy migrated documents often do not report a numerical range. This
    # is only reachable when the caller has explicitly accepted their pending
    # scientific-validation status.
    pressure = (-np.inf, np.inf) if pressure is None else tuple(pressure)
    temperature = (
        (np.nextafter(0.0, 1.0), np.inf) if temperature is None else tuple(temperature)
    )
    return ValidityRange(
        pressure_gpa=(float(pressure[0]), float(pressure[1])),
        temperature_k=(float(temperature[0]), float(temperature[1])),
        volume_ratio=(
            None
            if volume_ratio is None
            else (float(volume_ratio[0]), float(volume_ratio[1]))
        ),
        notes=tuple(str(note) for note in notes),
    )


def _eosmat_provenance(value: Any, *, thermal: bool) -> Mapping[str, str]:
    if not isinstance(value, Mapping):
        return MappingProxyType({})
    if any(
        key in value
        for key in ("reference_isotherm", "thermal_correction", "additional")
    ):
        flattened = _flatten_component_mapping(value, thermal=thermal)
    else:
        flattened = dict(value)
    return MappingProxyType({str(key): str(item) for key, item in flattened.items()})


def _material_from_eosmat(
    document: Mapping[str, Any],
    *,
    require_primary_validation: bool,
    record_identifiers: Iterable[str] | None = None,
) -> Material:
    from peritheos.eosmat import validate_eosmat_document

    validate_eosmat_document(document)
    records_data = document.get("eos_records")
    if not isinstance(records_data, list) or not records_data:
        raise MaterialError("An executable Material requires at least one EOS record")
    if record_identifiers is not None:
        selected = frozenset(str(identifier) for identifier in record_identifiers)
        available = {str(record.get("identifier")) for record in records_data}
        missing = selected - available
        if missing:
            raise MaterialLookupError(
                f"Unknown EOS record identifiers {sorted(missing)}; available: "
                f"{sorted(available)}"
            )
        records_data = [
            record
            for record in records_data
            if str(record.get("identifier")) in selected
        ]
        if not records_data:
            raise MaterialError(
                "record_identifiers must select at least one EOS record"
            )
    units = document.get("units", {})
    volume_unit = str(units.get("volume", "angstrom^3/conventional_unit_cell"))
    formula = str(document["formula"])
    phase = str(document.get("phase") or document.get("symmetry") or "unspecified")
    formula_units = document.get("formula_units_per_cell")
    cell_contents = str(
        document.get("cell_contents")
        or (
            f"{formula_units:g} formula units per conventional unit cell"
            if formula_units is not None
            else "unspecified conventional unit cell"
        )
    )

    eos_records: list[EOSRecord] = []
    for raw_record in records_data:
        try:
            validation = raw_record["scientific_validation"]
            validation_status = str(validation["status"])
            if (
                require_primary_validation
                and validation_status != "primary_source_validated"
            ):
                raise MaterialError(
                    f"record {raw_record.get('identifier')!r} is "
                    f"{validation_status!r}; pass require_primary_validation=False "
                    "only to accept unaudited parameters explicitly"
                )
            thermal_data = raw_record.get("thermal")
            thermal_model = (
                None if thermal_data is None else str(thermal_data.get("model"))
            )
            volume_data = raw_record.get("volume", {})
            if (
                isinstance(volume_data, Mapping)
                and "public_to_model_scale" in volume_data
            ):
                volume_scale = float(volume_data["public_to_model_scale"])
            elif thermal_model in _MOLAR_VOLUME_THERMAL_MODELS:
                if formula_units is None:
                    raise MaterialError(
                        "formula_units_per_cell or volume.public_to_model_scale is "
                        "required for molar-energy thermal EOS"
                    )
                volume_scale = float(Avogadro * 1.0e-25 / float(formula_units))
            else:
                volume_scale = 1.0

            eos_data = _plain_data(raw_record["eos"])
            if volume_scale != 1.0:
                eos_data["parameters"]["V0"] *= volume_scale
            reference_eos = _load_component(eos_data)
            if thermal_data is None:
                eos = reference_eos
            else:
                thermal_component = _plain_data(thermal_data)
                configuration = dict(thermal_component.get("configuration", {}))
                for name in (
                    "debye_temperature_law",
                    "thermal_expansion_law",
                    "reference_volume_law",
                ):
                    if name in thermal_component:
                        configuration[name] = thermal_component[name]
                thermal_component["configuration"] = configuration
                eos = _load_component(thermal_component, rt_eos=reference_eos)

            is_thermal = isinstance(eos, ThermalEOS)
            reference_errors = dict(raw_record.get("parameter_errors", {}))
            errors: dict[str, float] = {
                (f"rt_eos.{name}" if is_thermal else str(name)): float(value)
                for name, value in reference_errors.items()
                if value is not None
            }
            if isinstance(thermal_data, Mapping):
                errors.update(
                    {
                        str(name): float(value)
                        for name, value in thermal_data.get(
                            "parameter_errors", {}
                        ).items()
                        if value is not None
                    }
                )
            covariance_data = raw_record.get("parameter_covariance")
            covariance = None
            covariance_parameters = None
            if (
                isinstance(covariance_data, Mapping)
                and covariance_data.get("matrix") is not None
            ):
                covariance = tuple(
                    tuple(float(value) for value in row)
                    for row in covariance_data["matrix"]
                )
                covariance_parameters = tuple(
                    str(name) for name in covariance_data.get("parameter_order", ())
                )
            reference_temperature = float(
                raw_record.get(
                    "temperature_ref",
                    getattr(eos, "Tr", 300.0),
                )
            )
            notes_value = raw_record.get("notes", "")
            notes = (
                tuple(str(note) for note in notes_value)
                if isinstance(notes_value, list)
                else ((str(notes_value),) if notes_value else ())
            )
            eos_records.append(
                EOSRecord(
                    identifier=str(raw_record["identifier"]),
                    name=str(raw_record["label"]),
                    material=formula,
                    phase=phase,
                    cell_contents=cell_contents,
                    eos=eos,
                    reference_temperature=reference_temperature,
                    reference=_eosmat_reference(raw_record["reference"]),
                    validity=_eosmat_validity(raw_record),
                    parameter_provenance=_eosmat_provenance(
                        raw_record.get("parameter_provenance", {}),
                        thermal=is_thermal,
                    ),
                    parameter_errors=(None if not errors else MappingProxyType(errors)),
                    parameter_error_confidence=(
                        None
                        if raw_record.get("parameter_error_confidence") is None
                        else float(raw_record["parameter_error_confidence"])
                    ),
                    parameter_covariance=covariance,
                    covariance_parameters=covariance_parameters,
                    notes=notes,
                    volume_unit=volume_unit,
                    volume_scale=volume_scale,
                    scientific_validation_status=validation_status,
                    scientific_validation_note=str(validation.get("note", "")),
                    eosmat_metadata=_plain_data(raw_record),
                    aliases=tuple(raw_record.get("aliases", ())),
                )
            )
        except (KeyError, TypeError, ValueError) as error:
            identifier = (
                raw_record.get("identifier", "<unknown>")
                if isinstance(raw_record, Mapping)
                else "<unknown>"
            )
            raise MaterialError(
                f"Invalid EOS record {identifier!r}: {error}"
            ) from error

    known_keys = {
        "format",
        "format_version",
        "identifier",
        "name",
        "formula",
        "phase",
        "cell_contents",
        "units",
        "aliases",
        "symmetry",
        "lattice",
        "formula_units_per_cell",
        "space_group",
        "space_group_number",
        "atom_sites",
        "peaks",
        "eos_records",
    }
    return Material(
        identifier=str(document["identifier"]),
        formula=formula,
        phase=phase,
        cell_contents=cell_contents,
        eos_records=tuple(eos_records),
        volume_unit=volume_unit,
        name=str(document["name"]),
        symmetry=(
            None if document.get("symmetry") is None else str(document["symmetry"])
        ),
        lattice=(None if document.get("lattice") is None else document["lattice"]),
        formula_units_per_cell=(
            None if formula_units is None else float(formula_units)
        ),
        space_group=(
            None
            if document.get("space_group") is None
            else str(document["space_group"])
        ),
        space_group_number=(
            None
            if document.get("space_group_number") is None
            else int(document["space_group_number"])
        ),
        atom_sites=tuple(document.get("atom_sites", ())),
        peaks=tuple(tuple(peak) for peak in document.get("peaks", ())),
        eosmat_metadata={
            key: _plain_data(value)
            for key, value in document.items()
            if key not in known_keys
        },
        aliases=tuple(document.get("aliases", ())),
    )


def material_from_dict(document: Mapping[str, Any]) -> Material:
    """Load canonical eosmat format 3 or deprecated executable snapshot 2.

    Model identifiers are resolved through Peritheos's fixed registry. The
    informational ``implementation`` string is never imported or executed.
    """
    if (
        document.get("format") == "peritheos.material"
        or document.get("format_version") == 3
    ):
        return _material_from_eosmat(
            document,
            require_primary_validation=True,
        )
    if document.get("format") != "peritheos.material-snapshot":
        raise MaterialError("Not a Peritheos executable material snapshot")
    if document.get("format_version") != 2:
        raise MaterialError("Only material format version 2 is supported")
    try:
        material = document["material"]
        units = document["units"]
        records = document["eos_records"]
        material_identifier = str(material["identifier"])
        material_formula = str(material["formula"])
        phase = str(material["phase"])
        cell_contents = str(material["cell_contents"])
        volume_unit = str(units["volume"])
    except (KeyError, TypeError) as error:
        raise MaterialError("Invalid material document header") from error
    if units.get("pressure") != "GPa" or units.get("temperature") != "K":
        raise MaterialError("Pressure and temperature units must be GPa and K")
    if not isinstance(records, list) or not records:
        raise MaterialError("eos_records must be a non-empty list")

    eos_records = []
    for record in records:
        try:
            equation = record["equation"]
            combination = equation["combination"]
            if (
                combination.get("status") != "source_parameterization"
                or combination.get("validated_as_composed") is not True
            ):
                raise MaterialError(
                    "EOSRecord loading requires a validated source composition"
                )
            reference_eos = _load_component(equation["reference_isotherm"])
            thermal_component = equation.get("thermal_correction")
            eos = (
                reference_eos
                if thermal_component is None
                else _load_component(thermal_component, rt_eos=reference_eos)
            )
            is_thermal = isinstance(eos, ThermalEOS)
            reference_data = record["provenance"]["primary_reference"]
            reference = LiteratureReference(
                authors=str(reference_data["authors"]),
                year=int(reference_data["year"]),
                title=str(reference_data["title"]),
                doi=str(reference_data["doi"]),
                locations=tuple(reference_data["locations"]),
            )
            validity_data = record["validity"]
            validity = ValidityRange(
                pressure_gpa=tuple(validity_data["pressure_gpa"]),
                temperature_k=tuple(validity_data["temperature_k"]),
                volume_ratio=(
                    None
                    if validity_data.get("volume_ratio") is None
                    else tuple(validity_data["volume_ratio"])
                ),
                notes=tuple(validity_data.get("notes", ())),
            )
            provenance = _flatten_component_mapping(
                record["provenance"]["parameters"], thermal=is_thermal
            )
            uncertainty = record["uncertainty"]
            errors = _flatten_component_mapping(
                uncertainty["standard_errors"], thermal=is_thermal
            )
            covariance_data = uncertainty.get("covariance")
            covariance = (
                None
                if covariance_data is None
                else tuple(tuple(row) for row in covariance_data)
            )
            parameter_order = uncertainty.get("parameter_order")
            covariance_parameters = (
                None
                if parameter_order is None
                else tuple(
                    _unqualify_parameter_name(name, thermal=is_thermal)
                    for name in parameter_order
                )
            )
            volume_data = record["volume"]
            if volume_data["public_unit"] != volume_unit:
                raise MaterialError("Record and document volume units differ")
            if record["reference_temperature"].get("unit") != "K":
                raise MaterialError("Reference temperature unit must be K")
            volume_scale = float(volume_data["public_to_model_scale"])
            if not np.isfinite(volume_scale) or volume_scale <= 0.0:
                raise MaterialError("Volume scale must be finite and greater than zero")
            eos_record = EOSRecord(
                identifier=str(record["identifier"]),
                name=str(record["name"]),
                material=material_formula,
                phase=phase,
                cell_contents=cell_contents,
                eos=eos,
                reference_temperature=float(record["reference_temperature"]["value"]),
                reference=reference,
                validity=validity,
                parameter_provenance=MappingProxyType(provenance),
                parameter_errors=None if not errors else MappingProxyType(errors),
                parameter_error_confidence=(
                    None
                    if uncertainty.get("error_confidence") is None
                    else float(uncertainty["error_confidence"])
                ),
                parameter_covariance=covariance,
                covariance_parameters=covariance_parameters,
                notes=tuple(record.get("notes", ())),
                volume_unit=volume_unit,
                volume_scale=volume_scale,
            )
            if not np.isclose(
                eos_record.reference_volume,
                float(volume_data["reference_value"]),
                rtol=1.0e-12,
                atol=0.0,
            ):
                raise MaterialError(
                    "Reference volume is inconsistent with EOS parameters"
                )
        except (KeyError, TypeError, ValueError) as error:
            identifier = (
                record.get("identifier", "<unknown>")
                if isinstance(record, Mapping)
                else "<unknown>"
            )
            raise MaterialError(f"Invalid EOS record {identifier!r}") from error
        eos_records.append(eos_record)
    return Material(
        identifier=material_identifier,
        formula=material_formula,
        phase=phase,
        cell_contents=cell_contents,
        eos_records=tuple(eos_records),
        volume_unit=volume_unit,
    )


_TANGE_REFERENCE = LiteratureReference(
    authors="Y. Tange, Y. Nishihara, and T. Tsuchiya",
    year=2009,
    title=(
        "Unified analyses for P-V-T equation of state of MgO: A solution for "
        "pressure-scale problems in high P-T experiments"
    ),
    doi="10.1029/2008JB005813",
    locations=("equations 2, 4, 5, 15, and 16", "Tables 1, 2, 4, and 5"),
)

_DORFMAN_REFERENCE = LiteratureReference(
    authors="S. M. Dorfman, V. B. Prakapenka, Y. Meng, and T. S. Duffy",
    year=2012,
    title="Intercomparison of pressure standards (Au, Pt, Mo, MgO, NaCl and Ne) to 2.5 Mbar",
    doi="10.1029/2012JB009292",
    locations=("equation 2", "Tables 1 and 2", "sections 3.1, 3.2, and 4"),
)

_DEWAELE_2019_REFERENCE = LiteratureReference(
    authors="A. Dewaele",
    year=2019,
    title=(
        "Equations of State of Simple Solids (Including Pb, NaCl and LiF) "
        "Compressed in Helium or Neon in the Mbar Range"
    ),
    doi="10.3390/min9110684",
    locations=("equation 1", "Table 1", "Tables 3 and 4", "sections 3.2 and 3.3"),
)

_DEWAELE_SALTS_REFERENCE = LiteratureReference(
    authors="A. Dewaele, A. B. Belonoshko, G. Garbarino, F. Occelli, P. Bouvier, M. Hanfland, and M. Mezouar",
    year=2012,
    title=("High-pressure-high-temperature equation of state of KCl and KBr"),
    doi="10.1103/PhysRevB.85.214105",
    locations=("equation 2", "Tables I, II, III, and V", "section IV"),
)

_DATCHI_CBN_REFERENCE = LiteratureReference(
    authors="F. Datchi, A. Dewaele, Y. Le Godec, and P. Loubeyre",
    year=2007,
    title=(
        "Equation of state of cubic boron nitride at high pressures and temperatures"
    ),
    doi="10.1103/PhysRevB.75.214104",
    locations=("Table I", "Figure 1 caption", "sections II and III"),
)

_DEWAELE_DIAMOND_REFERENCE = LiteratureReference(
    authors="A. Dewaele, F. Datchi, P. Loubeyre, and M. Mezouar",
    year=2008,
    title="High pressure-high temperature equations of state of neon and diamond",
    doi="10.1103/PhysRevB.77.094106",
    locations=("equations 2, 3, and 6", "Tables I and III", "sections IV and V"),
)

_CORREA_DIAMOND_REFERENCE = LiteratureReference(
    authors=(
        "A. A. Correa, L. X. Benedict, D. A. Young, E. Schwegler, and S. A. Bonev"
    ),
    year=2008,
    title=(
        "First-principles multiphase equation of state of carbon under "
        "extreme conditions"
    ),
    doi="10.1103/PhysRevB.78.024101",
    locations=(
        "equations 2-7 and 13-18",
        "Table I, diamond row",
        "Figures 3 and 5-8",
        "section II.A",
    ),
)

_BENEDICT_DIAMOND_REFERENCE = LiteratureReference(
    authors=(
        "L. X. Benedict, K. P. Driver, S. Hamel, B. Militzer, T. Qi, "
        "A. A. Correa, and E. Schwegler"
    ),
    year=2014,
    title=(
        "A multiphase equation of state for carbon addressing high pressures "
        "and temperatures"
    ),
    doi="10.1103/PhysRevB.89.224109",
    locations=(
        "equations 3-7",
        "Table I, diamond column",
        "section III.A",
        "Figure 6",
    ),
)

_DEWAELE_METALS_REFERENCE = LiteratureReference(
    authors="A. Dewaele, M. Torrent, P. Loubeyre, and M. Mezouar",
    year=2008,
    title=(
        "Compression curves of transition metals in the Mbar range: "
        "Experiments and projector augmented-wave calculations"
    ),
    doi="10.1103/PhysRevB.78.104102",
    locations=("equation 1", "Tables I, II, and IV", "Figure 1 caption"),
)

_SOKOLOVA_2016_IMPLEMENTATION_REFERENCE = LiteratureReference(
    authors=(
        "T. S. Sokolova, P. I. Dorogokupets, A. M. Dymshits, "
        "B. S. Danilov, and K. D. Litasov"
    ),
    year=2016,
    title=(
        "Microsoft Excel spreadsheets for calculation of P-V-T relations and "
        "thermodynamic properties from equations of state of MgO, diamond and "
        "nine metals as pressure markers in high-pressure and high-temperature "
        "experiments"
    ),
    doi="10.1016/j.cageo.2016.06.002",
    locations=(
        "equations 1-12",
        "Table 1",
        "Tables 2 and 3",
        "sections 3 and 4",
        "Appendix A",
    ),
)

_SOKOLOVA_2013_REFERENCE = LiteratureReference(
    authors="T. S. Sokolova, P. I. Dorogokupets, and K. D. Litasov",
    year=2013,
    title=(
        "Self-consistent pressure scales based on the equations of state for "
        "ruby, diamond, MgO, B2-NaCl, as well as Au, Pt, and other metals to "
        "4 Mbar and 3000 K"
    ),
    doi="10.1016/j.rgg.2013.01.005",
    locations=(
        "equations 1-15",
        "Table 1 input reference values",
        "Table 4 optimized Holzapfel parameters",
        "optimized equations of state under quasi-hydrostatic conditions",
        "Figure 3",
    ),
)

_FEI_REFERENCE = LiteratureReference(
    authors=("Y. Fei, A. Ricolleau, M. Frank, K. Mibe, G. Shen, and V. Prakapenka"),
    year=2007,
    title="Toward an internally consistent pressure scale",
    doi="10.1073/pnas.0609013104",
    locations=("equations 2 and 3", "Table 1", "Figures 1-5"),
)

_ANZELLINI_RE_REFERENCE = LiteratureReference(
    authors=("S. Anzellini, A. Dewaele, F. Occelli, P. Loubeyre, and M. Mezouar"),
    year=2014,
    title=(
        "Equation of state of rhenium and application for ultra high pressure "
        "calibration"
    ),
    doi="10.1063/1.4863300",
    locations=("equation 6", "Tables I, III, and IV", "Figure 3"),
)


def _molar_scale(formula_units_per_cell: int) -> float:
    """Convert A^3/cell to J bar^-1 mol^-1 of formula units."""
    return float(Avogadro * 1.0e-25 / formula_units_per_cell)


_MGO_SCALE = _molar_scale(4)
_MGO_V0_INTERNAL = 74.698 * _MGO_SCALE
_MGO_TANGE_EOS = Tange2009Debye(
    Vinet(_MGO_V0_INTERNAL, 160.63, 4.367),
    Tr=300.0,
    theta0=761.0,
    gamma0=1.442,
    a=0.138,
    b=5.4,
    n=2.0,
)

MGO_TANGE_2009 = EOSRecord(
    identifier="mgo_b1_tange_2009_vinet",
    name="MgO B1 (Tange 2009 Fit3-Vinet)",
    material="MgO",
    phase="B1 (periclase/rock-salt), cubic Fm-3m",
    cell_contents="4 MgO formula units per conventional cubic cell",
    eos=_MGO_TANGE_EOS,
    reference_temperature=300.0,
    reference=_TANGE_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.0, 196.0),
        temperature_k=(300.0, 3700.0),
        volume_ratio=(0.652, 1.150),
        notes=(
            "Marginal envelope of pressure-scale-free data in Table 1; it is not a rectangular joint domain.",
            "Table 5 extends calculations to 4000 K, but values outside the data envelope are extrapolations.",
        ),
    ),
    parameter_provenance=MappingProxyType(
        {
            "rt_eos.V0": "Table 2 (fixed); 74.698(16) A^3/conventional cell",
            "rt_eos.K0": "Table 4, Fit3-Vinet; K_T0 = 160.63(18) GPa",
            "rt_eos.K0_prime": "Table 4, Fit3-Vinet; 4.367(13)",
            "theta0": "Table 4, Fit3-Vinet; 761(13) K",
            "gamma0": "Table 4, Fit3-Vinet; 1.442(15)",
            "a": "Table 4, Fit3-Vinet; 0.138(19)",
            "b": "Table 4, Fit3-Vinet; 5.4(11)",
            "Tr": "300 K reference isotherm throughout sections 2.2 and 3",
            "n": "Equation 5; n=2 atoms per MgO formula unit",
        }
    ),
    parameter_errors=MappingProxyType(
        {
            "rt_eos.V0": 0.016 * _MGO_SCALE,
            "rt_eos.K0": 0.18,
            "rt_eos.K0_prime": 0.013,
            "theta0": 13.0,
            "gamma0": 0.015,
            "a": 0.019,
            "b": 1.1,
        }
    ),
    notes=(
        "Fit3-Vinet is one of two models the authors retained; Fit3-BM3 was not statistically rejected.",
        "The 2010 correction (doi:10.1029/2010JB007959) changes Figure 11 only, not the EOS equations, parameters, or Table 5 values used here.",
        "No parameter covariance is published; reported standard errors are propagated as independent.",
        "The paper reports 0.8 GPa total RMS pressure residual, which is model scatter and is not added to parameter uncertainty.",
    ),
    volume_scale=_MGO_SCALE,
)


def _sokolova_eos_record(
    identifier: str,
    material: str,
    phase: str,
    cell_contents: str,
    cell_formula_units: int,
    V0: float,
    K0: float,
    K0_prime: float,
    theta1: float,
    multiplicity1: float,
    theta2: float,
    multiplicity2: float,
    delta: float,
    t: float,
    a_0: float | None,
    m: float | None,
    e_0: float | None,
    g: float | None,
    n: float,
    Z: float,
    z_provenance: str,
    *,
    scientific_reference: LiteratureReference = _SOKOLOVA_2013_REFERENCE,
    source_short_name: str = "Sokolova 2013; 2016 workbook",
    reference_value_table: str = "Sokolova et al. (2013), Table 1",
    optimized_parameter_table: str = "Sokolova et al. (2013), Table 4",
    earlier_fit_note: str = (
        "Dorogokupets et al. (2012), doi:10.5800/GT-2012-3-2-0067, "
        "publishes the preceding simultaneous optimization for diamond and "
        "the nine metals."
    ),
    a0_provenance: str | None = None,
    lineage_notes: tuple[str, ...] = (),
) -> EOSRecord:
    scale = _molar_scale(cell_formula_units)
    anharmonicity = 0.0 if a_0 is None else a_0
    anharmonic_exponent = 0.0 if m is None else m
    electronic = 0.0 if e_0 is None else e_0
    electronic_exponent = 0.0 if g is None else g
    eos = MultiOscillatorGruneisenThermalEOS(
        Holzapfel(V0, K0, K0_prime, n, Z),
        Tr=298.15,
        QE1o=theta1,
        mE1=multiplicity1,
        QE2o=theta2,
        mE2=multiplicity2,
        delta=delta,
        t=t,
        a_0=anharmonicity,
        m=anharmonic_exponent,
        g=electronic_exponent,
        e_0=electronic,
        n=n,
    )

    def table_value(value: float | None, unit: str = "") -> str:
        if value is None:
            return f"{optimized_parameter_table} dash; inactive term encoded as zero"
        suffix = f" {unit}" if unit else ""
        return f"{optimized_parameter_table}; {value:g}{suffix}"

    return EOSRecord(
        identifier=identifier,
        name=f"{material} ({phase}; {source_short_name})",
        material=material,
        phase=phase,
        cell_contents=cell_contents,
        eos=eos,
        reference_temperature=298.15,
        reference=scientific_reference,
        validity=ValidityRange(
            pressure_gpa=(0.0, 400.0),
            temperature_k=(298.15, 3000.0),
            notes=(
                "Sokolova et al. (2016), section 4, states the calculator range is at least 0-400 GPa and 298.15-3000 K; this is not a uniform experimental-data rectangle.",
                "Equation 11 is discussed through compression V/V0 = 0.6, but no joint per-material volume limit is published.",
            ),
        ),
        parameter_provenance=MappingProxyType(
            {
                "rt_eos.V0": f"{reference_value_table}; {V0:g} J bar^-1 mol^-1",
                "rt_eos.K0": f"{optimized_parameter_table}; {10.0 * K0:g} kbar, converted to {K0:g} GPa",
                "rt_eos.K0_prime": f"{optimized_parameter_table}; {K0_prime:g}",
                "rt_eos.n": f"2013 Table 1 and 2016 Table 2 definition; {n:g} atoms per chemical formula",
                "rt_eos.Z": f"Sokolova et al. (2016), {z_provenance}",
                "Tr": "Sokolova et al. (2016), Table 2 and section 2.1; T0 = 298.15 K",
                "QE1o": f"{optimized_parameter_table} theta_01; {theta1:g} K",
                "mE1": f"{optimized_parameter_table} m1; {multiplicity1:g}",
                "QE2o": f"{optimized_parameter_table} theta_02; {theta2:g} K",
                "mE2": f"{optimized_parameter_table} m2; {multiplicity2:g}",
                "delta": f"{optimized_parameter_table}; {delta:g}",
                "t": f"{optimized_parameter_table}; {t:g}",
                "a_0": a0_provenance or table_value(a_0, "10^-6 K^-1"),
                "m": table_value(m),
                "e_0": table_value(e_0, "10^-6 K^-1"),
                "g": table_value(g),
                "beta": "Not present in equations 7-12 or Table 1; disabled (zero)",
                "QBo": "No first generalized Bose mode in equations 7-12; unused default",
                "d": "No first generalized Bose mode in equations 7-12; unused default",
                "mb": "No first generalized Bose mode in equations 7-12; multiplicity disabled (zero)",
                "QB1o": "No second generalized Bose mode in equations 7-12; unused default",
                "d1": "No second generalized Bose mode in equations 7-12; unused default",
                "mb1": "No second generalized Bose mode in equations 7-12; multiplicity disabled (zero)",
                "n": f"Equation 12; {n:g} atoms per chemical formula",
            }
        ),
        notes=(
            "The scientific coefficient source is Sokolova et al. (2013): Table 1 gives the reference volume and composition inputs, and Table 4 gives the final cross-calibrated Holzapfel and thermal parameters.",
            "The 2013 optimization combines shock-wave, ultrasonic, X-ray diffraction, dilatometric, and thermochemical measurements; it is a self-consistent pressure-scale fit, not one experimental dataset.",
            earlier_fit_note,
            "Sokolova et al. (2016) republishes the coefficients and supplies the executable Excel/VBA calculator, reference-temperature convention, and corrected implementation equations; it is an implementation source rather than a new fit dataset.",
            "The source tables supply no individual parameter errors or covariance, and the complete point-by-point fit inputs and weights are not published; only measurement-state uncertainty can be propagated.",
            "Sokolova et al. (2016), section 4, estimates marker uncertainty at no more than 3-4% above 200 GPa and near 3000 K, without a confidence convention; it is recorded here but not treated as one-sigma parameter uncertainty.",
            "The *_sokolova_2013 identifier names the scientific fit year; the 2016 workbook remains explicit in the implementation lineage.",
            *lineage_notes,
        ),
        volume_scale=scale,
    )


MGO_SOKOLOVA_2013 = _sokolova_eos_record(
    "mgo_b1_sokolova_2013",
    "MgO",
    "B1 (periclase/rock-salt), cubic Fm-3m",
    "4 MgO formula units per conventional cubic cell",
    4,
    1.1248,
    160.3,
    4.10,
    748.0,
    3.0,
    401.0,
    3.0,
    -0.235,
    0.301,
    -17.4,
    4.95,
    None,
    None,
    2.0,
    10.34,
    "Equation 3 and Figure 1 spreadsheet input; effective atomic number 10.34 for MgO",
    earlier_fit_note=(
        "The MgO lineage runs through Dorogokupets (2010), "
        "doi:10.1007/s00269-010-0367-2, and the earlier Dorogokupets-Oganov "
        "formalism before its joint cross-calibration in 2013."
    ),
    a0_provenance=(
        "Sokolova et al. (2016), Table 1; -17.4 10^-6 K^-1. The 2016 paper "
        "corrects the earlier printed MgO a0 value, so this coefficient has "
        "explicit dual 2013/2016 provenance."
    ),
)

DIAMOND_SOKOLOVA_2013 = _sokolova_eos_record(
    "diamond_sokolova_2013",
    "C",
    "diamond, cubic Fd-3m",
    "8 C atoms per conventional cubic cell",
    8,
    0.3414,
    441.5,
    3.90,
    1561.0,
    2.436,
    684.0,
    0.564,
    -0.506,
    1.085,
    None,
    None,
    None,
    None,
    1.0,
    6.0,
    "Table 2 definition; atomic number Z = 6",
    lineage_notes=(
        "Temperatures above 3000 K, including 6000 K DAC calculations, are extrapolations beyond the source's stated calculation range.",
    ),
)

AL_SOKOLOVA_2013 = _sokolova_eos_record(
    "al_fcc_sokolova_2013",
    "Al",
    "fcc, cubic Fm-3m",
    "4 Al atoms per conventional cubic cell",
    4,
    0.998,
    72.8,
    4.51,
    381.0,
    1.5,
    202.0,
    1.5,
    -0.242,
    -0.958,
    None,
    None,
    64.1,
    0.33,
    1.0,
    13.0,
    "Table 2 definition; atomic number Z = 13",
)
CU_SOKOLOVA_2013 = _sokolova_eos_record(
    "cu_fcc_sokolova_2013",
    "Cu",
    "fcc, cubic Fm-3m",
    "4 Cu atoms per conventional cubic cell",
    4,
    0.7112,
    133.5,
    5.32,
    296.0,
    1.5,
    169.0,
    1.5,
    -0.07,
    1.401,
    None,
    None,
    27.7,
    2.18,
    1.0,
    29.0,
    "Table 2 definition; atomic number Z = 29",
)
AG_SOKOLOVA_2013 = _sokolova_eos_record(
    "ag_fcc_sokolova_2013",
    "Ag",
    "fcc, cubic Fm-3m",
    "4 Ag atoms per conventional cubic cell",
    4,
    1.025,
    100.0,
    6.15,
    199.0,
    1.5,
    115.0,
    1.5,
    0.178,
    2.210,
    None,
    None,
    22.1,
    0.19,
    1.0,
    47.0,
    "Table 2 definition; atomic number Z = 47",
)
AU_SOKOLOVA_2013 = _sokolova_eos_record(
    "au_fcc_sokolova_2013",
    "Au",
    "fcc, cubic Fm-3m",
    "4 Au atoms per conventional cubic cell",
    4,
    1.0215,
    167.0,
    5.90,
    179.5,
    1.5,
    83.0,
    1.5,
    0.134,
    0.087,
    None,
    None,
    None,
    None,
    1.0,
    79.0,
    "Table 2 definition; atomic number Z = 79",
)
PT_SOKOLOVA_2013 = _sokolova_eos_record(
    "pt_fcc_sokolova_2013",
    "Pt",
    "fcc, cubic Fm-3m",
    "4 Pt atoms per conventional cubic cell",
    4,
    0.9091,
    275.0,
    5.35,
    177.0,
    1.5,
    143.0,
    1.5,
    0.167,
    -0.343,
    None,
    None,
    80.6,
    0.06,
    1.0,
    78.0,
    "Table 2 definition; atomic number Z = 78",
)
NB_SOKOLOVA_2013 = _sokolova_eos_record(
    "nb_bcc_sokolova_2013",
    "Nb",
    "bcc, cubic Im-3m",
    "2 Nb atoms per conventional cubic cell",
    2,
    1.0828,
    170.5,
    3.65,
    302.0,
    1.5,
    134.0,
    1.5,
    -0.326,
    -0.763,
    None,
    None,
    115.9,
    0.90,
    1.0,
    41.0,
    "Table 2 definition; atomic number Z = 41",
)
TA_SOKOLOVA_2013 = _sokolova_eos_record(
    "ta_bcc_sokolova_2013",
    "Ta",
    "bcc, cubic Im-3m",
    "2 Ta atoms per conventional cubic cell",
    2,
    1.0861,
    191.0,
    3.83,
    254.0,
    1.5,
    101.0,
    1.5,
    -0.101,
    -0.148,
    None,
    None,
    82.3,
    0.12,
    1.0,
    73.0,
    "Table 2 definition; atomic number Z = 73",
)
MO_SOKOLOVA_2013 = _sokolova_eos_record(
    "mo_bcc_sokolova_2013",
    "Mo",
    "bcc, cubic Im-3m",
    "2 Mo atoms per conventional cubic cell",
    2,
    0.9369,
    260.0,
    4.20,
    353.0,
    1.5,
    222.0,
    1.5,
    -0.802,
    -0.791,
    None,
    None,
    143.2,
    2.66,
    1.0,
    42.0,
    "Table 2 definition; atomic number Z = 42",
)
W_SOKOLOVA_2013 = _sokolova_eos_record(
    "w_bcc_sokolova_2013",
    "W",
    "bcc, cubic Im-3m",
    "2 W atoms per conventional cubic cell",
    2,
    0.9552,
    308.0,
    4.12,
    309.0,
    1.5,
    172.0,
    1.5,
    -0.686,
    -0.591,
    None,
    None,
    100.1,
    2.77,
    1.0,
    74.0,
    "Table 2 definition; atomic number Z = 74",
)

SOKOLOVA_2013_EOS_RECORDS = (
    MGO_SOKOLOVA_2013,
    DIAMOND_SOKOLOVA_2013,
    AL_SOKOLOVA_2013,
    CU_SOKOLOVA_2013,
    AG_SOKOLOVA_2013,
    AU_SOKOLOVA_2013,
    PT_SOKOLOVA_2013,
    NB_SOKOLOVA_2013,
    TA_SOKOLOVA_2013,
    MO_SOKOLOVA_2013,
    W_SOKOLOVA_2013,
)


def _fei_eos_record(
    identifier: str,
    material: str,
    phase: str,
    cell_contents: str,
    cell_formula_units: int,
    V0: float,
    K0: float,
    K0_prime: float,
    theta0: float,
    gamma0: float,
    q: float,
    n: float,
    pressure_range: tuple[float, float],
    temperature_range: tuple[float, float],
    parameter_errors: Mapping[str, float],
    provenance_note: str,
) -> EOSRecord:
    scale = _molar_scale(cell_formula_units)
    eos = MieGruneisenDebye(
        Vinet(V0 * scale, K0, K0_prime),
        Tr=300.0,
        theta0=theta0,
        gamma0=gamma0,
        q=q,
        n=n,
        debye_temperature_law="variable_exponent",
    )
    scaled_errors = {
        name: (error * scale if name == "rt_eos.V0" else error)
        for name, error in parameter_errors.items()
    }
    return EOSRecord(
        identifier=identifier,
        name=f"{material} ({phase}; Fei 2007)",
        material=material,
        phase=phase,
        cell_contents=cell_contents,
        eos=eos,
        reference_temperature=300.0,
        reference=_FEI_REFERENCE,
        validity=ValidityRange(
            pressure_gpa=pressure_range,
            temperature_k=temperature_range,
            notes=(
                "Marginal envelope of the measurements and fitted curves in Figures 1-5; not every pressure-temperature combination was measured.",
            ),
        ),
        parameter_provenance=MappingProxyType(
            {
                "rt_eos.V0": f"Table 1; {V0:g} A^3/conventional cell",
                "rt_eos.K0": f"Table 1 Vinet fit; {K0:g} GPa",
                "rt_eos.K0_prime": f"Table 1 Vinet fit; {K0_prime:g}",
                "Tr": "Equation 3; 300 K reference isotherm",
                "theta0": f"Table 1; {theta0:g} K",
                "gamma0": f"Table 1; {gamma0:g}",
                "q": f"Table 1; {q:g}",
                "n": f"Equation 3 Debye energy; {n:g} atoms per formula unit",
                "source_fit": provenance_note,
            }
        ),
        parameter_errors=MappingProxyType(scaled_errors),
        notes=(
            "The Fei convention theta_D = theta0*(V/V0)^(-gamma(V)) is implemented explicitly; it is not replaced by the integrated constant-q Debye-temperature relation.",
            "Table 1 reports parenthetical/plus-minus parameter uncertainties without a confidence convention; they are propagated as one-standard-deviation errors and this assumption must be retained when reporting results.",
            "No parameter covariance is published; reported errors are propagated independently. Fixed/adopted parameters without errors are omitted rather than assigned zero uncertainty.",
        ),
        volume_scale=scale,
    )


AU_FEI_2007 = _fei_eos_record(
    "au_fcc_fei_2007",
    "Au",
    "fcc, cubic Fm-3m",
    "4 Au atoms per conventional cubic cell",
    4,
    67.850,
    167.0,
    6.00,
    170.0,
    2.97,
    0.6,
    1.0,
    (0.0, 125.0),
    (300.0, 2330.0),
    {"rt_eos.V0": 0.004, "rt_eos.K0_prime": 0.02, "gamma0": 0.03, "q": 0.3},
    "Table 1 footnote: optimized from Dewaele et al. (2004), Fei et al. (2004), and this study",
)
PT_FEI_2007 = _fei_eos_record(
    "pt_fcc_fei_2007",
    "Pt",
    "fcc, cubic Fm-3m",
    "4 Pt atoms per conventional cubic cell",
    4,
    60.38,
    277.0,
    5.08,
    230.0,
    2.72,
    0.5,
    1.0,
    (0.0, 94.0),
    (300.0, 1873.0),
    {"rt_eos.V0": 0.01, "rt_eos.K0_prime": 0.02, "gamma0": 0.03, "q": 0.5},
    "Table 1 footnote: optimized from Dewaele et al. (2004) and Fei et al. (2004)",
)
NACL_B2_FEI_2007 = _fei_eos_record(
    "nacl_b2_fei_2007",
    "NaCl",
    "B2 (CsCl-type), cubic Pm-3m",
    "1 NaCl formula unit per conventional cubic cell",
    1,
    41.35,
    26.86,
    5.25,
    290.0,
    1.70,
    0.5,
    2.0,
    (34.0, 107.0),
    (300.0, 1000.0),
    {"rt_eos.K0": 2.90, "rt_eos.K0_prime": 0.26, "q": 0.3},
    "This study; Table 1 footnote says V0, theta0, and gamma0 are adopted from Bukowinski and Aidun (1985)",
)
NE_FEI_2007 = _fei_eos_record(
    "ne_fcc_fei_2007",
    "Ne",
    "fcc, cubic Fm-3m",
    "4 Ne atoms per conventional cubic cell",
    4,
    88.967,
    1.16,
    8.23,
    75.1,
    2.05,
    0.6,
    1.0,
    (5.0, 115.0),
    (300.0, 1000.0),
    {"rt_eos.K0": 0.14, "rt_eos.K0_prime": 0.31, "q": 0.3},
    "This study; Table 1 footnote says V0, theta0, and gamma0 are adopted from Finger et al. (1981)",
)

FEI_2007_EOS_RECORDS = (AU_FEI_2007, PT_FEI_2007, NACL_B2_FEI_2007, NE_FEI_2007)


def _dorfman_eos_record(
    identifier: str,
    name: str,
    material: str,
    phase: str,
    cell_contents: str,
    V0: float,
    K0: float,
    K0_prime: float,
    pressure_range: tuple[float, float],
    formal_errors: str,
) -> EOSRecord:
    return EOSRecord(
        identifier=identifier,
        name=name,
        material=material,
        phase=phase,
        cell_contents=cell_contents,
        eos=Vinet(V0, K0, K0_prime),
        reference_temperature=300.0,
        reference=_DORFMAN_REFERENCE,
        validity=ValidityRange(
            pressure_gpa=pressure_range,
            temperature_k=(300.0, 300.0),
            notes=(
                "Pressure range is the experimental run envelope in Table 1, limited to 250 GPa where the paper claims internal 3% agreement.",
                "This is a 300 K isotherm; no thermal-pressure term was determined.",
            ),
        ),
        parameter_provenance=MappingProxyType(
            {
                "V0": f"Table 2, K0-fixed fit; {V0:g} A^3/conventional cell (fixed)",
                "K0": f"Table 2, K0-fixed fit; {K0:g} GPa",
                "K0_prime": f"Table 2, K0-fixed fit; {K0_prime:g}",
            }
        ),
        parameter_errors=MappingProxyType(
            {"K0": 0.02 * K0, "K0_prime": 0.02 * K0_prime}
        ),
        notes=(
            "The simultaneous relative scale is anchored to Tange et al. (2009) MgO, not an independent absolute calibration.",
            "The authors state constrained formal errors are artificially low and recommend about 2% realistic errors for K0 and K0_prime; those 2% errors are propagated here as independent because covariance was not published.",
            formal_errors,
            "The reported 2-3% non-hydrostatic stress contribution is not folded into parameter uncertainty.",
        ),
    )


AU_DORFMAN_2012 = _dorfman_eos_record(
    "au_fcc_dorfman_2012",
    "Au fcc (Dorfman 2012)",
    "Au",
    "fcc, cubic Fm-3m",
    "4 Au atoms per conventional cubic cell",
    67.85,
    167.0,
    5.88,
    (1.0, 250.0),
    "Table 2 formal K0_prime error: 0.02; V0 and K0 were fixed.",
)

PT_DORFMAN_2012 = _dorfman_eos_record(
    "pt_fcc_dorfman_2012",
    "Pt fcc (Dorfman 2012)",
    "Pt",
    "fcc, cubic Fm-3m",
    "4 Pt atoms per conventional cubic cell",
    60.38,
    277.0,
    5.43,
    (5.0, 228.0),
    "Table 2 formal K0_prime error: 0.02; V0 and K0 were fixed.",
)

MO_DORFMAN_2012 = _dorfman_eos_record(
    "mo_bcc_dorfman_2012",
    "Mo bcc (Dorfman 2012)",
    "Mo",
    "bcc, cubic Im-3m",
    "2 Mo atoms per conventional cubic cell",
    31.17,
    261.0,
    4.19,
    (43.0, 213.0),
    "Table 2 formal K0_prime error: 0.02; V0 and K0 were fixed.",
)

NACL_B2_DORFMAN_2012 = _dorfman_eos_record(
    "nacl_b2_dorfman_2012",
    "NaCl B2 (Dorfman 2012)",
    "NaCl",
    "B2 (CsCl-type), cubic Pm-3m",
    "1 NaCl formula unit per conventional cubic cell",
    41.35,
    24.2,
    5.76,
    (34.0, 250.0),
    "Table 2 gives K0 = 24.2(3) GPa and K0_prime = 5.76(3); Table 3 instead prints 5.76(4).",
)

NE_DORFMAN_2012 = _dorfman_eos_record(
    "ne_fcc_dorfman_2012",
    "Ne fcc (Dorfman 2012)",
    "Ne",
    "fcc, cubic Fm-3m",
    "4 Ne atoms per conventional cubic cell",
    88.967,
    1.04,
    8.48,
    (5.0, 250.0),
    "Table 2 formal errors: K0 = 0.02 GPa and K0_prime = 0.04; V0 was fixed.",
)


def _dewaele_2019_eos_record(
    identifier: str,
    name: str,
    material: str,
    phase: str,
    cell_contents: str,
    formula_units_per_cell: int,
    V0_per_formula_unit: float,
    K0: float,
    K0_prime: float,
    errors_95: tuple[float | None, float, float],
    pressure_range: tuple[float, float],
    fixed: tuple[str, ...] = (),
) -> EOSRecord:
    cell_v0 = formula_units_per_cell * V0_per_formula_unit
    parameter_errors = {
        name: error
        for name, error in {
            "V0": (
                None if errors_95[0] is None else formula_units_per_cell * errors_95[0]
            ),
            "K0": errors_95[1],
            "K0_prime": errors_95[2],
        }.items()
        if error is not None
    }
    return EOSRecord(
        identifier=identifier,
        name=name,
        material=material,
        phase=phase,
        cell_contents=cell_contents,
        eos=Vinet(cell_v0, K0, K0_prime),
        reference_temperature=300.0,
        reference=_DEWAELE_2019_REFERENCE,
        validity=ValidityRange(
            pressure_gpa=pressure_range,
            temperature_k=(300.0, 300.0),
            notes=(
                "Table 1 experimental fit domain; this is a 300 K isotherm.",
                "Pressures use the unified Dorogokupets ruby calibration identified as P_Dor in Table 1.",
            ),
        ),
        parameter_provenance=MappingProxyType(
            {
                "V0": f"Table 1 P_Dor column; {V0_per_formula_unit:g} A^3/formula unit"
                + (" (fixed)" if "V0" in fixed else ""),
                "K0": f"Table 1 P_Dor column; {K0:g} GPa"
                + (" (fixed)" if "K0" in fixed else ""),
                "K0_prime": f"Table 1 P_Dor column; {K0_prime:g}"
                + (" (fixed)" if "K0_prime" in fixed else ""),
            }
        ),
        parameter_errors=MappingProxyType(parameter_errors),
        parameter_error_confidence=0.95,
        notes=(
            "Table 1 states that the unprinted P_Dor errors equal the P_Mao-column 95% fit confidence intervals.",
            "Published errors are stored as 95% confidence-interval half-widths and converted to normal-equivalent standard errors only during uncertainty propagation.",
            "No covariance is published; retained parameter errors are treated as independent.",
            "The fit is calibration-dependent and is not an independent absolute pressure scale.",
        ),
    )


LIF_B1_DEWAELE_2019 = _dewaele_2019_eos_record(
    "lif_b1_dewaele_2019",
    "LiF B1 (Dewaele 2019, P_Dor)",
    "LiF",
    "B1 (rock-salt), cubic Fm-3m",
    "4 LiF formula units per conventional cubic cell",
    4,
    16.391,
    62.3,
    5.01,
    (0.030, 1.4, 0.60),
    (0.0, 109.0),
)

NACL_B1_DEWAELE_2019 = _dewaele_2019_eos_record(
    "nacl_b1_dewaele_2019",
    "NaCl B1 (Dewaele 2019, P_Dor)",
    "NaCl",
    "B1 (rock-salt), cubic Fm-3m",
    "4 NaCl formula units per conventional cubic cell",
    4,
    44.93,
    23.4,
    5.29,
    (0.12, 0.8, 0.06),
    (0.0, 35.0),
)

NACL_B2_DEWAELE_2019 = _dewaele_2019_eos_record(
    "nacl_b2_dewaele_2019",
    "NaCl B2 (Dewaele 2019, P_Dor)",
    "NaCl",
    "B2 (CsCl-type), cubic Pm-3m",
    "1 NaCl formula unit per conventional cubic cell",
    1,
    42.3,
    22.664,
    5.735,
    (None, 1.2, 0.20),
    (37.0, 155.0),
    fixed=("V0",),
)


def _dewaele_salt_eos_record(
    identifier: str,
    material: str,
    phase: str,
    cell_contents: str,
    formula_units_per_cell: int,
    V0_per_formula_unit: float,
    K0: float,
    K0_prime: float,
    pressure_range: tuple[float, float],
    *,
    alpha_KT: float | None = None,
    fixed: str,
) -> EOSRecord:
    V0 = formula_units_per_cell * V0_per_formula_unit
    rt_eos = Vinet(V0, K0, K0_prime)
    eos = (
        rt_eos
        if alpha_KT is None
        else LinearThermalPressure(rt_eos, Tr=300.0, alpha_KT=alpha_KT)
    )
    reference_temperature = 298.0 if alpha_KT is None else 300.0
    temperature_range = (
        (reference_temperature, reference_temperature)
        if alpha_KT is None
        else (300.0, 7000.0)
    )
    upper_pressure = pressure_range[1] if alpha_KT is None else 200.0
    volume_ratio = None if alpha_KT is None else (0.4, 1.0)
    return EOSRecord(
        identifier=identifier,
        name=f"{material} {phase.split()[0]} (Dewaele 2012)",
        material=material,
        phase=phase,
        cell_contents=cell_contents,
        eos=eos,
        reference_temperature=reference_temperature,
        reference=_DEWAELE_SALTS_REFERENCE,
        validity=ValidityRange(
            pressure_gpa=(pressure_range[0], upper_pressure),
            temperature_k=temperature_range,
            volume_ratio=volume_ratio,
            notes=(
                "Table III gives the 298 K experimental compression domain.",
                *(
                    ()
                    if alpha_KT is None
                    else (
                        "Equation 2 and Table V extend the B2 calibration using molecular-dynamics thermal pressure.",
                        "The published 0-200 GPa headline range includes states below the observed B1-B2 transition; this record enforces the phase-specific lower bound.",
                    )
                ),
            ),
        ),
        parameter_provenance=MappingProxyType(
            {
                "rt_eos.V0" if alpha_KT is not None else "V0": (
                    f"Table III, experimental Rydberg-Vinet fit; {V0_per_formula_unit:g} A^3/formula unit"
                    + (" (fixed)" if fixed == "V0" else "")
                ),
                "rt_eos.K0" if alpha_KT is not None else "K0": (
                    f"Table III; {K0:g} GPa"
                ),
                "rt_eos.K0_prime" if alpha_KT is not None else "K0_prime": (
                    f"Table III; {K0_prime:g}"
                    + (" (fixed to ultrasonic value)" if fixed == "K0_prime" else "")
                ),
                **(
                    {}
                    if alpha_KT is None
                    else {
                        "Tr": "Equation 2; 300 K reference isotherm",
                        "alpha_KT": f"Table V; {alpha_KT:g} GPa/K",
                    }
                ),
            }
        ),
        notes=(
            "The paper reports no fitted-parameter covariance or individual parameter errors.",
            "Measurement uncertainty can be propagated, but parameter uncertainty is therefore omitted rather than set to zero.",
            *(
                ()
                if alpha_KT is None
                else (
                    "The thermal term is computation-derived; the 298 K compression term is experimental.",
                    "The equation does not model melting or an independent high-temperature phase boundary.",
                )
            ),
        ),
    )


KCL_B1_DEWAELE_2012 = _dewaele_salt_eos_record(
    "kcl_b1_dewaele_2012",
    "KCl",
    "B1 (rock-salt), cubic Fm-3m",
    "4 KCl formula units per conventional cubic cell",
    4,
    62.36,
    17.1,
    5.5,
    (0.0, 2.6),
    fixed="K0_prime",
)

KCL_B2_DEWAELE_2012 = _dewaele_salt_eos_record(
    "kcl_b2_dewaele_2012",
    "KCl",
    "B2 (CsCl-type), cubic Pm-3m",
    "1 KCl formula unit per cubic cell",
    1,
    54.5,
    17.2,
    5.89,
    (2.6, 165.0),
    alpha_KT=0.00224,
    fixed="V0",
)

KBR_B1_DEWAELE_2012 = _dewaele_salt_eos_record(
    "kbr_b1_dewaele_2012",
    "KBr",
    "B1 (rock-salt), cubic Fm-3m",
    "4 KBr formula units per conventional cubic cell",
    4,
    71.89,
    14.2,
    5.5,
    (0.0, 2.3),
    fixed="K0_prime",
)

KBR_B2_DEWAELE_2012 = _dewaele_salt_eos_record(
    "kbr_b2_dewaele_2012",
    "KBr",
    "B2 (CsCl-type), cubic Pm-3m",
    "1 KBr formula unit per cubic cell",
    1,
    63.4,
    14.9,
    5.81,
    (2.3, 165.0),
    alpha_KT=0.00222,
    fixed="V0",
)


CBN_DATCHI_2007 = EOSRecord(
    identifier="cbn_zincblende_datchi_2007",
    name="c-BN zinc-blende (Datchi 2007)",
    material="BN",
    phase="cubic zinc-blende, F-43m",
    cell_contents="4 BN formula units (8 atoms) per conventional cubic cell",
    eos=Vinet(5.9062 * 8.0, 395.0, 3.62),
    reference_temperature=295.0,
    reference=_DATCHI_CBN_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.0, 162.0),
        temperature_k=(295.0, 295.0),
        volume_ratio=(0.77, 1.0),
        notes=("Section III and Figure 1; H2005 ruby-calibrated 295 K data.",),
    ),
    parameter_provenance=MappingProxyType(
        {
            "V0": "Section III and Figure 1 caption; 5.9062(6) A^3/atom, fixed in fit",
            "K0": "Table I, Vinet column; 395(2) GPa",
            "K0_prime": "Table I, Vinet column; 3.62(5)",
        }
    ),
    parameter_errors=MappingProxyType(
        {"V0": 0.0006 * 8.0, "K0": 2.0, "K0_prime": 0.05}
    ),
    notes=(
        "Parentheses are reported fitting standard deviations, not absolute uncertainties; the paper identifies ruby-scale accuracy as dominant.",
        "No covariance is published; errors are propagated independently.",
        "The high-temperature c-BN model is deferred because its cold-curve/zero-point convention needs separate representation.",
    ),
)


_DIAMOND_SCALE = _molar_scale(8)
_ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR = _molar_scale(1)
_EV_PER_ATOM_TO_J_PER_MOL = electron_volt * Avogadro

DIAMOND_CORREA_2008 = EOSRecord(
    identifier="diamond_correa_2008",
    name="Diamond (Correa 2008 logarithmic-moment double-Debye Helmholtz)",
    material="C",
    phase="diamond, cubic Fd-3m",
    cell_contents="8 C atoms per conventional cubic cell",
    eos=DoubleDebyeLogMomentHelmholtz(
        Vinet(
            5.785 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
            368.2,
            4.038,
        ),
        Vp=5.571 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        theta_a0=1887.8,
        a_a=-0.316 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_a=0.913,
        theta_b0=1887.8,
        a_b=0.168 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_b=0.429,
        theta_0_0=1887.8,
        a_0=0.131 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_0=0.202,
        n=1.0,
        anharmonic_a=3.8e-5,
        phi0=-155.059 * _EV_PER_ATOM_TO_J_PER_MOL,
    ),
    reference_temperature=300.0,
    reference=_CORREA_DIAMOND_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.0, 1075.0),
        temperature_k=(1.0, 10000.0),
        volume_ratio=(2.32 / 5.785, 1.0),
        notes=(
            "The lower volume is the compressed diamond phonon-DOS state in Figure 3; the cold curve itself was fitted over approximately 1-9 A^3/atom.",
            "The upper pressure is the calculated 0 K diamond-BC8 crossing. These are marginal model bounds, not a phase-stability surface at finite temperature.",
        ),
    ),
    parameter_provenance=MappingProxyType(
        {
            "rt_eos.V0": "Table I; 5.785 A^3/atom",
            "rt_eos.K0": "Table I; 368.2 GPa",
            "rt_eos.K0_prime": "Table I; 4.038",
            "Vp": "Table I caption; Vref(diamond)=5.571 A^3/atom",
            "theta_a0": "Table I; theta_A^(0)=1887.8 K",
            "a_a": "Table I; beta_A=-0.316 A^-3",
            "b_a": "Table I; alpha_A=0.913",
            "theta_b0": "Table I; theta_B^(0)=1887.8 K",
            "a_b": "Table I; beta_B=0.168 A^-3",
            "b_b": "Table I; alpha_B=0.429",
            "theta_0_0": "Table I; theta_0^(0)=1887.8 K",
            "a_0": "Table I; beta_0=0.131 A^-3",
            "b_0": "Table I; alpha_0=0.202",
            "n": "one atom per elemental-carbon formula unit",
            "anharmonic_a": "Equation 18 and Table I; 3.8e-5 K^-1/atom",
            "phi0": "Table I; -155.059 eV/atom",
        }
    ),
    notes=(
        "Complete diamond free-energy branch from equations 2-18; the double-Debye weights conserve the logarithmic phonon moment theta_0 (equation 13).",
        "Table I volumes are per atom; the public API uses A^3 per eight-atom conventional diamond cell and converts internally to J bar^-1 mol^-1 of atoms.",
        "The Vinet term is a motionless-ion DFT-GGA cold curve. The authors state that its theoretical V0 is about 3% too large after zero-point and thermal effects and may need an application-specific density shift; Peritheos preserves the published value.",
        "The constant anharmonic coefficient changes free energy, internal energy, and heat capacity but contributes no pressure. Electronic excitations are neglected for insulating diamond as in the source.",
        "Factor-of-two source ambiguity: Correa equation 18 defines F_anh=-a*T^2 with a=3.8e-5 K^-1, whereas Benedict 2014 says it carries over the same correction but writes F_anh=-alpha*T^2/2 and tabulates alpha=3.79e-5 K^-1. This record follows Correa literally; in Benedict's normalization the equivalent coefficient would be alpha=2*a.",
        "The source reports no coefficient uncertainties or covariance.",
    ),
    volume_scale=_DIAMOND_SCALE,
)

DIAMOND_BENEDICT_2014 = EOSRecord(
    identifier="diamond_benedict_2014",
    name="Diamond (Benedict 2014 double-Debye Helmholtz)",
    material="C",
    phase="diamond, cubic Fd-3m",
    cell_contents="8 C atoms per conventional cubic cell",
    eos=DoubleDebyeHelmholtz(
        Vinet(
            5.7034 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
            432.4,
            3.793,
        ),
        Vp=5.571 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        theta_a0=1887.8,
        a_a=-0.316 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_a=0.913,
        theta_b0=1887.8,
        a_b=0.168 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_b=0.429,
        theta_1_0=1887.8,
        a_1=0.0846 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_1=0.499,
        n=1.0,
        alpha0=3.79e-5,
        Ve=5.785 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        kappa=0.0,
        phi0=-9.066 * _EV_PER_ATOM_TO_J_PER_MOL,
    ),
    reference_temperature=300.0,
    reference=_BENEDICT_DIAMOND_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.0, 1000.0),
        temperature_k=(300.0, 9000.0),
        volume_ratio=(2.5 / 5.7034, 7.0 / 5.7034),
        notes=(
            "The paper constructs the diamond free-energy branch over 2.5-7.0 A^3/atom; Figure 6 directly checks DFT-MD states only over 3.0-5.6 A^3/atom and 2000-9000 K.",
            "These are marginal model bounds, not a phase-stability surface; check the diamond-BC8 boundary and melting separately, especially near 1 TPa.",
        ),
    ),
    parameter_provenance=MappingProxyType(
        {
            "rt_eos.V0": "Table I; 5.7034 A^3/atom",
            "rt_eos.K0": "Table I; 432.4 GPa",
            "rt_eos.K0_prime": "Table I; 3.793",
            "Vp": "Table I; 5.571 A^3/atom",
            "theta_a0": "Table I; 1887.8 K",
            "a_a": "Table I; -0.316 A^-3",
            "b_a": "Table I; 0.913",
            "theta_b0": "Table I; 1887.8 K",
            "a_b": "Table I; 0.168 A^-3",
            "b_b": "Table I; 0.429",
            "theta_1_0": "Table I; 1887.8 K",
            "a_1": "Table I; 0.0846 A^-3",
            "b_1": "Table I; 0.499",
            "n": "one atom per elemental-carbon formula unit",
            "alpha0": "Table I; 3.79e-5 K^-1",
            "Ve": "Table I; 5.785 A^3/atom",
            "kappa": "Table I; 0.0",
            "phi0": "Table I; -9.066 eV/atom",
        }
    ),
    notes=(
        "This is the complete equations 3-7 Helmholtz model, not a Vinet curve combined with a generic Mie-Gruneisen-Debye correction.",
        "Table I volumes are per atom. The public API uses A^3 per eight-atom conventional diamond cell and converts internally to J bar^-1 mol^-1 of atoms.",
        "The Vinet parameters describe the motionless-ion 0 K cold curve. V0=5.7034 A^3/atom is therefore not a 300 K, zero-total-pressure reference volume; zero-point and thermal pressure remain present.",
        "Factor-of-two source ambiguity: Benedict writes F_anh=-alpha*T^2/2 and tabulates alpha=3.79e-5 K^-1, while the Correa 2008 correction it claims to retain is F_anh=-a*T^2 with a=3.8e-5 K^-1. This record follows Benedict literally; matching Correa's anharmonic energy would require alpha approximately 2*a.",
        "The source reports fitted coefficients without parameter uncertainties or covariance, so uncertainty propagation is unavailable for this record.",
    ),
    volume_scale=_DIAMOND_SCALE,
)

DIAMOND_DEWAELE_2008 = EOSRecord(
    identifier="diamond_datchi_dewaele_2008",
    name="Diamond (Dewaele 2008 MGD-Vinet)",
    material="C",
    phase="diamond, cubic Fd-3m",
    cell_contents="8 C atoms per conventional cubic cell",
    eos=MieGruneisenDebye(
        Vinet(5.6693 * 8.0 * _DIAMOND_SCALE, 444.5, 4.18),
        Tr=298.0,
        theta0=1860.0,
        gamma0=0.85,
        q=3.6,
        n=1.0,
    ),
    reference_temperature=298.0,
    reference=_DEWAELE_DIAMOND_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.0, 80.0),
        temperature_k=(298.0, 900.0),
        volume_ratio=(0.85, 1.01),
        notes=(
            "298 K reference-isotherm data and the measured P-T envelope stated in section V; marginal bounds only.",
        ),
    ),
    parameter_provenance=MappingProxyType(
        {
            "rt_eos.V0": "Table III, 298 K fit; 5.6693(16) A^3/atom, fixed in all-data MGD fit",
            "rt_eos.K0": "Table III; 444.5 GPa, fixed from Brillouin data",
            "rt_eos.K0_prime": "Table III, 298 K fit; 4.18(15), fixed in all-data MGD fit",
            "Tr": "Equation 6; 298 K reference isotherm",
            "theta0": "Table III; 1860 K, fixed from heat-capacity data",
            "gamma0": "Table III; 0.85, fixed from ambient-pressure thermal expansion",
            "q": "Table III, all-data MGD fit; 3.6(1.5)",
            "n": "one atom per elemental-carbon formula unit",
        }
    ),
    parameter_errors=MappingProxyType(
        {
            "rt_eos.V0": 0.0016 * 8.0 * _DIAMOND_SCALE / 1.96,
            "rt_eos.K0_prime": 0.15 / 1.96,
            "q": 1.5 / 1.96,
        }
    ),
    notes=(
        "Table III reports 95% confidence intervals; they are divided by 1.96 for one-standard-deviation propagation.",
        "K0, theta0, and gamma0 were fixed and no errors are supplied; their uncertainty is omitted rather than treated as known zero.",
        "No covariance is published; retained parameter errors are propagated independently.",
        "The authors caution that the simple MGD form is not reliable for wide extrapolation beyond the measured range.",
    ),
    volume_scale=_DIAMOND_SCALE,
)

DIAMOND_CORREA_2008_DEWAELE_ANCHORED = EOSRecord(
    identifier="diamond_correa_2008_dewaele_anchored",
    name="Diamond (Correa thermal model, Dewaele 298 K Vinet anchor)",
    material="C",
    phase="diamond, cubic Fd-3m",
    cell_contents="8 C atoms per conventional cubic cell",
    eos=DoubleDebyeLogMomentHelmholtz(
        Vinet(
            5.6693 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
            444.5,
            4.18,
        ),
        Vp=5.571 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        theta_a0=1887.8,
        a_a=-0.316 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_a=0.913,
        theta_b0=1887.8,
        a_b=0.168 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_b=0.429,
        theta_0_0=1887.8,
        a_0=0.131 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_0=0.202,
        n=1.0,
        anharmonic_a=3.8e-5,
        phi0=0.0,
        Tr=298.0,
    ),
    reference_temperature=298.0,
    reference=_DEWAELE_DIAMOND_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.0, 1075.0),
        temperature_k=(1.0, 10000.0),
        volume_ratio=(2.32 / 5.6693, 5.785 / 5.6693),
        notes=(
            "The Dewaele 298 K reference isotherm was measured to 80 GPa; its Vinet representation remains evaluable beyond that calibration coverage.",
            "The Correa thermal contribution uses the simulation range recorded by the unanchored source record.",
        ),
    ),
    parameter_provenance=MappingProxyType(
        {
            "rt_eos.V0": "Dewaele et al. (2008), Table III; 5.6693(16) A^3/atom",
            "rt_eos.K0": "Dewaele et al. (2008), Table III; 444.5 GPa",
            "rt_eos.K0_prime": "Dewaele et al. (2008), Table III; 4.18(15)",
            "Vp": "Correa et al. (2008), Table I; 5.571 A^3/atom",
            "theta_a0": "Correa et al. (2008), Table I; 1887.8 K",
            "a_a": "Correa et al. (2008), Table I; beta_A=-0.316 A^-3",
            "b_a": "Correa et al. (2008), Table I; alpha_A=0.913",
            "theta_b0": "Correa et al. (2008), Table I; 1887.8 K",
            "a_b": "Correa et al. (2008), Table I; beta_B=0.168 A^-3",
            "b_b": "Correa et al. (2008), Table I; alpha_B=0.429",
            "theta_0_0": "Correa et al. (2008), Table I; 1887.8 K",
            "a_0": "Correa et al. (2008), Table I; beta_0=0.131 A^-3",
            "b_0": "Correa et al. (2008), Table I; alpha_0=0.202",
            "n": "one atom per elemental-carbon formula unit",
            "anharmonic_a": "Correa et al. (2008), equation 18 and Table I; 3.8e-5 K^-1/atom",
            "phi0": "Arbitrary zero for the integrated experimental reference isotherm",
            "Tr": "Dewaele et al. (2008), Table III; 298 K reference isotherm",
        }
    ),
    parameter_errors=MappingProxyType(
        {
            "rt_eos.V0": 0.0016 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR / 1.96,
            "rt_eos.K0_prime": 0.15 / 1.96,
        }
    ),
    notes=(
        "Derived reference-isotherm composition: P(V,T)=P_Dewaele(V,298 K)+P_Correa(V,T)-P_Correa(V,298 K).",
        "The record is exactly the experimental Dewaele Vinet isotherm at 298 K and retains Correa's simulated logarithmic-moment double-Debye thermal increment.",
        "Its Helmholtz energy uses an arbitrary additive reference constant; pressure, volume, thermal increments, heat capacities, and same-phase energy differences are unaffected.",
        "The Dewaele and Correa source records remain available independently and preserve their literal published formulations.",
    ),
    volume_scale=_DIAMOND_SCALE,
    scientific_validation_note=(
        "Derived composition of separately primary-source-validated Dewaele "
        "reference-isotherm and Correa thermal components."
    ),
)

DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED = EOSRecord(
    identifier="diamond_benedict_2014_dewaele_anchored",
    name="Diamond (Benedict thermal model, Dewaele 298 K Vinet anchor)",
    material="C",
    phase="diamond, cubic Fd-3m",
    cell_contents="8 C atoms per conventional cubic cell",
    eos=DoubleDebyeHelmholtz(
        Vinet(
            5.6693 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
            444.5,
            4.18,
        ),
        Vp=5.571 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        theta_a0=1887.8,
        a_a=-0.316 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_a=0.913,
        theta_b0=1887.8,
        a_b=0.168 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_b=0.429,
        theta_1_0=1887.8,
        a_1=0.0846 / _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        b_1=0.499,
        n=1.0,
        alpha0=3.79e-5,
        Ve=5.785 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR,
        kappa=0.0,
        phi0=0.0,
        Tr=298.0,
    ),
    reference_temperature=298.0,
    reference=_DEWAELE_DIAMOND_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.0, 1000.0),
        temperature_k=(298.0, 9000.0),
        volume_ratio=(2.5 / 5.6693, 7.0 / 5.6693),
        notes=(
            "The Dewaele 298 K reference isotherm was measured to 80 GPa; its Vinet representation remains evaluable beyond that calibration coverage.",
            "The Benedict thermal contribution uses the simulation range recorded by the unanchored source record.",
        ),
    ),
    parameter_provenance=MappingProxyType(
        {
            "rt_eos.V0": "Dewaele et al. (2008), Table III; 5.6693(16) A^3/atom",
            "rt_eos.K0": "Dewaele et al. (2008), Table III; 444.5 GPa",
            "rt_eos.K0_prime": "Dewaele et al. (2008), Table III; 4.18(15)",
            "Vp": "Benedict et al. (2014), Table I; 5.571 A^3/atom",
            "theta_a0": "Benedict et al. (2014), Table I; 1887.8 K",
            "a_a": "Benedict et al. (2014), Table I; -0.316 A^-3",
            "b_a": "Benedict et al. (2014), Table I; 0.913",
            "theta_b0": "Benedict et al. (2014), Table I; 1887.8 K",
            "a_b": "Benedict et al. (2014), Table I; 0.168 A^-3",
            "b_b": "Benedict et al. (2014), Table I; 0.429",
            "theta_1_0": "Benedict et al. (2014), Table I; 1887.8 K",
            "a_1": "Benedict et al. (2014), Table I; 0.0846 A^-3",
            "b_1": "Benedict et al. (2014), Table I; 0.499",
            "n": "one atom per elemental-carbon formula unit",
            "alpha0": "Benedict et al. (2014), Table I; 3.79e-5 K^-1",
            "Ve": "Benedict et al. (2014), Table I; 5.785 A^3/atom",
            "kappa": "Benedict et al. (2014), Table I; 0.0",
            "phi0": "Arbitrary zero for the integrated experimental reference isotherm",
            "Tr": "Dewaele et al. (2008), Table III; 298 K reference isotherm",
        }
    ),
    parameter_errors=MappingProxyType(
        {
            "rt_eos.V0": 0.0016 * _ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR / 1.96,
            "rt_eos.K0_prime": 0.15 / 1.96,
        }
    ),
    notes=(
        "Derived reference-isotherm composition: P(V,T)=P_Dewaele(V,298 K)+P_Benedict(V,T)-P_Benedict(V,298 K).",
        "The record is exactly the experimental Dewaele Vinet isotherm at 298 K and retains Benedict's simulated first-moment double-Debye thermal increment.",
        "Its Helmholtz energy uses an arbitrary additive reference constant; pressure, volume, thermal increments, heat capacities, and same-phase energy differences are unaffected.",
        "The Dewaele and Benedict source records remain available independently and preserve their literal published formulations.",
    ),
    volume_scale=_DIAMOND_SCALE,
    scientific_validation_note=(
        "Derived composition of separately primary-source-validated Dewaele "
        "reference-isotherm and Benedict thermal components."
    ),
)


def _dewaele_metal_eos_record(
    identifier: str,
    material: str,
    V0_per_atom: float,
    K0: float,
    K0_prime: float,
    K0_error_95: float,
    K0_prime_error_95: float,
    maximum_pressure: float,
) -> EOSRecord:
    return EOSRecord(
        identifier=identifier,
        name=f"{material} fcc (Dewaele 2008)",
        material=material,
        phase="fcc, cubic Fm-3m",
        cell_contents=f"4 {material} atoms per conventional cubic cell",
        eos=Vinet(V0_per_atom * 4.0, K0, K0_prime),
        reference_temperature=300.0,
        reference=_DEWAELE_METALS_REFERENCE,
        validity=ValidityRange(
            pressure_gpa=(0.0, maximum_pressure),
            temperature_k=(300.0, 300.0),
            notes=("Experimental DAC domain from Tables I and II.",),
        ),
        parameter_provenance=MappingProxyType(
            {
                "V0": f"Table IV; {V0_per_atom:g} A^3/atom (fixed)",
                "K0": f"Table IV, unconstrained DAC P'_R fit; {K0:g} GPa",
                "K0_prime": f"Table IV, unconstrained DAC P'_R fit; {K0_prime:g}",
            }
        ),
        parameter_errors=MappingProxyType(
            {
                "K0": K0_error_95 / 1.96,
                "K0_prime": K0_prime_error_95 / 1.96,
            }
        ),
        notes=(
            "Figure 1 uses the unconstrained DAC fit in Table IV; the separate ultrasonic-K0-constrained fit is not substituted.",
            "Published 95% confidence intervals are divided by 1.96 for one-standard-deviation propagation.",
            "No covariance is published; retained errors are propagated independently.",
            "Pressures use the revised P'_R ruby scale defined in the paper and are calibration-dependent.",
        ),
    )


NI_DEWAELE_2008 = _dewaele_metal_eos_record(
    "ni_fcc_dewaele_2008", "Ni", 10.942, 176.7, 5.23, 2.5, 0.09, 157.0
)

AG_DEWAELE_2008 = _dewaele_metal_eos_record(
    "ag_fcc_dewaele_2008", "Ag", 17.021, 103.9, 5.78, 1.4, 0.09, 124.0
)

_RE_V0_CELL = 8.8726e24 / Avogadro * 2.0
_RE_V0_ERROR_95_CELL = 0.01e24 / Avogadro * 2.0
RE_HCP_ANZELLINI_2014 = EOSRecord(
    identifier="re_hcp_anzellini_2014",
    name="Re hcp (Anzellini 2014 quasi-hydrostatic Vinet)",
    material="Re",
    phase="hcp, hexagonal P6_3/mmc",
    cell_contents="2 Re atoms per conventional hexagonal cell",
    eos=Vinet(_RE_V0_CELL, 352.6, 4.56),
    reference_temperature=300.0,
    reference=_ANZELLINI_RE_REFERENCE,
    validity=ValidityRange(
        pressure_gpa=(0.64, 144.0),
        temperature_k=(300.0, 300.0),
        notes=(
            "Quasi-hydrostatic room-temperature XRD range from Tables I and IV.",
            "This record describes the sample EOS; the paper separately tests gasket diffraction under non-hydrostatic conditions.",
        ),
    ),
    parameter_provenance=MappingProxyType(
        {
            "V0": (
                "Equation 6 fit and Table IV: 8.8726 +/- 0.01 cm^3/mol; "
                f"converted to {_RE_V0_CELL:.12g} A^3 per two-atom hcp cell"
            ),
            "K0": "Equation 6 fit and Table IV: 352.6 +/- 8 GPa",
            "K0_prime": "Equation 6 fit and Table IV: 4.56 +/- 0.17",
            "error_convention": (
                "Text immediately following equation 6: error bars are 95% "
                "confidence intervals of fitted values"
            ),
        }
    ),
    parameter_errors=MappingProxyType(
        {
            "V0": _RE_V0_ERROR_95_CELL,
            "K0": 8.0,
            "K0_prime": 0.17,
        }
    ),
    parameter_error_confidence=0.95,
    notes=(
        "Reported 95% confidence half-widths are converted to normal-equivalent standard errors only when uncertainty is propagated.",
        "No parameter covariance is published; the converted errors are therefore propagated independently.",
        "The paper reports pressures from ruby, W, and He gauges and a maximum macroscopic non-hydrostatic stress below its XRD detection limit in the sample geometry.",
    ),
)


_EOS_RECORD_CATALOG = MappingProxyType(
    {
        record.identifier: record
        for record in (
            MGO_TANGE_2009,
            *SOKOLOVA_2013_EOS_RECORDS,
            *FEI_2007_EOS_RECORDS,
            AU_DORFMAN_2012,
            PT_DORFMAN_2012,
            MO_DORFMAN_2012,
            NACL_B2_DORFMAN_2012,
            NE_DORFMAN_2012,
            LIF_B1_DEWAELE_2019,
            NACL_B1_DEWAELE_2019,
            NACL_B2_DEWAELE_2019,
            KCL_B1_DEWAELE_2012,
            KCL_B2_DEWAELE_2012,
            KBR_B1_DEWAELE_2012,
            KBR_B2_DEWAELE_2012,
            CBN_DATCHI_2007,
            DIAMOND_CORREA_2008,
            DIAMOND_BENEDICT_2014,
            DIAMOND_DEWAELE_2008,
            DIAMOND_CORREA_2008_DEWAELE_ANCHORED,
            DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED,
            NI_DEWAELE_2008,
            AG_DEWAELE_2008,
            RE_HCP_ANZELLINI_2014,
        )
    }
)


def _material_identifier(record: EOSRecord) -> str:
    """Return a compact stable identifier from formula and phase identity."""
    if record.material == "C" and record.phase.lower().startswith("diamond"):
        return "diamond"
    formula = re.sub(r"[^a-z0-9]+", "_", record.material.lower()).strip("_")
    phase = record.phase.lower()
    phase_match = re.search(r"\b(b1|b2|fcc|bcc|hcp|diamond)\b", phase)
    phase_key = phase_match.group(1) if phase_match else phase.split(",", 1)[0]
    phase_key = re.sub(r"[^a-z0-9]+", "_", phase_key).strip("_")
    return f"{formula}_{phase_key}"


def _build_material_catalog() -> Mapping[str, Material]:
    grouped: dict[tuple[str, str, str, str], list[EOSRecord]] = {}
    for record in _EOS_RECORD_CATALOG.values():
        key = (record.material, record.phase, record.cell_contents, record.volume_unit)
        grouped.setdefault(key, []).append(record)

    materials: dict[str, Material] = {}
    for records in grouped.values():
        first = records[0]
        identifier = _material_identifier(first)
        if identifier in materials:
            raise MaterialError(
                f"Duplicate generated material identifier {identifier!r}",
                code="material.duplicate_identifier",
                field="identifier",
            )
        materials[identifier] = Material(
            identifier=identifier,
            formula=first.material,
            phase=first.phase,
            cell_contents=first.cell_contents,
            eos_records=tuple(records),
            volume_unit=first.volume_unit,
        )
    return MappingProxyType(materials)


_MATERIAL_CATALOG = _build_material_catalog()

DEFERRED_EOS_RECORDS = (
    DeferredEOSRecord(
        material="Re",
        phase="hcp thermal/high-temperature",
        references=("10.1080/08957959.2018.1448082",),
        reason=(
            "The 300 K Anzellini Vinet scale is implemented separately. A thermal "
            "Re scale is deferred until its complete primary P-V-T formulation, "
            "reference state, covariance, and stated validity domain are audited."
        ),
    ),
)


def get_material(identifier: str) -> Material:
    """Return one bundled material by its stable identifier or alias."""
    from peritheos.catalog import get_material as catalog_get_material

    return catalog_get_material(identifier)


def list_materials(*, formula: str | None = None) -> tuple[Material, ...]:
    """List all bundled materials, optionally filtered by formula."""
    from peritheos.catalog import list_materials as catalog_list_materials

    return catalog_list_materials(formula=formula)


def get_eos_record(identifier: str) -> EOSRecord:
    """Return one bundled EOS record by its stable identifier or alias."""
    from peritheos.catalog import get_eos_record as catalog_get_eos_record

    return catalog_get_eos_record(identifier)


def list_eos_records(*, formula: str | None = None) -> tuple[EOSRecord, ...]:
    """List all bundled EOS records, optionally filtered by formula."""
    from peritheos.catalog import list_eos_records as catalog_list_eos_records

    return catalog_list_eos_records(formula=formula)


# Keep discovery beside the historical material lookup imports as well as at
# the package root. catalog_search imports catalog listing lazily, so this does
# not make construction of data classes depend on catalog initialization.
from peritheos.catalog_search import (  # noqa: E402, I001
    search_eos_records,
    search_materials,
)


__all__ = [
    "AG_DEWAELE_2008",
    "AG_SOKOLOVA_2013",
    "AL_SOKOLOVA_2013",
    "AU_DORFMAN_2012",
    "AU_FEI_2007",
    "AU_SOKOLOVA_2013",
    "CBN_DATCHI_2007",
    "DEFERRED_EOS_RECORDS",
    "DeferredEOSRecord",
    "DIAMOND_BENEDICT_2014",
    "DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED",
    "DIAMOND_CORREA_2008",
    "DIAMOND_CORREA_2008_DEWAELE_ANCHORED",
    "DIAMOND_DEWAELE_2008",
    "DIAMOND_SOKOLOVA_2013",
    "FEI_2007_EOS_RECORDS",
    "KBR_B1_DEWAELE_2012",
    "KBR_B2_DEWAELE_2012",
    "KCL_B1_DEWAELE_2012",
    "KCL_B2_DEWAELE_2012",
    "LIF_B1_DEWAELE_2019",
    "LiteratureReference",
    "Material",
    "MGO_TANGE_2009",
    "MGO_SOKOLOVA_2013",
    "MO_DORFMAN_2012",
    "MO_SOKOLOVA_2013",
    "NACL_B2_FEI_2007",
    "NACL_B1_DEWAELE_2019",
    "NACL_B2_DORFMAN_2012",
    "NACL_B2_DEWAELE_2019",
    "NE_DORFMAN_2012",
    "NE_FEI_2007",
    "NB_SOKOLOVA_2013",
    "NI_DEWAELE_2008",
    "PT_DORFMAN_2012",
    "PT_FEI_2007",
    "PT_SOKOLOVA_2013",
    "RE_HCP_ANZELLINI_2014",
    "EOSRecord",
    "SOKOLOVA_2013_EOS_RECORDS",
    "TA_SOKOLOVA_2013",
    "ValidityRange",
    "W_SOKOLOVA_2013",
    "CU_SOKOLOVA_2013",
    "get_eos_record",
    "get_material",
    "list_eos_records",
    "list_materials",
    "material_from_dict",
    "search_eos_records",
    "search_materials",
]
