"""
Equations of state module for Peritheos
"""

import inspect
from functools import cache
from typing import Any, Callable, Union

import numpy as np
from numpy.typing import NDArray
from scipy import optimize

from peritheos import _rust as _rust
from peritheos.errors import (
    ConfigurationError,
    EosNumericalError,
    EosValidationError,
    UnsupportedOperationError,
)

# Type alias for numeric values (scalar or array)
NumericType = Union[float, NDArray[np.float64]]


@cache
def _native_evaluation_types() -> tuple[type, ...]:
    """Return classes whose inherited behavior exactly matches a Rust model."""
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
        LinearThermalPressure,
        LogVolumeThermalPressure,
        MieGruneisenDebye,
        MieGruneisenEinstein,
        Sokolova2016,
        Tange2009Debye,
        ThermalModifiedTait,
        ThermalReferenceStateEOS,
    )

    return (
        BM2,
        BM3,
        BM4,
        Murnaghan,
        ModifiedTait,
        NaturalStrain2,
        NaturalStrain3,
        NaturalStrain4,
        Vinet,
        Holzapfel,
        LinearThermalPressure,
        LogVolumeThermalPressure,
        MieGruneisenDebye,
        MieGruneisenEinstein,
        Tange2009Debye,
        ThermalModifiedTait,
        ThermalReferenceStateEOS,
        Sokolova2016,
    )


def _native_for_exact_model(model):
    """Return a native handle only when no subclass behavior can be bypassed."""
    if type(model) not in _native_evaluation_types():
        return None
    return getattr(model, "_native", None)


def _native_rt_evaluate(native, quantity: str, values: NumericType) -> NumericType:
    """Evaluate a private native RT model while preserving NumPy shape semantics."""
    array = np.asarray(values, dtype=float)
    if array.ndim == 0:
        return float(getattr(native, f"{quantity}_scalar")(float(array)))
    return np.asarray(getattr(native, f"{quantity}_array")(array), dtype=float)


def _native_thermal_evaluate(
    native, quantity: str, first: NumericType, second: NumericType
) -> NumericType:
    """Evaluate a private native thermal model with NumPy broadcasting."""
    left = np.asarray(first, dtype=float)
    right = np.asarray(second, dtype=float)
    try:
        left, right = np.broadcast_arrays(left, right)
    except ValueError as error:
        raise EosValidationError(
            "V and T must have broadcast-compatible shapes"
        ) from error
    if left.ndim == 0:
        return float(native.evaluate_scalar(quantity, float(left), float(right)))
    return np.asarray(native.evaluate_array(quantity, left, right), dtype=float)


def validate_finite_scalar(value: float, name: str) -> float:
    """Return *value* as a finite float or raise a descriptive error."""
    value = float(value)
    if not np.isfinite(value):
        raise EosValidationError(f"{name} must be finite")
    return value


def validate_positive_scalar(value: float, name: str) -> float:
    """Return *value* as a positive finite float."""
    value = validate_finite_scalar(value, name)
    if value <= 0:
        raise EosValidationError(f"{name} must be greater than zero")
    return value


def validate_volume(V: NumericType) -> NumericType:
    """Validate volume input while preserving scalar or array behaviour."""
    values = np.asarray(V, dtype=float)
    if not np.all(np.isfinite(values)):
        raise EosValidationError("Volume must be finite")
    if np.any(values <= 0):
        raise EosValidationError("Volume must be greater than zero")
    if values.ndim == 0:
        return float(values)
    return values


def _scalar_pressure(
    pressure_function: Callable[[float], NumericType],
    value: float,
    *,
    variable_name: str = "V",
) -> float:
    """Evaluate a pressure callable and require a finite scalar result."""
    pressure = np.asarray(pressure_function(value), dtype=float)
    if pressure.ndim != 0:
        raise ConfigurationError(
            f"The pressure function must return a scalar for scalar {variable_name}"
        )
    result = float(pressure)
    if not np.isfinite(result):
        raise EosNumericalError(
            f"EOS returned a non-finite pressure at {variable_name}={value}"
        )
    return result


def solve_volume(
    pressure_function: Callable[[float], NumericType],
    pressure: float,
    reference_volume: float,
) -> float:
    """Solve ``pressure_function(V) == pressure`` on the branch nearest V0.

    The EOS is assumed to be locally monotonic around its reference volume:
    pressure increases on compression and decreases on expansion. Starting at
    ``V0``, the search geometrically decreases volume for a higher target
    pressure or increases it for a lower target pressure until the pressure
    residual changes sign. ``scipy.optimize.brentq`` then solves inside that
    positive-volume bracket. If no sign change is found, the requested state is
    reported as outside the model's invertible range.
    """
    target = validate_finite_scalar(pressure, "Pressure")
    V0 = validate_positive_scalar(reference_volume, "Reference volume")
    p0 = _scalar_pressure(pressure_function, V0)
    f0 = p0 - target
    pressure_tolerance = 1e-10 * max(1.0, abs(target), abs(p0))
    if abs(f0) <= pressure_tolerance:
        return V0

    if f0 < 0:
        # The target is above P(V0), so search towards compression.
        upper, f_upper = V0, f0
        lower = V0
        f_lower = f0
        for _ in range(160):
            lower *= 0.8
            if lower <= V0 * 1e-14:
                break
            f_lower = _scalar_pressure(pressure_function, lower) - target
            if f_lower >= 0:
                break
        else:
            f_lower = np.nan
        if not np.isfinite(f_lower) or f_lower < 0:
            raise EosValidationError(
                f"Could not bracket a positive volume for pressure {target}"
            )
    else:
        # The target is below P(V0), so search along the first expansion branch.
        lower, f_lower = V0, f0
        upper = V0
        f_upper = f0
        for _ in range(160):
            upper *= 1.05
            if upper >= V0 * 1e4:
                break
            f_upper = _scalar_pressure(pressure_function, upper) - target
            if f_upper <= 0:
                break
        else:
            f_upper = np.nan
        if not np.isfinite(f_upper) or f_upper > 0:
            raise EosValidationError(
                f"Pressure {target} is outside the invertible expansion range"
            )

    result = optimize.brentq(
        lambda volume: _scalar_pressure(pressure_function, volume) - target,
        lower,
        upper,
        xtol=max(np.finfo(float).eps * V0, 1e-14),
        rtol=1e-12,
    )
    residual = abs(_scalar_pressure(pressure_function, result) - target)
    if residual > 1e-8 * max(1.0, abs(target)):
        raise EosNumericalError(
            f"Volume inversion did not converge to the requested pressure; residual={residual}"
        )
    return float(result)


def solve_temperature(
    pressure_function: Callable[[float], NumericType],
    pressure: float,
    reference_temperature: float,
) -> float:
    """Solve ``pressure_function(T) == pressure`` near the reference isotherm.

    Starting at the reference temperature, the search expands geometrically
    towards both lower and higher positive temperatures. The first bracketed
    root is therefore on the branch nearest the reference isotherm. A target
    that cannot be bracketed over the positive-temperature search range is
    reported as outside the model's invertible range.
    """
    target = validate_finite_scalar(pressure, "Pressure")
    Tr = validate_positive_scalar(reference_temperature, "Reference temperature")
    p0 = _scalar_pressure(pressure_function, Tr, variable_name="T")
    f0 = p0 - target
    pressure_tolerance = 1e-10 * max(1.0, abs(target), abs(p0))
    if abs(f0) <= pressure_tolerance:
        return Tr

    lower = upper = Tr
    f_lower = f_upper = f0
    minimum_temperature = Tr * 1.0e-14
    maximum_temperature = Tr * 1.0e8
    lower_active = upper_active = True
    brackets: list[tuple[float, float]] = []

    for _ in range(160):
        if lower_active:
            next_lower = max(lower * 0.8, minimum_temperature)
            next_f_lower = (
                _scalar_pressure(pressure_function, next_lower, variable_name="T")
                - target
            )
            if next_f_lower * f_lower <= 0.0:
                brackets.append((next_lower, lower))
            lower, f_lower = next_lower, next_f_lower
            lower_active = lower > minimum_temperature

        if upper_active:
            next_upper = min(upper * 1.25, maximum_temperature)
            next_f_upper = (
                _scalar_pressure(pressure_function, next_upper, variable_name="T")
                - target
            )
            if next_f_upper * f_upper <= 0.0:
                brackets.append((upper, next_upper))
            upper, f_upper = next_upper, next_f_upper
            upper_active = upper < maximum_temperature

        if brackets:
            break
        if not lower_active and not upper_active:
            break

    if not brackets:
        raise EosValidationError(
            f"Pressure {target} is outside the invertible temperature range"
        )

    roots = [
        optimize.brentq(
            lambda temperature: (
                _scalar_pressure(pressure_function, temperature, variable_name="T")
                - target
            ),
            bracket_lower,
            bracket_upper,
            xtol=np.finfo(float).eps * max(1.0, Tr),
            rtol=1e-12,
        )
        for bracket_lower, bracket_upper in brackets
    ]
    result = min(roots, key=lambda temperature: abs(np.log(temperature / Tr)))
    residual = abs(
        _scalar_pressure(pressure_function, result, variable_name="T") - target
    )
    if residual > 1e-8 * max(1.0, abs(target)):
        raise EosNumericalError(
            "Temperature inversion did not converge to the requested pressure; "
            f"residual={residual}"
        )
    return float(result)


class EosBase:
    """
    Base class for equation of state implementations.

    This abstract class defines the interface that all equation of state
    implementations should follow.
    """

    _constructor_configuration_names: tuple[str, ...] = ()

    def _own_parameter_names(self) -> tuple[str, ...]:
        """Return numeric constructor parameters represented by public attributes."""
        signature = inspect.signature(type(self).__init__)
        names = []
        for name, parameter in signature.parameters.items():
            if name in {"self", "rt_eos"} or parameter.kind in {
                inspect.Parameter.VAR_POSITIONAL,
                inspect.Parameter.VAR_KEYWORD,
            }:
                continue
            if name in self._constructor_configuration_names:
                continue
            if not hasattr(self, name):
                raise UnsupportedOperationError(
                    f"{type(self).__name__} does not expose constructor parameter {name!r}"
                )
            names.append(name)
        return tuple(names)

    def configuration_values(self) -> dict[str, Any]:
        """Return fixed non-numeric constructor choices needed for reconstruction."""
        values = {}
        for name in self._constructor_configuration_names:
            if not hasattr(self, name):
                raise UnsupportedOperationError(
                    f"{type(self).__name__} does not expose constructor configuration "
                    f"{name!r}"
                )
            values[name] = getattr(self, name)
        return values

    def parameter_values(self, *, include_reference: bool = True) -> dict[str, float]:
        """Return reconstructable EOS parameters in constructor order.

        Thermal EOS reference-isotherm parameters use dotted names such as
        ``rt_eos.V0`` when ``include_reference`` is true.
        """
        values = {
            name: float(getattr(self, name)) for name in self._own_parameter_names()
        }
        reference = getattr(self, "rt_eos", None)
        if include_reference and isinstance(reference, EosBase):
            values.update(
                {
                    f"rt_eos.{name}": value
                    for name, value in reference.parameter_values(
                        include_reference=True
                    ).items()
                }
            )
        return values

    def with_parameters(self, **updates: float) -> "EosBase":
        """Reconstruct the EOS after replacing selected parameter values.

        Reconstruction, rather than attribute mutation, ensures that derived
        constants maintained by an EOS constructor remain consistent.
        """
        available = self.parameter_values(include_reference=True)
        unknown = set(updates) - set(available)
        if unknown:
            raise EosValidationError(
                f"Unknown parameters for {type(self).__name__}: {sorted(unknown)}"
            )

        own_values: dict[str, Any] = self.parameter_values(include_reference=False)
        own_values.update(self.configuration_values())
        own_values.update(
            {
                name: float(value)
                for name, value in updates.items()
                if not name.startswith("rt_eos.")
            }
        )
        reference = getattr(self, "rt_eos", None)
        if isinstance(reference, EosBase):
            reference_updates = {
                name.removeprefix("rt_eos."): float(value)
                for name, value in updates.items()
                if name.startswith("rt_eos.")
            }
            own_values["rt_eos"] = reference.with_parameters(**reference_updates)
        return type(self)(**own_values)

    def pressure(self, V: NumericType) -> NumericType:
        """
        Calculate pressure at a given volume.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in cubic angstroms or any consistent unit)

        Returns
        -------
        float or numpy.ndarray
            Pressure (in the same units as K0)
        """
        raise UnsupportedOperationError(
            "Subclasses must implement the pressure method."
        )

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Calculate the bulk modulus.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in cubic angstroms or any consistent unit)

        Returns
        -------
        float or numpy.ndarray
            Bulk modulus (in the same units as K0)
        """
        raise UnsupportedOperationError(
            "Subclasses must implement the bulk_modulus method."
        )

    def bulk_modulus_derivative(
        self, V: NumericType, relative_step: float = 1.0e-6
    ) -> NumericType:
        """Return the pressure derivative of bulk modulus, ``dK/dP``.

        The default implementation differentiates ``K(V)`` and ``P(V)`` with
        the same central volume step. Individual EOS classes may override this
        with an analytical expression or a more efficient equivalent.
        """
        relative_step = validate_positive_scalar(relative_step, "relative_step")
        volumes = np.asarray(validate_volume(V), dtype=float)
        steps = relative_step * volumes
        lower_volumes = volumes - steps
        upper_volumes = volumes + steps
        pressure_difference = np.asarray(
            self.pressure(upper_volumes), dtype=float
        ) - np.asarray(self.pressure(lower_volumes), dtype=float)
        if np.any(pressure_difference == 0.0):
            raise EosNumericalError("Cannot evaluate dK/dP where dP/dV is zero")
        result = (
            np.asarray(self.bulk_modulus(upper_volumes), dtype=float)
            - np.asarray(self.bulk_modulus(lower_volumes), dtype=float)
        ) / pressure_difference
        if not np.all(np.isfinite(result)):
            raise EosNumericalError("Bulk-modulus pressure derivative is not finite")
        if result.ndim == 0:
            return float(result)
        return result

    def calculate_volume(self, P: NumericType) -> NumericType:
        """
        Calculate volume at a given pressure using a bracketed root solver.

        Parameters
        ----------
        P : float or numpy.ndarray
            Pressure (in the same units as K0)

        Returns
        -------
        float or numpy.ndarray
            Volume (in the same units as V0)
        """
        pressures = np.asarray(P, dtype=float)
        if not np.all(np.isfinite(pressures)):
            raise EosValidationError("Pressure must be finite")
        native = _native_for_exact_model(self)
        if native is not None:
            return _native_rt_evaluate(native, "volume", pressures)
        if pressures.ndim == 0:
            return solve_volume(self.pressure, float(pressures), self.V0)
        return np.array(
            [solve_volume(self.pressure, value, self.V0) for value in pressures.flat]
        ).reshape(pressures.shape)

    def volume(self, P: NumericType) -> NumericType:
        """Alias for :meth:`calculate_volume` using the P-to-V terminology."""
        return self.calculate_volume(P)


class ThermalEOS(EosBase):
    """Base class for molar thermal equations of state.

    Thermal EOS implementations use temperature in K, pressure in GPa, and
    molar volume in J bar^-1 mol^-1. The latter is equal to cm^3 mol^-1 / 10.
    """

    def __init__(self, rt_eos: EosBase):
        self.rt_eos = rt_eos

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        raise UnsupportedOperationError(
            "This method should be implemented by the subclass"
        )

    def _thermal_pressure_function(self, V: float) -> Callable[[float], NumericType]:
        """Return a temperature-only thermal-pressure function at fixed volume."""
        return lambda temperature: self.thermal_pressure(V, temperature)

    @staticmethod
    def _validate_f_dac(f_dac: float) -> float:
        f_dac = validate_finite_scalar(f_dac, "f_dac")
        if f_dac < 0.0 or f_dac >= 1.0:
            raise EosValidationError(
                "f_dac must be greater than or equal to zero and less than one"
            )
        return f_dac

    def dac_thermal_pressure(
        self, V: NumericType, T: NumericType, f_dac: float
    ) -> NumericType:
        """Return the effective DAC pressure contribution in GPa.

        The contribution is ``f_dac * thermal_pressure(V, T)``. It models
        the portion of the EOS thermal pressure retained by confinement in a
        diamond-anvil cell.
        """
        f_dac = self._validate_f_dac(f_dac)
        result = f_dac * np.asarray(self.thermal_pressure(V, T), dtype=float)
        return self._scalar_or_array(result)

    def pressure(self, V: NumericType, T: NumericType) -> NumericType:
        native = _native_for_exact_model(self)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(native, "pressure", volumes, temperatures)
        return self.thermal_pressure(V, T) + self.rt_eos.pressure(V)

    @staticmethod
    def _broadcast_state(V: NumericType, T: NumericType):
        volumes = np.asarray(validate_volume(V), dtype=float)
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise EosValidationError("Temperature must be finite and greater than zero")
        try:
            return np.broadcast_arrays(volumes, temperatures)
        except ValueError as error:
            raise EosValidationError(
                "V and T must have broadcast-compatible shapes"
            ) from error

    @staticmethod
    def _scalar_or_array(values: NDArray[np.float64]) -> NumericType:
        if values.ndim == 0:
            return float(values)
        return values

    def bulk_modulus(
        self, V: NumericType, T: NumericType, relative_step: float = 1.0e-6
    ) -> NumericType:
        """Return the isothermal bulk modulus ``-V (dP/dV)_T`` in GPa."""
        relative_step = validate_positive_scalar(relative_step, "relative_step")
        native = _native_for_exact_model(self)
        if native is not None and relative_step == 1.0e-6:
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                native, "bulk_modulus", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        steps = relative_step * volumes
        derivative = (
            self.pressure(volumes + steps, temperatures)
            - self.pressure(volumes - steps, temperatures)
        ) / (2.0 * steps)
        result = -volumes * derivative
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def isothermal_compressibility(self, V: NumericType, T: NumericType) -> NumericType:
        """Return isothermal compressibility in GPa^-1."""
        native = _native_for_exact_model(self)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                native, "isothermal_compressibility", volumes, temperatures
            )
        result = 1.0 / np.asarray(self.bulk_modulus(V, T), dtype=float)
        return self._scalar_or_array(result)

    def thermal_expansivity(
        self, V: NumericType, T: NumericType, relative_step: float = 1.0e-5
    ) -> NumericType:
        """Return volumetric thermal expansivity in K^-1.

        This uses ``alpha = (dP/dT)_V / K_T``.
        """
        relative_step = validate_positive_scalar(relative_step, "relative_step")
        native = _native_for_exact_model(self)
        if native is not None and relative_step == 1.0e-5:
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                native, "thermal_expansivity", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        steps = np.minimum(relative_step * temperatures, 0.49 * temperatures)
        pressure_derivative = (
            self.pressure(volumes, temperatures + steps)
            - self.pressure(volumes, temperatures - steps)
        ) / (2.0 * steps)
        result = pressure_derivative / np.asarray(
            self.bulk_modulus(volumes, temperatures), dtype=float
        )
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Return constant-volume molar heat capacity in J mol^-1 K^-1."""
        raise UnsupportedOperationError(
            f"{type(self).__name__} does not provide a caloric model"
        )

    def molar_heat_capacity_p(self, V: NumericType, T: NumericType) -> NumericType:
        """Return constant-pressure molar heat capacity in J mol^-1 K^-1."""
        native = _native_for_exact_model(self)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                native, "molar_heat_capacity_p", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        cv = np.asarray(self.molar_heat_capacity_v(volumes, temperatures), dtype=float)
        alpha = np.asarray(self.thermal_expansivity(volumes, temperatures), dtype=float)
        bulk_modulus = np.asarray(self.bulk_modulus(volumes, temperatures), dtype=float)
        # GPa * J bar^-1 = 10^4 J.
        result = cv + alpha**2 * bulk_modulus * volumes * temperatures * 1.0e4
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Alias for :meth:`molar_heat_capacity_v`."""
        return self.molar_heat_capacity_v(V, T)

    def heat_capacity_p(self, V: NumericType, T: NumericType) -> NumericType:
        """Alias for :meth:`molar_heat_capacity_p`."""
        return self.molar_heat_capacity_p(V, T)

    def gruneisen_parameter(self, V: NumericType, T: NumericType) -> NumericType:
        """Return thermodynamic Gruneisen parameter ``alpha K_T V / C_V``."""
        volumes, temperatures = self._broadcast_state(V, T)
        alpha = np.asarray(self.thermal_expansivity(volumes, temperatures), dtype=float)
        bulk_modulus = np.asarray(self.bulk_modulus(volumes, temperatures), dtype=float)
        cv = np.asarray(self.molar_heat_capacity_v(volumes, temperatures), dtype=float)
        result = alpha * bulk_modulus * volumes * 1.0e4 / cv
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def adiabatic_bulk_modulus(self, V: NumericType, T: NumericType) -> NumericType:
        """Return adiabatic bulk modulus ``K_S = K_T C_P / C_V`` in GPa."""
        native = _native_for_exact_model(self)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                native, "adiabatic_bulk_modulus", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        kt = np.asarray(self.bulk_modulus(volumes, temperatures), dtype=float)
        cv = np.asarray(self.molar_heat_capacity_v(volumes, temperatures), dtype=float)
        cp = np.asarray(self.molar_heat_capacity_p(volumes, temperatures), dtype=float)
        result = kt * cp / cv
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def calculate_volume(self, P: NumericType, T: NumericType) -> NumericType:
        """Calculate volume at pressure and temperature using bracketed roots."""
        pressures, temperatures = np.broadcast_arrays(
            np.asarray(P, dtype=float), np.asarray(T, dtype=float)
        )
        if not np.all(np.isfinite(pressures)):
            raise EosValidationError("Pressure must be finite")
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise EosValidationError("Temperature must be finite and greater than zero")

        native = _native_for_exact_model(self)
        if native is not None:
            return _native_thermal_evaluate(native, "volume", pressures, temperatures)

        volumes = np.array(
            [
                solve_volume(
                    lambda volume, temperature=float(temperature): self.pressure(
                        volume, temperature
                    ),
                    float(pressure),
                    self.rt_eos.V0,
                )
                for pressure, temperature in zip(pressures.flat, temperatures.flat)
            ]
        ).reshape(pressures.shape)
        if volumes.ndim == 0:
            return float(volumes)
        return volumes

    def volume(self, P: NumericType, T: NumericType) -> NumericType:
        """Alias for :meth:`calculate_volume` for a thermal EOS."""
        return self.calculate_volume(P, T)

    def calculate_temperature(
        self,
        P: NumericType,
        V: NumericType,
    ) -> NumericType:
        """Calculate temperature at total pressure and volume using bracketed roots.

        The returned root is on the positive-temperature branch nearest the
        reference temperature.

        Parameters
        ----------
        P : float or numpy.ndarray
            Pressure in GPa.
        V : float or numpy.ndarray
            Molar volume in J bar^-1 mol^-1.

        Returns
        -------
        float or numpy.ndarray
            Temperature in K.
        """
        volumes = np.asarray(validate_volume(V), dtype=float)
        pressures = np.asarray(P, dtype=float)
        if not np.all(np.isfinite(pressures)):
            raise EosValidationError("Pressure must be finite")
        try:
            pressures, volumes = np.broadcast_arrays(pressures, volumes)
        except ValueError as error:
            raise EosValidationError(
                "P and V must have broadcast-compatible shapes"
            ) from error

        native = _native_for_exact_model(self)
        if native is not None:
            return _native_thermal_evaluate(native, "temperature", pressures, volumes)

        cold_pressures = np.asarray(self.rt_eos.pressure(volumes), dtype=float)
        target_thermal_pressures = pressures - cold_pressures
        temperatures = np.array(
            [
                solve_temperature(
                    self._thermal_pressure_function(float(volume)),
                    float(target_thermal_pressure),
                    self.Tr,
                )
                for target_thermal_pressure, volume in zip(
                    target_thermal_pressures.flat, volumes.flat
                )
            ]
        ).reshape(target_thermal_pressures.shape)
        if temperatures.ndim == 0:
            return float(temperatures)
        return temperatures

    def temperature(
        self,
        P: NumericType,
        V: NumericType,
    ) -> NumericType:
        """Alias for :meth:`calculate_temperature` for a thermal EOS."""
        return self.calculate_temperature(P, V)

    def temperature_from_volumes(
        self,
        V_ambient: NumericType,
        V_heated: NumericType,
        *,
        f_dac: float,
    ) -> NumericType:
        """Infer temperature from ambient and heated molar volumes.

        The ambient pressure is obtained from ``V_ambient`` on the reference
        isotherm. At the heated volume, this method solves

        ``pressure(V_heated, T) = P_ambient +``
        ``f_dac * thermal_pressure(V_heated, T)``.

        Since total EOS pressure is cold pressure plus thermal pressure, the
        solved expression is equivalently

        ``thermal_pressure(V_heated, T) =``
        ``(P_cold(V_ambient) - P_cold(V_heated)) / (1 - f_dac)``.

        Here ``f_dac`` is defined as ``(P_hot - P_ambient)`` divided by the
        EOS thermal pressure. It is not a fraction of the cold pressure and
        should be treated as an experimental boundary-condition assumption
        unless independently calibrated.

        Both volumes use J bar^-1 mol^-1 and may be broadcast-compatible
        scalars or arrays. ``f_dac`` must be greater than or equal to zero and
        less than one. Only heated states at or above the reference temperature are
        accepted. This empirical confinement model is distinct from ordinary
        total-pressure inversion with :meth:`temperature`. The result is in K.
        """
        f_dac = self._validate_f_dac(f_dac)
        ambient_volumes = np.asarray(validate_volume(V_ambient), dtype=float)
        heated_volumes = np.asarray(validate_volume(V_heated), dtype=float)
        try:
            ambient_volumes, heated_volumes = np.broadcast_arrays(
                ambient_volumes, heated_volumes
            )
        except ValueError as error:
            raise EosValidationError(
                "V_ambient and V_heated must have broadcast-compatible shapes"
            ) from error

        ambient_pressures = np.asarray(
            self.rt_eos.pressure(ambient_volumes), dtype=float
        )
        heated_cold_pressures = np.asarray(
            self.rt_eos.pressure(heated_volumes), dtype=float
        )
        target_thermal_pressures = (ambient_pressures - heated_cold_pressures) / (
            1.0 - f_dac
        )
        if np.any(target_thermal_pressures < 0.0):
            raise EosValidationError(
                "The volume pair implies a temperature below the reference "
                "temperature, not a heated state"
            )
        temperatures = np.array(
            [
                solve_temperature(
                    self._thermal_pressure_function(float(heated_volume)),
                    float(target_thermal_pressure),
                    self.Tr,
                )
                for target_thermal_pressure, heated_volume in zip(
                    target_thermal_pressures.flat, heated_volumes.flat
                )
            ]
        ).reshape(target_thermal_pressures.shape)
        temperature_tolerance = 1.0e-10 * max(1.0, self.Tr)
        if np.any(temperatures < self.Tr - temperature_tolerance):
            raise EosValidationError(
                "The volume pair implies a temperature below the reference "
                "temperature, not a heated state"
            )
        if temperatures.ndim == 0:
            return float(temperatures)
        return temperatures
