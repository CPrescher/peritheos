"""Mie-Gruneisen thermal equations of state."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy.constants import R
from scipy.integrate import quad

from peritheos.eos import (
    EosBase,
    NumericType,
    ThermalEOS,
    _native_for_exact_model,
    _native_thermal_evaluate,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class _MieGruneisenBase(ThermalEOS, ABC):
    """Shared quasi-harmonic Mie-Gruneisen thermal-pressure model."""

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        q: float,
        n: float,
    ) -> None:
        if not isinstance(rt_eos, EosBase):
            raise TypeError("rt_eos must be an equation of state")
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.theta0 = validate_positive_scalar(theta0, "theta0")
        self.gamma0 = validate_finite_scalar(gamma0, "gamma0")
        self.q = validate_finite_scalar(q, "q")
        self.n = validate_positive_scalar(n, "n")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) in (
            MieGruneisenDebye,
            MieGruneisenEinstein,
        ):
            from peritheos import _rust

            factory = (
                _rust.ThermalEos.mie_gruneisen_debye
                if type(self) is MieGruneisenDebye
                else _rust.ThermalEos.mie_gruneisen_einstein
            )
            self._native = factory(
                reference_native, self.Tr, self.theta0, self.gamma0, self.q, self.n
            )

    def _native_evaluate(
        self, quantity: str, first: NumericType, second: NumericType
    ) -> NumericType:
        return _native_thermal_evaluate(self._native, quantity, first, second)

    def gruneisen_parameter(
        self, V: NumericType, T: NumericType | None = None
    ) -> NumericType:
        """Return ``gamma(V) = gamma0 * (V / V0)**q``."""
        native = getattr(self, "_native", None)
        if native is not None:
            temperatures = self.Tr if T is None else T
            volumes, temperatures = self._broadcast_state(V, temperatures)
            return self._native_evaluate("gruneisen_parameter", volumes, temperatures)
        V = validate_volume(V)
        return self.gamma0 * np.exp(self.q * np.log(V / self.rt_eos.V0))

    def characteristic_temperature(self, V: NumericType) -> NumericType:
        """Return the characteristic lattice temperature at volume *V*.

        Its volume dependence is thermodynamically consistent with
        ``gamma = -d(log(theta)) / d(log(V))``.
        """
        native = getattr(self, "_native", None)
        if native is not None:
            volumes = validate_volume(V)
            return self._native_evaluate("characteristic_temperature", volumes, self.Tr)
        V = validate_volume(V)
        logarithmic_volume = np.log(V / self.rt_eos.V0)
        if self.q == 0.0:
            exponent = -self.gamma0 * logarithmic_volume
        else:
            exponent = -self.gamma0 * np.expm1(self.q * logarithmic_volume) / self.q
        theta = self.theta0 * np.exp(exponent)
        if not np.all(np.isfinite(theta)):
            raise ArithmeticError("Characteristic temperature is not finite")
        return theta

    @abstractmethod
    def thermal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return molar vibrational thermal energy in J mol^-1."""

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return thermal pressure relative to ``Tr`` in GPa.

        Volumes must be molar volumes in J bar^-1 mol^-1. This is equivalent
        to cm^3 mol^-1 divided by ten.
        """
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("thermal_pressure", volumes, temperatures)
        V = validate_volume(V)
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise ValueError("Temperature must be finite and greater than zero")
        try:
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), temperatures
            )
        except ValueError as error:
            raise ValueError("V and T must have broadcast-compatible shapes") from error

        energy_difference = self.thermal_energy(
            volumes, temperatures
        ) - self.thermal_energy(volumes, self.Tr)
        # gamma * E / V is in bar for E [J/mol] and V [J/bar/mol].
        pressure = (
            self.gruneisen_parameter(volumes) * energy_difference / volumes / 10000.0
        )
        if pressure.ndim == 0:
            return float(pressure)
        return pressure

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Return vibrational ``C_V`` in J mol^-1 K^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("molar_heat_capacity_v", volumes, temperatures)
        volumes, temperatures = self._broadcast_state(V, T)
        steps = 1.0e-5 * temperatures
        result = (
            self.thermal_energy(volumes, temperatures + steps)
            - self.thermal_energy(volumes, temperatures - steps)
        ) / (2.0 * steps)
        return self._scalar_or_array(np.asarray(result, dtype=float))

    @abstractmethod
    def thermal_entropy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return vibrational entropy in J mol^-1 K^-1."""

    def thermal_internal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return vibrational internal energy in J mol^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate(
                "thermal_internal_energy", volumes, temperatures
            )
        return self.thermal_energy(V, T)

    def vibrational_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the unreferenced vibrational pressure in GPa."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("vibrational_pressure", volumes, temperatures)
        volumes, temperatures = self._broadcast_state(V, T)
        pressure = (
            self.gruneisen_parameter(volumes)
            * np.asarray(self.thermal_energy(volumes, temperatures), dtype=float)
            / volumes
            / 1.0e4
        )
        return self._scalar_or_array(np.asarray(pressure, dtype=float))

    def thermal_helmholtz_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return vibrational Helmholtz energy in J mol^-1.

        Zero-point and static reference energies are omitted.
        """
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate(
                "thermal_helmholtz_free_energy", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        result = np.asarray(
            self.thermal_energy(volumes, temperatures), dtype=float
        ) - temperatures * np.asarray(
            self.thermal_entropy(volumes, temperatures), dtype=float
        )
        return self._scalar_or_array(result)

    def thermal_enthalpy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the vibrational enthalpy contribution in J mol^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("thermal_enthalpy", volumes, temperatures)
        volumes, temperatures = self._broadcast_state(V, T)
        result = (
            np.asarray(self.thermal_energy(volumes, temperatures), dtype=float)
            + np.asarray(self.vibrational_pressure(volumes, temperatures), dtype=float)
            * volumes
            * 1.0e4
        )
        return self._scalar_or_array(result)

    def thermal_gibbs_free_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the vibrational Gibbs-energy contribution in J mol^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate(
                "thermal_gibbs_free_energy", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        result = (
            np.asarray(
                self.thermal_helmholtz_free_energy(volumes, temperatures), dtype=float
            )
            + np.asarray(self.vibrational_pressure(volumes, temperatures), dtype=float)
            * volumes
            * 1.0e4
        )
        return self._scalar_or_array(result)


class MieGruneisenDebye(_MieGruneisenBase):
    """Mie-Gruneisen-Debye thermal equation of state.

    The thermal pressure is

    ``Delta P = gamma(V) / V * (E_D(V, T) - E_D(V, Tr))``,

    where ``E_D`` is the Debye vibrational energy, ``gamma(V) = gamma0 *
    (V/V0)**q``, and the Debye temperature follows from
    ``gamma = -d(log(theta))/d(log(V))``.

    Reference
    ---------
    Jackson, I. & Rigden, S. M. (1996). Analysis of P-V-T data: constraints
    on the thermoelastic properties of high-pressure minerals. Physics of the
    Earth and Planetary Interiors, 96, 85-112.
    doi:10.1016/0031-9201(96)03143-3
    """

    def thermal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye vibrational thermal energy in J mol^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("thermal_energy", volumes, temperatures)
        V = validate_volume(V)
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise ValueError("Temperature must be finite and greater than zero")
        try:
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), temperatures
            )
        except ValueError as error:
            raise ValueError("V and T must have broadcast-compatible shapes") from error

        ratio = self.characteristic_temperature(volumes) / temperatures
        energy = 3.0 * self.n * R * temperatures * _debye_function_3(ratio)
        if energy.ndim == 0:
            return float(energy)
        return energy

    def thermal_entropy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye vibrational entropy in J mol^-1 K^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("thermal_entropy", volumes, temperatures)
        volumes, temperatures = self._broadcast_state(V, T)
        ratio = self.characteristic_temperature(volumes) / temperatures
        log_term = np.log(-np.expm1(-ratio))
        entropy = self.n * R * (4.0 * _debye_function_3(ratio) - 3.0 * log_term)
        return self._scalar_or_array(np.asarray(entropy, dtype=float))


class MieGruneisenEinstein(_MieGruneisenBase):
    """Mie-Gruneisen-Einstein thermal equation of state.

    This uses the same Mie-Gruneisen thermal-pressure relation and volume
    dependence as :class:`MieGruneisenDebye`, with a single-frequency
    Einstein model for the vibrational energy.

    Reference
    ---------
    Dorogokupets, P. I. (2010). P-V-T equations of state of MgO and
    thermodynamics. Physics and Chemistry of Minerals, 37, 677-684.
    doi:10.1007/s00269-010-0367-2
    """

    def thermal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Einstein vibrational thermal energy in J mol^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("thermal_energy", volumes, temperatures)
        V = validate_volume(V)
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise ValueError("Temperature must be finite and greater than zero")
        try:
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), temperatures
            )
        except ValueError as error:
            raise ValueError("V and T must have broadcast-compatible shapes") from error

        theta = self.characteristic_temperature(volumes)
        ratio = theta / temperatures
        occupation = np.exp(-ratio) / (-np.expm1(-ratio))
        energy = 3.0 * self.n * R * theta * occupation
        if energy.ndim == 0:
            return float(energy)
        return energy

    def thermal_entropy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Einstein vibrational entropy in J mol^-1 K^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("thermal_entropy", volumes, temperatures)
        volumes, temperatures = self._broadcast_state(V, T)
        ratio = self.characteristic_temperature(volumes) / temperatures
        occupation = np.exp(-ratio) / (-np.expm1(-ratio))
        entropy = 3.0 * self.n * R * (ratio * occupation - np.log(-np.expm1(-ratio)))
        return self._scalar_or_array(np.asarray(entropy, dtype=float))


def _debye_function_3(x: NumericType) -> NumericType:
    """Return the third-order Debye function with stable limiting forms."""
    values = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Debye-function arguments must be finite and positive")

    def evaluate(value: float) -> float:
        if value < 1.0e-3:
            return 1.0 - 3.0 * value / 8.0 + value**2 / 20.0 - value**4 / 1680.0
        if value > 150.0:
            return np.pi**4 / (5.0 * value**3)
        integral = quad(lambda y: y**3 / np.expm1(y), 0.0, value)[0]
        return 3.0 * integral / value**3

    result = np.array([evaluate(float(value)) for value in values.flat]).reshape(
        values.shape
    )
    if result.ndim == 0:
        return float(result)
    return result
