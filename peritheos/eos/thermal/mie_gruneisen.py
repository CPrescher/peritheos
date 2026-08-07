"""Mie-Gruneisen thermal equations of state."""

from abc import ABC, abstractmethod

import numpy as np
from scipy.constants import R
from scipy.integrate import quad

from peritheos.eos import (
    EosBase,
    NumericType,
    ThermalEOS,
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

    def gruneisen_parameter(self, V: NumericType) -> NumericType:
        """Return ``gamma(V) = gamma0 * (V / V0)**q``."""
        V = validate_volume(V)
        return self.gamma0 * np.exp(self.q * np.log(V / self.rt_eos.V0))

    def characteristic_temperature(self, V: NumericType) -> NumericType:
        """Return the characteristic lattice temperature at volume *V*.

        Its volume dependence is thermodynamically consistent with
        ``gamma = -d(log(theta)) / d(log(V))``.
        """
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
            self.gruneisen_parameter(volumes)
            * energy_difference
            / volumes
            / 10000.0
        )
        if pressure.ndim == 0:
            return float(pressure)
        return pressure


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
