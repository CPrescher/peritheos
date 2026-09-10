"""Dewaele et al. (2006) hcp-Fe thermal pressure scale."""

from __future__ import annotations

import numpy as np
from scipy.constants import R

from peritheos.errors import ConfigurationError, EosNumericalError, EosValidationError

from .. import (
    EosBase,
    NumericType,
    ThermalEOS,
    _native_for_exact_model,
    _native_thermal_evaluate,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)
from .mie_gruneisen import _debye_function_3


class Dewaele2006(ThermalEOS):
    """Single-Debye hcp-Fe pressure scale of Dewaele et al. (2006).

    Equations (1)--(2) combine a complete reference-temperature Vinet
    isotherm with a single-Debye quasiharmonic pressure and quadratic
    intrinsic-anharmonic and electronic contributions. The Gruneisen and
    Debye-temperature laws are the simplified Dorogokupets--Oganov form used
    by the source.

    Volume uses the Peritheos molar convention,
    ``J bar^-1 mol^-1 == cm^3 mol^-1 / 10``. The ``anharmonic_a`` and
    ``electronic_e`` coefficients are supplied in K^-1, matching Dewaele's
    printed equation rather than the ``10^-6 K^-1`` table convention used by
    :class:`DorogokupetsOganov2007`.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        gamma_inf: float,
        beta: float,
        anharmonic_a: float,
        anharmonic_m: float,
        electronic_e: float,
        electronic_g: float,
        n: float,
    ) -> None:
        if not isinstance(rt_eos, EosBase):
            raise ConfigurationError("rt_eos must be an equation of state")
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.theta0 = validate_positive_scalar(theta0, "theta0")
        self.gamma0 = validate_positive_scalar(gamma0, "gamma0")
        self.gamma_inf = validate_positive_scalar(gamma_inf, "gamma_inf")
        if self.gamma_inf > self.gamma0:
            raise EosValidationError("gamma_inf must not exceed gamma0")
        self.beta = validate_positive_scalar(beta, "beta")
        self.anharmonic_a = validate_finite_scalar(anharmonic_a, "anharmonic_a")
        if self.anharmonic_a < 0:
            raise EosValidationError("anharmonic_a must not be negative")
        self.anharmonic_m = validate_positive_scalar(anharmonic_m, "anharmonic_m")
        self.electronic_e = validate_finite_scalar(electronic_e, "electronic_e")
        if self.electronic_e < 0:
            raise EosValidationError("electronic_e must not be negative")
        self.electronic_g = validate_positive_scalar(electronic_g, "electronic_g")
        self.n = validate_positive_scalar(n, "n")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is Dewaele2006:
            from peritheos import _rust

            self._native = _rust.ThermalEos.dewaele_2006(
                reference_native,
                self.Tr,
                self.theta0,
                self.gamma0,
                self.gamma_inf,
                self.beta,
                self.anharmonic_a,
                self.anharmonic_m,
                self.electronic_e,
                self.electronic_g,
                self.n,
            )

    def gruneisen_parameter(
        self, V: NumericType, T: NumericType | None = None
    ) -> NumericType:
        """Return ``gamma_inf + (gamma0-gamma_inf)*(V/V0)**beta``."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        if hasattr(self, "_native"):
            temperatures = self.Tr if T is None else T
            volumes, temperatures = self._broadcast_state(volumes, temperatures)
            return _native_thermal_evaluate(
                self._native, "gruneisen_parameter", volumes, temperatures
            )
        ratio = volumes / self.rt_eos.V0
        result = self.gamma_inf + (self.gamma0 - self.gamma_inf) * ratio**self.beta
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def characteristic_temperature(self, V: NumericType) -> NumericType:
        """Return the integrated Debye-temperature law in kelvin."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native,
                "characteristic_temperature",
                volumes,
                np.full_like(volumes, self.Tr),
            )
        ratio = volumes / self.rt_eos.V0
        result = (
            self.theta0
            * ratio ** (-self.gamma_inf)
            * np.exp(
                (self.gamma0 - self.gamma_inf) / self.beta * (1.0 - ratio**self.beta)
            )
        )
        if not np.all(np.isfinite(result)):
            raise EosNumericalError("Characteristic temperature is not finite")
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def _vibrational_energy(self, V: NumericType, T: NumericType) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        theta = np.asarray(self.characteristic_temperature(volumes), dtype=float)
        result = (
            3.0 * self.n * R * temperatures * _debye_function_3(theta / temperatures)
        )
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def vibrational_pressure_increment(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return the reference-relative single-Debye pressure in GPa."""
        volumes, temperatures = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native,
                "vibrational_pressure_increment",
                volumes,
                temperatures,
            )
        energy_difference = np.asarray(
            self._vibrational_energy(volumes, temperatures), dtype=float
        ) - np.asarray(self._vibrational_energy(volumes, self.Tr), dtype=float)
        result = self.gruneisen_parameter(volumes) * energy_difference / volumes / 1.0e4
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def _quadratic_pressure_increment(
        self,
        V: NumericType,
        T: NumericType,
        coefficient: float,
        exponent: float,
    ) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        ratio = volumes / self.rt_eos.V0
        result = (
            1.5
            * self.n
            * R
            * coefficient
            * exponent
            * ratio**exponent
            * (temperatures**2 - self.Tr**2)
            / volumes
            / 1.0e4
        )
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def anharmonic_pressure_increment(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return the reference-relative intrinsic-anharmonic pressure in GPa."""
        volumes, temperatures = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native,
                "anharmonic_pressure_increment",
                volumes,
                temperatures,
            )
        return self._quadratic_pressure_increment(
            volumes, temperatures, self.anharmonic_a, self.anharmonic_m
        )

    def electronic_pressure_increment(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return the reference-relative electronic pressure in GPa."""
        volumes, temperatures = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native,
                "electronic_pressure_increment",
                volumes,
                temperatures,
            )
        return self._quadratic_pressure_increment(
            volumes, temperatures, self.electronic_e, self.electronic_g
        )

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return total pressure above the reference-temperature isotherm."""
        volumes, temperatures = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        result = (
            np.asarray(
                self.vibrational_pressure_increment(volumes, temperatures), dtype=float
            )
            + np.asarray(
                self.anharmonic_pressure_increment(volumes, temperatures), dtype=float
            )
            + np.asarray(
                self.electronic_pressure_increment(volumes, temperatures), dtype=float
            )
        )
        return self._scalar_or_array(result)
