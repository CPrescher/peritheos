"""Dorogokupets--Oganov (2007) semiempirical thermal equation of state."""

from __future__ import annotations

import numpy as np
from scipy.constants import R

from peritheos.errors import EosValidationError

from .. import (
    EosBase,
    NumericType,
    ThermalEOS,
    _native_for_exact_model,
    _native_thermal_evaluate,
    validate_finite_scalar,
    validate_positive_scalar,
)


class DorogokupetsOganov2007(ThermalEOS):
    """Four-oscillator Helmholtz model from Dorogokupets and Oganov (2007).

    The supplied reference EOS is the published 298.15 K isotherm. Thermal
    pressure is evaluated from equations (7)--(14) and explicitly subtracts
    the contribution at ``Tr``. Volume uses the Peritheos molar convention,
    ``J bar^-1 mol^-1 == cm^3 mol^-1 / 10``.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta_B1: float,
        d_B1: float,
        m_B1: float,
        theta_B2: float,
        d_B2: float,
        m_B2: float,
        theta_E1: float,
        m_E1: float,
        theta_E2: float,
        m_E2: float,
        gamma0: float,
        gamma_inf: float,
        beta: float,
        anharmonic_a: float,
        anharmonic_m: float,
        electronic_e: float,
        electronic_g: float,
        defect_H: float,
        defect_S: float,
        n: float,
    ):
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.theta_B1 = validate_positive_scalar(theta_B1, "theta_B1")
        self.d_B1 = validate_positive_scalar(d_B1, "d_B1")
        self.m_B1 = validate_positive_scalar(m_B1, "m_B1")
        self.theta_B2 = validate_positive_scalar(theta_B2, "theta_B2")
        self.d_B2 = validate_positive_scalar(d_B2, "d_B2")
        self.m_B2 = validate_positive_scalar(m_B2, "m_B2")
        self.theta_E1 = validate_positive_scalar(theta_E1, "theta_E1")
        self.m_E1 = validate_positive_scalar(m_E1, "m_E1")
        self.theta_E2 = validate_positive_scalar(theta_E2, "theta_E2")
        self.m_E2 = validate_positive_scalar(m_E2, "m_E2")
        self.gamma0 = validate_positive_scalar(gamma0, "gamma0")
        self.gamma_inf = validate_positive_scalar(gamma_inf, "gamma_inf")
        self.beta = validate_positive_scalar(beta, "beta")
        self.anharmonic_a = validate_finite_scalar(
            anharmonic_a, "anharmonic_a"
        )
        self.anharmonic_m = validate_finite_scalar(
            anharmonic_m, "anharmonic_m"
        )
        self.electronic_e = validate_finite_scalar(electronic_e, "electronic_e")
        self.electronic_g = validate_finite_scalar(electronic_g, "electronic_g")
        self.defect_H = validate_positive_scalar(defect_H, "defect_H")
        self.defect_S = validate_finite_scalar(defect_S, "defect_S")
        self.n = validate_positive_scalar(n, "n")
        multiplicity = self.m_B1 + self.m_B2 + self.m_E1 + self.m_E2
        if not np.isclose(multiplicity, 3.0 * self.n, rtol=0.0, atol=5.0e-3):
            raise EosValidationError(
                "Oscillator multiplicities must sum to three times the atom count"
            )
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is DorogokupetsOganov2007:
            from peritheos import _rust

            self._native = _rust.ThermalEos.dorogokupets_oganov_2007(
                reference_native,
                self.Tr,
                self.theta_B1,
                self.d_B1,
                self.m_B1,
                self.theta_B2,
                self.d_B2,
                self.m_B2,
                self.theta_E1,
                self.m_E1,
                self.theta_E2,
                self.m_E2,
                self.gamma0,
                self.gamma_inf,
                self.beta,
                self.anharmonic_a,
                self.anharmonic_m,
                self.electronic_e,
                self.electronic_g,
                self.defect_H,
                self.defect_S,
                self.n,
            )

    def _gamma(self, ratio: np.ndarray) -> np.ndarray:
        return self.gamma_inf + (self.gamma0 - self.gamma_inf) * ratio**self.beta

    def _theta(self, theta0: float, ratio: np.ndarray) -> np.ndarray:
        return theta0 * ratio ** (-self.gamma_inf) * np.exp(
            (self.gamma0 - self.gamma_inf)
            / self.beta
            * (1.0 - ratio**self.beta)
        )

    @staticmethod
    def _einstein_occupation(theta: np.ndarray, temperature: np.ndarray) -> np.ndarray:
        exponent = theta / temperature
        decay = np.exp(-exponent)
        return decay / (-np.expm1(-exponent))

    @classmethod
    def _einstein_energy(
        cls, theta: np.ndarray, temperature: np.ndarray
    ) -> np.ndarray:
        return theta * (0.5 + cls._einstein_occupation(theta, temperature))

    @staticmethod
    def _bose_energy(
        theta: np.ndarray, temperature: np.ndarray, dispersion: float
    ) -> np.ndarray:
        exponent = dispersion * np.log1p(theta / (temperature * dispersion))
        decay = np.exp(-exponent)
        occupation = decay / (-np.expm1(-exponent))
        return theta * (dispersion - 1.0) / (2.0 * dispersion) + (
            temperature
            * theta
            * dispersion
            * occupation
            / (temperature * dispersion + theta)
        )

    @classmethod
    def _anharmonic_bracket(
        cls, theta: np.ndarray, temperature: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        occupation = cls._einstein_occupation(theta, temperature)
        fluctuation = occupation * (occupation + 1.0)
        energy = theta * (0.5 + occupation)
        energy_derivative = 0.5 + occupation - theta / temperature * fluctuation
        bracket = energy**2 + 2.0 * theta**2 * fluctuation
        derivative = (
            2.0 * energy * energy_derivative
            + 4.0 * theta * fluctuation
            - 2.0
            * theta**2
            / temperature
            * fluctuation
            * (2.0 * occupation + 1.0)
        )
        return bracket, derivative

    def _absolute_nonreference_pressure(
        self, volume: np.ndarray, temperature: np.ndarray
    ) -> np.ndarray:
        ratio = volume / self.rt_eos.V0
        gamma = self._gamma(ratio)
        mode_parameters = (
            (self.theta_B1, self.m_B1, self.d_B1),
            (self.theta_B2, self.m_B2, self.d_B2),
            (self.theta_E1, self.m_E1, None),
            (self.theta_E2, self.m_E2, None),
        )
        quasiharmonic = np.zeros_like(volume)
        anharmonic_derivative = np.zeros_like(volume)
        for theta0, multiplicity, dispersion in mode_parameters:
            theta = self._theta(theta0, ratio)
            if dispersion is None:
                energy = self._einstein_energy(theta, temperature)
            else:
                energy = self._bose_energy(theta, temperature, dispersion)
            quasiharmonic += multiplicity * R * energy * gamma / volume

            bracket, bracket_derivative = self._anharmonic_bracket(
                theta, temperature
            )
            anharmonic_derivative += (
                multiplicity
                * R
                * self.anharmonic_a
                * 1.0e-6
                * ratio**self.anharmonic_m
                / (6.0 * volume)
                * (
                    self.anharmonic_m * bracket
                    - gamma * theta * bracket_derivative
                )
            )

        electronic = (
            1.5
            * self.n
            * R
            * self.electronic_e
            * 1.0e-6
            * self.electronic_g
            * ratio**self.electronic_g
            * temperature**2
            / volume
        )
        defect_exponent = self.defect_S / ratio - self.defect_H / (
            temperature * ratio**2
        )
        defect = (
            1.5
            * self.n
            * R
            * temperature
            * np.exp(defect_exponent)
            * (-self.defect_S / ratio + 2.0 * self.defect_H / (temperature * ratio**2))
            / volume
        )
        pressure = quasiharmonic - anharmonic_derivative + electronic + defect
        if not np.all(np.isfinite(pressure)):
            raise EosValidationError("Thermal pressure is not finite")
        return pressure / 1.0e4

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return pressure relative to the published 298.15 K isotherm."""
        volumes, temperatures = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        pressure = self._absolute_nonreference_pressure(volumes, temperatures)
        reference = self._absolute_nonreference_pressure(
            volumes, np.full_like(temperatures, self.Tr)
        )
        return self._scalar_or_array(np.asarray(pressure - reference, dtype=float))
