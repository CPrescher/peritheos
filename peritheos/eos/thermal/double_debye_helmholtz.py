"""Vinet plus double-Debye Helmholtz equation of state."""

from __future__ import annotations

import numpy as np
from scipy.constants import R

from peritheos.eos import (
    NumericType,
    ThermalEOS,
    solve_temperature,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)
from peritheos.eos.rt import Vinet
from peritheos.errors import (
    ConfigurationError,
    EosNumericalError,
    EosValidationError,
)

from .mie_gruneisen import _debye_function_3


class DoubleDebyeHelmholtz(ThermalEOS):
    """Full Vinet/double-Debye Helmholtz equation of state.

    This is a generic implementation of

    ``F(V,T) = E_cold(V) + F_ion(V,T) + F_anh(V,T)``.

    ``rt_eos`` must be a :class:`~peritheos.eos.rt.Vinet` object.  Unlike the
    reference-isotherm objects used by most :class:`ThermalEOS` subclasses, it
    represents the classical, motionless-ion 0 K cold curve.  The ionic term
    is unreferenced and includes zero-point energy and pressure.

    The three characteristic temperatures (the two Debye cutoffs and their
    first phonon moment) independently follow

    ``theta(V) = theta0 * (V / Vp)**(-b) * exp(a * (Vp - V))``.

    Their Gruneisen parameters are consequently ``gamma(V) = a*V + b``.
    The double-Debye weights are chosen so that the weighted first moment is
    ``theta_1``.  Their volume derivatives are retained in the pressure.

    The anharmonic term is

    ``F_anh = -n*R*alpha0*(V/Ve)**kappa*T**2/2``.

    Here ``alpha0`` uses the Benedict normalization, for which the associated
    heat-capacity term is ``n*R*alpha(V)*T``.  The factor ``1/2`` is therefore
    a coefficient convention, not additional physics: a free-energy
    coefficient defined directly by ``F_anh = -n*R*a(V)*T**2`` satisfies
    ``alpha(V) = 2*a(V)``.

    Parameters use Peritheos thermal units: ``V``, ``Vp``, and ``Ve`` are in
    J bar^-1 mol^-1, ``a_*`` in (J bar^-1 mol^-1)^-1, temperatures in K,
    energy in J mol^-1, and pressure in GPa.  ``phi0`` is the cold energy at
    ``rt_eos.V0`` in J mol^-1.  ``n`` is the number of atoms per formula unit.

    Notes
    -----
    This formulation is the solid free energy used by Benedict et al. for
    diamond, BC8, and simple-cubic carbon.  Their Table I coefficients are a
    parameter set for this generic class, not defaults of the model.

    Reference
    ---------
    Benedict, L. X. et al. (2014), Physical Review B 89, 224109, equations
    (3)--(7) and Table I. doi:10.1103/PhysRevB.89.224109
    """

    # Used only to choose the nearest branch in the inherited numerical
    # temperature inversion.  It is not a thermodynamic reference isotherm.
    Tr = 300.0

    def __init__(
        self,
        rt_eos: Vinet,
        Vp: float,
        theta_a0: float,
        a_a: float,
        b_a: float,
        theta_b0: float,
        a_b: float,
        b_b: float,
        theta_1_0: float,
        a_1: float,
        b_1: float,
        n: float = 1.0,
        alpha0: float = 0.0,
        Ve: float = 1.0,
        kappa: float = 0.0,
        phi0: float = 0.0,
    ) -> None:
        if not isinstance(rt_eos, Vinet):
            raise ConfigurationError("rt_eos must be a Vinet cold curve")
        super().__init__(rt_eos)
        self.Vp = validate_positive_scalar(Vp, "Vp")
        self.theta_a0 = validate_positive_scalar(theta_a0, "theta_a0")
        self.a_a = validate_finite_scalar(a_a, "a_a")
        self.b_a = validate_finite_scalar(b_a, "b_a")
        self.theta_b0 = validate_positive_scalar(theta_b0, "theta_b0")
        self.a_b = validate_finite_scalar(a_b, "a_b")
        self.b_b = validate_finite_scalar(b_b, "b_b")
        self.theta_1_0 = validate_positive_scalar(theta_1_0, "theta_1_0")
        self.a_1 = validate_finite_scalar(a_1, "a_1")
        self.b_1 = validate_finite_scalar(b_1, "b_1")
        self.n = validate_positive_scalar(n, "n")
        self.alpha0 = validate_finite_scalar(alpha0, "alpha0")
        if self.alpha0 < 0.0:
            raise EosValidationError("alpha0 must not be negative")
        self.Ve = validate_positive_scalar(Ve, "Ve")
        self.kappa = validate_finite_scalar(kappa, "kappa")
        self.phi0 = validate_finite_scalar(phi0, "phi0")

    @staticmethod
    def _state(V: NumericType, T: NumericType) -> tuple[np.ndarray, np.ndarray]:
        volumes = np.asarray(validate_volume(V), dtype=float)
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures < 0.0):
            raise EosValidationError("Temperature must be finite and non-negative")
        try:
            broadcast_volumes, broadcast_temperatures = np.broadcast_arrays(
                volumes, temperatures
            )
            return broadcast_volumes, broadcast_temperatures
        except ValueError as error:
            raise EosValidationError(
                "V and T must have broadcast-compatible shapes"
            ) from error

    @staticmethod
    def _result(values: NumericType) -> NumericType:
        result = np.asarray(values, dtype=float)
        if result.ndim == 0:
            return float(result)
        return result

    def _temperature_law(
        self, V: NumericType, theta0: float, a: float, b: float
    ) -> tuple[np.ndarray, np.ndarray]:
        volumes = np.asarray(validate_volume(V), dtype=float)
        logarithmic_ratio = np.log(volumes / self.Vp)
        theta = theta0 * np.exp(-b * logarithmic_ratio + a * (self.Vp - volumes))
        gamma = a * volumes + b
        if not np.all(np.isfinite(theta)) or np.any(theta <= 0.0):
            raise EosNumericalError("Debye temperature is not finite and positive")
        return theta, gamma

    def debye_temperatures(
        self, V: NumericType
    ) -> tuple[NumericType, NumericType, NumericType]:
        """Return ``(theta_a, theta_b, theta_1)`` in K."""
        theta_a, _ = self._temperature_law(V, self.theta_a0, self.a_a, self.b_a)
        theta_b, _ = self._temperature_law(V, self.theta_b0, self.a_b, self.b_b)
        theta_1, _ = self._temperature_law(V, self.theta_1_0, self.a_1, self.b_1)
        return (
            self._result(theta_a),
            self._result(theta_b),
            self._result(theta_1),
        )

    def _mode_terms(self, V: NumericType) -> tuple[np.ndarray, ...]:
        volumes = np.asarray(validate_volume(V), dtype=float)
        theta_a, gamma_a = self._temperature_law(
            volumes, self.theta_a0, self.a_a, self.b_a
        )
        theta_b, gamma_b = self._temperature_law(
            volumes, self.theta_b0, self.a_b, self.b_b
        )
        theta_1, gamma_1 = self._temperature_law(
            volumes, self.theta_1_0, self.a_1, self.b_1
        )

        denominator = theta_b - theta_a
        numerator = theta_b - theta_1
        scale = np.maximum.reduce((theta_a, theta_b, theta_1))
        regular = np.abs(denominator) > 1.0e-10 * scale

        safe_denominator = np.where(regular, denominator, 1.0)
        weight_a = numerator / safe_denominator
        theta_a_prime = -gamma_a * theta_a / volumes
        theta_b_prime = -gamma_b * theta_b / volumes
        theta_1_prime = -gamma_1 * theta_1 / volumes
        weight_a_prime = (
            (theta_b_prime - theta_1_prime) * safe_denominator
            - numerator * (theta_b_prime - theta_a_prime)
        ) / safe_denominator**2

        gamma_denominator = gamma_b - gamma_a
        gamma_regular = np.abs(gamma_denominator) > 1.0e-12
        limiting_weight = np.divide(
            gamma_b - gamma_1,
            gamma_denominator,
            out=np.full_like(gamma_denominator, 0.5),
            where=gamma_regular,
        )
        weight_a = np.where(regular, weight_a, limiting_weight)
        weight_a_prime = np.where(regular, weight_a_prime, 0.0)
        return (
            volumes,
            theta_a,
            gamma_a,
            theta_b,
            gamma_b,
            theta_1,
            weight_a,
            weight_a_prime,
        )

    def double_debye_weights(self, V: NumericType) -> tuple[NumericType, NumericType]:
        """Return the volume-dependent weights of modes A and B."""
        weight_a = self._mode_terms(V)[6]
        return self._result(weight_a), self._result(1.0 - weight_a)

    @staticmethod
    def _single_debye_free_energy(theta: np.ndarray, T: np.ndarray) -> np.ndarray:
        """Return one mole-of-atoms Debye free energy, including zero point."""
        theta, temperatures = np.broadcast_arrays(theta, T)
        result = 9.0 * R * theta / 8.0
        positive = temperatures > 0.0
        if np.any(positive):
            ratio = theta[positive] / temperatures[positive]
            thermal = (
                R
                * temperatures[positive]
                * (3.0 * np.log(-np.expm1(-ratio)) - _debye_function_3(ratio))
            )
            result = np.array(result, copy=True)
            result[positive] += thermal
        return result

    @staticmethod
    def _single_debye_internal_energy(theta: np.ndarray, T: np.ndarray) -> np.ndarray:
        """Return one mole-of-atoms Debye internal energy, including zero point."""
        theta, temperatures = np.broadcast_arrays(theta, T)
        result = 9.0 * R * theta / 8.0
        positive = temperatures > 0.0
        if np.any(positive):
            ratio = theta[positive] / temperatures[positive]
            result = np.array(result, copy=True)
            result[positive] += (
                3.0 * R * temperatures[positive] * _debye_function_3(ratio)
            )
        return result

    @staticmethod
    def _single_debye_heat_capacity(theta: np.ndarray, T: np.ndarray) -> np.ndarray:
        """Return one mole-of-atoms Debye constant-volume heat capacity."""
        theta, temperatures = np.broadcast_arrays(theta, T)
        result = np.zeros_like(temperatures)
        positive = temperatures > 0.0
        if np.any(positive):
            ratio = theta[positive] / temperatures[positive]
            occupation_term = ratio * np.exp(-ratio) / (-np.expm1(-ratio))
            result[positive] = (
                3.0 * R * (4.0 * _debye_function_3(ratio) - 3.0 * occupation_term)
            )
        return result

    def cold_energy(self, V: NumericType) -> NumericType:
        """Return the Vinet cold-curve energy in J mol^-1."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        delta = self.rt_eos.K0_prime - 1.0
        x = np.cbrt(volumes / self.rt_eos.V0)
        if abs(delta) < 1.0e-7:
            y = x - 1.0
            reduced = 1.125 * y**2 - 1.125 * delta * y**3
        else:
            X = 1.5 * delta * (x - 1.0)
            reduced = (-np.expm1(-X) - X * np.exp(-X)) / delta**2
        energy = self.phi0 + 4.0 * self.rt_eos.V0 * self.rt_eos.K0 * 1.0e4 * reduced
        return self._result(energy)

    def zero_point_energy(self, V: NumericType) -> NumericType:
        """Return the weighted double-Debye zero-point energy in J mol^-1."""
        terms = self._mode_terms(V)
        theta_a, theta_b, weight_a = terms[1], terms[3], terms[6]
        energy = (
            self.n * 9.0 * R / 8.0 * (weight_a * theta_a + (1.0 - weight_a) * theta_b)
        )
        return self._result(energy)

    def ion_helmholtz_free_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the double-Debye ionic free energy, including zero point."""
        volumes, temperatures = self._state(V, T)
        terms = self._mode_terms(volumes)
        theta_a, theta_b, weight_a = terms[1], terms[3], terms[6]
        free_a = self._single_debye_free_energy(theta_a, temperatures)
        free_b = self._single_debye_free_energy(theta_b, temperatures)
        result = self.n * (weight_a * free_a + (1.0 - weight_a) * free_b)
        return self._result(result)

    def anharmonic_coefficient(self, V: NumericType) -> NumericType:
        """Return ``alpha(V)`` in K^-1."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        result = self.alpha0 * np.exp(self.kappa * np.log(volumes / self.Ve))
        return self._result(result)

    def anharmonic_helmholtz_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return the anharmonic/electronic ``T^2`` free energy in J mol^-1."""
        volumes, temperatures = self._state(V, T)
        alpha = np.asarray(self.anharmonic_coefficient(volumes), dtype=float)
        return self._result(-0.5 * self.n * R * alpha * temperatures**2)

    def helmholtz_free_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the complete Helmholtz free energy in J mol^-1."""
        volumes, temperatures = self._state(V, T)
        result = (
            np.asarray(self.cold_energy(volumes), dtype=float)
            + np.asarray(
                self.ion_helmholtz_free_energy(volumes, temperatures), dtype=float
            )
            + np.asarray(
                self.anharmonic_helmholtz_free_energy(volumes, temperatures),
                dtype=float,
            )
        )
        return self._result(result)

    def ion_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return ionic pressure, including zero-point pressure, in GPa."""
        volumes, temperatures = self._state(V, T)
        terms = self._mode_terms(volumes)
        (
            _,
            theta_a,
            gamma_a,
            theta_b,
            gamma_b,
            _,
            weight_a,
            weight_a_prime,
        ) = terms
        free_a = self._single_debye_free_energy(theta_a, temperatures)
        free_b = self._single_debye_free_energy(theta_b, temperatures)
        energy_a = self._single_debye_internal_energy(theta_a, temperatures)
        energy_b = self._single_debye_internal_energy(theta_b, temperatures)
        pressure_bar = self.n * (
            (weight_a * gamma_a * energy_a + (1.0 - weight_a) * gamma_b * energy_b)
            / volumes
            - weight_a_prime * (free_a - free_b)
        )
        return self._result(pressure_bar / 1.0e4)

    def anharmonic_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return pressure from the anharmonic/electronic term in GPa."""
        volumes, temperatures = self._state(V, T)
        alpha = np.asarray(self.anharmonic_coefficient(volumes), dtype=float)
        pressure_bar = 0.5 * self.n * R * self.kappa * alpha * temperatures**2 / volumes
        return self._result(pressure_bar / 1.0e4)

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return all non-cold pressure, including zero point, in GPa.

        This is an absolute contribution, not a difference from ``Tr``.
        """
        result = np.asarray(self.ion_pressure(V, T), dtype=float) + np.asarray(
            self.anharmonic_pressure(V, T), dtype=float
        )
        return self._result(result)

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Return total non-cold constant-volume heat capacity in J mol^-1 K^-1."""
        volumes, temperatures = self._state(V, T)
        terms = self._mode_terms(volumes)
        theta_a, theta_b, weight_a = terms[1], terms[3], terms[6]
        result = self.n * (
            weight_a * self._single_debye_heat_capacity(theta_a, temperatures)
            + (1.0 - weight_a) * self._single_debye_heat_capacity(theta_b, temperatures)
            + R
            * np.asarray(self.anharmonic_coefficient(volumes), dtype=float)
            * temperatures
        )
        return self._result(result)

    def dac_thermal_pressure(
        self, V: NumericType, T: NumericType, f_dac: float
    ) -> NumericType:
        """Return the confined pressure increment above the ``Tr`` isotherm."""
        f_dac = self._validate_f_dac(f_dac)
        increment = np.asarray(self.pressure(V, T), dtype=float) - np.asarray(
            self.pressure(V, self.Tr), dtype=float
        )
        return self._result(f_dac * increment)

    def temperature_from_volumes(
        self,
        V_ambient: NumericType,
        V_heated: NumericType,
        *,
        f_dac: float,
    ) -> NumericType:
        """Infer temperature from ambient and heated volumes in a DAC.

        For this absolute Helmholtz EOS, the pressure increment caused by
        heating at fixed volume is

        ``Delta P_th(V,T) = P(V,T) - P(V,Tr)``.

        This removes the zero-point and finite reference-temperature pressure
        without discarding them from either total pressure.  The empirical DAC
        boundary condition then reduces to

        ``Delta P_th(V_heated,T) =``
        ``(P(V_ambient,Tr) - P(V_heated,Tr)) / (1 - f_dac)``.

        ``f_dac`` is the fraction of that reference-relative thermal-pressure
        increment appearing as ``P_hot - P_ambient``.  It must lie in
        ``[0, 1)``.  Only roots at or above ``Tr`` are accepted.
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

        ambient_reference_pressures = np.asarray(
            self.pressure(ambient_volumes, self.Tr), dtype=float
        )
        heated_reference_pressures = np.asarray(
            self.pressure(heated_volumes, self.Tr), dtype=float
        )
        target_increments = (
            ambient_reference_pressures - heated_reference_pressures
        ) / (1.0 - f_dac)
        if np.any(target_increments < 0.0):
            raise EosValidationError(
                "The volume pair implies a temperature below the reference "
                "temperature, not a heated state"
            )

        temperatures = np.array(
            [
                solve_temperature(
                    lambda temperature, volume=float(heated_volume): (
                        self.pressure(volume, temperature)
                        - self.pressure(volume, self.Tr)
                    ),
                    float(target_increment),
                    self.Tr,
                )
                for target_increment, heated_volume in zip(
                    target_increments.flat, heated_volumes.flat
                )
            ]
        ).reshape(target_increments.shape)
        temperature_tolerance = 1.0e-10 * max(1.0, self.Tr)
        if np.any(temperatures < self.Tr - temperature_tolerance):
            raise EosValidationError(
                "The volume pair implies a temperature below the reference "
                "temperature, not a heated state"
            )
        return self._result(temperatures)


class DoubleDebyeLogMomentHelmholtz(DoubleDebyeHelmholtz):
    """Vinet/double-Debye Helmholtz EOS constrained by ``theta_0``.

    The two Debye-mode weights satisfy the logarithmic phonon-moment
    constraint

    ``xi_a = log(theta_b / theta_0) / log(theta_b / theta_a)``.

    Each characteristic temperature follows

    ``theta(V) = theta_ref * (V / Vp)**(-b) * exp(a * (Vp - V))``.

    The volume derivative of the weights is included in the pressure.  The
    optional Correa anharmonic term is volume independent and has the form
    ``F_anh = -n*R*anharmonic_a*T**2``; it therefore changes energy and heat
    capacity but not pressure.  Correa's coefficient is the direct
    free-energy coefficient; in the later Benedict ``-alpha*T**2/2``
    normalization, the equivalent value is ``alpha = 2*anharmonic_a``.

    Parameters use the same molar units as :class:`DoubleDebyeHelmholtz`.

    Reference
    ---------
    Correa, A. A. et al. (2008), Physical Review B 78, 024101, equations
    (2)--(7) and (13)--(18), and Table I.
    doi:10.1103/PhysRevB.78.024101
    """

    def __init__(
        self,
        rt_eos: Vinet,
        Vp: float,
        theta_a0: float,
        a_a: float,
        b_a: float,
        theta_b0: float,
        a_b: float,
        b_b: float,
        theta_0_0: float,
        a_0: float,
        b_0: float,
        n: float = 1.0,
        anharmonic_a: float = 0.0,
        phi0: float = 0.0,
    ) -> None:
        if not isinstance(rt_eos, Vinet):
            raise ConfigurationError("rt_eos must be a Vinet cold curve")
        ThermalEOS.__init__(self, rt_eos)
        self.Vp = validate_positive_scalar(Vp, "Vp")
        self.theta_a0 = validate_positive_scalar(theta_a0, "theta_a0")
        self.a_a = validate_finite_scalar(a_a, "a_a")
        self.b_a = validate_finite_scalar(b_a, "b_a")
        self.theta_b0 = validate_positive_scalar(theta_b0, "theta_b0")
        self.a_b = validate_finite_scalar(a_b, "a_b")
        self.b_b = validate_finite_scalar(b_b, "b_b")
        self.theta_0_0 = validate_positive_scalar(theta_0_0, "theta_0_0")
        self.a_0 = validate_finite_scalar(a_0, "a_0")
        self.b_0 = validate_finite_scalar(b_0, "b_0")
        self.n = validate_positive_scalar(n, "n")
        self.anharmonic_a = validate_finite_scalar(anharmonic_a, "anharmonic_a")
        if self.anharmonic_a < 0.0:
            raise EosValidationError("anharmonic_a must not be negative")
        self.phi0 = validate_finite_scalar(phi0, "phi0")

    def debye_temperatures(
        self, V: NumericType
    ) -> tuple[NumericType, NumericType, NumericType]:
        """Return ``(theta_a, theta_b, theta_0)`` in K."""
        theta_a, _ = self._temperature_law(V, self.theta_a0, self.a_a, self.b_a)
        theta_b, _ = self._temperature_law(V, self.theta_b0, self.a_b, self.b_b)
        theta_0, _ = self._temperature_law(V, self.theta_0_0, self.a_0, self.b_0)
        return (
            self._result(theta_a),
            self._result(theta_b),
            self._result(theta_0),
        )

    def _mode_terms(self, V: NumericType) -> tuple[np.ndarray, ...]:
        volumes = np.asarray(validate_volume(V), dtype=float)
        theta_a, gamma_a = self._temperature_law(
            volumes, self.theta_a0, self.a_a, self.b_a
        )
        theta_b, gamma_b = self._temperature_law(
            volumes, self.theta_b0, self.a_b, self.b_b
        )
        theta_0, gamma_0 = self._temperature_law(
            volumes, self.theta_0_0, self.a_0, self.b_0
        )

        denominator = np.log(theta_b / theta_a)
        numerator = np.log(theta_b / theta_0)
        regular = np.abs(denominator) > 1.0e-10
        safe_denominator = np.where(regular, denominator, 1.0)
        weight_a = numerator / safe_denominator

        numerator_prime = (gamma_0 - gamma_b) / volumes
        denominator_prime = (gamma_a - gamma_b) / volumes
        weight_a_prime = (
            numerator_prime * safe_denominator - numerator * denominator_prime
        ) / safe_denominator**2

        gamma_denominator = gamma_b - gamma_a
        gamma_regular = np.abs(gamma_denominator) > 1.0e-12
        limiting_weight = np.divide(
            gamma_b - gamma_0,
            gamma_denominator,
            out=np.full_like(gamma_denominator, 0.5),
            where=gamma_regular,
        )
        weight_a = np.where(regular, weight_a, limiting_weight)
        weight_a_prime = np.where(regular, weight_a_prime, 0.0)
        return (
            volumes,
            theta_a,
            gamma_a,
            theta_b,
            gamma_b,
            theta_0,
            weight_a,
            weight_a_prime,
        )

    def anharmonic_coefficient(self, V: NumericType) -> NumericType:
        """Return Correa's volume-independent ``a`` coefficient in K^-1."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        result = np.full_like(volumes, self.anharmonic_a, dtype=float)
        return self._result(result)

    def anharmonic_helmholtz_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return ``-n*R*a*T^2`` in J mol^-1."""
        _, temperatures = self._state(V, T)
        return self._result(-self.n * R * self.anharmonic_a * temperatures**2)

    def anharmonic_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return zero because Correa's anharmonic coefficient is constant."""
        volumes, _ = self._state(V, T)
        return self._result(np.zeros_like(volumes))

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Return non-cold constant-volume heat capacity in J mol^-1 K^-1."""
        volumes, temperatures = self._state(V, T)
        terms = self._mode_terms(volumes)
        theta_a, theta_b, weight_a = terms[1], terms[3], terms[6]
        result = self.n * (
            weight_a * self._single_debye_heat_capacity(theta_a, temperatures)
            + (1.0 - weight_a) * self._single_debye_heat_capacity(theta_b, temperatures)
            + 2.0 * R * self.anharmonic_a * temperatures
        )
        return self._result(result)
