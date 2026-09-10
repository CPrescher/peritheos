"""Mie-Gruneisen thermal equations of state."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy.constants import R
from scipy.integrate import quad
from scipy.optimize import brentq

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
from peritheos.errors import (
    ConfigurationError,
    EosNumericalError,
    EosValidationError,
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
            raise ConfigurationError("rt_eos must be an equation of state")
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.theta0 = validate_positive_scalar(theta0, "theta0")
        self.gamma0 = validate_finite_scalar(gamma0, "gamma0")
        self.q = validate_finite_scalar(q, "q")
        self.n = validate_positive_scalar(n, "n")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is MieGruneisenEinstein:
            from peritheos import _rust

            self._native = _rust.ThermalEos.mie_gruneisen_einstein(
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
            raise EosNumericalError("Characteristic temperature is not finite")
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
            raise EosValidationError("Temperature must be finite and greater than zero")
        try:
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), temperatures
            )
        except ValueError as error:
            raise EosValidationError(
                "V and T must have broadcast-compatible shapes"
            ) from error

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

    where ``E_D`` is the Debye vibrational energy and ``gamma(V) = gamma0 *
    (V/V0)**q``. ``thermal_pressure_reference="absolute_zero"`` instead uses
    ``P_th = gamma(V) E_D(V, T) / V`` so that the supplied isothermal EOS is a
    0 K cold curve. ``thermal_pressure_reference="reference_isentrope"``
    subtracts the Debye energy at ``T_S(V) = Tr * theta(V) / theta0``. This
    maps publications whose BM3 term is a reference isentrope rather than an
    isotherm. ``Cvmax`` optionally replaces the Dulong--Petit limit ``3 n R``.
    ``debye_temperature_law`` selects either the conventional
    thermodynamically integrated relation (the default) or the direct
    variable-exponent relation printed by Fei et al. (2007).

    Reference
    ---------
    Jackson, I. & Rigden, S. M. (1996). Analysis of P-V-T data: constraints
    on the thermoelastic properties of high-pressure minerals. Physics of the
    Earth and Planetary Interiors, 96, 85-112.
    doi:10.1016/0031-9201(96)03143-3

    Fei, Y. et al. (2007). Toward an internally consistent pressure scale.
    Proceedings of the National Academy of Sciences, 104, 9182-9186.
    Equation (3), the definition following it, and Table 1.
    doi:10.1073/pnas.0609013104
    """

    _constructor_configuration_names: tuple[str, ...] = (
        "Cvmax",
        "debye_temperature_law",
        "thermal_pressure_reference",
    )
    _DEBYE_TEMPERATURE_LAWS = {"integrated_gruneisen", "variable_exponent"}
    _THERMAL_PRESSURE_REFERENCES = {
        "reference_temperature",
        "reference_isentrope",
        "absolute_zero",
    }

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        q: float,
        n: float,
        debye_temperature_law: str = "integrated_gruneisen",
        thermal_pressure_reference: str = "reference_temperature",
        Cvmax: float | None = None,
    ) -> None:
        super().__init__(rt_eos, Tr, theta0, gamma0, q, n)
        self._cvmax_explicit = Cvmax is not None
        self.Cvmax = (
            3.0 * self.n * R
            if Cvmax is None
            else validate_positive_scalar(Cvmax, "Cvmax")
        )
        if (
            not isinstance(debye_temperature_law, str)
            or debye_temperature_law not in self._DEBYE_TEMPERATURE_LAWS
        ):
            raise EosValidationError(
                "debye_temperature_law must be 'integrated_gruneisen' or "
                "'variable_exponent'"
            )
        self.debye_temperature_law = debye_temperature_law
        if (
            not isinstance(thermal_pressure_reference, str)
            or thermal_pressure_reference not in self._THERMAL_PRESSURE_REFERENCES
        ):
            raise EosValidationError(
                "thermal_pressure_reference must be 'reference_temperature', "
                "'reference_isentrope', or 'absolute_zero'"
            )
        self.thermal_pressure_reference = thermal_pressure_reference
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is MieGruneisenDebye:
            from peritheos import _rust

            self._native = _rust.ThermalEos.mie_gruneisen_debye(
                reference_native,
                self.Tr,
                self.theta0,
                self.gamma0,
                self.q,
                self.n,
                self.debye_temperature_law,
                self.thermal_pressure_reference,
                self.Cvmax,
            )

    def _own_parameter_names(self) -> tuple[str, ...]:
        """Expose ``Cvmax`` as a fit parameter only when explicitly supplied."""
        names = super()._own_parameter_names()
        if self._cvmax_explicit:
            return (*names, "Cvmax")
        return names

    def configuration_values(self) -> dict[str, str | float]:
        """Return non-numeric choices, omitting the default pressure baseline."""
        if type(self) is not MieGruneisenDebye:
            return super().configuration_values()
        configuration = {"debye_temperature_law": self.debye_temperature_law}
        if self.thermal_pressure_reference != "reference_temperature":
            configuration["thermal_pressure_reference"] = (
                self.thermal_pressure_reference
            )
        return configuration

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return referenced or absolute-zero Debye thermal pressure in GPa."""
        if self.thermal_pressure_reference == "reference_temperature":
            return super().thermal_pressure(V, T)
        if self.thermal_pressure_reference == "reference_isentrope":
            volumes, temperatures = self._broadcast_state(V, T)
            reference_temperature = (
                self.Tr * self.characteristic_temperature(volumes) / self.theta0
            )
            energy_difference = np.asarray(
                self.thermal_energy(volumes, temperatures), dtype=float
            ) - np.asarray(
                self.thermal_energy(volumes, reference_temperature), dtype=float
            )
            pressure = (
                self.gruneisen_parameter(volumes) * energy_difference / volumes / 1.0e4
            )
            return self._scalar_or_array(np.asarray(pressure, dtype=float))
        return self.vibrational_pressure(V, T)

    def thermal_pressure_increment(self, V: NumericType, T: NumericType) -> NumericType:
        """Return pressure above the configured ``Tr`` isotherm in GPa."""
        if self.thermal_pressure_reference == "reference_temperature":
            return self.thermal_pressure(V, T)
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate(
                "thermal_pressure_increment", volumes, temperatures
            )
        increment = np.asarray(
            self.vibrational_pressure(V, T), dtype=float
        ) - np.asarray(self.vibrational_pressure(V, self.Tr), dtype=float)
        return self._scalar_or_array(increment)

    def characteristic_temperature(self, V: NumericType) -> NumericType:
        """Return Debye temperature using the selected volume relation."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes = validate_volume(V)
            return self._native_evaluate("characteristic_temperature", volumes, self.Tr)
        if self.debye_temperature_law == "integrated_gruneisen":
            return super().characteristic_temperature(V)

        volumes = np.asarray(validate_volume(V), dtype=float)
        ratio = volumes / self.rt_eos.V0
        gamma = np.asarray(self.gruneisen_parameter(volumes), dtype=float)
        result = self.theta0 * np.exp(-gamma * np.log(ratio))
        if not np.all(np.isfinite(result)):
            raise EosNumericalError("Characteristic temperature is not finite")
        if result.ndim == 0:
            return float(result)
        return result

    def thermal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye vibrational thermal energy in J mol^-1."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes, temperatures = self._broadcast_state(V, T)
            return self._native_evaluate("thermal_energy", volumes, temperatures)
        V = validate_volume(V)
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise EosValidationError("Temperature must be finite and greater than zero")
        try:
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), temperatures
            )
        except ValueError as error:
            raise EosValidationError(
                "V and T must have broadcast-compatible shapes"
            ) from error

        ratio = self.characteristic_temperature(volumes) / temperatures
        energy = self.Cvmax * temperatures * _debye_function_3(ratio)
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
        entropy = (self.Cvmax / 3.0) * (4.0 * _debye_function_3(ratio) - 3.0 * log_term)
        return self._scalar_or_array(np.asarray(entropy, dtype=float))


class Tange2009Debye(MieGruneisenDebye):
    """Fit3 Mie-Gruneisen-Debye thermal model of Tange et al. (2009).

    This model replaces the constant-``q`` power law with the authors'
    volume dependence

    ``gamma(V) = gamma0 * (1 + a * ((V/V0)**b - 1))``.

    The characteristic temperature is the analytic integral required by
    ``gamma = -d(log(theta))/d(log(V))``. Pressure remains relative to the
    reference isotherm supplied as ``rt_eos``.

    Reference
    ---------
    Tange, Y., Nishihara, Y. & Tsuchiya, T. (2009). Unified analyses for
    P-V-T equation of state of MgO: A solution for pressure-scale problems
    in high P-T experiments. Journal of Geophysical Research, 114, B03208.
    Equations (4), (5), (15), and (16), and Table 4.
    doi:10.1029/2008JB005813
    """

    _constructor_configuration_names: tuple[str, ...] = ()

    @staticmethod
    def constrained_reference_parameters(
        *,
        gamma0: float,
        reference_temperature: float,
        adiabatic_bulk_modulus: float,
        thermal_expansivity: float,
        molar_heat_capacity_p: float,
        n: float,
    ) -> dict[str, float]:
        """Derive the dependent reference parameters used by Tange Fit 3.

        Tange et al. fixed ``K_S0``, ``alpha0``, and ``C_P0`` at 300 K and
        optimized ``gamma0``.  Their thermodynamic constraints therefore make
        ``K_T0`` and ``theta0`` dependent parameters, rather than two additional
        least-squares variables::

            K_T0 = K_S0 / (1 + alpha0 * gamma0 * T0)
            C_V0 = C_P0 / (1 + alpha0 * gamma0 * T0)

        ``theta0`` is the Debye temperature whose heat capacity at ``T0`` is
        ``C_V0``.  Exposing this convention here prevents reproductions from
        accidentally fitting the six printed Table 4 coefficients independently.
        """
        gamma0 = validate_positive_scalar(gamma0, "gamma0")
        temperature = validate_positive_scalar(
            reference_temperature, "reference_temperature"
        )
        ks0 = validate_positive_scalar(adiabatic_bulk_modulus, "adiabatic_bulk_modulus")
        alpha0 = validate_positive_scalar(thermal_expansivity, "thermal_expansivity")
        cp0 = validate_positive_scalar(molar_heat_capacity_p, "molar_heat_capacity_p")
        atoms = validate_positive_scalar(n, "n")
        ratio = 1.0 + alpha0 * gamma0 * temperature
        cv0 = cp0 / ratio
        classical_limit = 3.0 * atoms * R
        if cv0 >= classical_limit:
            raise EosValidationError(
                "The constrained C_V0 must be below the Debye classical limit 3 n R"
            )

        def debye_heat_capacity(theta: float) -> float:
            x = theta / temperature
            occupation_term = 0.0 if x > 700.0 else x / np.expm1(x)
            return classical_limit * (
                4.0 * float(_debye_function_3(x)) - 3.0 * occupation_term
            )

        theta0 = brentq(
            lambda theta: debye_heat_capacity(theta) - cv0,
            1.0e-9 * temperature,
            1.0e4 * temperature,
        )
        return {
            "K0": ks0 / ratio,
            "theta0": float(theta0),
            "Cv0": cv0,
        }

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        a: float,
        b: float,
        n: float,
    ) -> None:
        super().__init__(rt_eos, Tr, theta0, gamma0, q=0.0, n=n)
        self.a = validate_finite_scalar(a, "a")
        self.b = validate_finite_scalar(b, "b")
        if not 0.0 <= self.a <= 1.0:
            raise EosValidationError("a must lie between zero and one")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is Tange2009Debye:
            from peritheos import _rust

            self._native = _rust.ThermalEos.asymptotic_power_law_mie_gruneisen_debye(
                reference_native,
                self.Tr,
                self.theta0,
                self.gamma0,
                self.a,
                self.b,
                self.n,
            )

    def gruneisen_parameter(
        self, V: NumericType, T: NumericType | None = None
    ) -> NumericType:
        """Return the Tange et al. equation (15) Gruneisen parameter."""
        native = getattr(self, "_native", None)
        if native is not None:
            temperatures = self.Tr if T is None else T
            volumes, temperatures = self._broadcast_state(V, temperatures)
            return self._native_evaluate("gruneisen_parameter", volumes, temperatures)
        volumes = np.asarray(validate_volume(V), dtype=float)
        ratio = volumes / self.rt_eos.V0
        result = self.gamma0 * (1.0 + self.a * (ratio**self.b - 1.0))
        if not np.all(np.isfinite(result)):
            raise EosNumericalError("Gruneisen parameter is not finite")
        if result.ndim == 0:
            return float(result)
        return result

    def characteristic_temperature(self, V: NumericType) -> NumericType:
        """Return Debye temperature from Tange et al. equation (16)."""
        native = getattr(self, "_native", None)
        if native is not None:
            volumes = validate_volume(V)
            return self._native_evaluate("characteristic_temperature", volumes, self.Tr)
        volumes = np.asarray(validate_volume(V), dtype=float)
        logarithmic_ratio = np.log(volumes / self.rt_eos.V0)
        if self.b == 0.0:
            exponent = -self.gamma0 * logarithmic_ratio
        else:
            power_minus_one = np.expm1(self.b * logarithmic_ratio)
            exponent = -self.gamma0 * (
                (1.0 - self.a) * logarithmic_ratio + self.a * power_minus_one / self.b
            )
        result = self.theta0 * np.exp(exponent)
        if not np.all(np.isfinite(result)):
            raise EosNumericalError("Characteristic temperature is not finite")
        if result.ndim == 0:
            return float(result)
        return result


class AsymptoticPowerLawMieGruneisenDebyeExcess(Tange2009Debye):
    """Asymptotic-power-law Debye EOS with a volume-dependent ``T^2`` term.

    The quasi-harmonic contribution is the same as :class:`Tange2009Debye`.
    The additional Helmholtz-energy and pressure terms are

    ``F_ex = -0.5 * beta0 * (V / V0)**m * T**2``

    and

    ``P_ex = 0.5 * beta0 * m / V0 * (V / V0)**(m - 1) * T**2``.

    ``beta0`` is in J mol^-1 K^-2 and model volumes are in J bar^-1 mol^-1,
    so pressures are converted from bar to GPa.  Total pressure is referenced
    to ``Tr`` by subtracting both the Debye and excess terms at that
    temperature.  This is the model used by Zhu et al. (2025), equations
    (2)-(8), and by their version-3 Au, Pt, and MgO pressure calculators.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        a: float,
        b: float,
        n: float,
        beta0: float,
        m: float,
    ) -> None:
        super().__init__(rt_eos, Tr, theta0, gamma0, a, b, n)
        self.beta0 = validate_finite_scalar(beta0, "beta0")
        self.m = validate_finite_scalar(m, "m")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None:
            from peritheos import _rust

            # Reuse the exact native quasi-harmonic kernel; this subclass adds
            # the inexpensive excess free-energy terms in Python.
            self._native = _rust.ThermalEos.asymptotic_power_law_mie_gruneisen_debye(
                reference_native,
                self.Tr,
                self.theta0,
                self.gamma0,
                self.a,
                self.b,
                self.n,
            )

    def _state(
        self, V: NumericType, T: NumericType
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        volumes, temperatures = self._broadcast_state(V, T)
        ratio = volumes / self.rt_eos.V0
        return volumes, temperatures, ratio

    def excess_helmholtz_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return the unreferenced ``T^2`` Helmholtz term in J mol^-1."""
        _, temperatures, ratio = self._state(V, T)
        result = -0.5 * self.beta0 * ratio**self.m * temperatures**2
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def excess_internal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the unreferenced ``T^2`` internal-energy term in J mol^-1."""
        _, temperatures, ratio = self._state(V, T)
        result = 0.5 * self.beta0 * ratio**self.m * temperatures**2
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def excess_pressure(
        self, V: NumericType, T: NumericType, *, referenced: bool = False
    ) -> NumericType:
        """Return the unreferenced or ``Tr``-referenced excess pressure in GPa."""
        volumes, temperatures, ratio = self._state(V, T)
        temperature_term = temperatures**2
        if referenced:
            temperature_term = temperature_term - self.Tr**2
        result = (
            0.5
            * self.beta0
            * self.m
            / self.rt_eos.V0
            * ratio ** (self.m - 1.0)
            * temperature_term
            / 1.0e4
        )
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye plus excess pressure relative to ``Tr`` in GPa."""
        volumes, temperatures, _ = self._state(V, T)
        debye_energy = np.asarray(
            MieGruneisenDebye.thermal_energy(self, volumes, temperatures), dtype=float
        )
        reference_energy = np.asarray(
            MieGruneisenDebye.thermal_energy(self, volumes, self.Tr), dtype=float
        )
        debye_pressure = (
            np.asarray(self.gruneisen_parameter(volumes), dtype=float)
            * (debye_energy - reference_energy)
            / volumes
            / 1.0e4
        )
        result = debye_pressure + np.asarray(
            self.excess_pressure(volumes, temperatures, referenced=True), dtype=float
        )
        return self._scalar_or_array(result)

    def thermal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye plus excess internal energy in J mol^-1."""
        debye = np.asarray(
            MieGruneisenDebye.thermal_energy(self, V, T), dtype=float
        )
        result = debye + np.asarray(self.excess_internal_energy(V, T), dtype=float)
        return self._scalar_or_array(result)

    def thermal_entropy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye plus excess entropy in J mol^-1 K^-1."""
        _, temperatures, ratio = self._state(V, T)
        debye = np.asarray(
            MieGruneisenDebye.thermal_entropy(self, V, T), dtype=float
        )
        result = debye + self.beta0 * ratio**self.m * temperatures
        return self._scalar_or_array(result)

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye plus excess ``C_V`` in J mol^-1 K^-1."""
        _, temperatures, ratio = self._state(V, T)
        steps = 1.0e-5 * temperatures
        debye = (
            np.asarray(
                MieGruneisenDebye.thermal_energy(
                    self, V, temperatures + steps
                ),
                dtype=float,
            )
            - np.asarray(
                MieGruneisenDebye.thermal_energy(
                    self, V, temperatures - steps
                ),
                dtype=float,
            )
        ) / (2.0 * steps)
        result = debye + self.beta0 * ratio**self.m * temperatures
        return self._scalar_or_array(result)

    def thermal_helmholtz_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return Debye plus excess Helmholtz energy in J mol^-1."""
        volumes, temperatures, _ = self._state(V, T)
        debye_energy = np.asarray(
            MieGruneisenDebye.thermal_energy(self, volumes, temperatures), dtype=float
        )
        debye_entropy = np.asarray(
            MieGruneisenDebye.thermal_entropy(self, volumes, temperatures), dtype=float
        )
        result = (
            debye_energy
            - temperatures * debye_entropy
            + np.asarray(
                self.excess_helmholtz_free_energy(volumes, temperatures), dtype=float
            )
        )
        return self._scalar_or_array(result)

    def vibrational_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the unreferenced quasi-harmonic Debye pressure in GPa."""
        volumes, temperatures, _ = self._state(V, T)
        debye_energy = np.asarray(
            MieGruneisenDebye.thermal_energy(self, volumes, temperatures), dtype=float
        )
        result = (
            np.asarray(self.gruneisen_parameter(volumes), dtype=float)
            * debye_energy
            / volumes
            / 1.0e4
        )
        return self._scalar_or_array(result)

    def thermal_enthalpy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the full thermal enthalpy contribution in J mol^-1."""
        volumes, temperatures, _ = self._state(V, T)
        pressure = np.asarray(
            self.vibrational_pressure(volumes, temperatures), dtype=float
        ) + np.asarray(self.excess_pressure(volumes, temperatures), dtype=float)
        result = (
            np.asarray(self.thermal_energy(volumes, temperatures), dtype=float)
            + pressure * volumes * 1.0e4
        )
        return self._scalar_or_array(result)

    def thermal_gibbs_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return the full thermal Gibbs-energy contribution in J mol^-1."""
        volumes, temperatures, _ = self._state(V, T)
        pressure = np.asarray(
            self.vibrational_pressure(volumes, temperatures), dtype=float
        ) + np.asarray(self.excess_pressure(volumes, temperatures), dtype=float)
        result = (
            np.asarray(
                self.thermal_helmholtz_free_energy(volumes, temperatures), dtype=float
            )
            + pressure * volumes * 1.0e4
        )
        return self._scalar_or_array(result)


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
            raise EosValidationError("Temperature must be finite and greater than zero")
        try:
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), temperatures
            )
        except ValueError as error:
            raise EosValidationError(
                "V and T must have broadcast-compatible shapes"
            ) from error

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
        raise EosValidationError("Debye-function arguments must be finite and positive")

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
