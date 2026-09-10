"""Sound-velocity-constrained quasi-Debye Helmholtz equations of state."""

from __future__ import annotations

import numpy as np
from scipy.constants import Avogadro, Boltzmann, R, hbar
from scipy.integrate import quad
from scipy.optimize import brentq

from peritheos.eos import (
    EosBase,
    NumericType,
    ThermalEOS,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)
from peritheos.eos.thermal.mie_gruneisen import _debye_function_3
from peritheos.errors import ConfigurationError, EosNumericalError, EosValidationError


class SoundVelocityDebyeHelmholtz(ThermalEOS):
    """Quasi-Debye EOS whose Debye temperature follows elastic velocities.

    The model implements Appendix A of Luo et al. (2023). Linear, independent
    room-temperature fits of longitudinal and shear velocity against density
    provide only their ratio and hence Poisson's ratio. The cold-curve bulk
    modulus and that Poisson ratio then determine the model velocities,
    effective Debye velocity, and characteristic temperature. The vibrational
    energy excludes zero-point energy, matching Luo et al.'s integral from
    zero temperature.

    Volumes use the :class:`~peritheos.eos.ThermalEOS` molar convention of
    J bar^-1 mol^-1. Velocity intercepts are in km/s and slopes in
    km s^-1 per g cm^-3. ``n`` is atoms per formula unit and
    ``molar_mass_g_mol`` is the formula-unit molar mass.

    References
    ----------
    Luo, Y. et al. (2023), Physical Review B 107, 134116, equations (1)--(6)
    and Appendix A. doi:10.1103/PhysRevB.107.134116

    Xian, Y. T. et al. (2022), AIP Advances 12, 055313, equations (1)--(15).
    doi:10.1063/5.0089292
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        molar_mass_g_mol: float,
        n: float,
        longitudinal_intercept: float,
        longitudinal_slope: float,
        shear_intercept: float,
        shear_slope: float,
    ) -> None:
        if not isinstance(rt_eos, EosBase):
            raise ConfigurationError("rt_eos must be an equation of state")
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.molar_mass_g_mol = validate_positive_scalar(
            molar_mass_g_mol, "molar_mass_g_mol"
        )
        self.n = validate_positive_scalar(n, "n")
        self.longitudinal_intercept = validate_finite_scalar(
            longitudinal_intercept, "longitudinal_intercept"
        )
        self.longitudinal_slope = validate_finite_scalar(
            longitudinal_slope, "longitudinal_slope"
        )
        self.shear_intercept = validate_finite_scalar(
            shear_intercept, "shear_intercept"
        )
        self.shear_slope = validate_finite_scalar(shear_slope, "shear_slope")

    @staticmethod
    def _result(values: np.ndarray) -> NumericType:
        if values.ndim == 0:
            return float(values)
        return values

    def density(self, V: NumericType) -> NumericType:
        """Return density in g cm^-3 for molar volume *V*."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        # 1 J bar^-1 = 10 cm^3.
        density = self.molar_mass_g_mol / (10.0 * volumes)
        return self._result(density)

    def velocity_fit_values(self, V: NumericType) -> tuple[NumericType, NumericType]:
        """Return the two empirical room-temperature velocity regressions."""
        density = np.asarray(self.density(V), dtype=float)
        longitudinal = self.longitudinal_intercept + self.longitudinal_slope * density
        shear = self.shear_intercept + self.shear_slope * density
        if (
            not np.all(np.isfinite(longitudinal))
            or not np.all(np.isfinite(shear))
            or np.any(longitudinal <= shear)
            or np.any(shear <= 0.0)
        ):
            raise EosValidationError(
                "Velocity regressions must satisfy longitudinal > shear > 0"
            )
        return self._result(longitudinal), self._result(shear)

    def poisson_ratio(self, V: NumericType) -> NumericType:
        """Return Poisson's ratio inferred from the empirical velocity ratio."""
        longitudinal, shear = self.velocity_fit_values(V)
        tau_squared = (
            np.asarray(longitudinal, dtype=float) / np.asarray(shear, dtype=float)
        ) ** 2
        poisson = (tau_squared - 2.0) / (2.0 * (tau_squared - 1.0))
        if not np.all(np.isfinite(poisson)) or np.any(
            (poisson <= -1.0) | (poisson >= 0.5)
        ):
            raise EosValidationError("Velocity ratio implies an invalid Poisson ratio")
        return self._result(poisson)

    def shear_modulus(self, V: NumericType) -> NumericType:
        """Return the Appendix-A shear modulus in GPa."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        poisson = np.asarray(self.poisson_ratio(volumes), dtype=float)
        bulk = np.asarray(self.rt_eos.bulk_modulus(volumes), dtype=float)
        shear = 3.0 * (1.0 - 2.0 * poisson) * bulk / (2.0 * (1.0 + poisson))
        if not np.all(np.isfinite(shear)) or np.any(shear <= 0.0):
            raise EosNumericalError("Calculated shear modulus is not positive and finite")
        return self._result(shear)

    def longitudinal_velocity(self, V: NumericType) -> NumericType:
        """Return model longitudinal velocity in km/s."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        bulk = np.asarray(self.rt_eos.bulk_modulus(volumes), dtype=float)
        shear = np.asarray(self.shear_modulus(volumes), dtype=float)
        density_kg_m3 = np.asarray(self.density(volumes), dtype=float) * 1000.0
        velocity = np.sqrt((bulk + 4.0 * shear / 3.0) * 1.0e9 / density_kg_m3)
        return self._result(velocity / 1000.0)

    def shear_velocity(self, V: NumericType) -> NumericType:
        """Return model shear velocity in km/s."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        shear = np.asarray(self.shear_modulus(volumes), dtype=float)
        density_kg_m3 = np.asarray(self.density(volumes), dtype=float) * 1000.0
        velocity = np.sqrt(shear * 1.0e9 / density_kg_m3)
        return self._result(velocity / 1000.0)

    def effective_sound_velocity(self, V: NumericType) -> NumericType:
        """Return the Debye-average sound velocity in km/s."""
        longitudinal = np.asarray(self.longitudinal_velocity(V), dtype=float)
        shear = np.asarray(self.shear_velocity(V), dtype=float)
        effective = (3.0 / (longitudinal**-3 + 2.0 * shear**-3)) ** (1.0 / 3.0)
        return self._result(effective)

    def characteristic_temperature(self, V: NumericType) -> NumericType:
        """Return the sound-velocity-derived Debye temperature in K."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        volume_m3_mol = volumes * 1.0e-5
        number_density = self.n * Avogadro / volume_m3_mol
        velocity_m_s = np.asarray(self.effective_sound_velocity(volumes)) * 1000.0
        theta = (
            hbar
            / Boltzmann
            * np.cbrt(6.0 * np.pi**2 * number_density)
            * velocity_m_s
        )
        if not np.all(np.isfinite(theta)) or np.any(theta <= 0.0):
            raise EosNumericalError("Characteristic temperature is not positive and finite")
        return self._result(theta)

    def gruneisen_parameter(
        self, V: NumericType, T: NumericType | None = None
    ) -> NumericType:
        """Return ``-d(log(theta))/d(log(V))`` by centered differentiation."""
        volumes = np.asarray(validate_volume(V), dtype=float)
        step = 1.0e-6
        theta_low = np.asarray(
            self.characteristic_temperature(volumes * np.exp(-step)), dtype=float
        )
        theta_high = np.asarray(
            self.characteristic_temperature(volumes * np.exp(step)), dtype=float
        )
        gamma = -(np.log(theta_high) - np.log(theta_low)) / (2.0 * step)
        if not np.all(np.isfinite(gamma)):
            raise EosNumericalError("Gruneisen parameter is not finite")
        return self._result(gamma)

    def thermal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return vibrational internal energy without zero-point energy."""
        volumes, temperatures = self._broadcast_state(V, T)
        ratio = np.asarray(self.characteristic_temperature(volumes)) / temperatures
        energy = 3.0 * self.n * R * temperatures * _debye_function_3(ratio)
        return self._result(np.asarray(energy, dtype=float))

    def thermal_entropy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return vibrational entropy in J mol^-1 K^-1."""
        volumes, temperatures = self._broadcast_state(V, T)
        ratio = np.asarray(self.characteristic_temperature(volumes)) / temperatures
        entropy = self.n * R * (
            4.0 * _debye_function_3(ratio) - 3.0 * np.log(-np.expm1(-ratio))
        )
        return self._result(np.asarray(entropy, dtype=float))

    def thermal_helmholtz_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Return vibrational Helmholtz energy in J mol^-1."""
        volumes, temperatures = self._broadcast_state(V, T)
        free_energy = np.asarray(self.thermal_energy(volumes, temperatures)) - (
            temperatures * np.asarray(self.thermal_entropy(volumes, temperatures))
        )
        return self._result(np.asarray(free_energy, dtype=float))

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Debye constant-volume heat capacity in J mol^-1 K^-1."""
        volumes, temperatures = self._broadcast_state(V, T)
        ratio = np.asarray(self.characteristic_temperature(volumes)) / temperatures
        heat_capacity = 3.0 * self.n * R * (
            4.0 * _debye_function_3(ratio)
            - 3.0 * ratio / np.expm1(ratio)
        )
        return self._result(np.asarray(heat_capacity, dtype=float))

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the absolute vibrational pressure in GPa."""
        volumes, temperatures = self._broadcast_state(V, T)
        pressure = (
            np.asarray(self.gruneisen_parameter(volumes))
            * np.asarray(self.thermal_energy(volumes, temperatures))
            / volumes
            / 1.0e4
        )
        return self._result(np.asarray(pressure, dtype=float))

    def thermal_pressure_increment(self, V: NumericType, T: NumericType) -> NumericType:
        """Return vibrational pressure above the configured reference temperature."""
        increment = np.asarray(self.thermal_pressure(V, T), dtype=float) - np.asarray(
            self.thermal_pressure(V, self.Tr), dtype=float
        )
        return self._result(increment)

    def cold_energy(self, V: NumericType) -> NumericType:
        """Return cold compression energy relative to ``rt_eos.V0`` in J/mol."""
        volumes = np.asarray(validate_volume(V), dtype=float)

        def one(volume: float) -> float:
            integral, _ = quad(
                lambda value: float(self.rt_eos.pressure(value)) * 1.0e4,
                float(self.rt_eos.V0),
                volume,
                epsabs=1.0e-7,
                epsrel=1.0e-10,
            )
            return -integral

        values = np.asarray([one(float(value)) for value in volumes.flat]).reshape(
            volumes.shape
        )
        return self._result(values)

    def internal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Return cold plus vibrational internal energy in J mol^-1."""
        volumes, temperatures = self._broadcast_state(V, T)
        energy = np.asarray(self.cold_energy(volumes)) + np.asarray(
            self.thermal_energy(volumes, temperatures)
        )
        return self._result(np.asarray(energy, dtype=float))

    def hugoniot_temperature_from_state(
        self,
        V: float,
        P: float,
        *,
        initial_volume: float,
        initial_temperature: float = 300.0,
        initial_pressure: float = 0.0,
    ) -> float:
        """Solve shock temperature from an observed pressure-volume state."""
        volume = validate_positive_scalar(V, "V")
        pressure = validate_finite_scalar(P, "P")
        initial_volume = validate_positive_scalar(initial_volume, "initial_volume")
        initial_temperature = validate_positive_scalar(
            initial_temperature, "initial_temperature"
        )
        initial_pressure = validate_finite_scalar(initial_pressure, "initial_pressure")
        initial_energy = float(self.internal_energy(initial_volume, initial_temperature))
        target_energy = (
            initial_energy
            + 0.5
            * (pressure + initial_pressure)
            * (initial_volume - volume)
            * 1.0e4
        )

        def residual(temperature: float) -> float:
            return float(self.internal_energy(volume, temperature)) - target_energy

        lower = 1.0e-9
        upper = max(initial_temperature, 1000.0)
        while residual(upper) < 0.0 and upper < 1.0e8:
            upper *= 2.0
        if residual(lower) * residual(upper) > 0.0:
            raise EosValidationError("Could not bracket the Hugoniot temperature")
        return float(brentq(residual, lower, upper, xtol=1.0e-8, rtol=1.0e-12))
