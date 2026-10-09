"""Bounded absolute Debye and electronic pressure with a separate cold volume."""

from __future__ import annotations

import numpy as np
from scipy.constants import R
from scipy.optimize import brentq

from peritheos.eos import (
    EosBase,
    NumericType,
    ThermalEOS,
    validate_positive_scalar,
    validate_volume,
)
from peritheos.errors import EosValidationError

from .mie_gruneisen import _debye_function_3


class AsymptoticDebyeTabulatedPressure(ThermalEOS):
    """Pc(V,0)+Pph(V,T)+Pel(T), in GPa and J/bar/mol volume units.

    ``rt_eos.V0`` is the cold 0 K volume Vc. The phonon reference V0 is
    ``Vc/cold_volume_ratio`` at ambient 300 K. Electronic pressure is absolute
    relative to 0 K; it is not subtracted at Tr. This pressure-only model has
    no electronic caloric potential. Numerical bounds do not certify phase
    stability. The source phase annotations belong to record provenance.
    """

    allows_zero_temperature = True
    _constructor_configuration_names = (
        "electronic_temperature_k",
        "electronic_pressure_gpa",
        "interpolation",
        "volume_ratio_range",
        "temperature_range_k",
        "residual_temperature_k",
        "residual_pressure_gpa",
    )

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        a: float,
        b: float,
        n: float,
        cold_volume_ratio: float,
        *,
        electronic_temperature_k,
        electronic_pressure_gpa,
        interpolation: str,
        volume_ratio_range,
        temperature_range_k,
        residual_temperature_k,
        residual_pressure_gpa,
    ) -> None:
        super().__init__(rt_eos)
        for name, value in dict(
            Tr=Tr,
            theta0=theta0,
            gamma0=gamma0,
            b=b,
            n=n,
            cold_volume_ratio=cold_volume_ratio,
        ).items():
            setattr(self, name, validate_positive_scalar(value, name))
        self.a = float(a)
        if not np.isfinite(self.a) or not 0 <= self.a <= 1:
            raise EosValidationError("a must lie between zero and one")
        t = np.asarray(electronic_temperature_k, dtype=float)
        p = np.asarray(electronic_pressure_gpa, dtype=float)
        if (
            t.ndim != 1
            or t.shape != p.shape
            or len(t) < 2
            or not np.all(np.isfinite(t))
            or not np.all(np.isfinite(p))
            or t[0] < 0
            or np.any(np.diff(t) <= 0)
            or np.any(np.diff(p) < 0)
        ):
            raise EosValidationError(
                "Electronic table requires increasing temperatures and nondecreasing finite pressures"
            )
        if interpolation != "linear":
            raise EosValidationError(
                "Electronic interpolation must be explicitly 'linear'"
            )
        self.electronic_temperature_k = tuple(t.tolist())
        self.electronic_pressure_gpa = tuple(p.tolist())
        self.interpolation = interpolation
        residual_t = np.asarray(residual_temperature_k, dtype=float)
        residual_p = np.asarray(residual_pressure_gpa, dtype=float)
        if (
            residual_t.ndim != 1
            or residual_t.shape != residual_p.shape
            or len(residual_t) < 2
            or not np.all(np.isfinite(residual_t))
            or not np.all(np.isfinite(residual_p))
            or residual_t[0] < 0
            or np.any(np.diff(residual_t) <= 0)
        ):
            raise EosValidationError(
                "Residual pressure table requires matching finite arrays and increasing temperatures"
            )
        self.residual_temperature_k = tuple(residual_t.tolist())
        self.residual_pressure_gpa = tuple(residual_p.tolist())
        self.volume_ratio_range = self._range(
            volume_ratio_range, "Volume ratio", positive=True
        )
        self.temperature_range_k = self._range(
            temperature_range_k, "Temperature", positive=False
        )
        lo, hi = self.temperature_range_k
        if lo < t[0] or hi > t[-1] or not lo <= self.Tr <= hi:
            raise EosValidationError(
                "Temperature range and Tr must lie within the electronic table"
            )
        if lo < residual_t[0] or hi > residual_t[-1]:
            raise EosValidationError(
                "Temperature range must lie within the residual pressure table"
            )

    @staticmethod
    def _range(values, name, *, positive):
        array = np.asarray(values, dtype=float)
        if (
            array.shape != (2,)
            or not np.all(np.isfinite(array))
            or array[0] >= array[1]
            or array[0] < 0
            or (positive and array[0] == 0)
        ):
            raise EosValidationError(
                f"{name} range must contain two ordered finite bounds"
            )
        return tuple(array.tolist())

    @property
    def ambient_reference_volume(self):
        """Ambient 300 K molar-volume normalization, distinct from cold Vc."""
        return self.rt_eos.V0 / self.cold_volume_ratio

    def _volumes(self, V):
        volumes = np.asarray(validate_volume(V), dtype=float)
        ratio = volumes / self.ambient_reference_volume
        lo, hi = self.volume_ratio_range
        tolerance = 16 * np.finfo(float).eps
        if np.any((ratio < lo - tolerance) | (ratio > hi + tolerance)):
            raise EosValidationError("Volume lies outside the reconstruction domain")
        return volumes

    def _temperatures(self, T):
        temperatures = np.asarray(T, dtype=float)
        lo, hi = self.temperature_range_k
        if np.any(~np.isfinite(temperatures)) or np.any(
            (temperatures < lo) | (temperatures > hi)
        ):
            raise EosValidationError(
                "Temperature lies outside the reconstruction domain"
            )
        return temperatures

    def _state(self, V, T):
        volumes, temperatures = self._volumes(V), self._temperatures(T)
        try:
            return np.broadcast_arrays(volumes, temperatures)
        except ValueError as error:
            raise EosValidationError(
                "V and T must have broadcast-compatible shapes"
            ) from error

    def gruneisen_parameter(self, V, T=None):
        ratio = self._volumes(V) / self.ambient_reference_volume
        return self._scalar_or_array(
            np.asarray(self.gamma0 * (1 + self.a * (ratio**self.b - 1)))
        )

    def characteristic_temperature(self, V):
        ratio = self._volumes(V) / self.ambient_reference_volume
        gamma = self.gruneisen_parameter(V)
        theta = (
            self.theta0
            * ratio ** (-(1 - self.a) * self.gamma0)
            * np.exp(-(gamma - self.gamma0) / self.b)
        )
        return self._scalar_or_array(np.asarray(theta))

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Absolute phonon plus electronic pressure; zero phonon energy at 0 K."""
        volumes, temperatures = self._state(V, T)
        theta = self.characteristic_temperature(volumes)
        energy = np.zeros(volumes.shape)
        positive = temperatures > 0
        energy[positive] = (
            3
            * self.n
            * R
            * temperatures[positive]
            * _debye_function_3(np.asarray(theta)[positive] / temperatures[positive])
        )
        phonon = self.gruneisen_parameter(volumes) * energy / volumes / 1e4
        electronic = np.interp(
            temperatures, self.electronic_temperature_k, self.electronic_pressure_gpa
        )
        residual = np.interp(
            temperatures, self.residual_temperature_k, self.residual_pressure_gpa
        )
        return self._scalar_or_array(np.asarray(phonon + electronic + residual))

    def thermal_pressure_increment(self, V, T):
        return self.thermal_pressure(V, T) - self.thermal_pressure(V, self.Tr)

    @staticmethod
    def _root(function, lower, upper):
        left, right = function(lower), function(upper)
        tolerance = 1e-10
        if abs(left) <= tolerance:
            return lower
        if abs(right) <= tolerance:
            return upper
        if left * right > 0:
            raise EosValidationError("No root within the reconstruction domain")
        return brentq(function, lower, upper, xtol=1e-12, rtol=1e-12)

    def _volume_roots(self, P, T, f_dac=0):
        pressures = np.asarray(P, dtype=float)
        if np.any(~np.isfinite(pressures)):
            raise EosValidationError("Pressure must be finite")
        try:
            pressures, temperatures = np.broadcast_arrays(
                pressures, self._temperatures(T)
            )
        except ValueError as error:
            raise EosValidationError(
                "P and T must have broadcast-compatible shapes"
            ) from error
        lower, upper = np.array(self.volume_ratio_range) * self.ambient_reference_volume
        roots = [
            self._root(
                lambda v: float(
                    self.pressure(v, t)
                    - f_dac * self.thermal_pressure_increment(v, t)
                    - p
                ),
                lower,
                upper,
            )
            for p, t in zip(pressures.flat, temperatures.flat)
        ]
        return self._scalar_or_array(np.array(roots).reshape(pressures.shape))

    def calculate_volume(self, P, T):
        return self._volume_roots(P, T)

    def calculate_temperature(self, P, V):
        pressures = np.asarray(P, dtype=float)
        if np.any(~np.isfinite(pressures)):
            raise EosValidationError("Pressure must be finite")
        try:
            pressures, volumes = np.broadcast_arrays(pressures, self._volumes(V))
        except ValueError as error:
            raise EosValidationError(
                "P and V must have broadcast-compatible shapes"
            ) from error
        lower, upper = self.temperature_range_k
        roots = [
            self._root(lambda t: float(self.pressure(v, t) - p), lower, upper)
            for p, v in zip(pressures.flat, volumes.flat)
        ]
        return self._scalar_or_array(np.array(roots).reshape(pressures.shape))

    def volume_with_dac_confinement(self, P_cold, T, *, f_dac):
        return self._volume_roots(P_cold, T, self._validate_f_dac(f_dac))

    def temperature_from_volumes(self, V_ambient, V_heated, *, f_dac):
        f_dac = self._validate_f_dac(f_dac)
        ambient, heated = np.broadcast_arrays(
            self._volumes(V_ambient), self._volumes(V_heated)
        )
        reference_hot = self.pressure(heated, self.Tr)
        increment = (self.pressure(ambient, self.Tr) - reference_hot) / (1 - f_dac)
        if np.any(increment < -1e-10):
            raise EosValidationError("The volume pair implies a temperature below Tr")
        return self.calculate_temperature(reference_hot + increment, heated)
