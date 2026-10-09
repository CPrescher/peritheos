"""Debye pressure plus a bounded, volume-independent tabulated correction."""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from peritheos.eos import EosBase, NumericType, ThermalEOS, validate_volume
from peritheos.errors import EosValidationError

from .mie_gruneisen import MieGruneisenDebye


class DebyeTabulatedThermalPressure(ThermalEOS):
    """Add MGD pressure and ``P_el(T) - P_el(Tr)`` in GPa.

    The supplied raw electronic pressures are fixed configuration, never
    scalar fitting parameters. Linear interpolation is an explicit numerical
    choice; no extrapolation outside the supplied temperature table is allowed.
    The integrated Gruneisen law and molar J/bar/mol volume convention follow
    :class:`MieGruneisenDebye`. This pressure surface supplies no electronic
    caloric potential or electronic heat capacity.
    """

    _constructor_configuration_names = (
        "electronic_temperature_k",
        "electronic_pressure_gpa",
        "interpolation",
    )

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        q: float,
        n: float,
        *,
        electronic_temperature_k,
        electronic_pressure_gpa,
        interpolation: str,
    ) -> None:
        super().__init__(rt_eos)
        self._debye = MieGruneisenDebye(rt_eos, Tr, theta0, gamma0, q, n)
        for name in ("Tr", "theta0", "gamma0", "q", "n"):
            setattr(self, name, getattr(self._debye, name))
        t = np.asarray(electronic_temperature_k, dtype=float)
        p = np.asarray(electronic_pressure_gpa, dtype=float)
        if (
            t.ndim != 1
            or p.shape != t.shape
            or len(t) < 2
            or not np.all(np.isfinite(t))
            or not np.all(np.isfinite(p))
            or t[0] < 0
            or np.any(np.diff(t) <= 0)
            or np.any(np.diff(p) < 0)
        ):
            raise EosValidationError(
                "Electronic table requires matching finite one-dimensional arrays, "
                "increasing nonnegative temperatures and nondecreasing pressures"
            )
        if not t[0] <= self.Tr <= t[-1]:
            raise EosValidationError(
                "Tr must lie within the electronic temperature table"
            )
        if interpolation != "linear":
            raise EosValidationError(
                "Electronic interpolation must be explicitly 'linear'"
            )
        self.electronic_temperature_k = tuple(t.tolist())
        self.electronic_pressure_gpa = tuple(p.tolist())
        self.interpolation = interpolation
        self._reference_electronic = float(np.interp(self.Tr, t, p))

    def electronic_pressure_increment(self, T: NumericType) -> NumericType:
        """Interpolate raw tabulated pressure and subtract its Tr value once."""
        t = np.asarray(T, dtype=float)
        lo, hi = self.electronic_temperature_k[0], self.electronic_temperature_k[-1]
        if not np.all(np.isfinite(t)) or np.any(t <= 0) or np.any((t < lo) | (t > hi)):
            raise EosValidationError(
                f"Temperature must be positive and within the electronic table [{lo}, {hi}] K"
            )
        result = np.interp(
            t, self.electronic_temperature_k, self.electronic_pressure_gpa
        )
        return self._scalar_or_array(np.asarray(result - self._reference_electronic))

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        electronic = self.electronic_pressure_increment(temperatures)
        result = self._debye.thermal_pressure(volumes, temperatures) + electronic
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def calculate_temperature(self, P: NumericType, V: NumericType) -> NumericType:
        """Invert only within the electronic table's positive-temperature domain."""
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
        # A 0 K source node defines the limiting correction. ThermalEOS accepts
        # positive T; this small lower bound avoids underflow in Debye kernels.
        lower = max(1e-8, self.electronic_temperature_k[0])
        upper = self.electronic_temperature_k[-1]
        result = []
        for pressure, volume in zip(pressures.flat, volumes.flat):

            def function(t):
                return float(self.pressure(float(volume), t) - pressure)

            a, b = function(lower), function(upper)
            if a == 0:
                result.append(lower)
            elif b == 0:
                result.append(upper)
            elif a * b > 0:
                raise EosValidationError(
                    "No temperature root within the electronic table"
                )
            else:
                result.append(brentq(function, lower, upper, xtol=1e-8))
        return self._scalar_or_array(np.asarray(result).reshape(pressures.shape))

    def temperature_from_volumes(self, V_ambient, V_heated, *, f_dac):
        """Use the bounded inversion for the standard DAC confinement equation."""
        f_dac = self._validate_f_dac(f_dac)
        ambient = np.asarray(validate_volume(V_ambient), dtype=float)
        heated = np.asarray(validate_volume(V_heated), dtype=float)
        try:
            ambient, heated = np.broadcast_arrays(ambient, heated)
        except ValueError as error:
            raise EosValidationError(
                "Volumes must have broadcast-compatible shapes"
            ) from error
        cold = np.asarray(self.rt_eos.pressure(heated), dtype=float)
        target = (self.rt_eos.pressure(ambient) - cold) / (1 - f_dac)
        if np.any(target < 0):
            raise EosValidationError("The volume pair implies a temperature below Tr")
        return self.calculate_temperature(cold + target, heated)
