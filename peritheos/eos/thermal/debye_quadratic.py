"""Debye thermal pressure with an empirical quadratic-temperature term."""

from __future__ import annotations

import numpy as np

from peritheos.eos import (
    EosBase,
    NumericType,
    ThermalEOS,
    _native_for_exact_model,
    _native_thermal_evaluate,
    validate_finite_scalar,
)

from .mie_gruneisen import MieGruneisenDebye


class DebyeQuadraticThermalPressure(ThermalEOS):
    r"""Add MGD pressure and ``A * (V/V0)**m * (T**2-Tr**2)``.

    ``A`` is in GPa/K², ``m`` is dimensionless, and the Debye parameters
    ``Tr``, ``theta0``, ``gamma0``, ``q``, ``n`` follow
    :class:`MieGruneisenDebye` with the integrated Gruneisen law. Volumes
    are molar, in J/bar/mol. The model implements Fei et al. (2016),
    doi:10.1002/2016GL069456, equation (2), with
    ``A = gamma_e * beta0 * rho0 / 2e9`` for SI mass-specific beta0 and
    density, and ``m = k - 1``.

    This is an empirical pressure surface. It does not assert a caloric
    potential: the source's independent electronic Gruneisen coefficient
    need not equal the volume exponent of its electronic heat capacity.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        q: float,
        n: float,
        A: float,
        m: float,
    ) -> None:
        super().__init__(rt_eos)
        self._debye = MieGruneisenDebye(rt_eos, Tr, theta0, gamma0, q, n)
        self.Tr = self._debye.Tr
        self.theta0 = self._debye.theta0
        self.gamma0 = self._debye.gamma0
        self.q = self._debye.q
        self.n = self._debye.n
        self.A = validate_finite_scalar(A, "A")
        self.m = validate_finite_scalar(m, "m")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is DebyeQuadraticThermalPressure:
            from peritheos import _rust

            self._native = _rust.ThermalEos.debye_quadratic_thermal_pressure(
                reference_native,
                self.Tr,
                self.theta0,
                self.gamma0,
                self.q,
                self.n,
                self.A,
                self.m,
            )

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        result = self._debye.thermal_pressure(volumes, temperatures) + (
            self.A
            * (volumes / self.rt_eos.V0) ** self.m
            * (temperatures**2 - self.Tr**2)
        )
        return self._scalar_or_array(np.asarray(result, dtype=float))
