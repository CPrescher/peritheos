"""Monatomic Debye plus anharmonic Helmholtz energy of Xiao et al. (2025)."""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from peritheos.eos import (
    EosBase,
    NumericType,
    ThermalEOS,
    _native_for_exact_model,
    _native_thermal_evaluate,
    validate_finite_scalar,
)
from peritheos.errors import EosValidationError

from .mie_gruneisen import MieGruneisenDebye


class DebyeAnharmonicHelmholtz(ThermalEOS):
    r"""Published Xiao et al. solid-argon thermal Helmholtz model.

    DOI: 10.1007/s10765-024-03469-2, equations 19–25. The monatomic Debye
    energy excludes zero-point energy and uses the integrated power-law gamma.
    The anharmonic energy is ``b1 R theta0 t**4/(1+b2*t**2)
    * exp(b3*(V/V0-1))``, with ``t=T/theta0``. The reference EOS must be
    the zero-temperature cold curve. ``Tr`` is only a positive-temperature
    reporting/increment reference; full pressure never subtracts its energy.

    Pressure is GPa; volume is J/bar/mol (cm³/mol divided by 10); energies
    are J/mol and entropy/heat capacity J/mol/K. The gas constant is 8.31451
    J/mol/K, matching the authors’ supplementary workbook. Full energies use zero cold
    energy at V0 and do not include a gas-phase reference-state offset.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        q: float,
        b1: float,
        b2: float,
        b3: float,
    ) -> None:
        super().__init__(rt_eos)
        self._debye = MieGruneisenDebye(
            rt_eos,
            Tr,
            theta0,
            gamma0,
            q,
            1.0,
            thermal_pressure_reference="absolute_zero",
            Cvmax=3.0 * 8.31451,
        )
        self.Tr, self.theta0 = self._debye.Tr, self._debye.theta0
        self.gamma0, self.q = self._debye.gamma0, self._debye.q
        self.b1 = validate_finite_scalar(b1, "b1")
        self.b2 = validate_finite_scalar(b2, "b2")
        self.b3 = validate_finite_scalar(b3, "b3")
        if self.b2 < 0:
            raise EosValidationError("b2 must be nonnegative")
        reference = _native_for_exact_model(rt_eos)
        if reference is not None and type(self) is DebyeAnharmonicHelmholtz:
            from peritheos import _rust

            self._native = _rust.ThermalEos.debye_anharmonic_helmholtz(
                reference,
                self.Tr,
                self.theta0,
                self.gamma0,
                self.q,
                self.b1,
                self.b2,
                self.b3,
            )

    def _anharmonic(self, V: NumericType, T: NumericType) -> NumericType:
        R = 8.31451  # Authors’ supplementary VBA gas-constant convention.

        t = np.asarray(T) / self.theta0
        return (
            self.b1
            * R
            * self.theta0
            * t**4
            / (1 + self.b2 * t**2)
            * np.exp(self.b3 * (np.asarray(V) / self.rt_eos.V0 - 1))
        )

    def _evaluate(self, quantity: str, V: NumericType, T: NumericType) -> NumericType:
        v, t = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(self._native, quantity, v, t)
        a = self._anharmonic(v, t)
        y = self.b2 * (t / self.theta0) ** 2
        if quantity == "thermal_pressure":
            result = (
                self._debye.thermal_pressure(v, t) - self.b3 / self.rt_eos.V0 / 1e4 * a
            )
        elif quantity == "entropy":
            result = self._debye.thermal_entropy(v, t) - a / t * (4 + 2 * y) / (1 + y)
        elif quantity == "molar_heat_capacity_v":
            result = (
                self._debye.molar_heat_capacity_v(v, t)
                - a / t * (12 + 6 * y + 2 * y * y) / (1 + y) ** 2
            )
        elif quantity == "thermal_helmholtz_free_energy":
            result = self._debye.thermal_helmholtz_free_energy(v, t) + a
        elif quantity == "helmholtz_free_energy":
            cold = np.vectorize(
                lambda volume: (
                    -1e4 * quad(self.rt_eos.pressure, self.rt_eos.V0, volume)[0]
                )
            )(v)
            result = cold + self._evaluate("thermal_helmholtz_free_energy", v, t)
        elif quantity == "internal_energy":
            result = self.helmholtz_free_energy(v, t) + t * self.entropy(v, t)
        else:
            raise ValueError(f"Unsupported quantity: {quantity}")
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Debye plus anharmonic pressure relative to absolute zero, in GPa."""
        return self._evaluate("thermal_pressure", V, T)

    def thermal_pressure_increment(self, V: NumericType, T: NumericType) -> NumericType:
        """Thermal pressure difference from the reporting temperature Tr."""
        return self.thermal_pressure(V, T) - self.thermal_pressure(V, self.Tr)

    def entropy(self, V: NumericType, T: NumericType) -> NumericType:
        """Molar entropy in J/mol/K."""
        return self._evaluate("entropy", V, T)

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Constant-volume molar heat capacity in J/mol/K."""
        return self._evaluate("molar_heat_capacity_v", V, T)

    def thermal_helmholtz_free_energy(
        self, V: NumericType, T: NumericType
    ) -> NumericType:
        """Debye and anharmonic energy without the cold term, in J/mol."""
        return self._evaluate("thermal_helmholtz_free_energy", V, T)

    def helmholtz_free_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Full Helmholtz energy with zero cold energy at V0, in J/mol."""
        return self._evaluate("helmholtz_free_energy", V, T)

    def internal_energy(self, V: NumericType, T: NumericType) -> NumericType:
        """Full internal energy with zero cold energy at V0, in J/mol."""
        return self._evaluate("internal_energy", V, T)
