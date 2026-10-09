"""Thermodynamically integrated variable-q Mie-Gruneisen-Debye EOS."""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_for_exact_model,
    validate_finite_scalar,
    validate_volume,
)
from peritheos.errors import EosNumericalError

from .mie_gruneisen import MieGruneisenDebye


class VariableQMieGruneisenDebye(MieGruneisenDebye):
    """Debye thermal pressure with ``q(V) = q0 * (V/V0)**q1``.

    Speziale et al. (2001), Equations (10)--(11), doi:10.1029/2000JB900318:
    ``gamma = gamma0 * exp(q0/q1 * ((V/V0)**q1 - 1))``.
    The Debye temperature follows ``d ln(theta) = -gamma d ln(V)`` with
    ``theta(V0) = theta0``. The author confirmed this integration prescription
    in correspondence supplied on 9 October 2026. Adaptive quadrature evaluates
    the relation here; it does not reproduce the author's trapezoidal grid.

    Volumes use J bar^-1 mol^-1, temperature K, pressure GPa, and ``n`` atoms
    per formula unit. Thermal energy is subtracted at ``Tr``. ``q1=0`` recovers
    the constant-q model; ``q0=0`` recovers constant gamma.
    """

    _constructor_configuration_names: tuple[str, ...] = ()

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        theta0: float,
        gamma0: float,
        q0: float,
        q1: float,
        n: float,
    ) -> None:
        self.q0 = validate_finite_scalar(q0, "q0")
        self.q1 = validate_finite_scalar(q1, "q1")
        super().__init__(rt_eos, Tr, theta0, gamma0, q=self.q0, n=n)
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is VariableQMieGruneisenDebye:
            from peritheos import _rust

            self._native = _rust.ThermalEos.variable_q_mie_gruneisen_debye(
                reference_native,
                self.Tr,
                self.theta0,
                self.gamma0,
                self.q0,
                self.q1,
                self.n,
            )

    def _gamma_at_log_ratio(self, logarithmic_ratio: NumericType) -> NumericType:
        with np.errstate(over="ignore", invalid="ignore"):
            exponent = (
                self.q0 * logarithmic_ratio
                if self.q1 == 0.0
                else self.q0 * np.expm1(self.q1 * logarithmic_ratio) / self.q1
            )
            result = self.gamma0 * np.exp(exponent)
        if not np.all(np.isfinite(result)):
            raise EosNumericalError("Gruneisen parameter is not finite")
        return result

    def gruneisen_parameter(
        self, V: NumericType, T: NumericType | None = None
    ) -> NumericType:
        """Return the Equation (11) Gruneisen parameter."""
        if getattr(self, "_native", None) is not None:
            volumes, temperatures = self._broadcast_state(
                V, self.Tr if T is None else T
            )
            return self._native_evaluate("gruneisen_parameter", volumes, temperatures)
        volumes = np.asarray(validate_volume(V), dtype=float)
        result = np.asarray(
            self._gamma_at_log_ratio(np.log(volumes / self.rt_eos.V0)), dtype=float
        )
        return self._scalar_or_array(result)

    def characteristic_temperature(self, V: NumericType) -> NumericType:
        """Integrate the defining Gruneisen relation from the reference volume."""
        if getattr(self, "_native", None) is not None:
            volumes = validate_volume(V)
            return self._native_evaluate("characteristic_temperature", volumes, self.Tr)
        volumes = np.asarray(validate_volume(V), dtype=float)
        log_ratios = np.log(volumes / self.rt_eos.V0)
        if self.q0 == 0.0:
            exponent = -self.gamma0 * log_ratios
        elif self.q1 == 0.0:
            exponent = -self.gamma0 * np.expm1(self.q0 * log_ratios) / self.q0
        else:
            unique, inverse = np.unique(log_ratios, return_inverse=True)
            integrals = np.array(
                [
                    quad(
                        self._gamma_at_log_ratio,
                        0.0,
                        float(upper),
                        epsabs=1e-11,
                        epsrel=1e-11,
                    )[0]
                    for upper in unique
                ]
            )
            exponent = -integrals[inverse].reshape(log_ratios.shape)
        with np.errstate(over="ignore", invalid="ignore"):
            result = self.theta0 * np.exp(exponent)
        if not np.all(np.isfinite(result)) or np.any(result <= 0):
            raise EosNumericalError(
                "Characteristic temperature must be positive and finite"
            )
        return self._scalar_or_array(np.asarray(result, dtype=float))
