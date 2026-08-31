"""Simple thermal-pressure models used by published material EOS records."""

from __future__ import annotations

import numpy as np

from .. import (
    EosBase,
    NumericType,
    ThermalEOS,
    validate_finite_scalar,
    validate_positive_scalar,
)


class LinearThermalPressure(ThermalEOS):
    r"""Add a constant-:math:`\alpha K_T` thermal-pressure correction.

    The model is

    .. math::

        P(V,T) = P_{\mathrm{ref}}(V) + \alpha K_T (T-T_r).

    It is the form used for B2 KCl and KBr by Dewaele et al. (2012),
    doi:10.1103/PhysRevB.85.214105, equation (2); for B1/B2 KCl by Walker et
    al. (2002), doi:10.2138/am-2002-0701, equation BE1; and for the approximate
    finite-temperature platinum scale of Holmes et al. (1989),
    doi:10.1063/1.344177, equation (12). ``alpha_KT`` is in GPa/K and ``Tr``
    in K. Because the correction is independent of volume, the reference EOS
    may use any internally consistent volume unit.
    """

    def __init__(self, rt_eos: EosBase, Tr: float, alpha_KT: float):
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.alpha_KT = validate_finite_scalar(alpha_KT, "alpha_KT")

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        result = self.alpha_KT * (temperatures - self.Tr)
        # Broadcast against volume even though the equation is V-independent.
        result = np.broadcast_to(result, volumes.shape)
        return self._scalar_or_array(np.asarray(result, dtype=float))


class ThermalReferenceStateEOS(ThermalEOS):
    r"""Vary the reference volume and bulk modulus with temperature.

    The model uses ``V0(T) = V0(Tr) * exp(alpha0 * (T - Tr))`` and
    ``K0(T) = K0(Tr) + dK_dT * (T - Tr)`` before evaluating the reference
    isotherm. This is the mechanism behind the ``AlphaKT`` Dioptas interchange
    type and equations (1)--(3) of Bezacier et al. (2014),
    doi:10.1063/1.4894421. The reference EOS must expose reconstructable
    ``V0`` and ``K0`` parameters.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        alpha0: float,
        dK_dT: float,
    ) -> None:
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.alpha0 = validate_finite_scalar(alpha0, "alpha0")
        self.dK_dT = validate_finite_scalar(dK_dT, "dK_dT")
        parameters = rt_eos.parameter_values(include_reference=False)
        if "V0" not in parameters or "K0" not in parameters:
            raise ValueError("rt_eos must expose reconstructable V0 and K0")

    def _state_eos(self, temperature: float) -> EosBase:
        delta_temperature = temperature - self.Tr
        with np.errstate(over="ignore", under="ignore"):
            V0 = self.rt_eos.V0 * np.exp(self.alpha0 * delta_temperature)
        K0 = self.rt_eos.K0 + self.dK_dT * delta_temperature
        if not np.isfinite(V0) or V0 <= 0.0:
            raise ValueError("Temperature produces a non-positive reference volume")
        if not np.isfinite(K0) or K0 <= 0.0:
            raise ValueError("Temperature produces a non-positive bulk modulus")
        return self.rt_eos.with_parameters(V0=float(V0), K0=float(K0))

    def pressure(self, V: NumericType, T: NumericType) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        result = np.fromiter(
            (
                float(self._state_eos(float(temperature)).pressure(float(volume)))
                for volume, temperature in zip(volumes.flat, temperatures.flat)
            ),
            dtype=float,
            count=volumes.size,
        ).reshape(volumes.shape)
        return self._scalar_or_array(result)

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        result = np.asarray(
            self.pressure(volumes, temperatures), dtype=float
        ) - np.asarray(self.rt_eos.pressure(volumes), dtype=float)
        return self._scalar_or_array(result)

    def bulk_modulus(
        self, V: NumericType, T: NumericType, relative_step: float = 1.0e-6
    ) -> NumericType:
        # Keep the common ThermalEOS method contract even though this model can
        # delegate the derivative analytically to its temperature-shifted curve.
        validate_positive_scalar(relative_step, "relative_step")
        volumes, temperatures = self._broadcast_state(V, T)
        result = np.fromiter(
            (
                float(self._state_eos(float(temperature)).bulk_modulus(float(volume)))
                for volume, temperature in zip(volumes.flat, temperatures.flat)
            ),
            dtype=float,
            count=volumes.size,
        ).reshape(volumes.shape)
        return self._scalar_or_array(result)


__all__ = ["LinearThermalPressure", "ThermalReferenceStateEOS"]
