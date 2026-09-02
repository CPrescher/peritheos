"""Simple thermal-pressure models used by published material EOS records."""

from __future__ import annotations

import numpy as np

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
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is LinearThermalPressure:
            from peritheos import _rust

            self._native = _rust.ThermalEos.linear_thermal_pressure(
                reference_native, self.Tr, self.alpha_KT
            )

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        if hasattr(self, "_native"):
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        result = self.alpha_KT * (temperatures - self.Tr)
        # Broadcast against volume even though the equation is V-independent.
        result = np.broadcast_to(result, volumes.shape)
        return self._scalar_or_array(np.asarray(result, dtype=float))


class LogVolumeThermalPressure(ThermalEOS):
    r"""Add a linear-in-temperature pressure with a logarithmic volume slope.

    The model is

    .. math::

        P(V,T) = P_{\mathrm{ref}}(V) +
        \left[\alpha K_{T,r} +
        \left(\frac{\partial K_T}{\partial T}\right)_V
        \ln\left(\frac{V_0}{V}\right)\right](T-T_r).

    It is the mechanism in Anderson, Isaak, and Yamamoto (1989),
    doi:10.1063/1.342969, equations (26)--(29). ``alpha_KT_ref`` and
    ``dK_dT_V`` are in GPa/K, ``Tr`` is in K, and the composed reference EOS
    must expose ``V0`` in the same volume convention as ``V``.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        alpha_KT_ref: float,
        dK_dT_V: float,
    ) -> None:
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.alpha_KT_ref = validate_finite_scalar(alpha_KT_ref, "alpha_KT_ref")
        self.dK_dT_V = validate_finite_scalar(dK_dT_V, "dK_dT_V")
        if not hasattr(rt_eos, "V0"):
            raise EosValidationError("rt_eos must expose V0")
        validate_positive_scalar(rt_eos.V0, "rt_eos.V0")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is LogVolumeThermalPressure:
            from peritheos import _rust

            self._native = _rust.ThermalEos.log_volume_thermal_pressure(
                reference_native,
                self.Tr,
                self.alpha_KT_ref,
                self.dK_dT_V,
            )

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        if hasattr(self, "_native"):
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        slope = self.alpha_KT_ref + self.dK_dT_V * np.log(self.rt_eos.V0 / volumes)
        result = slope * (temperatures - self.Tr)
        return self._scalar_or_array(np.asarray(result, dtype=float))


class ThermalReferenceStateEOS(ThermalEOS):
    r"""Vary the reference volume and bulk modulus with temperature.

    By default, the reference volume is obtained by integrating a selected
    volumetric thermal-expansion law,

    .. math::

        V_0(T) = V_0(T_r)\exp\left[\int_{T_r}^{T}\alpha(T')\,dT'\right],

    while ``K0(T) = K0(Tr) + dK_dT * (T - Tr)``. The default ``constant``
    law uses :math:`\alpha(T)=\alpha_0`. ``linear_temperature`` uses
    :math:`\alpha(T)=\alpha_0+\alpha_1T`, following Martinez et al. (1996),
    doi:10.2138/am-1996-5-608, equations (2), (4), and (5).

    ``reference_volume_law="linear_temperature"`` instead applies the direct
    relation

    .. math::

        V_0(T) = V_0(T_r)[1 + \alpha_0(T-T_r)],

    where :math:`\alpha_0` is a mean expansion coefficient over the represented
    interval rather than a constant instantaneous expansivity. This is
    Martinez et al. (1996), equation (3). The default
    ``integrated_expansivity`` reference-volume law is the mechanism behind
    equations (1)--(3) of Bezacier et al. (2014),
    doi:10.1063/1.4894421. The reference EOS must expose reconstructable
    ``V0`` and ``K0`` parameters.
    """

    _constructor_configuration_names = (
        "thermal_expansion_law",
        "reference_volume_law",
    )

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        alpha0: float,
        dK_dT: float,
        alpha1: float = 0.0,
        thermal_expansion_law: str = "constant",
        reference_volume_law: str = "integrated_expansivity",
    ) -> None:
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.alpha0 = validate_finite_scalar(alpha0, "alpha0")
        self.dK_dT = validate_finite_scalar(dK_dT, "dK_dT")
        self.alpha1 = validate_finite_scalar(alpha1, "alpha1")
        if thermal_expansion_law not in {"constant", "linear_temperature"}:
            raise EosValidationError(
                "thermal_expansion_law must be 'constant' or 'linear_temperature'"
            )
        self.thermal_expansion_law = thermal_expansion_law
        if reference_volume_law not in {
            "integrated_expansivity",
            "linear_temperature",
        }:
            raise EosValidationError(
                "reference_volume_law must be 'integrated_expansivity' or "
                "'linear_temperature'"
            )
        self.reference_volume_law = reference_volume_law
        if thermal_expansion_law == "constant" and self.alpha1 != 0.0:
            raise EosValidationError("alpha1 must be zero for constant thermal expansion")
        if reference_volume_law == "linear_temperature" and (
            thermal_expansion_law != "constant" or self.alpha1 != 0.0
        ):
            raise EosValidationError(
                "linear_temperature reference volume requires constant thermal "
                "expansion configuration and alpha1=0"
            )
        parameters = rt_eos.parameter_values(include_reference=False)
        if "V0" not in parameters or "K0" not in parameters:
            raise EosValidationError("rt_eos must expose reconstructable V0 and K0")
        reference_native = _native_for_exact_model(rt_eos)
        if reference_native is not None and type(self) is ThermalReferenceStateEOS:
            from peritheos import _rust

            self._native = _rust.ThermalEos.thermal_reference_state(
                reference_native,
                self.Tr,
                self.alpha0,
                self.dK_dT,
                self.alpha1,
                self.thermal_expansion_law,
                self.reference_volume_law,
            )

    def _state_eos(self, temperature: float) -> EosBase:
        delta_temperature = temperature - self.Tr
        if self.reference_volume_law == "linear_temperature":
            V0 = self.rt_eos.V0 * (1.0 + self.alpha0 * delta_temperature)
        else:
            exponent = self.alpha0 * delta_temperature
            if self.thermal_expansion_law == "linear_temperature":
                exponent += 0.5 * self.alpha1 * (temperature**2 - self.Tr**2)
            with np.errstate(over="ignore", under="ignore"):
                V0 = self.rt_eos.V0 * np.exp(exponent)
        K0 = self.rt_eos.K0 + self.dK_dT * delta_temperature
        if not np.isfinite(V0) or V0 <= 0.0:
            raise EosValidationError("Temperature produces a non-positive reference volume")
        if not np.isfinite(K0) or K0 <= 0.0:
            raise EosValidationError("Temperature produces a non-positive bulk modulus")
        return self.rt_eos.with_parameters(V0=float(V0), K0=float(K0))

    def pressure(self, V: NumericType, T: NumericType) -> NumericType:
        if hasattr(self, "_native"):
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native, "pressure", volumes, temperatures
            )
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
        if hasattr(self, "_native"):
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
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
        if hasattr(self, "_native") and relative_step == 1.0e-6:
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native, "bulk_modulus", volumes, temperatures
            )
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


__all__ = [
    "LinearThermalPressure",
    "LogVolumeThermalPressure",
    "ThermalReferenceStateEOS",
]
