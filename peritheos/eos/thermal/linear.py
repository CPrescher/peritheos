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


class SecondOrderTaylorThermalPressure(ThermalEOS):
    r"""Add an absolute bivariate second-order thermal pressure.

    The model is

    .. math::

        P(V,T) = P_c(V) + c_0 + c_1\delta\eta + c_2\delta T
        + \frac{c_3}{2}\delta\eta^2 + \frac{c_4}{2}\delta T^2
        + \frac{c_5}{2}\delta\eta\delta T,

    where ``eta = 1 - V/V0``, ``delta_eta = eta - eta0``, and
    ``delta_T = T - Tr``. Unlike the reference-isotherm thermal models, this
    pressure is absolute: ``rt_eos`` is a cold curve and the thermal term is
    generally nonzero at ``Tr``. :meth:`thermal_pressure_increment` subtracts
    the value at ``Tr`` for DAC and ambient-isotherm operations.

    This is the pressure-scale form published for MgO by Luo et al. (2023),
    doi:10.1103/PhysRevB.107.134116, Appendix B, equations (B1)--(B2).
    ``c0``, ``c1``, and ``c3`` use GPa; ``c2`` and ``c5`` use GPa/K; and
    ``c4`` uses GPa/K^2. The volume convention is inherited from the cold
    curve because only ``V/V0`` enters the thermal expression.
    """

    def __init__(
        self,
        rt_eos: EosBase,
        Tr: float,
        eta0: float,
        c0: float,
        c1: float,
        c2: float,
        c3: float,
        c4: float,
        c5: float,
    ) -> None:
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.eta0 = validate_finite_scalar(eta0, "eta0")
        self.c0 = validate_finite_scalar(c0, "c0")
        self.c1 = validate_finite_scalar(c1, "c1")
        self.c2 = validate_finite_scalar(c2, "c2")
        self.c3 = validate_finite_scalar(c3, "c3")
        self.c4 = validate_finite_scalar(c4, "c4")
        self.c5 = validate_finite_scalar(c5, "c5")
        if not hasattr(rt_eos, "V0"):
            raise EosValidationError("rt_eos must expose V0")
        validate_positive_scalar(rt_eos.V0, "rt_eos.V0")
        reference_native = _native_for_exact_model(rt_eos)
        if (
            reference_native is not None
            and type(self) is SecondOrderTaylorThermalPressure
        ):
            from peritheos import _rust

            self._native = _rust.ThermalEos.second_order_taylor_thermal_pressure(
                reference_native,
                self.Tr,
                self.eta0,
                self.c0,
                self.c1,
                self.c2,
                self.c3,
                self.c4,
                self.c5,
            )

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return the absolute non-cold pressure contribution in GPa."""
        if hasattr(self, "_native"):
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        volumes, temperatures = self._broadcast_state(V, T)
        delta_eta = 1.0 - volumes / self.rt_eos.V0 - self.eta0
        delta_temperature = temperatures - self.Tr
        result = (
            self.c0
            + self.c1 * delta_eta
            + self.c2 * delta_temperature
            + 0.5 * self.c3 * delta_eta**2
            + 0.5 * self.c4 * delta_temperature**2
            + 0.5 * self.c5 * delta_eta * delta_temperature
        )
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def thermal_pressure_increment(self, V: NumericType, T: NumericType) -> NumericType:
        """Return pressure gained above the complete ``Tr`` isotherm."""
        if hasattr(self, "_native"):
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native,
                "thermal_pressure_increment",
                volumes,
                temperatures,
            )
        increment = np.asarray(self.thermal_pressure(V, T), dtype=float) - np.asarray(
            self.thermal_pressure(V, self.Tr), dtype=float
        )
        return self._scalar_or_array(increment)


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
    ``linear_reference_temperature`` instead defines the intercept at the
    reference temperature,
    :math:`\alpha(T)=\alpha_0+\alpha_1(T-T_r)`, as in Suzuki (2016),
    doi:10.2465/jmps.160719c.

    ``reference_volume_law="berman"`` applies the truncated quadratic form

    .. math::

        V_0(T) = V_0(T_r)\left[1 + \alpha_0(T-T_r)
        + \frac{1}{2}\alpha_1(T-T_r)^2\right],

    used by Berman (1988) and EosFit7 (Angel et al. 2014,
    doi:10.1515/zkri-2013-1711). Here :math:`\alpha_0` is the expansion
    coefficient at :math:`T_r`; the derivative of the truncated polynomial
    only approximately equals :math:`\alpha_0+\alpha_1(T-T_r)`.

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

    ``bulk_modulus_law="reciprocal_cubic"`` instead uses
    ``1/K0(T) = 1/K0(Tr) + beta1*(T-Tr) + beta2*(T**2-Tr**2)
    + beta3*(T**3-Tr**3)`` and requires ``dK_dT=0``. The reference
    pressure derivative can vary by ``kprime_log_coefficient*(T-Tr)*ln(T/Tr)``
    for supported three-parameter reference families. These are the laws of
    Hirose et al. (2008), Table 2, fit #2, with the intercept constrained by
    the stated reference bulk modulus. They define a mechanical P-V-T
    surface, not a caloric free-energy model.
    """

    _constructor_configuration_names = (
        "thermal_expansion_law",
        "reference_volume_law",
        "bulk_modulus_law",
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
        bulk_modulus_law: str = "linear_temperature",
        beta1: float = 0.0,
        beta2: float = 0.0,
        beta3: float = 0.0,
        kprime_log_coefficient: float = 0.0,
    ) -> None:
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.alpha0 = validate_finite_scalar(alpha0, "alpha0")
        self.dK_dT = validate_finite_scalar(dK_dT, "dK_dT")
        self.alpha1 = validate_finite_scalar(alpha1, "alpha1")
        self.beta1 = validate_finite_scalar(beta1, "beta1")
        self.beta2 = validate_finite_scalar(beta2, "beta2")
        self.beta3 = validate_finite_scalar(beta3, "beta3")
        self.kprime_log_coefficient = validate_finite_scalar(
            kprime_log_coefficient, "kprime_log_coefficient"
        )
        if bulk_modulus_law not in {"linear_temperature", "reciprocal_cubic"}:
            raise EosValidationError("Invalid bulk_modulus_law")
        self.bulk_modulus_law = bulk_modulus_law
        if bulk_modulus_law == "reciprocal_cubic" and self.dK_dT != 0.0:
            raise EosValidationError("dK_dT must be zero with cubic compressibility")
        if bulk_modulus_law == "linear_temperature" and any((beta1, beta2, beta3)):
            raise EosValidationError("beta coefficients require reciprocal_cubic")
        if self.kprime_log_coefficient and type(rt_eos).__name__ not in {
            "BM3",
            "Baonza",
            "Murnaghan",
            "Morse3",
            "NaturalStrain3",
            "SunMorse3",
            "SunMorse4",
            "Vinet",
        }:
            raise EosValidationError("Unsupported reference EOS for K0_prime shift")
        if thermal_expansion_law not in {
            "constant",
            "linear_temperature",
            "linear_reference_temperature",
        }:
            raise EosValidationError(
                "thermal_expansion_law must be 'constant', 'linear_temperature', "
                "or 'linear_reference_temperature'"
            )
        self.thermal_expansion_law = thermal_expansion_law
        if reference_volume_law not in {
            "integrated_expansivity",
            "linear_temperature",
            "berman",
        }:
            raise EosValidationError(
                "reference_volume_law must be 'integrated_expansivity' or "
                "'linear_temperature' or 'berman'"
            )
        self.reference_volume_law = reference_volume_law
        if thermal_expansion_law == "constant" and self.alpha1 != 0.0:
            raise EosValidationError(
                "alpha1 must be zero for constant thermal expansion"
            )
        if reference_volume_law == "linear_temperature" and (
            thermal_expansion_law != "constant" or self.alpha1 != 0.0
        ):
            raise EosValidationError(
                "linear_temperature reference volume requires constant thermal "
                "expansion configuration and alpha1=0"
            )
        if (
            reference_volume_law == "berman"
            and thermal_expansion_law != "linear_temperature"
        ):
            raise EosValidationError(
                "berman reference volume requires linear_temperature thermal expansion"
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
                self.bulk_modulus_law,
                self.beta1,
                self.beta2,
                self.beta3,
                self.kprime_log_coefficient,
            )

    def configuration_values(self) -> dict[str, str]:
        """Preserve the established configuration for the default modulus law."""
        values = super().configuration_values()
        if self.bulk_modulus_law == "linear_temperature":
            values.pop("bulk_modulus_law")
        return values

    def _state_eos(self, temperature: float) -> EosBase:
        delta_temperature = temperature - self.Tr
        if self.reference_volume_law == "linear_temperature":
            V0 = self.rt_eos.V0 * (1.0 + self.alpha0 * delta_temperature)
        elif self.reference_volume_law == "berman":
            V0 = self.rt_eos.V0 * (
                1.0
                + self.alpha0 * delta_temperature
                + 0.5 * self.alpha1 * delta_temperature**2
            )
        else:
            exponent = self.alpha0 * delta_temperature
            if self.thermal_expansion_law == "linear_temperature":
                exponent += 0.5 * self.alpha1 * (temperature**2 - self.Tr**2)
            elif self.thermal_expansion_law == "linear_reference_temperature":
                exponent += 0.5 * self.alpha1 * delta_temperature**2
            with np.errstate(over="ignore", under="ignore"):
                V0 = self.rt_eos.V0 * np.exp(exponent)
        if self.bulk_modulus_law == "reciprocal_cubic":
            compressibility = (
                1.0 / self.rt_eos.K0
                + self.beta1 * delta_temperature
                + self.beta2 * (temperature**2 - self.Tr**2)
                + self.beta3 * (temperature**3 - self.Tr**3)
            )
            if not np.isfinite(compressibility) or compressibility <= 0.0:
                raise EosValidationError(
                    "Temperature produces non-positive compressibility"
                )
            K0 = 1.0 / compressibility
        else:
            K0 = self.rt_eos.K0 + self.dK_dT * delta_temperature
        if not np.isfinite(V0) or V0 <= 0.0:
            raise EosValidationError(
                "Temperature produces a non-positive reference volume"
            )
        if not np.isfinite(K0) or K0 <= 0.0:
            raise EosValidationError("Temperature produces a non-positive bulk modulus")
        parameters = {"V0": float(V0), "K0": float(K0)}
        if self.kprime_log_coefficient:
            parameters["K0_prime"] = (
                self.rt_eos.K0_prime
                + self.kprime_log_coefficient
                * delta_temperature
                * np.log(temperature / self.Tr)
            )
        return self.rt_eos.with_parameters(**parameters)

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
