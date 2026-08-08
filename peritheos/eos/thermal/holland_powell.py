"""Holland-Powell thermal modified Tait equation of state."""

from __future__ import annotations

import numpy as np
from scipy.constants import R

from peritheos.eos import (
    NumericType,
    ThermalEOS,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)
from peritheos.eos.rt import ModifiedTait


class ThermalModifiedTait(ThermalEOS):
    """Holland-Powell thermal modified Tait EOS.

    The model combines a modified Tait reference isotherm with an Einstein
    thermal pressure for which ``alpha * K / C_V`` is constant.

    Reference: Holland, T. J. B. & Powell, R. (2011), Journal of Metamorphic
    Geology 29, 333-383, doi:10.1111/j.1525-1314.2010.00923.x.
    """

    def __init__(
        self,
        rt_eos: ModifiedTait,
        Tr: float,
        theta: float,
        alpha0: float,
        n: float,
    ) -> None:
        if not isinstance(rt_eos, ModifiedTait):
            raise TypeError("ThermalModifiedTait requires a ModifiedTait EOS")
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.theta = validate_positive_scalar(theta, "theta")
        self.alpha0 = validate_finite_scalar(alpha0, "alpha0")
        self.n = validate_positive_scalar(n, "n")
        self._cv0 = float(self._einstein_heat_capacity(self.Tr))
        self._pressure_factor = self.alpha0 * self.rt_eos.K0 / self._cv0

    def _einstein_energy(self, T: NumericType) -> NumericType:
        temperatures = np.asarray(T, dtype=float)
        ratio = self.theta / temperatures
        return 3.0 * self.n * R * self.theta * np.exp(-ratio) / (-np.expm1(-ratio))

    def _einstein_heat_capacity(self, T: NumericType) -> NumericType:
        temperatures = np.asarray(T, dtype=float)
        ratio = self.theta / temperatures
        decay = np.exp(-ratio)
        return 3.0 * self.n * R * ratio**2 * decay / (1.0 - decay) ** 2

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Return thermal pressure relative to ``Tr`` in GPa."""
        volumes, temperatures = self._broadcast_state(V, T)
        pressure = self._pressure_factor * (
            self._einstein_energy(temperatures) - self._einstein_energy(self.Tr)
        )
        pressure = np.broadcast_to(pressure, volumes.shape)
        return self._scalar_or_array(np.asarray(pressure, dtype=float))

    def molar_heat_capacity_v(self, V: NumericType, T: NumericType) -> NumericType:
        """Return Einstein constant-volume heat capacity in J mol^-1 K^-1."""
        volumes, temperatures = self._broadcast_state(V, T)
        result = np.broadcast_to(
            self._einstein_heat_capacity(temperatures), volumes.shape
        )
        return self._scalar_or_array(np.asarray(result, dtype=float))

    def gruneisen_parameter(
        self, V: NumericType, T: NumericType | None = None
    ) -> NumericType:
        """Return gamma implied by the constant ``alpha K / C_V`` model."""
        if T is None:
            volumes = np.asarray(validate_volume(V), dtype=float)
        else:
            volumes, _ = self._broadcast_state(V, T)
        result = volumes * self._pressure_factor * 1.0e4
        return self._scalar_or_array(np.asarray(result, dtype=float))


HollandPowell2011 = ThermalModifiedTait
