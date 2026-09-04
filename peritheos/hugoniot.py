"""Shock Hugoniot equations of state.

This compatibility module backs the canonical ``peritheos.eos.hugoniot``
namespace introduced with the common :class:`peritheos.eos.EquationOfState`
hierarchy.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Union

import numpy as np
from numpy.typing import NDArray

from peritheos import _rust
from peritheos.eos import (
    EquationOfState,
    validate_finite_scalar,
    validate_positive_scalar,
)
from peritheos.errors import EosValidationError

NumericType = Union[float, NDArray[np.float64]]


def _native_evaluate(native: Any, quantity: str, values: Any) -> NumericType:
    array = np.asarray(values, dtype=float)
    if array.ndim == 0:
        return float(native.evaluate_scalar(quantity, float(array)))
    return np.asarray(native.evaluate_array(quantity, array), dtype=float)


@dataclass(frozen=True)
class HugoniotState:
    """Coupled mechanical state on a shock Hugoniot."""

    volume: NumericType
    density: NumericType
    pressure: NumericType
    particle_velocity: NumericType
    shock_velocity: NumericType
    specific_internal_energy_change: NumericType


class HugoniotBase(EquationOfState):
    """Base interface for a one-dimensional pressure-volume shock path."""

    V0: float
    rho0: float
    P0: float

    def parameter_values(self) -> dict[str, float]:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def pressure(self, V: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def volume(self, P: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def particle_velocity(self, V: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def shock_velocity(self, V: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def density(self, V: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def specific_internal_energy_change(self, V: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def tangent_modulus(self, V: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def shock_velocity_from_particle_velocity(self, up: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def pressure_from_particle_velocity(self, up: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def volume_from_particle_velocity(self, up: Any) -> NumericType:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def state_from_particle_velocity(self, up: Any) -> HugoniotState:
        raise NotImplementedError  # pragma: no cover - interface declaration

    def with_parameters(self, **updates: float) -> HugoniotBase:
        parameters = self.parameter_values()
        unknown = set(updates) - set(parameters)
        if unknown:
            raise EosValidationError(
                f"Unknown parameters for {type(self).__name__}: {sorted(unknown)}"
            )
        parameters.update({name: float(value) for name, value in updates.items()})
        return type(self)(**parameters)


class LinearUsUpHugoniot(HugoniotBase):
    """Shock Hugoniot defined by ``Us = c0 + s * up``.

    ``rho0`` is in g/cm^3, ``c0`` and particle/shock velocities are in km/s,
    and ``P0`` and returned pressures are in GPa. ``V0`` and evaluated volumes
    may use any consistent volume unit. Specific internal-energy changes are
    returned in MJ/kg.
    """

    def __init__(
        self,
        V0: float,
        rho0: float,
        c0: float,
        s: float,
        P0: float = 0.0,
    ) -> None:
        self.V0 = validate_positive_scalar(V0, "V0")
        self.rho0 = validate_positive_scalar(rho0, "rho0")
        self.c0 = validate_positive_scalar(c0, "c0")
        self.s = validate_positive_scalar(s, "s")
        self.P0 = validate_finite_scalar(P0, "P0")
        self._native = _rust.HugoniotEos.linear_us_up(
            self.V0, self.rho0, self.c0, self.s, self.P0
        )

    def parameter_values(self) -> dict[str, float]:
        return {
            "V0": self.V0,
            "rho0": self.rho0,
            "c0": self.c0,
            "s": self.s,
            "P0": self.P0,
        }

    def pressure(self, V: Any) -> NumericType:
        """Return Hugoniot pressure in GPa at compressed volume *V*."""
        return _native_evaluate(self._native, "pressure", V)

    def volume(self, P: Any) -> NumericType:
        """Return volume on the compression Hugoniot at pressure *P* in GPa."""
        return _native_evaluate(self._native, "volume", P)

    def particle_velocity(self, V: Any) -> NumericType:
        """Return particle velocity in km/s at compressed volume *V*."""
        return _native_evaluate(self._native, "particle_velocity", V)

    def shock_velocity(self, V: Any) -> NumericType:
        """Return shock velocity in km/s at compressed volume *V*."""
        return _native_evaluate(self._native, "shock_velocity", V)

    def density(self, V: Any) -> NumericType:
        """Return density in g/cm^3 at compressed volume *V*."""
        return _native_evaluate(self._native, "density", V)

    def specific_internal_energy_change(self, V: Any) -> NumericType:
        """Return the Rankine--Hugoniot internal-energy increase in MJ/kg."""
        return _native_evaluate(self._native, "specific_internal_energy_change", V)

    def tangent_modulus(self, V: Any) -> NumericType:
        """Return ``-V dP_H/dV`` along the Hugoniot in GPa."""
        return _native_evaluate(self._native, "tangent_modulus", V)

    def shock_velocity_from_particle_velocity(self, up: Any) -> NumericType:
        """Return shock velocity in km/s from particle velocity in km/s."""
        return _native_evaluate(
            self._native, "shock_velocity_from_particle_velocity", up
        )

    def pressure_from_particle_velocity(self, up: Any) -> NumericType:
        """Return pressure in GPa from particle velocity in km/s."""
        return _native_evaluate(self._native, "pressure_from_particle_velocity", up)

    def volume_from_particle_velocity(self, up: Any) -> NumericType:
        """Return volume from particle velocity in km/s."""
        return _native_evaluate(self._native, "volume_from_particle_velocity", up)

    def state(self, V: Any) -> HugoniotState:
        """Return all mechanical Hugoniot quantities at volume *V*."""
        return HugoniotState(
            volume=(
                float(V) if np.asarray(V).ndim == 0 else np.asarray(V, dtype=float)
            ),
            density=self.density(V),
            pressure=self.pressure(V),
            particle_velocity=self.particle_velocity(V),
            shock_velocity=self.shock_velocity(V),
            specific_internal_energy_change=self.specific_internal_energy_change(V),
        )

    def state_from_particle_velocity(self, up: Any) -> HugoniotState:
        """Return the coupled state at particle velocity *up* in km/s."""
        return self.state(self.volume_from_particle_velocity(up))


__all__ = ["HugoniotBase", "HugoniotState", "LinearUsUpHugoniot"]
