"""Finite-strain models for jointly measured density and acoustic velocities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from peritheos.eos import validate_finite_scalar, validate_positive_scalar
from peritheos.errors import EosNumericalError, EosValidationError


def _validated_density(density: Any) -> NDArray[np.float64]:
    values = np.asarray(density, dtype=float)
    if not np.all(np.isfinite(values)):
        raise EosValidationError("Density must be finite")
    if np.any(values <= 0.0):
        raise EosValidationError("Density must be greater than zero")
    return values


@dataclass(frozen=True)
class EulerianFiniteStrainAcoustic:
    """Third-order acoustic finite-strain model of Davies and Dziewonski.

    Density is expressed in g/cm^3 and velocity in km/s, so ``rho * velocity^2``
    is numerically a modulus in GPa. ``K_S0`` and ``G0`` are the reference-state
    adiabatic bulk and shear moduli; the primed quantities are their pressure
    derivatives.
    """

    rho0: float
    K_S0: float
    K_S0_prime: float
    G0: float
    G0_prime: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "rho0", validate_positive_scalar(self.rho0, "rho0"))
        object.__setattr__(self, "K_S0", validate_positive_scalar(self.K_S0, "K_S0"))
        object.__setattr__(
            self, "K_S0_prime", validate_finite_scalar(self.K_S0_prime, "K_S0_prime")
        )
        object.__setattr__(self, "G0", validate_positive_scalar(self.G0, "G0"))
        object.__setattr__(
            self, "G0_prime", validate_finite_scalar(self.G0_prime, "G0_prime")
        )

    def strain(self, density: Any) -> NDArray[np.float64]:
        """Return compression-positive Eulerian strain from density."""
        values = _validated_density(density)
        return 0.5 * ((values / self.rho0) ** (2.0 / 3.0) - 1.0)

    def moduli(self, density: Any) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return adiabatic bulk and shear moduli in GPa."""
        strain = self.strain(density)
        factor = (1.0 + 2.0 * strain) ** 2.5
        bulk = factor * (
            self.K_S0 + (3.0 * self.K_S0 * self.K_S0_prime - 5.0 * self.K_S0) * strain
        )
        shear = factor * (
            self.G0 + (3.0 * self.K_S0 * self.G0_prime - 5.0 * self.G0) * strain
        )
        return bulk, shear

    def velocities(
        self, density: Any
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return compressional and shear velocities in km/s."""
        values = _validated_density(density)
        bulk, shear = self.moduli(values)
        longitudinal = bulk + 4.0 * shear / 3.0
        if np.any(longitudinal <= 0.0) or np.any(shear <= 0.0):
            raise EosNumericalError(
                "Acoustic finite-strain model produced a non-positive modulus"
            )
        return np.sqrt(longitudinal / values), np.sqrt(shear / values)

    def parameter_values(self) -> dict[str, float]:
        """Return model parameters in stable constructor order."""
        return {
            "rho0": self.rho0,
            "K_S0": self.K_S0,
            "K_S0_prime": self.K_S0_prime,
            "G0": self.G0,
            "G0_prime": self.G0_prime,
        }


__all__ = ["EulerianFiniteStrainAcoustic"]
