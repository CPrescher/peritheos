"""Finger et al. (1981) argon EOS, DOI 10.1063/1.92597, Eqs. (1)-(6).

This study-specific evaluator is independent of the generic eosmat factories.
Volumes are four-atom fcc conventional-cell angstrom^3, pressures GPa and
absolute temperatures K. The cold curve is second-order *Murnaghan*, not BM2.
See docs/literature-reproductions/argon-fcc-finger-1981.md for source caveats.
"""

from __future__ import annotations

import numpy as np
from scipy.constants import Avogadro, R
from scipy.integrate import quad
from scipy.optimize import brentq

from peritheos.eos import EquationOfState
from peritheos.errors import EosValidationError

CELL_TO_MOLAR = Avogadro * 1e-24 / 4  # A^3/cell -> cm^3/mol of Ar


def _scalar_or_array(value):
    value = np.asarray(value, dtype=float)
    return float(value) if value.ndim == 0 else value


class Finger1981Argon(EquationOfState):
    """Published Table II coefficients, including the printed static offset.

    ``K0_prime`` and ``K0_double_prime`` may be varied for diagnostic refits;
    the latter is in GPa^-1 (Table II reports -0.040 +/- 0.010 kbar^-1).
    V0 and K0 are adopted 0 K low-temperature constraints, not 293 K values.
    Only the 293 +/- 1 K observations were measured in this study; computed
    isotherms at other temperatures are predictions. No coefficient refit is
    silently substituted for the published values.
    """

    V0_molar = 22.557
    V0 = V0_molar / CELL_TO_MOLAR
    K0 = 2.3701
    P_static_offset = -0.10289
    theta0 = 93.3
    gamma1 = 2.20
    temperature_ref = 293.0

    def __init__(self, K0_prime=6.97, K0_double_prime=-0.40):
        self.K0_prime = float(K0_prime)
        self.K0_double_prime = float(K0_double_prime)
        discriminant = self.K0_prime**2 - 2 * self.K0 * self.K0_double_prime
        if not np.isfinite(discriminant) or discriminant <= 0:
            raise EosValidationError("Second-order Murnaghan requires real positive q")
        self.q = float(np.sqrt(discriminant))

    def parameter_values(self):
        return {
            name: float(getattr(self, name))
            for name in (
                "V0",
                "K0",
                "K0_prime",
                "K0_double_prime",
                "P_static_offset",
                "theta0",
                "gamma1",
            )
        }

    def pressure_components(self, V, T=293.0):
        """Return static, zero-point and thermal contributions separately.

        Zero K is accepted as the analytic Debye limit. The printed Eq. (1)
        has an apparent extra equals sign before P_T; the additive energy
        decomposition described in the text and Eqs. (2)-(6) is used.
        """
        try:
            v, t = np.broadcast_arrays(np.asarray(V, float), np.asarray(T, float))
        except ValueError as error:
            raise EosValidationError("V and T must broadcast") from error
        if np.any(~np.isfinite(v)) or np.any(v <= 0):
            raise EosValidationError("Volume must be finite and positive")
        if np.any(~np.isfinite(t)) or np.any(t < 0):
            raise EosValidationError("Temperature must be finite and nonnegative")
        vm = v * CELL_TO_MOLAR
        x = vm / self.V0_molar
        y = np.expm1(-self.q * np.log(x))
        denominator = self.q * (y + 2) - self.K0_prime * y
        if np.any(denominator <= 0) or np.any(~np.isfinite(denominator)):
            raise EosValidationError("Volume is outside the regular Murnaghan branch")
        static = self.P_static_offset + 2 * self.K0 * y / denominator
        gamma = 0.5 + self.gamma1 * x
        theta = self.theta0 / np.sqrt(x) * np.exp(self.gamma1 * (1 - x))
        zero_point = 9 / 8 * gamma * R * theta / vm * 1e-3
        d3 = np.zeros_like(t)
        for index in np.ndindex(t.shape):
            if t[index] > 0:
                z = float(theta[index] / t[index])
                # Integral tail beyond 700 is below floating-point relevance.
                integral = quad(
                    lambda u: u**3 / np.expm1(u) if u else 0.0,
                    0,
                    min(z, 700),
                    epsabs=1e-11,
                    epsrel=1e-11,
                )[0]
                d3[index] = 3 * integral / z**3
        thermal = 3 * gamma * R * t * d3 / vm * 1e-3
        return {
            key: _scalar_or_array(value)
            for key, value in (
                ("static_gpa", static),
                ("zero_point_gpa", zero_point),
                ("thermal_gpa", thermal),
            )
        }

    def pressure(self, V, T=293.0):
        """Total pressure (GPa); extrapolation is not experimental validation."""
        return sum(self.pressure_components(V, T).values())

    def volume(self, P, T=293.0):
        """Invert on the compressed branch V0/2 <= V <= V0, in A^3/cell.

        This covers the observations and avoids selecting another mathematical
        branch when extrapolating the rational cold curve. Targets outside this
        bracket raise EosValidationError instead of returning an arbitrary root.
        """
        try:
            p, t = np.broadcast_arrays(np.asarray(P, float), np.asarray(T, float))
        except ValueError as error:
            raise EosValidationError("P and T must broadcast") from error
        if np.any(~np.isfinite(p)):
            raise EosValidationError("Pressure must be finite")
        result = np.empty_like(p)
        for index in np.ndindex(p.shape):
            try:
                result[index] = brentq(
                    lambda v: self.pressure(v, float(t[index])) - float(p[index]),
                    self.V0 / 2,
                    self.V0,
                    xtol=1e-11,
                )
            except ValueError as error:
                raise EosValidationError("No volume root in [V0/2, V0]") from error
        return _scalar_or_array(result)
