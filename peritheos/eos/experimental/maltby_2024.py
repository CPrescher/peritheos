"""Literal Maltby et al. (2024) equations; NOT a validated parameterization.

DOI 10.1063/5.0237497, equations 4 and 7--19, Tables 1 and 3.
The numerical CSM cutoff is not specified by the article, so the caller must
supply it. Table SI.1 extends to squared nearest-neighbor distance 64.
See docs/literature-reproductions/argon-maltby-2024.md for the unresolved
Table 8 discrepancy. This class is deliberately absent from EOSMAT dispatch.

Public volume: J/bar/mol (10 cm^3/mol); pressure: GPa; energy: J/mol.
The energy is unshifted: the fluid-reference adjustments in equations 21--23
are not included. They do not affect pressure or heat capacities.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

R = 8.314462618
N_A_ANGSTROM = 0.602214076  # number density in A^-3 for v=1 cm^3/mol
V_REF = 22.56
WEIGHTS = np.array([0.0261, 0.03784, 0.04512])
EINSTEIN_THETA = np.array([77.81, 550.0, 45.36])
EINSTEIN_GAMMA = np.array([6.221, 1.617e-6, 3.127])


@lru_cache(maxsize=32)
def fcc_shells(cutoff_squared: int) -> tuple[np.ndarray, np.ndarray]:
    """Exact fcc coordination shells, using half-cell integer coordinates."""
    if isinstance(cutoff_squared, bool) or not isinstance(cutoff_squared, int):
        raise ValueError("shell_cutoff_squared must be an integer")
    if not 1 <= cutoff_squared <= 256:
        raise ValueError("shell_cutoff_squared must lie in [1, 256]")
    extent = int(np.ceil(np.sqrt(2 * cutoff_squared)))
    counts: Counter[int] = Counter()
    for i in range(-extent, extent + 1):
        for j in range(-extent, extent + 1):
            for k in range(-extent, extent + 1):
                squared = i * i + j * j + k * k
                if (i + j + k) % 2 == 0 and 0 < squared <= 2 * cutoff_squared:
                    counts[squared // 2] += 1
    if cutoff_squared not in counts:
        raise ValueError("cutoff must coincide with an occupied fcc shell")
    distances = np.array(sorted(counts), dtype=float)
    populations = np.array([counts[int(m)] for m in distances], dtype=float)
    distances.setflags(write=False)
    populations.setflags(write=False)
    return distances, populations


def _debye_integral(x: float) -> float:
    """Source D3: one THIRD of the common normalized Debye function."""
    if x < 1e-3:
        return (1 - 3 * x / 8 + x * x / 20 - x**4 / 1680) / 3
    if x > 150:
        return np.pi**4 / (15 * x**3)
    return quad(lambda u: u**3 / np.expm1(u), 0, x, epsabs=1e-10)[0] / x**3


class Maltby2024Published:
    """Unvalidated literal published coefficients with an explicit CSM cutoff.

    ``shell_cutoff_squared`` is (r_cut/r_NN)^2, not the number of shells.
    Volume inversion is restricted to 8--26 cm^3/mol and positive bulk modulus;
    this is a numerical search interval, not a claim of phase stability.
    """

    def __init__(self, shell_cutoff_squared: int):
        self.shell_cutoff_squared = shell_cutoff_squared
        self._shells, self._populations = fcc_shells(shell_cutoff_squared)
        # Outer physical zero of Eq. (15); the Buckingham catastrophe's inner
        # zero is excluded. This is a derived quantity, not a fitted coefficient.
        self.sigma_angstrom = 3.802 * brentq(
            lambda x: 6 * np.exp(14.19 * (1 - x)) - 14.19 / x**6,
            0.8,
            1.0,
            xtol=1e-14,
        )

    @staticmethod
    def _state(volume: float, temperature: float) -> tuple[float, float]:
        if not np.isfinite(volume) or volume <= 0:
            raise ValueError("volume must be positive and finite")
        if not np.isfinite(temperature) or temperature < 0:
            raise ValueError("temperature must be nonnegative and finite")
        return volume * 10, temperature

    def _energy_and_dv(self, volume: float, temperature: float) -> tuple[float, float]:
        """Return A and dA/dv for v in cm^3/mol; derivative includes density."""
        v, t = self._state(volume, temperature)
        rnn = (np.sqrt(2) * v / N_A_ANGSTROM) ** (1 / 3)
        radii = np.sqrt(self._shells) * rnn
        rc = np.sqrt(self.shell_cutoff_squared) * rnn
        eps, steepness, rmin = 134.7, 14.19, 3.802
        decay = steepness / rmin
        rep = eps * 6 / (steepness - 6) * np.exp(steepness - decay * radii)
        att = eps * steepness / (steepness - 6) * (rmin / radii) ** 6
        pair = np.dot(self._populations, rep - att) / 2
        pair_dv = np.dot(self._populations, -decay * radii * rep / 3 + 2 * att) / (
            2 * v
        )
        # Eq. (13) integrated analytically, retaining moving r_cut(v).
        rep_rc = eps * 6 / (steepness - 6) * np.exp(steepness - decay * rc)
        b = eps * steepness / (steepness - 6) * rmin**6
        tail = (
            2
            * np.pi
            * N_A_ANGSTROM
            / v
            * (
                rep_rc * (rc**2 / decay + 2 * rc / decay**2 + 2 / decay**3)
                - b / (3 * rc**3)
            )
        )
        tail_dv = (
            -tail / v
            - 2 * np.pi * N_A_ANGSTROM / v**2 * (rep_rc - b / rc**6) * rc**3 / 3
        )
        pair += tail
        pair_dv += tail_dv
        correction = 3.202e5 * N_A_ANGSTROM / (eps * self.sigma_angstrom**6)
        energy = pair * (1 - correction / v)
        derivative = pair_dv * (1 - correction / v) + pair * correction / v**2
        zpv = 140 * np.exp(2.34 / 0.683 * (1 - (v / 19.7) ** 0.683))
        energy += zpv
        derivative -= zpv * 2.34 * (v / 19.7) ** 0.683 / v
        cold1 = 0.7204 * np.exp(-1.614 * (1 - v / V_REF))
        cold2 = -0.01943 * V_REF / v * np.exp(-27.64 * (1 - v / V_REF))
        energy += cold1 + cold2
        derivative += cold1 * 1.614 / V_REF + cold2 * (-1 / v + 27.64 / V_REF)
        if t > 0:
            theta_t = 92 + 10.4 * np.expm1(-0.02503 * t**2 - 0.0001568 * t**3)
            gamma_d = 2.563 * (v / V_REF) ** 0.2874
            theta = theta_t * np.exp(2.563 / 0.2874 * (1 - (v / V_REF) ** 0.2874))
            x = theta / t
            d3 = _debye_integral(x)
            weight = 1 - WEIGHTS.sum()
            energy += 3 * weight * t * (np.log(-np.expm1(-x)) - d3)
            derivative -= 9 * weight * t * gamma_d * d3 / v
            for w, theta0, gamma0 in zip(WEIGHTS, EINSTEIN_THETA, EINSTEIN_GAMMA):
                theta_i = theta0 * np.exp(gamma0 * (1 - v / V_REF))
                y = theta_i / t
                energy += 3 * w * t * np.log(-np.expm1(-y))
                occupation = 0 if y > 700 else 1 / np.expm1(y)
                derivative -= 3 * w * theta_i * occupation * gamma0 / V_REF
            anh = (
                -0.0004475
                * theta_t
                * (t / 92) ** 4
                / (1 + 2.041e-6 * (t / 92) ** 2)
                * np.exp(5.75e-7 * (v / V_REF - 1))
            )
            energy += anh
            derivative += anh * 5.75e-7 / V_REF
        result = R * energy, R * derivative
        if not all(np.isfinite(result)):
            raise ValueError("nonfinite Maltby model result")
        return result

    def molar_helmholtz_energy(self, volume: float, temperature: float) -> float:
        return self._energy_and_dv(volume, temperature)[0]

    def pressure(self, volume: float, temperature: float) -> float:
        return -self._energy_and_dv(volume, temperature)[1] * 1e-3

    def thermal_pressure_increment(self, volume: float, temperature: float) -> float:
        return self.pressure(volume, temperature) - self.pressure(volume, 300.0)

    def bulk_modulus(
        self, volume: float, temperature: float, relative_step: float = 1e-5
    ) -> float:
        if not np.isfinite(relative_step) or not 0 < relative_step < 1:
            raise ValueError("relative_step must lie in (0, 1)")
        h = volume * relative_step
        return (
            -volume
            * (
                self.pressure(volume + h, temperature)
                - self.pressure(volume - h, temperature)
            )
            / (2 * h)
        )

    def volume(self, pressure: float, temperature: float) -> float:
        if not np.isfinite(pressure):
            raise ValueError("pressure must be finite")
        self._state(1.0, temperature)
        grid = np.linspace(0.8, 2.6, 181)
        residual = [self.pressure(v, temperature) - pressure for v in grid]
        roots = []
        for lo, hi, flo, fhi in zip(grid[:-1], grid[1:], residual[:-1], residual[1:]):
            if flo >= 0 and fhi <= 0:
                root = brentq(
                    lambda v: self.pressure(v, temperature) - pressure,
                    lo,
                    hi,
                    xtol=1e-13,
                )
                if self.bulk_modulus(root, temperature) > 0:
                    roots.append(root)
        if not roots:
            raise ValueError(
                "no stable root in numerical 8--26 cm^3/mol search interval"
            )
        return min(roots, key=lambda v: abs(np.log(v / (V_REF / 10))))
