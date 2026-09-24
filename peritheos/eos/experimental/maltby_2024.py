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
from dataclasses import dataclass, fields
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

R = 8.314462618
N_A_ANGSTROM = 0.602214076  # number density in A^-3 for v=1 cm^3/mol


@dataclass(frozen=True)
class Maltby2024Parameters:
    """Tables 1 and 3; explicit trial changes do not modify published defaults.

    No parameter covariance is available. Construction validates numerical
    domains, not physical accuracy or phase stability. Fluid-reference offsets
    are excluded because they do not affect the evaluated derivatives.
    """

    epsilon_k: float = 134.7
    alpha_r: float = 14.19
    rmin_a: float = 3.802
    lambda_k_a9: float = 3.202e5
    z1_k: float = 140.0
    z2: float = 2.34
    z3: float = 0.683
    z4_cm3: float = 19.7
    vref_cm3: float = 22.56
    theta_d0_k: float = 92.0
    a_d_k: float = 10.4
    b_d: float = 0.02503
    c_d: float = 0.0001568
    gamma_d0: float = 2.563
    q_d: float = 0.2874
    b1: float = -0.0004475
    b2: float = 2.041e-6
    b3: float = 5.75e-7
    c1_k: float = 0.7204
    c2: float = -1.614
    c3_k: float = -0.01943
    c4: float = -27.64
    a0: float = 0.0261
    a1: float = 0.03784
    a2: float = 0.04512
    theta0_k: float = 77.81
    theta1_k: float = 550.0
    theta2_k: float = 45.36
    gamma0: float = 6.221
    gamma1: float = 1.617e-6
    gamma2: float = 3.127

    def __post_init__(self):
        if not all(np.isfinite(getattr(self, f.name)) for f in fields(self)):
            raise ValueError("all trial coefficients must be finite")
        positive = (
            "epsilon_k",
            "rmin_a",
            "z1_k",
            "z2",
            "z3",
            "z4_cm3",
            "vref_cm3",
            "theta_d0_k",
            "gamma_d0",
            "q_d",
            "theta0_k",
            "theta1_k",
            "theta2_k",
            "gamma0",
            "gamma1",
            "gamma2",
        )
        if any(getattr(self, name) <= 0 for name in positive):
            raise ValueError("length, temperature and exponent scales must be positive")
        if self.alpha_r <= 6 or self.lambda_k_a9 < 0:
            raise ValueError(
                "trial Buckingham potential requires alpha_r > 6 and lambda >= 0"
            )
        if self.alpha_r - 6 + 7 * np.log(6 / self.alpha_r) <= 0:
            raise ValueError("trial Buckingham potential has no outer zero")
        if (
            not 0 <= self.a_d_k < self.theta_d0_k
            or min(self.b_d, self.c_d, self.b2) < 0
        ):
            raise ValueError("invalid Debye temperature or anharmonic denominator")
        if min(self.a0, self.a1, self.a2) < 0 or self.a0 + self.a1 + self.a2 >= 1:
            raise ValueError(
                "vibrational weights must be nonnegative and sum to less than one"
            )


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

    _parameters = Maltby2024Parameters()

    @property
    def parameters(self) -> Maltby2024Parameters:
        return self._parameters

    def __init__(self, shell_cutoff_squared: int):
        self.shell_cutoff_squared = shell_cutoff_squared
        self._shells, self._populations = fcc_shells(shell_cutoff_squared)
        # Outer physical zero of Eq. (15); the Buckingham catastrophe's inner
        # zero is excluded. This is a derived quantity, not a fitted coefficient.
        p = self.parameters
        # The maximum of log(repulsion/attraction) is at x=6/alpha.
        # This bracket isolates the outer zero for every allowed alpha.
        self.sigma_angstrom = p.rmin_a * brentq(
            lambda x: np.log(6 / p.alpha_r) + p.alpha_r * (1 - x) + 6 * np.log(x),
            6 / p.alpha_r,
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
        p = self.parameters
        rnn = (np.sqrt(2) * v / N_A_ANGSTROM) ** (1 / 3)
        radii = np.sqrt(self._shells) * rnn
        rc = np.sqrt(self.shell_cutoff_squared) * rnn
        eps, steepness, rmin = p.epsilon_k, p.alpha_r, p.rmin_a
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
        correction = p.lambda_k_a9 * N_A_ANGSTROM / (eps * self.sigma_angstrom**6)
        energy = pair * (1 - correction / v)
        derivative = pair_dv * (1 - correction / v) + pair * correction / v**2
        zpv = p.z1_k * np.exp(p.z2 / p.z3 * (1 - (v / p.z4_cm3) ** p.z3))
        energy += zpv
        derivative -= zpv * p.z2 * (v / p.z4_cm3) ** p.z3 / v
        cold1 = p.c1_k * np.exp(p.c2 * (1 - v / p.vref_cm3))
        cold2 = p.c3_k * p.vref_cm3 / v * np.exp(p.c4 * (1 - v / p.vref_cm3))
        energy += cold1 + cold2
        derivative += -cold1 * p.c2 / p.vref_cm3 + cold2 * (-1 / v - p.c4 / p.vref_cm3)
        if t > 0:
            theta_t = p.theta_d0_k + p.a_d_k * np.expm1(-p.b_d * t**2 - p.c_d * t**3)
            gamma_d = p.gamma_d0 * (v / p.vref_cm3) ** p.q_d
            theta = theta_t * np.exp(
                p.gamma_d0 / p.q_d * (1 - (v / p.vref_cm3) ** p.q_d)
            )
            x = theta / t
            d3 = _debye_integral(x)
            weight = 1 - p.a0 - p.a1 - p.a2
            energy += 3 * weight * t * (np.log(-np.expm1(-x)) - d3)
            derivative -= 9 * weight * t * gamma_d * d3 / v
            for w, theta0, gamma0 in [
                (p.a0, p.theta0_k, p.gamma0),
                (p.a1, p.theta1_k, p.gamma1),
                (p.a2, p.theta2_k, p.gamma2),
            ]:
                theta_i = theta0 * np.exp(gamma0 * (1 - v / p.vref_cm3))
                y = theta_i / t
                energy += 3 * w * t * np.log(-np.expm1(-y))
                occupation = 0 if y > 700 else 1 / np.expm1(y)
                derivative -= 3 * w * theta_i * occupation * gamma0 / p.vref_cm3
            anh = (
                p.b1
                * theta_t
                * (t / p.theta_d0_k) ** 4
                / (1 + p.b2 * (t / p.theta_d0_k) ** 2)
                * np.exp(p.b3 * (v / p.vref_cm3 - 1))
            )
            energy += anh
            derivative += anh * p.b3 / p.vref_cm3
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
        return min(
            roots, key=lambda v: abs(np.log(v / (self.parameters.vref_cm3 / 10)))
        )


class Maltby2024Trial(Maltby2024Published):
    """Explicit research-only parameter trial; not a published or validated EOS.

    Use separate provenance for every fitted parameter set. This class is not
    exposed through EOSMAT and has no native dispatcher or default fitted values.
    """

    def __init__(self, shell_cutoff_squared: int, parameters: Maltby2024Parameters):
        if not isinstance(parameters, Maltby2024Parameters):
            raise TypeError("parameters must be Maltby2024Parameters")
        self._parameters = parameters
        super().__init__(shell_cutoff_squared)
