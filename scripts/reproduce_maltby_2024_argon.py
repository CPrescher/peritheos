#!/usr/bin/env python3
"""Audit literal Maltby coefficients without silently correcting the publication."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from peritheos.eos.experimental.maltby_2024 import Maltby2024Published

ROOT = Path(__file__).resolve().parents[1]


def independent_energy(v_cm3: float, t: float, cutoff: int = 64) -> float:
    """Separate direct SI quadrature, using the source potential and fcc sites.

    Do not call the production energy or pressure function. The potential tail
    is integrated in dimensionless nearest-neighbor coordinates so that its
    magnitude is large enough for the quadrature absolute tolerance.
    """
    na, kb = 6.02214076e23, 1.380649e-23
    v = v_cm3 * 1e-6
    rnn = (np.sqrt(2) * v / na) ** (1 / 3)
    sigma = 3.802e-10 * brentq(
        lambda x: 6 * np.exp(14.19 * (1 - x)) - 14.19 / x**6, 0.8, 1.0
    )
    extent = int(np.ceil(np.sqrt(2 * cutoff)))
    shell_counts = {}
    for i in range(-extent, extent + 1):
        for j in range(-extent, extent + 1):
            for k in range(-extent, extent + 1):
                m = (i * i + j * j + k * k) / 2
                if (i + j + k) % 2 == 0 and 0 < m <= cutoff:
                    shell_counts[m] = shell_counts.get(m, 0) + 1

    def potential_k(r):
        return (
            134.7
            / (14.19 - 6)
            * (6 * np.exp(14.19 * (1 - r / 3.802e-10)) - 14.19 * (3.802e-10 / r) ** 6)
            * (1 - 3.202e5 * 1e-90 * na / (134.7 * v * sigma**6))
        )

    cold = sum(n * potential_k(np.sqrt(m) * rnn) for m, n in shell_counts.items()) / 2
    cold += (
        2
        * np.pi
        * na
        / v
        * rnn**3
        * quad(
            lambda x: potential_k(x * rnn) * x * x,
            np.sqrt(cutoff),
            np.inf,
            epsabs=1e-10,
        )[0]
    )
    zpv = 140 * np.exp(2.34 / 0.683 * (1 - (v_cm3 / 19.7) ** 0.683))
    ac = 0.7204 * np.exp(
        -1.614 * (1 - v_cm3 / 22.56)
    ) - 0.01943 * 22.56 / v_cm3 * np.exp(-27.64 * (1 - v_cm3 / 22.56))
    if t == 0:
        return (cold + zpv + ac) * na * kb
    theta_t = 92 + 10.4 * np.expm1(-0.02503 * t * t - 0.0001568 * t**3)
    theta_d = theta_t * np.exp(2.563 / 0.2874 * (1 - (v_cm3 / 22.56) ** 0.2874))

    # Direct density-of-states integral, rather than the integrated Eq. (7).
    def log_occupation(x):
        return np.log(-np.expm1(-theta_d * x / t)) if x else 0

    ev = (
        9
        * (1 - 0.0261 - 0.03784 - 0.04512)
        * t
        * quad(lambda x: x * x * log_occupation(x), 0, 1, epsabs=1e-10)[0]
    )
    for w, th, g in zip(
        [0.0261, 0.03784, 0.04512], [77.81, 550, 45.36], [6.221, 1.617e-6, 3.127]
    ):
        ev += 3 * w * t * np.log(-np.expm1(-th * np.exp(g * (1 - v_cm3 / 22.56)) / t))
    anh = (
        -0.0004475
        * theta_t
        * (t / 92) ** 4
        / (1 + 2.041e-6 * (t / 92) ** 2)
        * np.exp(5.75e-7 * (v_cm3 / 22.56 - 1))
    )
    return (cold + zpv + ac + ev + anh) * na * kb


def properties(model, v, t):
    h = 0.01
    cv = (
        -t
        * (
            model.molar_helmholtz_energy(v, t + h)
            + model.molar_helmholtz_energy(v, t - h)
            - 2 * model.molar_helmholtz_energy(v, t)
        )
        / h**2
    )
    kt = model.bulk_modulus(v, t)
    alpha = (model.pressure(v, t + h) - model.pressure(v, t - h)) / (2 * h * kt)
    cp = cv + alpha**2 * kt * v * t * 1e4
    return dict(
        pressure_gpa=model.pressure(v, t),
        volume_at_1mpa_cm3_mol=10 * model.volume(0.001, t),
        unshifted_helmholtz_j_mol=model.molar_helmholtz_energy(v, t),
        cv_j_mol_k=cv,
        cp_j_mol_k=cp,
        alpha_mk_inverse=alpha * 1000,
        gamma=alpha * kt * v * 1e4 / cv,
        chi_t_gpa_inverse=1 / kt,
        chi_s_gpa_inverse=cv / cp / kt,
    )


def main():
    model = Maltby2024Published(64)
    candidate = properties(model, 2.397, 70)
    published = dict(
        pressure_gpa=0.001,
        volume_at_1mpa_cm3_mol=23.97,
        cv_j_mol_k=22.93,
        cp_j_mol_k=30.35,
        alpha_mk_inverse=1.684,
        gamma=2.745,
        chi_t_gpa_inverse=0.6412,
        chi_s_gpa_inverse=0.4844,
    )
    independent = []
    for v, t in [(23.97, 70), (18, 300), (12, 0)]:
        h = v * 1e-5
        pressure = (
            -(independent_energy(v + h, t) - independent_energy(v - h, t))
            / (2 * h)
            * 0.001
        )
        independent.append(
            dict(
                volume_cm3_mol=v,
                temperature_k=t,
                pressure_gpa=pressure,
                analytic_minus_quadrature_gpa=model.pressure(v / 10, t) - pressure,
            )
        )
    observations = []
    path = ROOT / "peritheos/data/datasets/argon-fcc-dewaele-2021-supplement.csv"
    with path.open() as f:
        for row in csv.DictReader(f):
            if not row["pressure_gpa"]:
                continue
            p, t = float(row["pressure_gpa"]), float(row["temperature_k"])
            # Comparison subset is intentionally NOT claimed as Maltby's fit.
            if 0 < p <= 16 and 0 <= t <= 300:
                cell = float(row["argon_a_angstrom"]) ** 3
                volume = cell * 0.602214076 / 4 / 10
                observations.append((p, model.pressure(volume, t) - p))
    residual = np.array(observations)
    report = dict(
        status="not_reproduced",
        equations="Literal Eqs.4,7-19; Tables1,3; explicit CSM cutoff squared64",
        sigma_angstrom=model.sigma_angstrom,
        published_table8=published,
        candidate_table8=candidate,
        relative_volume_discrepancy_percent=100
        * (candidate["volume_at_1mpa_cm3_mol"] / 23.97 - 1),
        independent_si_quadrature=independent,
        cutoff_sensitivity=[
            dict(cutoff_squared=c, **properties(Maltby2024Published(c), 2.397, 70))
            for c in [16, 36, 64, 100, 144, 256]
        ],
        experimental_comparison=dict(
            dataset="argon_fcc_dewaele_2021_supplement",
            selection="Every finite pressure in (0,16] GPa and T in [0,300] K; all runs, repeated rows retained. Not the published fit selection.",
            count=len(observations),
            rms_pressure_gpa=float(np.sqrt(np.mean(residual[:, 1] ** 2))),
            mean_absolute_relative_pressure_percent=float(
                np.mean(abs(residual[:, 1] / residual[:, 0])) * 100
            ),
        ),
        limitations=[
            "No independent refit: full weighted primary property observations not supplied in this paper/supplement.",
            "Numerical lattice cutoff unspecified; explicit sensitivity cases do not repair Table8.",
            "Fluid-reference Gibbs/entropy offsets not evaluated; source A,g,s cannot directly validate unshifted energy.",
            "No active EOSMAT or Studio curve until the sample-volume discrepancy is resolved.",
        ],
    )
    output = ROOT / "docs/data/argon-maltby-2024-reproduction.json"
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
