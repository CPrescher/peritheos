"""Test published Yokoo thermal increments without fitting pressure corrections."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from scripts.fit_yokoo_2009_gold_thermal import check_reconstruction, quadrature

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/yokoo-2009-pressure-conventions.json"
ELECTRONIC = "tsuchiya-kawamura-2002-table1-electronic-pressure.csv"
PYTHEOS_COMMIT = "b524b3439dfacfee8fa46a3e0718d50520ac37c5"
PYTHEOS_POLYNOMIALS = {
    "gold": [-0.00021606, -4.3795e-6, 1.4526e-8, 7.8072e-14],
    "platinum": [0.011316, -5.6486e-7, 2.67e-7, -2.8531e-11],
}
SOURCES = {
    "gold": dict(
        table="gold-yokoo-2009-table3-isochores.csv",
        v0=67.72,
        gamma0=2.96,
        a=0.45,
        b=4.2,
        theta0=170.0,
        cold_ratio=0.9901781199449332,
        density=19.32,
    ),
    "platinum": dict(
        table="platinum-yokoo-2009-table5-isochores.csv",
        v0=60.55,
        gamma0=2.63,
        a=0.39,
        b=5.2,
        theta0=230.0,
        cold_ratio=0.9938200084195474,
        density=21.40,
    ),
}


def phonon(ratio, temperature, source, *, parameters=None, reference_ratio=1.0):
    """SI quadrature of Eqs. 4-12; volume denominator is always physical V."""
    r, t = np.broadcast_arrays(np.asarray(ratio, float), np.asarray(temperature, float))
    gamma0, a, b, theta0 = (
        (source[k] for k in ["gamma0", "a", "b", "theta0"])
        if parameters is None
        else parameters
    )
    normalized = r / reference_ratio
    gamma = gamma0 * (1 + a * (normalized**b - 1))
    theta = theta0 * normalized ** (-(1 - a) * gamma0) * np.exp(-(gamma - gamma0) / b)
    result = np.zeros(r.shape)
    positive = t > 0
    x = theta[positive] / t[positive]
    nodes, weights = quadrature(64)
    z = x[:, None] * nodes
    integral = x * np.sum(weights * z**3 / np.expm1(z), axis=1)
    energy = 9 * R * t[positive] * integral / x**3
    result[positive] = (
        gamma[positive] * energy / (source["v0"] * r[positive] * N_A * 1e-30 / 4) / 1e9
    )
    return result


def metrics(values):
    return {
        "states": len(values),
        "rmse_gpa": float(np.sqrt(np.mean(values**2))),
        "max_abs_difference_gpa": float(np.max(abs(values))),
    }


def audit():
    with (DATA / ELECTRONIC).open() as stream:
        electronic = list(csv.DictReader(stream))
    ts = np.array([float(x["temperature_k"]) for x in electronic])
    metals = {}
    for metal, s in SOURCES.items():
        with (DATA / s["table"]).open() as stream:
            rows = list(csv.DictReader(stream))
        r, t, p = [
            np.array([float(x[k]) for x in rows])
            for k in ["volume_ratio", "temperature_k", "pressure_gpa"]
        ]
        cold = {
            float(x["volume_ratio"]): float(x["pressure_gpa"])
            for x in rows
            if float(x["temperature_k"]) == 0
        }
        delta = p - np.array([cold[v] for v in r])
        column = (
            "gold_corrected_electronic_pressure_gpa"
            if metal == "gold"
            else "platinum_electronic_pressure_gpa"
        )
        ps = np.array([float(x[column]) for x in electronic])
        # Every source temperature is a recovered electronic table node. There
        # is no between-node interpolation choice in this source comparison.
        assert all(float(temp) in ts for temp in t)
        pel = np.interp(t, ts, ps)
        _, linear, quadratic, cubic = PYTHEOS_POLYNOMIALS[metal]
        polynomial_increment = linear * t + quadratic * t**2 + cubic * t**3
        ph = phonon(r, t, s)
        residual = ph + pel - delta
        warm = (t > 0) & np.array(
            [x["source_phase_annotation"] == "unmarked" for x in rows]
        )
        parameters = np.array([s[k] for k in ["gamma0", "a", "b", "theta0"]])
        half_steps = np.array([0.005, 0.005, 0.05, 0.5])
        lower, upper = parameters - half_steps, parameters + half_steps
        fit = least_squares(
            lambda x: (phonon(r, t, s, parameters=x) + pel - delta)[warm],
            parameters,
            bounds=(lower, upper),
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        assert fit.success
        intervals = []
        for i in np.flatnonzero(warm & (r == 1)):
            # At r=1, a/b disappear. Debye energy decreases with theta.
            lo_parameters = [lower[0], s["a"], s["b"], upper[3]]
            hi_parameters = [upper[0], s["a"], s["b"], lower[3]]
            # Include last-digit rounding of electronic pressure and density.
            lo = (
                float(phonon(1, t[i], s, parameters=lo_parameters))
                * (s["density"] - 0.005)
                / s["density"]
                + pel[i]
                - 0.005
            )
            hi = (
                float(phonon(1, t[i], s, parameters=hi_parameters))
                * (s["density"] + 0.005)
                / s["density"]
                + pel[i]
                + 0.005
            )
            source_lo, source_hi = delta[i] - 0.01, delta[i] + 0.01
            intervals.append(
                dict(
                    temperature_k=float(t[i]),
                    source_increment_gpa=float(delta[i]),
                    source_rounding_interval_gpa=[float(source_lo), float(source_hi)],
                    printed_model_rounding_interval_gpa=[lo, hi],
                    interval_gap_gpa=float(max(0, source_lo - hi, lo - source_hi)),
                    required_electronic_pressure_if_printed_phonon_exact_gpa=float(
                        delta[i] - ph[i]
                    ),
                    source_electronic_pressure_gpa=float(pel[i]),
                )
            )
        adaptive = []
        for ratio in [0.6, 0.8, 1.0]:
            gamma = s["gamma0"] * (1 + s["a"] * (ratio ** s["b"] - 1))
            theta = (
                s["theta0"]
                * ratio ** (-(1 - s["a"]) * s["gamma0"])
                * np.exp(-(gamma - s["gamma0"]) / s["b"])
            )
            for temp in [300, 1000, 3000]:
                x = theta / temp
                integral = quad(lambda z: z**3 / np.expm1(z), 0, x, epsabs=1e-11)[0]
                pressure = (
                    gamma
                    * 9
                    * R
                    * temp
                    * integral
                    / x**3
                    / (s["v0"] * ratio * N_A * 1e-30 / 4)
                    / 1e9
                )
                adaptive.append(abs(float(phonon(ratio, temp, s)) - pressure))
        metals[metal] = {
            "published_phonon_parameters": dict(
                zip(["gamma0", "a", "b", "theta0_k"], parameters.tolist())
            ),
            "ambient_volume_cell_a3": s["v0"],
            "printed_increment_unmarked_warm": metrics(residual[warm]),
            "best_within_assumed_last_digit_rounding": {
                "parameters": fit.x.tolist(),
                "lower_bounds": lower.tolist(),
                "upper_bounds": upper.tolist(),
                "residuals": metrics(
                    (phonon(r, t, s, parameters=fit.x) + pel - delta)[warm]
                ),
            },
            "ambient_interval_checks": intervals,
            "cold_normalized_phonon_alternative": metrics(
                (phonon(r, t, s, reference_ratio=s["cold_ratio"]) + pel - delta)[warm]
            ),
            "pytheos_electronic_polynomial_increment_alternative": {
                "coefficients_gpa_temperature_powers_0_to_3": PYTHEOS_POLYNOMIALS[
                    metal
                ],
                "residuals": metrics((ph + polynomial_increment - delta)[warm]),
                "qualification": "Electronic polynomial only; this is not a replay of Pytheos' separate 300 K BM3 reference curve. Constant polynomial and 300 K reference offsets cancel when forming increments from 0 K. These coefficients are an independent software convention, not recovered Yokoo author inputs.",
            },
            "adaptive_quad_max_difference_gpa": max(adaptive),
            "cold_curve_cancels_in_increment": True,
            "volume_fixed_zero_point_energy_cancels_in_increment": True,
            "constant_electronic_reference_offset_cancels_in_increment": True,
            "all_source_temperatures_are_electronic_nodes": True,
            "missing_source_cells_added": 0,
            "states": [
                {
                    **row,
                    "source_increment_from_0k_gpa": float(dp),
                    "printed_phonon_pressure_gpa": float(pph),
                    "source_electronic_pressure_gpa": float(pe),
                    "printed_increment_minus_source_gpa": float(err),
                }
                for row, dp, pph, pe, err in zip(rows, delta, ph, pel, residual)
            ],
        }
    return {
        "scope": "Published-equation thermal-increment audit; no fitted empirical pressure correction",
        "primary_source_pdf_sha256": {
            "10.1103/PhysRevB.80.104114": "aff3f39fee51031532f8251a3f8d5af5e08bcd58856d9fd09b22bcc68bef9907",
            "10.1103/PhysRevB.66.094115": "2f35128ef0a86d500896d7930a3b13092b1d33085d64fa628a3e3271088211e8",
        },
        "published_analytical_pvt_reproduced": False,
        "conclusion": "The supplied printed phonon parameters plus recovered electronic-pressure nodes do not reproduce the source thermal increments within the stated last-digit rounding test. Cold-volume choice, volume-fixed zero-point terms and constant pressure-reference offsets cancel and cannot account for this discrepancy. The adopted source phonon/electronic pressure implementation or unrounded coefficients are needed to resolve it; the audit does not identify a unique physical cause.",
        "rounding_policy": "Conditional half-last-digit intervals, not parameter uncertainties: gamma0/a +/-0.005, b +/-0.05, theta0 +/-0.5 K, electronic pressure +/-0.005 GPa, density +/-0.005 Mg/m3; source warm-minus-cold pressure +/-0.01 GPa. Source a/b uncertainties are not substituted for numerical rounding.",
        "cold_phonon_alternative_policy": "Diagnostic substitution of inferred cold-volume ratios, contradicting the ambient normalization stated for Eq. 10; it is not an accepted source convention.",
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [ELECTRONIC, *[s["table"] for s in SOURCES.values()]]
        },
        "metals": metals,
        "independent_software_source": {
            "repository": "https://github.com/SHDShim/pytheos",
            "commit": PYTHEOS_COMMIT,
            "source_files": [
                "pytheos/scales/gold.py",
                "pytheos/scales/platinum.py",
                "pytheos/eqn_electronic.py",
            ],
            "qualification": "Independent implementation evidence only; not author-supplied code or confirmation. Au's class explicitly offers a changed room-temperature K0-prime to improve table agreement; it does not recover the original cold-plus-phonon implementation.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = audit()
    if args.check:
        check_reconstruction(json.loads(OUTPUT.read_text()), report)
    else:
        OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                m: d["printed_increment_unmarked_warm"]
                for m, d in report["metals"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
