"""Independent, diagnostic replay of Litasov et al. (2013) Fe3C.

Repository maintenance tool; no catalog registration or source-data mutations.
Direct cell-volume equations, SciPy quadrature/least squares, and independently
printed KE checkpoints. See docs/literature-reproductions/litasov-2013-fe3c.md.
"""

from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import Avogadro, Boltzmann, R
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/cementite-litasov-2013-table1.csv"
TOKENS = DATA.with_name("cementite-litasov-2013-table1-tokens.json")
PDF = ROOT / "output/cementite-eos-audit-2026-10-25/sources/litasov-2013.pdf"
OUTPUT = ROOT / "docs/data/litasov-2013-fe3c-reproduction.json"
TR = 298.15
NODES, WEIGHTS = leggauss(64)
NODES, WEIGHTS = (NODES + 1) / 2, WEIGHTS / 2
# All values are independent literal transcriptions, not loaded from a catalog.
CASES = {
    "fm_table2": [155.2, 175.0, 5.0],
    "fm_text": [155.27, 175.0, 5.0],
    "ht_table2": [154.56, 190.0, 4.8, 3.90e-5, 1.22e-8, -0.029],
    "ht_text": [154.5, 190.0, 4.8, 3.90e-5, 1.22e-8, -0.029],
    "mgd_body": [154.42, 194.0, 4.6, 2.15, -0.03, 314.0],
    "mgd_abstract": [154.6, 192.0, 4.5, 2.09, -0.1, 490.0],
    # The fixed-500 K paragraph does not reprint the cold coefficients.
    # Reuse of the body cold curve is an explicit diagnostic assumption.
    "mgd_fixed500_conditional": [154.42, 194.0, 4.6, 2.13, -0.36, 500.0],
    "ke_text": [154.56, 190.8, 4.68, 2.10, 0.3, 305.0, 9.2e-5, 0.66],
    "ke_table2": [154.56, 191.0, 4.68, 2.10, 0.3, 305.0, 9.2e-5, 0.66],
}
ERRORS = {
    "fm_text": [0.05, 2.0, None],
    "fm_table2": [None, None, None],
    "ht_table2": [None, 2, 0.1, 0.15e-5, 0.17e-8, 0.002],
    "ht_text": [0.1, 3, 0.2, 0.15e-5, 0.17e-8, 0.002],
    "mgd_body": [None, 1, 0.1, 0.03, 0.27, 120],
    "mgd_abstract": [0.1, 3, 0.1, 0.04, 0.3, 120],
    "mgd_fixed500_conditional": [None, 1, 0.1, 0.02, 0.23, None],
    "ke_text": [None] * 8,
    "ke_table2": [None, 2, 0.08, None, None, None, None, None],
}
NAMES = {
    "fm": ["V0", "K0", "K0_prime"],
    "ht": ["V0", "K0", "K0_prime", "alpha0", "alpha1", "dK_dT"],
    "mgd": ["V0", "K0", "K0_prime", "gamma0", "q", "theta0"],
    "ke": ["V0", "K0", "K0_prime", "gamma0", "beta", "thetaE", "e0", "g"],
}
# Table 2 comparison rows are Litasov's reported/recalibrated comparisons,
# not authorized replacements for records of the original publications.
COMPARISON_ROWS = [
    (
        "li_2002",
        "unspecified",
        "Vinet",
        298,
        [155.3, 168, 5.8],
        [None, 3, 0.3],
        "0-31 GPa",
    ),
    (
        "ono_mibe_2010_experiment",
        "FM",
        "Vinet",
        298,
        [155.4, 166, 6.7],
        [None, 6, 0.3],
        "14-35 GPa",
    ),
    (
        "prescher_2012",
        "PM",
        "BM_order_unspecified",
        298,
        [154.2, 145, 8.5],
        [None, 3, 0.7],
        "8-22 GPa",
    ),
    ("sata_2010", "NM", "Vinet", 298, [148.5, 316, 3.5], [None, 6, 0.1], "50-184 GPa"),
    (
        "vocadlo_2002_fm",
        "FM",
        "BM_order_unspecified",
        0,
        [153.3, 173, 5.8],
        [None] * 3,
        "not printed",
    ),
    (
        "huang_2005_fm",
        "FM",
        "not_printed",
        0,
        [152.0, 212, 4.5],
        [None] * 3,
        "not printed",
    ),
    (
        "ono_mibe_2010_calculation_fm",
        "FM",
        "BM_order_unspecified",
        0,
        [151.9, 216, 4.15],
        [None] * 3,
        "0-35 GPa",
    ),
    (
        "mookherjee_2011_fm",
        "FM",
        "Vinet",
        0,
        [151.6, 183, 6.0],
        [None] * 3,
        "not printed",
    ),
    (
        "vocadlo_2002_nm",
        "NM",
        "BM_order_unspecified",
        0,
        [148.9, 317, 4.3],
        [None] * 3,
        "not printed",
    ),
    (
        "huang_2005_nm",
        "NM",
        "not_printed",
        0,
        [149.5, 322, 3.7],
        [None] * 3,
        "not printed",
    ),
    (
        "ono_mibe_2010_calculation_nm",
        "NM",
        "BM_order_unspecified",
        0,
        [142.5, 316, 4.4],
        [None] * 3,
        "<400 GPa",
    ),
    (
        "mookherjee_2011_nm",
        "NM",
        "Vinet",
        0,
        [143.2, 297, 4.9],
        [None] * 3,
        "not printed",
    ),
]


def comparison_inventory():
    out = []
    for name, phase, family, temp, parameters, errors, domain in COMPARISON_ROWS:
        item = {
            "id": name,
            "phase": phase,
            "family": family,
            "T_K": temp,
            "V0_K0_Kprime": parameters,
            "printed_error_widths": errors,
            "printed_pressure_domain": domain,
            "workflow_reproduced": False,
            "catalog_integration": False,
        }
        if family == "Vinet":
            item["equation_evaluated_at_V_over_V0"] = [0.95, 0.85]
            item["pressure_gpa"] = vinet(
                np.array([0.95, 0.85]) * parameters[0], *parameters
            ).tolist()
            item["verdict"] = (
                "Literal Vinet calculation evaluated; no independent benchmark or recalibrated fit inputs supplied by this paper."
            )
        else:
            item["verdict"] = (
                "Exact equation cannot be recovered from this comparison row: BM order or entire EOS family not printed. Original-publication audit required."
            )
        out.append(item)
    return out


def rows():
    with DATA.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_audit():
    """Freshly extract both Table 1 pages and compare every token/CSV value."""
    text = subprocess.check_output(
        ["pdftotext", "-f", "3", "-l", "4", "-layout", str(PDF), "-"],
        text=True,
    )
    recovered = []
    run = None
    for line in text.splitlines():
        match = re.fullmatch(r"\s*Run #(\d)\s*", line)
        if match:
            run = int(match[1])
        if "Run #7 includes" in line:
            break
        tokens = line.split()
        if (
            run is not None
            and len(tokens) in (6, 7)
            and all(
                re.fullmatch(r"\d+(?:\.\d+)?(?:\(\d+\))?", token) for token in tokens
            )
        ):
            if len(tokens) == 6:
                tokens.insert(1, None)
            recovered.append({"row": len(recovered) + 1, "run": run, "tokens": tokens})
    preserved = json.loads(TOKENS.read_text())
    assert recovered == preserved, "PDF and preserved Table 1 tokens differ"

    def number(token):
        if not token:
            return None, None
        match = re.fullmatch(r"([\d.]+)(?:\((\d+)\))?", token)
        central, error = match.groups()
        decimals = len(central.partition(".")[2])
        return float(central), None if error is None else int(error) * 10.0**-decimals

    for raw, row in zip(recovered, rows()):
        assert int(row["source_row"]) == raw["row"]
        assert int(row["run"]) == raw["run"]
        for token, key in zip(
            raw["tokens"],
            [
                "temperature_k",
                "mgo_volume_a3",
                "pressure_gpa",
                "a_angstrom",
                "b_angstrom",
                "c_angstrom",
                "volume_a3",
            ],
        ):
            value, error = number(token)
            assert (None if row[key] == "" else float(row[key])) == value
            sigma_key = key.replace("_a3", "_sigma_a3").replace("_gpa", "_sigma_gpa")
            if sigma_key in row and sigma_key != key:
                actual = None if row[sigma_key] == "" else float(row[sigma_key])
                assert actual == error or (
                    actual is not None
                    and error is not None
                    and abs(actual - error) < 1e-12
                )
    data = rows()
    return {
        "pdf_sha256": sha256(PDF),
        "csv_sha256": sha256(DATA),
        "tokens_sha256": sha256(TOKENS),
        "all_pdf_tokens_and_csv_values_match": True,
        "observations": len(data),
        "run_counts": dict(Counter(r["run"] for r in data)),
        "missing_calibrant_rows": [
            int(r["source_row"]) for r in data if not r["mgo_volume_a3"]
        ],
        "range": {
            key: [min(float(r[key]) for r in data), max(float(r[key]) for r in data)]
            for key in ["pressure_gpa", "temperature_k", "volume_a3"]
        },
        "volume_molar_cm3_per_formula": CASES["ke_text"][0] * Avogadro * 1e-24 / 4,
        "max_abs_abc_minus_printed_volume_a3": max(
            abs(
                float(r["a_angstrom"]) * float(r["b_angstrom"]) * float(r["c_angstrom"])
                - float(r["volume_a3"])
            )
            for r in data
        ),
    }


def vinet(v, v0, k0, kp):
    y = np.cbrt(np.asarray(v) / v0)
    return 3 * k0 * (1 - y) / y**2 * np.exp(1.5 * (kp - 1) * (1 - y))


def debye_energy(theta, t, nodes=NODES, weights=WEIGHTS):
    """Physical 16-atom conventional cell; J, excluding zero point energy."""
    y = np.asarray(theta) / np.asarray(t)
    z = y[..., None] * nodes
    return (
        9
        * 16
        * Boltzmann
        * np.asarray(t)
        * y
        * np.sum(weights * nodes**3 / np.expm1(z), axis=-1)
    )


def pressure(case, v, t, parameters=None, tr=TR, literal_debye_sign=False):
    p = np.asarray(CASES[case] if parameters is None else parameters)
    v, t = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(t, dtype=float))
    v0, k0, kp = p[:3]
    if case.startswith("fm"):
        return vinet(v, v0, k0, kp)
    if case.startswith("ht"):
        a0, a1, dk = p[3:]
        vt = v0 * np.exp(a0 * (t - tr) + a1 / 2 * (t**2 - tr**2))
        return vinet(v, vt, k0 + dk * (t - tr), kp)
    if case.startswith("mgd"):
        g0, q, theta0 = p[3:]
        x = v / v0
        gamma = g0 * x**q
        exponent = (
            -g0 * np.log(x) if abs(q) < 1e-9 else -g0 * np.expm1(q * np.log(x)) / q
        )
        theta = theta0 * np.exp(exponent)
        if literal_debye_sign:
            theta = theta0 * np.exp(-exponent)
        return vinet(v, v0, k0, kp) + gamma * (
            debye_energy(theta, t) - debye_energy(theta, tr)
        ) / (v * 1e-21)
    g0, beta, theta0, e0, g = p[3:]
    x = v / v0
    gamma = g0 * x**beta  # gamma_infinity = 0 fixed by source.
    theta = theta0 * np.exp((g0 - gamma) / beta)
    thermal_cell = 48 * Boltzmann * gamma * theta * (
        1 / np.expm1(theta / t) - 1 / np.expm1(theta / tr)
    ) + 24 * Boltzmann * e0 * x**g * g * (t**2 - tr**2)
    return vinet(v, v0, k0, kp) + thermal_cell / (v * 1e-21)


def ke_state(p, t, case="ke_text", parameters=None):
    """Thermodynamics from Helmholtz energy and numerical derivatives."""
    par = np.asarray(CASES[case] if parameters is None else parameters)
    v0, k0, kp, g0, beta, theta0, e0, g = par
    v = brentq(
        lambda vv: pressure(case, vv, t, par) - p, v0 * 0.45, v0 * 1.5, xtol=1e-11
    )
    x = v / v0
    vm = v * Avogadro * 1e-30 / 4  # m3 / mol Fe3C
    theta = theta0 * np.exp(g0 * (1 - x**beta) / beta)
    y = theta / t
    cv = 12 * R * y**2 * np.exp(-y) / (-np.expm1(-y)) ** 2 + 12 * R * e0 * x**g * t
    entropy = (
        12 * R * (y / np.expm1(y) - np.log(-np.expm1(-y))) + 12 * R * e0 * x**g * t
    )
    kt = (
        -v
        * (pressure(case, v + 0.001, t, par) - pressure(case, v - 0.001, t, par))
        / 0.002
    )
    dpdt = (pressure(case, v, t + 0.01, par) - pressure(case, v, t - 0.01, par)) / 0.02
    alpha = dpdt / kt
    cp = cv + alpha**2 * vm * kt * 1e9 * t
    ks = kt * cp / cv

    def free_energy(temp):
        return (
            12 * R * temp * np.log(-np.expm1(-theta / temp))
            - 6 * R * e0 * x**g * temp**2
        )

    cold_energy = -quad(lambda vv: vinet(vv, v0, k0, kp), v0, v, epsabs=1e-9)[0]
    cold_energy *= Avogadro * 1e-21 / 4
    gibbs = cold_energy + free_energy(t) - free_energy(TR) + p * 1e9 * vm
    # Table 3's gamma follows the oscillator law, not alpha*V*Kt/Cv in
    # the presence of an electronic heat capacity with its own exponent g.
    return np.array(
        [x, alpha * 1e5, entropy, cp, cv, kt, ks, g0 * x**beta, gibbs / 1000]
    )


# Table 3, p. 5280. Calculated outputs, NEVER measured fitting observations.
TABLE3 = np.array(
    [
        [
            0.0001,
            298.15,
            1.0000,
            4.36,
            104.48,
            96.76,
            94.25,
            190.8,
            195.9,
            2.096,
            0.002,
        ],
        [
            0.0001,
            500,
            1.0093,
            4.80,
            157.13,
            106.42,
            101.47,
            183.5,
            192.4,
            2.102,
            -26.775,
        ],
        [
            0.0001,
            1000,
            1.0355,
            5.44,
            235.32,
            120.26,
            108.50,
            165.0,
            182.8,
            2.118,
            -126.736,
        ],
        [
            0.0001,
            1500,
            1.0659,
            6.14,
            286.67,
            134.32,
            113.87,
            145.9,
            172.1,
            2.137,
            -257.861,
        ],
        [
            0.0001,
            2000,
            1.1015,
            7.05,
            327.55,
            151.36,
            119.21,
            126.1,
            160.1,
            2.158,
            -411.702,
        ],
        [10, 298.15, 0.9541, 3.58, 95.53, 94.49, 92.49, 236.0, 241.1, 2.067, 227.129],
        [10, 500, 0.9614, 3.95, 147.19, 104.67, 100.68, 228.9, 238.0, 2.072, 202.272],
        [10, 1000, 0.9818, 4.40, 223.90, 117.34, 108.01, 211.4, 229.7, 2.085, 107.662],
        [10, 1500, 1.0046, 4.82, 273.65, 128.99, 113.25, 193.5, 220.4, 2.099, -17.367],
        [10, 2000, 1.0303, 5.31, 312.50, 142.02, 118.33, 175.2, 210.2, 2.115, -164.212],
        [20, 298.15, 0.9177, 3.06, 88.34, 92.48, 90.82, 278.6, 283.6, 2.043, 444.765],
        [20, 500, 0.9237, 3.40, 139.17, 103.31, 99.94, 271.7, 280.9, 2.047, 421.456],
        [20, 1000, 0.9405, 3.75, 214.84, 115.44, 107.60, 254.8, 273.4, 2.058, 331.128],
        [20, 1500, 0.9590, 4.05, 263.59, 125.78, 112.76, 237.6, 265.0, 2.070, 210.873],
        [20, 2000, 0.9793, 4.37, 301.27, 136.86, 117.67, 220.1, 256.0, 2.083, 69.339],
        [30, 298.15, 0.8874, 2.68, 82.33, 90.63, 89.21, 319.3, 324.4, 2.023, 654.690],
        [30, 500, 0.8926, 3.01, 132.42, 102.17, 99.23, 312.6, 321.9, 2.026, 632.681],
        [30, 1000, 0.9069, 3.30, 207.29, 114.05, 107.23, 296.1, 314.9, 2.036, 545.944],
        [30, 1500, 0.9225, 3.53, 255.33, 123.57, 112.35, 279.3, 307.2, 2.046, 429.640],
        [30, 2000, 0.9395, 3.77, 292.22, 133.47, 117.14, 262.4, 298.9, 2.057, 292.429],
        [40, 298.15, 0.8616, 2.40, 77.17, 88.90, 87.66, 358.7, 363.7, 2.005, 858.122],
        [40, 500, 0.8661, 2.71, 126.59, 101.15, 98.55, 352.1, 361.4, 2.008, 837.234],
        [40, 1000, 0.8786, 2.97, 200.78, 112.96, 106.89, 335.8, 354.9, 2.016, 753.598],
        [40, 1500, 0.8922, 3.16, 248.29, 121.92, 112.00, 319.4, 347.7, 2.026, 640.680],
        [40, 2000, 0.9068, 3.35, 284.60, 131.01, 116.69, 302.8, 340.0, 2.036, 507.133],
    ]
)
TABLE4_T = [298.15, 500, 750, 1000, 1250, 1500, 1750, 2000]
TABLE4_X = np.array(
    [
        [1, 1.0093, 1.0220, 1.0355, 1.0501, 1.0659, 1.0829, 1.1015],
        [0.9756, 0.9838, 0.9948, 1.0066, 1.0192, 1.0326, 1.0470, 1.0624],
        [0.9541, 0.9614, 0.9713, 0.9818, 0.9929, 1.0046, 1.0171, 1.0303],
        [0.9350, 0.9416, 0.9505, 0.9599, 0.9699, 0.9804, 0.9914, 1.0031],
        [0.9177, 0.9237, 0.9319, 0.9405, 0.9495, 0.9590, 0.9689, 0.9793],
        [0.9019, 0.9075, 0.9150, 0.9229, 0.9312, 0.9398, 0.9489, 0.9583],
        [0.8874, 0.8926, 0.8996, 0.9069, 0.9145, 0.9225, 0.9308, 0.9395],
        [0.8741, 0.8789, 0.8854, 0.8922, 0.8993, 0.9067, 0.9144, 0.9224],
        [0.8616, 0.8661, 0.8722, 0.8786, 0.8853, 0.8922, 0.8993, 0.9068],
    ]
)


def checkpoints():
    names = [
        "x",
        "alpha_times_1e5",
        "S_J_mol_K",
        "Cp_J_mol_K",
        "Cv_J_mol_K",
        "Kt_GPa",
        "Ks_GPa",
        "gamma_oscillator",
        "deltaG_kJ_mol",
    ]
    result = {}
    for case in ("ke_text", "ke_table2"):
        calculated = np.array([ke_state(p, t, case) for p, t in TABLE3[:, :2]])
        difference = calculated - TABLE3[:, 2:]
        result[case] = {
            "table3_source": TABLE3.tolist(),
            "table3_calculated": calculated.tolist(),
            "max_absolute_difference": dict(
                zip(names, np.max(abs(difference), axis=0).tolist())
            ),
        }
        # Linear worst-direction sensitivity to half of the printed last digit.
        # Fixed beta=0.3 is held exact, as are a0=0 and gamma_inf=0.
        widths = [
            0.005,
            0.05 if case == "ke_text" else 0.5,
            0.005,
            0.005,
            0,
            0.5,
            0.05e-5,
            0.005,
        ]
        envelope = np.zeros_like(calculated)
        for j, width in enumerate(widths):
            if width:
                plus, minus = np.array(CASES[case]), np.array(CASES[case])
                plus[j] += width
                minus[j] -= width
                envelope += (
                    np.abs(
                        np.array(
                            [
                                ke_state(p, t, case, plus) - ke_state(p, t, case, minus)
                                for p, t in TABLE3[:, :2]
                            ]
                        )
                    )
                    / 2
                )
        output_halfdigits = np.array(
            [0.00005, 0.005, 0.005, 0.005, 0.005, 0.05, 0.05, 0.0005, 0.0005]
        )
        result[case]["table3_within_linear_printing_envelope"] = bool(
            np.all(abs(difference) <= envelope + output_halfdigits)
        )
        result[case]["table3_printing_envelope_max_fraction_used"] = float(
            np.max(abs(difference) / (envelope + output_halfdigits))
        )
        result[case]["input_half_last_digit_widths"] = widths
    xx = np.array(
        [
            [ke_state(0.0001 if i == 0 else i * 5, t)[0] for t in TABLE4_T]
            for i in range(9)
        ]
    )
    result["table4"] = {
        "source": TABLE4_X.tolist(),
        "calculated": xx.tolist(),
        "max_abs_x_difference": float(np.max(abs(xx - TABLE4_X))),
    }
    return result


def selection(name):
    data = rows()
    if name == "fm":
        return [
            r
            for r in data
            if float(r["temperature_k"]) == 300 and float(r["pressure_gpa"]) < 8
        ]
    # No invented Curie line. Two bounded alternatives bracket 473 K states.
    cutoff = 523 if name == "pm_conservative" else 473
    return [
        r
        for r in data
        if float(r["temperature_k"]) >= cutoff or float(r["pressure_gpa"]) >= 9.3
    ]


def metrics(case, data, parameters=None, tr=TR, literal_debye_sign=False):
    v, t, p = (
        np.array([float(r[key]) for r in data])
        for key in ("volume_a3", "temperature_k", "pressure_gpa")
    )
    residual = pressure(case, v, t, parameters, tr, literal_debye_sign) - p
    return {
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "mean_gpa": float(np.mean(residual)),
        "max_abs_gpa": float(max(abs(residual))),
        "row_ids": [int(r["source_row"]) for r in data],
        "residuals_gpa": residual.tolist(),
    }


def fit(
    case,
    mask,
    weighted=False,
    fixed_cold=False,
    tr=TR,
    fixed_theta=None,
    literal_debye_sign=False,
):
    data = selection(mask)
    v, t, obs = (
        np.array([float(r[key]) for r in data])
        for key in ("volume_a3", "temperature_k", "pressure_gpa")
    )
    source = np.array(CASES[case])
    initial = source.copy()
    if fixed_theta is not None:
        initial[-1] = fixed_theta
    free = np.ones(len(initial), dtype=bool)
    if case.startswith("fm"):
        free[2] = False
    if "fixed500" in case or fixed_theta is not None:
        free[-1] = False
    if fixed_cold:
        free[:3] = False
    lower, upper = [150, 100, 1], [160, 300, 10]
    if case.startswith("ht"):
        lower += [0, -1e-7, -0.1]
        upper += [1e-4, 1e-7, 0]
    elif case.startswith("mgd"):
        lower += [0.1, -5, 1]
        upper += [5, 5, 2000]
    # KE optimization requires external thermochemical inputs. Only the cold
    # coefficients are refit here with all thermal terms held to source values.
    elif case.startswith("ke"):
        free[3:] = False
        lower += [0.1, 0.01, 1, 0, 0]
        upper += [5, 3, 2000, 0.001, 5]
    scales = {
        "fm": [1, 10, 1],
        "ht": [1, 10, 1, 1e-5, 1e-8, 0.01],
        "mgd": [1, 10, 1, 1, 1, 300],
        "ke": [1, 10, 1, 1, 1, 300, 1e-4, 1],
    }
    scale = np.array(scales[case.split("_")[0]])[free]
    sigma = np.ones(len(data))
    if weighted:
        dv = (
            pressure(case, v + 0.001, t, initial, tr, literal_debye_sign)
            - pressure(case, v - 0.001, t, initial, tr, literal_debye_sign)
        ) / 0.002
        sp = np.array(
            [
                float(r["pressure_sigma_gpa"]) if r["pressure_sigma_gpa"] else 0
                for r in data
            ]
        )
        sv = np.array([float(r["volume_sigma_a3"]) for r in data])
        sigma = np.sqrt(sp**2 + (dv * sv) ** 2)

    def unpack(x):
        par = initial.copy()
        par[free] = x * scale
        return par

    def objective(x):
        return (pressure(case, v, t, unpack(x), tr, literal_debye_sign) - obs) / sigma

    starts = [initial.copy()]
    if case.startswith("mgd") and free[-1]:
        for theta in [500, 800]:
            candidate = initial.copy()
            candidate[-1] = theta
            starts.append(candidate)
    solutions = [
        least_squares(
            objective,
            start[free] / scale,
            bounds=(np.array(lower)[free] / scale, np.array(upper)[free] / scale),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-10,
            max_nfev=3000,
        )
        for start in starts
    ]
    sol = min(solutions, key=lambda s: float(s.fun @ s.fun))
    par = unpack(sol.x)
    dof = len(data) - int(sum(free))
    covariance_scaled = np.linalg.pinv(sol.jac.T @ sol.jac) * (sol.fun @ sol.fun) / dof
    covariance = covariance_scaled * np.outer(scale, scale)
    se = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(se, se)
    names = NAMES[case.split("_")[0]]
    comparison = {}
    for i, (key, a, b) in enumerate(zip(names, source, par)):
        width = ERRORS.get(case, [None] * len(initial))[i]
        comparison[key] = {
            "source": a,
            "refit": b,
            "delta": b - a,
            "source_error_width": width,
            "free": bool(free[i]),
            "within_source_error": None if width is None else bool(abs(b - a) <= width),
        }
    source_p = pressure(case, v, t, source, tr, literal_debye_sign)
    refit_p = pressure(case, v, t, par, tr, literal_debye_sign)
    return {
        "case": case,
        "literal_cited_2007_debye_sign": literal_debye_sign,
        "mask": mask,
        "n": len(data),
        "Tr": tr,
        "fixed_cold": fixed_cold,
        "profile_fixed_theta_K": fixed_theta,
        "weighting": "fixed_effective_sigma_at_source"
        if weighted
        else "equal_pressure",
        "parameters": dict(zip(names, par.tolist())),
        "comparison": comparison,
        "free_parameters": [name for name, flag in zip(names, free) if flag],
        "standard_errors_residual_scaled": se.tolist(),
        "covariance": covariance.tolist(),
        "correlation": correlation.tolist(),
        "degrees_of_freedom": dof,
        "singular_values_scaled_jacobian": np.linalg.svd(
            sol.jac, compute_uv=False
        ).tolist(),
        "source_metrics": metrics(
            case, data, tr=tr, literal_debye_sign=literal_debye_sign
        ),
        "refit_metrics": metrics(case, data, par, tr, literal_debye_sign),
        "max_source_refit_pressure_difference_gpa": float(max(abs(source_p - refit_p))),
        "converged": bool(sol.success),
        "active_bounds": sol.active_mask.tolist(),
        "multistart_parameters": [unpack(s.x).tolist() for s in solutions],
        "objective_sum_squares": float(sol.fun @ sol.fun),
        "author_workflow_reproduced": False,
    }


def calibration_replay():
    """Check BOTH original Sokolova 2013 MgO variants, without recalibration.

    Tables 1-3 and Eqs. 2, 4, 6, 9, 14. The Fe3C paper does not identify
    which cold curve variant or private optimized parameter precision was used.
    """
    output = {}
    vm0 = 11.248e-6
    cell0 = vm0 / Avogadro * 4e30
    for kind, kp, th1, delta, tt, a0, m in [
        ("holzapfel_table2", 4.23, 749, -0.226, 0.555, 14.6e-6, 5.3),
        ("vinet_table3", 4.38, 747, -0.165, 1.064, 16.5e-6, 4.96),
    ]:

        def cold(x):
            if kind.startswith("vinet"):
                return vinet(x, 1, 160.3, kp)
            y = x ** (1 / 3)
            c0 = -np.log(3 * 160.3 / (1003.6 * (2 * 10.34 / 11.248) ** (5 / 3)))
            c2 = 1.5 * (kp - 3) - c0
            return (
                3
                * 160.3
                * y**-5
                * (1 - y)
                * np.exp(c0 * (1 - y))
                * (1 + c2 * y * (1 - y))
            )

        def gamma(x):
            h = 1e-4
            value = cold(x)
            first = (cold(x + h) - cold(x - h)) / (2 * h)
            second = (cold(x + h) - 2 * value + cold(x - h)) / h**2
            modulus = -x * first
            derivative = -(first + x * second) / first
            return (derivative / 2 - 1 / 6 - tt / 3 * (1 - value / (2 * modulus))) / (
                1 - 2 * tt * value / (3 * modulus)
            ) + delta

        predicted, actual, ids = [], [], []
        for row in rows():
            if not row["mgo_volume_a3"]:
                continue
            x = float(row["mgo_volume_a3"]) / cell0
            temp = float(row["temperature_k"])
            multiplier = np.exp(quad(lambda u: gamma(u) / u, x, 1, epsabs=1e-7)[0])
            gam, anh = gamma(x), a0 * x**m

            def thermal(t):
                energy = sum(
                    3
                    * R
                    * theta
                    * multiplier
                    * np.exp(anh * t / 2)
                    / np.expm1(theta * multiplier * np.exp(anh * t / 2) / t)
                    for theta in [th1, 401]
                )
                return energy * (gam - m * anh * t / 2) / (vm0 * x) * 1e-9

            predicted.append(float(cold(x) + thermal(temp) - thermal(TR)))
            actual.append(float(row["pressure_gpa"]))
            ids.append(int(row["source_row"]))
        residual = np.array(predicted) - actual
        output[kind] = {
            "row_ids": ids,
            "predicted_gpa": predicted,
            "residual_gpa": residual.tolist(),
            "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "max_abs_gpa": float(max(abs(residual))),
            "cell0_a3": cell0,
        }
    return output


def axial():
    data = [
        r
        for r in rows()
        if float(r["temperature_k"]) == 300 and float(r["pressure_gpa"]) >= 9.3
    ]
    p = np.array([float(r["pressure_gpa"]) for r in data])
    out = {}
    for axis, zero, modulus in [
        ("a", 5.085, 195),
        ("b", 6.743, 186),
        ("c", 4.523, 163),
    ]:
        length = np.array([float(r[axis + "_angstrom"]) for r in data])

        # Conventional linearized EOS: insert (a/a0)^3, K is volume-equivalent.
        def calc(par):
            return vinet(length**3, par[0] ** 3, par[1], 4.8)

        sol = least_squares(
            lambda par: calc(par) - p,
            [zero, modulus],
            x_scale=[zero, modulus],
            ftol=1e-13,
            xtol=1e-13,
            gtol=1e-10,
        )
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * np.sum(sol.fun**2) / (len(p) - 2)
        residual = calc([zero, modulus]) - p
        out[axis] = {
            "source": {"length0": zero, "K0": modulus, "K0_prime_fixed": 4.8},
            "n": len(data),
            "row_ids": [int(r["source_row"]) for r in data],
            "refit": sol.x.tolist(),
            "covariance": cov.tolist(),
            "se": np.sqrt(np.diag(cov)).tolist(),
            "source_residual_gpa": residual.tolist(),
            "source_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "refit_rmse_gpa": float(np.sqrt(np.mean(sol.fun**2))),
        }
    return out


def library_check():
    """Optional installed Python/Rust cross-check; never used as fit engine."""
    from peritheos.eos.rt import Vinet
    from peritheos.eos.thermal import MieGruneisenDebye, ThermalReferenceStateEOS

    v = np.array([142.27, 150.49, 148.77, 161.66])
    t = np.array([1473.0, 873.0, 473.0, 1273.0])
    factor = Avogadro * 1e-25 / 4  # A3/cell -> J/bar/mol formula
    out = {}
    for case in ["mgd_body", "mgd_abstract", "mgd_fixed500_conditional"]:
        v0, k0, kp, g0, q, theta = CASES[case]
        model = MieGruneisenDebye(Vinet(v0 * factor, k0, kp), TR, theta, g0, q, 4)
        direct = pressure(case, v, t)
        public = np.asarray(model.pressure(v * factor, t))
        assert hasattr(model, "_native")
        del model._native
        python = np.asarray(model.pressure(v * factor, t))
        out[case] = {
            "max_abs_difference_gpa": float(max(abs(direct - public))),
            "direct_gpa": direct.tolist(),
            "public_gpa": public.tolist(),
            "max_abs_python_native_difference_gpa": float(max(abs(python - public))),
        }
    for case in ["ht_table2", "ht_text"]:
        v0, k0, kp, a0, a1, dk = CASES[case]
        model = ThermalReferenceStateEOS(
            Vinet(v0 * factor, k0, kp),
            TR,
            a0,
            dk,
            a1,
            thermal_expansion_law="linear_temperature",
        )
        direct = pressure(case, v, t)
        public = np.asarray(model.pressure(v * factor, t))
        out[case] = {"max_abs_difference_gpa": float(max(abs(direct - public)))}
    for case in ["fm_text", "fm_table2"]:
        v0, k0, kp = CASES[case]
        vv = np.array([155.25, 154.0, 153.43, 149.57])
        out[case] = {
            "max_abs_difference_gpa": float(
                max(abs(vinet(vv, v0, k0, kp) - Vinet(v0, k0, kp).pressure(vv)))
            )
        }
    return out


def nm_extrapolation():
    """Literal Sata Table 2 cold curve + each named PM thermal alternative.

    These are conditional interpretations of Figure 13's unspecified thermal
    transfer. No Table 1 refit or assertion of NM physical validity is made.
    """
    mass = 3 * 55.845 + 12.011
    out = {}
    for case in ["ke_text", "mgd_body", "ht_table2"]:
        par = np.array(CASES[case])
        par[:3] = [148.5, 316, 3.5]
        states = []
        for p in [240, 300, 330, 360]:
            for t in [300, 3000, 5000, 7000]:
                v = brentq(lambda vv: pressure(case, vv, t, par) - p, 50, 200)
                states.append(
                    {
                        "P_GPa": p,
                        "T_K": t,
                        "volume_a3": v,
                        "density_g_cm3": mass / (v * Avogadro * 1e-24 / 4),
                    }
                )
        out[case] = states
    return out


def graph_benchmarks():
    """Read original embedded figure rasters, preserving pixel calibration.

    Graphical checkpoints are source calculations, not experimental rows.
    Fixed columns avoid legends and visible marker centres/error bars.
    """
    from PIL import Image

    def raster(page, width):
        with tempfile.TemporaryDirectory() as folder:
            dest = Path(folder) / "source.svg"
            subprocess.run(
                [
                    "pdftocairo",
                    "-f",
                    str(page),
                    "-l",
                    str(page),
                    "-svg",
                    str(PDF),
                    str(dest),
                ],
                check=True,
            )
            for element in ET.parse(dest).getroot().iter():
                if element.tag.endswith("image") and int(element.get("width")) == width:
                    blob = base64.b64decode(
                        element.get("{http://www.w3.org/1999/xlink}href").split(",")[1]
                    )
                    return np.array(Image.open(io.BytesIO(blob))), hashlib.sha256(
                        blob
                    ).hexdigest()
        raise RuntimeError("Source figure raster not found")

    im, digest = raster(7, 884)
    points = []
    for column in [180, 200, 220, 280, 310, 350, 400, 430, 480, 700, 750, 800, 840]:
        strip = im[:650, column - 1 : column + 2, :3]
        yy, _ = np.where(
            (strip[:, :, 0] > 180) & (strip[:, :, 1] < 90) & (strip[:, :, 2] < 90)
        )
        yy = yy[(yy > 45) & (yy < 618)]
        ypos = float(np.median(yy))
        p = (column - 123) * 40 / (867 - 123)
        x = 0.85 + (617 - ypos) * 0.15 / (617 - 43)
        points.append({"pixel_x": column, "pixel_y": ypos, "P_GPa": p, "V_over_V0": x})
    comparisons = {}
    for case in ["ht_table2", "mgd_body", "mgd_abstract", "ke_text"]:
        residual = [
            brentq(lambda vv: pressure(case, vv, TR) - point["P_GPa"], 120, 165)
            / CASES[case][0]
            - point["V_over_V0"]
            for point in points
        ]
        comparisons[case] = {
            "max_abs_V_over_V0_difference": float(max(abs(np.array(residual)))),
            "residual": residual,
            "within_graphical_bound": bool(max(abs(np.array(residual))) < 0.0015),
        }
    fig8 = {
        "raster_sha256": digest,
        "axes_pixels": {"P0": 123, "P40": 867, "x1": 43, "x085": 617},
        "points": points,
        "graphical_absolute_bound_V_over_V0": 0.0015,
        "comparisons": comparisons,
    }
    im, digest = raster(10, 912)
    points = []
    for column in [192, 400, 445, 650, 750, 800]:
        strip = im[15:537, column - 1 : column + 2, :3].astype(int)
        yy, _ = np.where(
            (strip[:, :, 0] < 110)
            & (strip[:, :, 1] < 85)
            & (strip[:, :, 2] < 85)
            & (strip[:, :, 0] > strip[:, :, 1] + 4)
        )
        yy = np.unique(yy + 15)
        yy = yy[yy > 190]
        groups = np.split(yy, np.where(np.diff(yy) > 3)[0] + 1)
        assert len(groups) == 3
        p = 240 + (column - 129) * 120 / (888 - 129)
        for temp, group in zip([300, 3000, 5000], groups):
            ypos = float(np.median(group))
            density = 11 + (539 - ypos) * 3.5 / (539 - 11)
            points.append(
                {
                    "pixel_x": column,
                    "pixel_y": ypos,
                    "P_GPa": p,
                    "T_K": temp,
                    "density_g_cm3": density,
                }
            )
    comparisons = {}
    for case in ["ke_text", "mgd_body", "ht_table2"]:
        par = np.array(CASES[case])
        par[:3] = [148.5, 316, 3.5]
        residual = []
        for point in points:
            v = brentq(
                lambda vv: pressure(case, vv, point["T_K"], par) - point["P_GPa"],
                50,
                200,
            )
            rho = (3 * 55.845 + 12.011) / (v * Avogadro * 1e-24 / 4)
            residual.append(rho - point["density_g_cm3"])
        comparisons[case] = {
            "density_residual_g_cm3": residual,
            "max_abs_g_cm3": float(max(abs(np.array(residual)))),
        }
    return {
        "figure8": fig8,
        "figure13": {
            "raster_sha256": digest,
            "axes_pixels": {"P240": 129, "P360": 888, "rho11": 539, "rho145": 11},
            "points": points,
            "graphical_reading_bound_g_cm3": 0.025,
            "comparisons": comparisons,
            "qualification": "Curve reading bound only; input coefficient precision and thermal transfer remain separate.",
        },
    }


def plot_report(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    pm = selection("pm_473_included")
    pp = np.linspace(0, 32, 100)
    for temp, color in [
        (300, "#285c95"),
        (673, "#427a64"),
        (1073, "#a57820"),
        (1473, "#ab354b"),
    ]:
        vv = [brentq(lambda v: pressure("ke_text", v, temp) - p, 120, 170) for p in pp]
        axes[0].plot(pp, vv, color=color, label=f"{temp} K, KE")
        rr = [r for r in pm if float(r["temperature_k"]) == temp]
        axes[0].scatter(
            [float(r["pressure_gpa"]) for r in rr],
            [float(r["volume_a3"]) for r in rr],
            color=color,
            edgecolors="white",
            linewidths=0.4,
            s=26,
            zorder=3,
        )
    axes[0].set(
        xlabel="Pressure (GPa)",
        ylabel="Cell volume (Å³)",
        title="Published KE curves and measured states",
    )
    axes[0].legend(fontsize=7)
    obs = np.array([float(r["pressure_gpa"]) for r in pm])
    for case, color in [
        ("ht_table2", "#427a64"),
        ("mgd_body", "#285c95"),
        ("ke_text", "#ab354b"),
    ]:
        residual = report["published_residuals"][case]["pm_473_included"][
            "residuals_gpa"
        ]
        axes[1].scatter(
            obs, residual, s=15, alpha=0.6, label=case.replace("_", " "), color=color
        )
    axes[1].axhline(0, color="black", linewidth=0.6)
    axes[1].set(
        xlabel="Measured pressure (GPa)",
        ylabel="P(model) − P(table) (GPa)",
        title="Literal source coefficients; 68-state mask",
    )
    axes[1].legend(fontsize=7)
    profile = report["mgd_theta_profile"]
    axes[2].plot(
        [p["parameters"]["theta0"] for p in profile],
        [p["objective_sum_squares"] for p in profile],
        "o-",
        color="#285c95",
    )
    for value, label in [(314, "body"), (490, "abstract")]:
        axes[2].axvline(
            value, linestyle="--", linewidth=0.8, label=f"{label}: {value} K"
        )
    axes[2].set(
        xlabel="Fixed Debye temperature (K)",
        ylabel="Profile pressure SSE (GPa²)",
        title="MGD: weak Debye-temperature constraint",
    )
    axes[2].legend(fontsize=7)
    for axis in axes:
        axis.grid(alpha=0.16)
    path = ROOT / "docs/literature-reproductions/figures/litasov-2013-fe3c-replay.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    audit = source_audit()
    published = {
        case: {
            mask: metrics(case, selection(mask))
            for mask in (
                ["fm"]
                if case.startswith("fm")
                else ["pm_conservative", "pm_473_included"]
            )
        }
        for case in CASES
    }
    fits = []
    for case in [
        "fm_text",
        "ht_table2",
        "ht_text",
        "mgd_body",
        "mgd_abstract",
        "mgd_fixed500_conditional",
        "ke_text",
    ]:
        for mask in (
            ["fm"] if case.startswith("fm") else ["pm_conservative", "pm_473_included"]
        ):
            for weighted in [False, True]:
                fits.append(fit(case, mask, weighted))
    # Source-supported cold/thermal sequencing sensitivity, not a replacement law.
    for case in ["mgd_body", "mgd_abstract", "mgd_fixed500_conditional"]:
        fits.append(fit(case, "pm_473_included", fixed_cold=True))
    sensitivity = {
        case: {
            str(tr): metrics(case, selection("pm_473_included"), tr=tr)
            for tr in [298.0, 298.15, 300.0]
        }
        for case in ["ht_table2", "mgd_body", "ke_text"]
    }
    report = {
        "source": "Litasov et al. 2013, doi:10.1002/2013JB010270",
        "diagnostic_only": True,
        "source_audit": audit,
        "source_coefficients": CASES,
        "source_errors": ERRORS,
        "table2_comparison_inventory": comparison_inventory(),
        "published_residuals": published,
        "ke_checkpoints": checkpoints(),
        "fits": fits,
        "axial": axial(),
        "reference_temperature_sensitivity": sensitivity,
        "library_crosscheck": library_check(),
        "calibration_diagnostics": calibration_replay(),
        "mgd_theta_profile": [
            fit("mgd_body", "pm_473_included", fixed_theta=theta)
            for theta in [1, 100, 200, 277.33, 314, 490, 500, 800, 1000]
        ],
        "nm_conditional_extrapolations": nm_extrapolation(),
        "cited_2007_literal_sign_controls": [
            fit(case, "pm_473_included", literal_debye_sign=True)
            for case in ["mgd_body", "mgd_abstract", "mgd_fixed500_conditional"]
        ],
        "graphical_source_calculations": graph_benchmarks(),
        "ke_fixed_or_inactive_parameters": {
            "gamma_infinity": 0,
            "a0": 0,
            "m": 4.3,
            "mE1": 12,
            "n": 4,
            "Z": 4,
            "k_Kunc": 2,
        },
    }
    OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    plot_report(report)
    print(
        f"Wrote {OUTPUT.relative_to(ROOT)}: {len(fits)} fits, 25 thermodynamic + 72 volume checkpoints"
    )


if __name__ == "__main__":
    main()
