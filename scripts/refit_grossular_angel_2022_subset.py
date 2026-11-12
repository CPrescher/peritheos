"""The predeclared 95-channel diagnostic; never writes material records.

Repository-only: python -m scripts.refit_grossular_angel_2022_subset.
All assumptions and exclusions are in grossular/refit-plan.json.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import least_squares

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import (
    DATA,
    GAS_CONSTANT,
    observations,
)
from scripts.refit_pyrope_angel_2022_expanded import recalibrate

V0, TR, N = 125.35, 298.0, 20
NODES, WEIGHTS = np.polynomial.legendre.leggauss(64)
NAMES = [
    "K0_gpa",
    "K0_prime",
    "theta0_k",
    "gamma0",
    "scale_milani_vt",
    "scale_du_vt",
    "scale_kono_ks",
]
SCALING = np.array([170, 5, 750, 1.15, 1, 1, 1])
LOWER = np.array([50, 2, 50, 0.05, 0.95, 0.95, 0.95])
UPPER = np.array([350, 10, 3000, 4, 1.05, 1.05, 1.05])


def thermal(t, c):
    t = np.atleast_1d(t)

    def integral(temp, cv=False):
        y = c[2] / temp
        x = 0.5 * y[:, None] * (NODES + 1)
        f = x**4 * np.exp(x) / np.expm1(x) ** 2 if cv else x**3 / np.expm1(x)
        return 0.5 * np.sum(f * WEIGHTS, axis=1) / y**2

    energy = 9 * N * GAS_CONSTANT * t * integral(t)
    reference = 9 * N * GAS_CONSTANT * TR * integral(np.array([TR]))[0]
    return c[3] * (energy - reference) / (1000 * V0), 9 * N * GAS_CONSTANT * integral(
        t, cv=True
    )


def cold(v, c):
    eta = (V0 / np.asarray(v)) ** (1 / 3)
    a = 0.75 * (c[1] - 4)
    pressure = 1.5 * c[0] * (eta**7 - eta**5) * (1 + a * (eta**2 - 1))
    kt = (
        0.5
        * c[0]
        * eta
        * (
            (7 * eta**6 - 5 * eta**4) * (1 + a * (eta**2 - 1))
            + (eta**7 - eta**5) * 2 * a * eta
        )
    )
    return pressure, kt


def pressure(c, d, *, y=None, t=None):
    y = d["Y"] if y is None else y
    t = d["T"] if t is None else t
    pth, cv = thermal(t, c)
    scale = np.array([c[i] if i is not None else 1 for i in d["scale"]])
    v = y / scale
    k = d["KS"]
    if np.any(k):
        ks = v[k].copy()
        vv = V0 * (c[0] / ks) ** (1 / c[1])
        slope = c[3] ** 2 * cv[k] * t[k] / (1000 * V0**2)
        for _ in range(30):
            derivative = (
                cold(vv + 0.0001, c)[1] - cold(vv - 0.0001, c)[1]
            ) / 0.0002 + slope
            step = (cold(vv, c)[1] + slope * vv - ks) / derivative
            vv = np.clip(vv - step, 0.6 * V0, 1.12 * V0)
            if np.max(np.abs(step)) < 1e-10:
                break
        if np.max(np.abs(cold(vv, c)[1] + slope * vv - ks)) > 1e-7:
            raise ValueError("Grossular KS inversion failed")
        v[k] = vv
    return cold(v, c)[0] + pth


def prepared_rows():
    rows = []

    def add(dataset, r, p, t, y, sy, *, ks=False, sp=None, st=None, scale=None):
        rows.append(
            dict(
                dataset=dataset,
                row=int(r["source_row"]),
                P=p,
                T=t,
                Y=y,
                sY=sy,
                sP=sp,
                sT=st,
                KS=ks,
                scale=scale,
            )
        )

    for name, filename, scale in [
        ("milani_vt", "milani-2017-table-2.csv", 4),
        ("du_vt", "du-2015-preprint-table-2.csv", 5),
    ]:
        for r in observations(filename):
            st = float(
                r["temperature_sd"] if name == "milani_vt" else r["temperature_error_k"]
            )
            add(
                name,
                r,
                0.0001,
                float(r["temperature_k"]),
                float(r["cell_volume_a3"]) * Avogadro / 8e24,
                float(r["cell_volume_sd_a3"]) * Avogadro / 8e24,
                st=st,
                scale=scale,
            )
    pv = observations("milani-2017-table-3.csv")
    anchor = float(pv[0]["cell_volume_a3"])
    for r in pv[1:]:
        p, sp, _ = recalibrate(float(r["pressure_gpa"]), float(r["pressure_sd"]))
        add(
            "milani_pv",
            r,
            p,
            TR,
            V0 * float(r["cell_volume_a3"]) / anchor,
            V0 * float(r["cell_volume_sd_a3"]) / anchor,
            sp=sp,
        )
    for r in observations("kono-2010-table-1.csv"):
        add(
            "kono_ks",
            r,
            float(r["pressure_gpa"]),
            float(r["temperature_k"]),
            float(r["ks_gpa"]),
            float(r["ks_gpa_error"]),
            ks=True,
            sp=float(r["pressure_gpa_error"]),
            scale=6,
        )
    for r in observations("isaak-1992-table-3b.csv"):
        if r["ks_gpa_reported_error"]:
            add(
                "isaak_ks",
                r,
                0.0001,
                float(r["temperature_k"]),
                float(r["ks_gpa"]),
                float(r["ks_gpa_reported_error"]),
                ks=True,
            )
    for r in observations("bass-1989-ambient.csv"):
        add(
            "bass_ks",
            r,
            0.0001,
            TR,
            float(r["ks_gpa"]),
            float(r["ks_gpa_reported_error"]),
            ks=True,
        )
    return rows


def arrays(rows):
    return {
        key: np.array([r[key] if r[key] is not None else 0 for r in rows])
        for key in ["T", "P", "Y", "sY", "sP", "sT", "KS"]
    } | {"scale": [r["scale"] for r in rows]}


def effective_sigma(c, d):
    py = (
        pressure(c, d, y=d["Y"] + 0.0001) - pressure(c, d, y=d["Y"] - 0.0001)
    ) / 0.0002
    pt = (pressure(c, d, t=d["T"] + 0.01) - pressure(c, d, t=d["T"] - 0.01)) / 0.02
    return np.sqrt(d["sP"] ** 2 + (py * d["sY"]) ** 2 + (pt * d["sT"]) ** 2)


def solve(rows, initial):
    d, c, trace = arrays(rows), np.array(initial, dtype=float), []
    converged = False
    for iteration in range(30):
        sigma = effective_sigma(c, d)
        fit = least_squares(
            lambda cc: (pressure(cc, d) - d["P"]) / sigma,
            c,
            bounds=(LOWER, UPPER),
            jac="3-point",
            x_scale=SCALING,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-9,
            max_nfev=400,
        )
        shift = float(np.max(np.abs((fit.x - c) / SCALING)))
        c = fit.x
        trace.append(
            dict(
                iteration=iteration + 1,
                max_scaled_shift=shift,
                success=bool(fit.success),
                nfev=fit.nfev,
            )
        )
        if shift < 1e-8 and fit.success:
            converged = True
            break
    sigma = effective_sigma(c, d)
    residual = pressure(c, d) - d["P"]
    j = []
    for i, scale in enumerate(SCALING):
        h = scale * 1e-5
        plus, minus = c.copy(), c.copy()
        plus[i] += h
        minus[i] -= h
        j.append((pressure(plus, d) - pressure(minus, d)) / (2 * h * sigma))
    _, singular, vectors = np.linalg.svd(np.array(j).T * SCALING, full_matrices=False)
    if np.sum(singular > singular[0] * 1e-12) != len(c):
        raise ValueError("Conditional fit covariance is rank deficient")
    inverse = vectors.T / singular
    covariance = np.einsum("ik,jk->ij", inverse, inverse) * np.outer(SCALING, SCALING)
    chisq = float(np.sum((residual / sigma) ** 2))
    return dict(
        parameters=dict(zip(NAMES, c.tolist())),
        numerically_converged=converged,
        at_parameter_bound=[
            NAMES[i]
            for i in range(len(c))
            if min(c[i] - LOWER[i], UPPER[i] - c[i]) < 1e-5 * SCALING[i]
        ],
        irls_trace=trace,
        objective_sum_squares=chisq,
        degrees_of_freedom=len(rows) - len(c),
        conditional_reduced_chi_square=chisq / (len(rows) - len(c)),
        covariance_order=NAMES,
        conditional_unscaled_covariance=covariance.tolist(),
        conditional_unscaled_widths=dict(
            zip(NAMES, np.sqrt(np.diag(covariance)).tolist())
        ),
        scaled_jacobian_singular_values=singular.tolist(),
        rows=[
            r
            | dict(
                calculated_pressure_gpa=float(r["P"] + v),
                pressure_residual_gpa=float(v),
                effective_error_gpa=float(s),
            )
            for r, v, s in zip(rows, residual, sigma)
        ],
    )


if __name__ == "__main__":
    p, _ = source("grossular")
    rows = prepared_rows()
    starts = [[p[2], p[3], p[11], p[18], 1, 1, p[53]], [170, 4.8, 800, 1.15, 1, 1, 1]]
    results = [solve(rows, start) for start in starts]
    report = dict(
        scope="conditional_95_channel_subset_refit",
        author_workflow_reproduced=False,
        source_coefficients_changed=False,
        plan="refit-plan.json",
        channels=len(rows),
        fits=results,
        plan_sha256=hashlib.sha256((DATA / "refit-plan.json").read_bytes()).hexdigest(),
        input_sha256={
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                "milani-2017-table-2.csv",
                "milani-2017-table-3.csv",
                "du-2015-preprint-table-2.csv",
                "kono-2010-table-1.csv",
                "isaak-1992-table-3b.csv",
                "bass-1989-ambient.csv",
            ]
        },
    )
    (DATA / "subset-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
