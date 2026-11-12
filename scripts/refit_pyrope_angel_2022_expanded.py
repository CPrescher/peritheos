"""Expanded, conditional recovered-input pyrope fits; no catalog changes.

Repository-only: python -m scripts.refit_pyrope_angel_2022_expanded.
Read expanded-refit-plan.json for selection, uncertainty and provenance limits.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re

import numpy as np
from scipy.optimize import brentq

from scripts.audit_garnet_angel_2022 import ROOT, source
from scripts.audit_pyrope_angel_2022_observations import DATA, observations
from scripts.refit_pyrope_angel_2022_subset import (
    LOWER,
    NAMES,
    SCALING,
    UPPER,
    V0,
    prepared_rows,
    solve,
)

EXPANDED_NAMES = NAMES + ["scale_sinogeikin_2000_ks"]


def token(value):
    """Source parentheses: last digits, or absolute uncertainty if decimal."""
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(?:\((\d+(?:\.\d+)?)\))?", value)
    if match is None:
        raise ValueError(f"Invalid source token: {value}")
    center, error = match.groups()
    if error is None:
        return float(center), None
    power = -len(center.split(".")[1]) if "." in center else 0
    return float(center), float(error) * (1 if "." in error else 10**power)


def quartz_pressure(v, *, updated):
    """Factual coefficients of author-distributed 1997/2016 quartz EOS files."""
    if updated:
        v0, k0, kp, kpp = 112.967644, 37.023453, 6.654295, -0.729518
    else:
        v0, k0, kp, kpp = 112.980888, 37.125893, 5.988268, None
    f = ((v0 / v) ** (2 / 3) - 1) / 2
    series = 1 + 1.5 * (kp - 4) * f
    if kpp is not None:
        series += 1.5 * (k0 * kpp + kp * (kp - 7) + 143 / 9) * f**2
    return 3 * k0 * f * (1 + 2 * f) ** 2.5 * series


def recalibrate(p, error):
    v = brentq(lambda x: quartz_pressure(x, updated=False) - p, 85, 114)
    h = 0.0001
    derivative = (
        quartz_pressure(v + h, updated=True) - quartz_pressure(v - h, updated=True)
    ) / (quartz_pressure(v + h, updated=False) - quartz_pressure(v - h, updated=False))
    return quartz_pressure(v, updated=True), error * abs(derivative), v


def row(dataset, source_row, p, t, y, sy, *, kind="V", sp=None, st=None, scale=None):
    return dict(
        dataset=dataset,
        row=source_row,
        kind=kind,
        P=p,
        T=t,
        Y=y,
        sY=sy,
        sP=sp,
        sT=st,
        scale_index=scale,
    )


def expanded_rows(include_temperature=False):
    rows = prepared_rows()
    for r in observations("zhang-1998-table-1.csv"):
        if r["medium"] == "ambient" or r["source_excluded"] == "true":
            continue
        p, sp = token(r["pressure_gpa_token"])
        v, sv = token(r["cell_volume_a3_token"])
        rows.append(
            row(
                "zhang_" + r["medium"],
                int(r["source_row"]),
                p,
                298,
                V0 * v / 1502.9,
                V0 * sv / 1502.9,
                sp=sp,
            )
        )
    for r in observations("du-2015-preprint-table-4.csv"):
        p = float(r["pressure_gpa"])
        if p == 0:
            continue
        v, sv = token(r["cell_volume_a3_token"])
        rows.append(
            row(
                "du_pv",
                int(r["source_row"]),
                p,
                298,
                V0 * v / 1503.2,
                V0 * sv / 1503.2,
                sp=float(r["pressure_error_gpa"]),
            )
        )
    path = ROOT / "peritheos/data/datasets/pyrope-milani-2015-table-a2-compression.csv"
    with path.open() as f:
        milani = list(csv.DictReader(f))
    for i, r in enumerate(milani[1:], 2):
        p, sp, quartz_v = recalibrate(
            float(r["pressure_gpa"]), float(r["pressure_sd_gpa"])
        )
        rows.append(
            row(
                "milani_pv",
                i,
                p,
                298,
                V0 * float(r["unit_cell_volume_a3"]) / 1505.987,
                V0 * float(r["unit_cell_volume_sd_a3"]) / 1505.987,
                sp=sp,
            )
            | {
                "original_pressure_gpa": float(r["pressure_gpa"]),
                "calculated_equivalent_quartz_volume_a3": quartz_v,
            }
        )
    for r in observations("chantel-2016-table-1.csv"):
        if int(r["source_row"]) == 1 or float(r["temperature_k"]) != 298:
            continue
        rho, srho = token(r["density_token"])
        vp, svp = token(r["vp_token"])
        vs, svs = token(r["vs_token"])
        p, i = float(r["calibrant_pressure_gpa"]), int(r["source_row"])
        v = V0 * 3.565 / rho
        k = rho * (vp**2 - 4 / 3 * vs**2)
        sk = np.sqrt(
            (k / rho * srho) ** 2
            + (2 * rho * vp * svp) ** 2
            + (8 / 3 * rho * vs * svs) ** 2
        )
        rows.append(row("chantel_pv", i, p, 298, v, v * srho / rho))
        rows.append(row("chantel_ks", i, p, 298, k, float(sk), kind="KS"))
    for r in observations("sinogeikin-2000-digitized-angel-figure-8.csv"):
        # Keep the extracted pressure, including the 0.002 GPa axis residual,
        # rather than rounding inferred source values. Bound is 0.01 GPa.
        rows.append(
            row(
                "sinogeikin_2000_ks",
                int(r["plot_row"]),
                float(r["pressure_gpa"]),
                298,
                float(r["ks_original_gpa"]),
                float(r["ks_error_bar_gpa"]),
                kind="KS",
                scale=6,
            )
        )
    if include_temperature:
        for r in observations("sinogeikin-2002-figure-2-pixels.csv"):
            # Pixel coordinates in the documented 4x PDF crop, clip rectangle
            # (300,390,550,640) pt. Readable pairs only; no fitted-line values.
            x = float(r["x_pixel"])
            t = (x - 122) * 800 / (703 - 122) + 273.15
            c11 = 300 + (70 - float(r["c11_center_y_pixel"])) * 30 / 328
            c12 = 110 + (452 - float(r["c12_center_y_pixel"])) * 10 / 109
            e11 = (
                (float(r["c11_bottom_y_pixel"]) - float(r["c11_top_y_pixel"]))
                * 15
                / 328
            )
            e12 = (
                (float(r["c12_bottom_y_pixel"]) - float(r["c12_top_y_pixel"])) * 5 / 109
            )
            k, sk = (c11 + 2 * c12) / 3, np.sqrt(e11**2 + 4 * e12**2) / 3
            rows.append(
                row(
                    "sinogeikin_2002_ks",
                    int(r["digitized_row"]),
                    0.0001,
                    t,
                    k,
                    float(sk),
                    kind="KS",
                    st=5,
                )
                | {"digitized_c11_gpa": c11, "digitized_c12_gpa": c12}
            )
    return rows


def run():
    p, _ = source("pyrope")
    initial = np.array([p[2], p[3], p[11], p[18], 1.001, 1.0003, p[53]])
    cases = {}
    for include_t in (False, True):
        rows = expanded_rows(include_t)
        fit = solve(
            rows,
            initial,
            list(range(7)),
            names=EXPANDED_NAMES,
            scaling=np.r_[SCALING, 1],
            lower=np.r_[LOWER, 0.9],
            upper=np.r_[UPPER, 1.1],
        )
        profiled = solve(
            rows,
            initial,
            [4, 5, 6],
            names=EXPANDED_NAMES,
            scaling=np.r_[SCALING, 1],
            lower=np.r_[LOWER, 0.9],
            upper=np.r_[UPPER, 1.1],
        )
        cases[str(len(rows))] = {
            "n_channels": len(rows),
            "refit": fit,
            "published_eos_scales_only": profiled,
            "counts": {
                name: sum(r["dataset"] == name for r in rows)
                for name in sorted({r["dataset"] for r in rows})
            },
        }
        print(
            len(rows),
            fit["parameters"],
            fit["conditional_reduced_chi_square"],
            fit["numerically_converged"],
            flush=True,
        )
    report = {
        "scope": "Expanded conditional recovered-input fits, not full 194-channel or author-workflow reproduction",
        "author_workflow_reproduced": False,
        "missing_selected_channels": 3,
        "published_parameters": dict(zip(EXPANDED_NAMES, initial.tolist())),
        "input_sha256": {
            f.name: hashlib.sha256(f.read_bytes()).hexdigest()
            for f in sorted(DATA.glob("*.csv"))
        },
        "plan_sha256": hashlib.sha256(
            (DATA / "expanded-refit-plan.json").read_bytes()
        ).hexdigest(),
        "cases": cases,
    }
    (DATA / "expanded-refits.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    run()
