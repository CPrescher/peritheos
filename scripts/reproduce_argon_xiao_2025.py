"""Independent Xiao (2025) checks against publisher workbook and neutron data.

No source VBA executes. Equation 18 and SciPy quadrature implement the
manuscript independently. Workbook cached values are frozen primary checks.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
R = 8.31451  # Author workbook VBA, not the modern SI value.
V00 = 22.555  # cm3/mol


def properties(v, t):
    """Return A,U,S,P,Cv,K; J/mol, K, MPa, cm3/mol conventions."""
    x = v / V00
    lz = -np.log(x)
    gamma = 2.68 * x**0.0024
    theta = 86.44 * np.exp(2.68 / 0.0024 * (1 - x**0.0024))
    cold = V00 * (2656.5 * lz**2 / 2 + 7298 * lz**3 / 3 + 10 * lz**4 / 4)
    pc = (2656.5 * lz + 7298 * lz**2 + 10 * lz**3) / x
    kc = (2656.5 + (2656.5 + 2 * 7298) * lz + (7298 + 30) * lz**2 + 10 * lz**3) / x
    if t == 0:
        return dict(A=cold, U=cold, S=0.0, P=pc, Cv=0.0, K=kc)
    y = theta / t
    d3 = 3 / y**3 * quad(lambda z: z**3 / np.expm1(z), 0, min(y, 700), epsabs=1e-12)[0]
    ad = R * t * (3 * np.log(-np.expm1(-y)) - d3)
    ud = 3 * R * t * d3
    cd = 3 * R * (4 * d3 - 3 * y / np.expm1(y))
    s = t / 86.44
    den = 1 + 0.388 * s * s
    f = s**4 / den
    f1 = 2 * s**3 * (2 + 0.388 * s * s) / den**2
    f2 = 2 * s * s * (6 + 3 * 0.388 * s * s + 0.388**2 * s**4) / den**3
    fv = np.exp(7.85 * (x - 1))
    aa = 0.0128 * R * 86.44 * f * fv
    ua = 0.0128 * R * 86.44 * (f - s * f1) * fv
    ca = -0.0128 * R * s * f2 * fv
    pa = -7.85 / V00 * aa
    kd = gamma / v * ((1 - 0.0024 + gamma) * ud - gamma * t * cd)
    ka = -v * 7.85 / V00 * pa
    return dict(
        A=cold + ad + aa,
        U=cold + ud + ua,
        S=(ud - ad) / t - 0.0128 * R * f1 * fv,
        P=pc + gamma * ud / v + pa,
        Cv=cd + ca,
        K=kc + kd + ka,
    )


def reproduce():
    # Source: official XLSM sample calculation D27:D38, D13, D43, D51.
    expected = dict(
        A=-998.2296606304144,
        U=1098.2008906760243,
        S=29.949007875806267,
        P=78.3517375843629,
        Cv=22.854917560882143,
        K=1 / 0.0004448515520965651,
    )
    actual = properties(23.0, 70.0)
    cases = []
    for t, p, v in [
        (70.0, 10.0, 23.827471251113984),
        (83.806, 0.068891, 24.604805616824237),
        (60.0, 0.0006965135426981063, 23.608512889347836),
    ]:
        calc = brentq(lambda vol: properties(vol, t)["P"] - p, 18, 26)
        cases.append(
            dict(
                temperature_k=t,
                pressure_mpa=p,
                source_volume_cm3_mol=v,
                calculated_volume_cm3_mol=calc,
                relative_difference=(calc - v) / v,
            )
        )
    with (ROOT / "peritheos/data/datasets/argon-fcc-xiao-2025-table5.csv").open() as f:
        data = list(csv.DictReader(f))
    # Source says pressure is slightly above sublimation; zero is only the
    # stated negligible-pressure approximation, never a measured coordinate.
    diffs = []
    for row in data:
        t = float(row["temperature_k"])
        v = float(row["molar_volume_cm3_mol"])
        calc = brentq(lambda vol: properties(vol, t)["P"], 20, 25)
        diffs.append((calc - v) / v)
    out = dict(
        reference="10.1007/s10765-024-03469-2",
        workbook_checkpoint=dict(
            temperature_k=70,
            volume_cm3_mol=23,
            expected=expected,
            independent=actual,
            relative_differences={k: (actual[k] - v) / v for k, v in expected.items()},
        ),
        workbook_inversions=cases,
        neutron_table5=dict(
            rows=len(data),
            pressure_assumption="P=0 approximation only; source pressure unreported and slightly above sublimation",
            relative_volume_aad_percent=float(np.mean(np.abs(diffs)) * 100),
            relative_volume_rms_percent=float(np.sqrt(np.mean(np.square(diffs))) * 100),
        ),
        equal_weight_diagnostic=dict(
            script="scripts/refit_argon_xiao_2025.py",
            report="docs/data/argon-xiao-2025-equal-weight-refit.json",
            objective="Equal weight per relative property residual across all entries in each recovered-data case",
            published_parameters_replaced=False,
        ),
        global_refit=dict(
            status="not_reproduced",
            reason="The original empirical property weights are published in Table2. Complete legacy row selection and quantitative penalties/iteration protocol remain unavailable. A separate user-directed equal-relative-weight recovered-data diagnostic is performed; it is not the original regression.",
        ),
    )
    return out


if __name__ == "__main__":
    output = ROOT / "docs/data/argon-xiao-2025-reproduction.json"
    output.write_text(json.dumps(reproduce(), indent=2, allow_nan=False) + "\n")
    print(output)
