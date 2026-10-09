"""Sensitivity of the Fei 2004 Pt refit to constraints and diagnostic weights.

These weights are explicit reconstruction choices, not recovered author weights.
Only printed source measurements/calibration pressures enter as fit targets.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts.reproduce_fei_2004_pressure_scales import load_rows, source_pressure

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/fei-2004-platinum-fit-weight-sensitivity.json"
NAMES = ("K0", "K0_prime", "gamma0", "q")
BASE = np.array([273.0, 4.8, 2.69, 0.5])


def main():
    rows = load_rows("platinum")
    v, t, p, sv, av, sav = [
        np.array([float(row[key]) for row in rows])
        for key in (
            "volume_a3",
            "temperature_k",
            "pressure_gpa",
            "volume_uncertainty_a3",
            "gold_volume_a3",
            "gold_volume_uncertainty_a3",
        )
    ]
    # Symmetric finite differences at published EOS coefficients. Weights stay
    # fixed during fitting. They omit T errors and Au/Pt calibration covariances.
    sigma_pt = (
        abs(
            source_pressure(v + sv, t, "platinum")
            - source_pressure(v - sv, t, "platinum")
        )
        / 2
    )
    sigma_au = (
        abs(source_pressure(av + sav, t, "gold") - source_pressure(av - sav, t, "gold"))
        / 2
    )
    combined = np.hypot(sigma_pt, sigma_au)
    all_rows = np.ones(len(rows), bool)
    hot = t > 300
    variants = (
        ("joint_unweighted", (0, 1, 2, 3), all_rows, np.ones(len(rows))),
        ("fixed_cold_unweighted", (2, 3), hot, np.ones(len(rows))),
        ("fixed_cold_gamma_unweighted", (3,), hot, np.ones(len(rows))),
        ("joint_pt_volume_weighted", (0, 1, 2, 3), all_rows, sigma_pt),
        ("joint_au_volume_weighted", (0, 1, 2, 3), all_rows, sigma_au),
        ("joint_combined_volume_weighted", (0, 1, 2, 3), all_rows, combined),
        ("fixed_cold_combined_volume_weighted", (2, 3), hot, combined),
    )
    result = {}
    for label, free, mask, sigma in variants:
        indices = list(free)

        def unpack(x):
            c = BASE.copy()
            c[indices] = x
            return c

        fits = []
        for q_start in (0.05, 0.5, 2.0):
            start = BASE.copy()
            start[3] = q_start
            fit = least_squares(
                lambda x: (
                    (source_pressure(v[mask], t[mask], "platinum", unpack(x)) - p[mask])
                    / sigma[mask]
                ),
                start[indices],
                x_scale="jac",
                max_nfev=2000,
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
            assert fit.success
            fits.append(fit)
        assert np.ptp([f.cost for f in fits]) < 1e-7
        fit = min(fits, key=lambda f: f.cost)
        c = unpack(fit.x)
        residual = source_pressure(v[mask], t[mask], "platinum", c) - p[mask]
        result[label] = {
            "rows": int(mask.sum()),
            "free_parameters": [NAMES[i] for i in indices],
            "fixed_parameters": [NAMES[i] for i in range(4) if i not in indices]
            + ["V0=60.38", "theta0=230"],
            "parameters": dict(zip(NAMES, c)),
            "unweighted_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "weighted_squared_residual_sum": float(2 * fit.cost),
            "published_weighted_squared_residual_sum": float(
                np.sum(
                    (
                        (source_pressure(v[mask], t[mask], "platinum", BASE) - p[mask])
                        / sigma[mask]
                    )
                    ** 2
                )
            ),
            "q_by_start": [float(unpack(f.x)[3]) for f in fits],
        }
    payload = {
        "scope": "Fei 2004 Pt, original 42 Table 2 observations and original pressure_gpa calibration column. No Dewaele data or Fei 2007 Au recalibration. No EOS record changed.",
        "published_parameters": dict(zip(NAMES, BASE)),
        "published_unweighted_rmse_all42_gpa": float(
            np.sqrt(np.mean((source_pressure(v, t, "platinum", BASE) - p) ** 2))
        ),
        "weight_qualification": "Approximate pressure-equivalent errors from symmetric +/- reported Au/Pt volume errors evaluated with published Fei 2004 EOSs; independent errors assumed and quadrature sum where combined. Weights fixed during fitting. Temperature uncertainty, calibration parameter uncertainty and correlated systematic errors excluded. Original author objective/weights not recovered; these are sensitivity diagnostics, not author-fit parity.",
        "variants": result,
        "source_sha256": hashlib.sha256(
            (ROOT / "peritheos/data/datasets/platinum-fei-2004-table2.csv").read_bytes()
        ).hexdigest(),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
