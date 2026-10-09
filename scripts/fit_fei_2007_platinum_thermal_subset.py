"""Fit Fei (2004) hot Pt rows with a frozen Fei (2007) cold EOS."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts.reproduce_fei_2007_platinum import (
    DATA,
    PARAMETERS,
    observations,
    source_pressure,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/fei-2007-platinum-fixed-rt-thermal-subset.json"


def fit_thermal(volume, temperature, target, cold, q_fixed=None, q_min=None):
    def unpack(x):
        coefficients = cold.copy()
        coefficients[3] = x[0]
        coefficients[4] = x[1] if q_fixed is None else q_fixed
        return coefficients

    starts = (0.05, 0.5, 2.0) if q_fixed is None else (None,)
    solutions = []
    for start in starts:
        initial = [cold[3], start] if q_fixed is None else [cold[3]]
        bounds = (
            ([-np.inf, q_min], [np.inf, np.inf])
            if q_min is not None
            else (-np.inf, np.inf)
        )
        result = least_squares(
            lambda x: (
                source_pressure(volume, temperature, coefficients=unpack(x)) - target
            ),
            initial,
            bounds=bounds,
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
            max_nfev=2000,
        )
        assert result.success
        solutions.append(result)
    rmses = [float(np.sqrt(np.mean(r.fun**2))) for r in solutions]
    assert max(rmses) - min(rmses) < 1e-8
    result = solutions[1] if len(solutions) > 1 else solutions[0]
    coefficients = unpack(result.x)
    assert np.array_equal(coefficients[[0, 1, 2, 5]], cold[[0, 1, 2, 5]])
    if q_min is not None:
        assert coefficients[4] >= q_min
    return {
        "gamma0": float(coefficients[3]),
        "q": float(coefficients[4]),
        "q_fixed": q_fixed,
        "q_lower_bound": q_min,
        "q_at_lower_bound": bool(q_min is not None and coefficients[4] - q_min < 1e-7),
        "hot_rmse_gpa": float(np.sqrt(np.mean(result.fun**2))),
        "coefficients": coefficients.tolist(),
        "solver_success": bool(result.success),
        "q_start_sensitivity": [
            dict(q_start=s, hot_rmse_gpa=r) for s, r in zip(starts, rmses)
        ],
    }


def main():
    volume, temperature, target, paired, reduced = observations()
    rt = temperature == 300
    hot = temperature > 300
    base = np.array(PARAMETERS["platinum"])
    variants = {}
    for label, cold_mask in (
        ("all_43_rt_rows", rt),
        ("dewaele_36_rt_rows", np.arange(len(volume)) < 36),
        ("published_cold_coefficients", None),
    ):
        cold = base.copy()
        if cold_mask is not None:

            def residual(x):
                c = base.copy()
                c[2] = x[0]
                return (
                    source_pressure(
                        volume[cold_mask], temperature[cold_mask], coefficients=c
                    )
                    - target[cold_mask]
                )

            result = least_squares(
                residual, [base[2]], xtol=1e-11, ftol=1e-11, gtol=1e-11
            )
            assert result.success
            cold[2] = result.x[0]
        frozen_rt = source_pressure(volume[rt], temperature[rt], coefficients=cold)
        fits = {
            "unconstrained_q": fit_thermal(
                volume[hot], temperature[hot], target[hot], cold
            ),
            "q_nonnegative": fit_thermal(
                volume[hot], temperature[hot], target[hot], cold, q_min=0.0
            ),
            "q_fixed_at_published": fit_thermal(
                volume[hot], temperature[hot], target[hot], cold, q_fixed=0.5
            ),
        }
        for fit in fits.values():
            coefficients = np.array(fit["coefficients"])
            delta = float(
                np.max(
                    abs(
                        source_pressure(
                            volume[rt], temperature[rt], coefficients=coefficients
                        )
                        - frozen_rt
                    )
                )
            )
            assert delta < 1e-10
            fit["rt_pressure_max_change_gpa"] = delta
        variants[label] = {
            "cold_parameters_fixed_in_thermal_stage": dict(
                zip(("V0", "K0", "K0_prime", "theta0"), cold[[0, 1, 2, 5]].tolist())
            ),
            "rt_rows_used_in_cold_fit": int(sum(cold_mask))
            if cold_mask is not None
            else None,
            "all_43_rt_rmse_gpa": float(
                np.sqrt(np.mean((frozen_rt - target[rt]) ** 2))
            ),
            "hot_rmse_with_published_gamma_q_gpa": float(
                np.sqrt(
                    np.mean(
                        (
                            source_pressure(
                                volume[hot], temperature[hot], coefficients=cold
                            )
                            - target[hot]
                        )
                        ** 2
                    )
                )
            ),
            "thermal_fits": fits,
        }
    report = {
        "objective": "Equal-weight pressure residuals. Fit K0-prime to RT observations first with V0 and K0 fixed, then freeze the complete cold EOS and fit gamma0/q only to the 35 hot Fei (2004) rows. No RT residuals enter the thermal objective.",
        "pressure_calibration": "Hot pressures independently recalculated from measured Au volumes using published Fei (2007) Au coefficients. Original source data and canonical EOS records unchanged.",
        "counts": {"rt": int(sum(rt)), "hot": int(sum(hot))},
        "hot_rows_per_temperature": {
            str(int(t)): int(sum(temperature[hot] == t))
            for t in sorted(set(temperature[hot]))
        },
        "hot_volume_a3_range": [float(volume[hot].min()), float(volume[hot].max())],
        "hot_pressure_gpa_range": [float(target[hot].min()), float(target[hot].max())],
        "hot_inputs": [
            dict(
                run=row["run"],
                temperature_k=float(row["temperature_k"]),
                platinum_volume_a3=float(row["volume_a3"]),
                gold_volume_a3=float(row["gold_volume_a3"]),
                recalculated_au_pressure_gpa=float(p),
            )
            for row, p in zip(paired, reduced)
            if float(row["temperature_k"]) > 300
        ],
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in (
                "platinum-dewaele-2004-table1-compression.csv",
                "platinum-fei-2004-table2.csv",
            )
        },
        "variants": variants,
        "qualification": "Diagnostic subset fits with explicit frozen cold coefficients. Original author weights, fitting constraints and calibration covariance remain unrecovered. RMS values in this report refer to the hot subset unless explicitly labeled RT. No symmetric boundary parameter errors assigned.",
    }
    OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for label, variant in variants.items():
        print(label, variant["cold_parameters_fixed_in_thermal_stage"])
        for name, fit in variant["thermal_fits"].items():
            print(
                name,
                "gamma0",
                fit["gamma0"],
                "q",
                fit["q"],
                "hot RMS",
                fit["hot_rmse_gpa"],
            )


if __name__ == "__main__":
    main()
