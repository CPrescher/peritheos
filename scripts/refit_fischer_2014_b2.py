#!/usr/bin/env python3
"""Independent, physically normalized B2 FeSi fits and uncertainty diagnostics.

Run with --output docs/data/fischer-2014-b2-refits.json to regenerate evidence.
The production objective is unweighted pressure residuals on all 114 exact B2
rows. The reported pressure uncertainties already depend on temperature;
their unknown joint covariance precludes inventing a full measurement-error
likelihood. Alternative weights, reference assumptions and pressure blocks
are sensitivity diagnostics, not additional production fits or confidence
intervals. This script does not use the Peritheos evaluator or optimizer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import least_squares

from scripts.reproduce_fischer_2014_fesi import ROOT, load_data, parameters, pressure

PARAMETER_ORDER = ("K0", "gamma0", "q")
DATA_PATH = ROOT / "peritheos/data/datasets/fischer-2014-table-s2-pvt.csv"
REPORT_PATH = ROOT / "docs/data/fischer-2014-b2-refits.json"


def reproduce():
    """Return full covariance, convergence and explicitly conditional diagnostics."""
    rows = load_data("fesi_b2")
    v = np.array([float(r["b2_volume_a3"]) for r in rows]) * Avogadro / 1e24 / 2
    t = np.array([float(r["temperature_k"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    sigma = np.array([float(r["pressure_uncertainty_gpa"]) for r in rows])
    all_rows = np.ones(len(p), dtype=bool)
    output = {
        "format": "peritheos.fischer-2014-b2-refit-evidence",
        "format_version": 1,
        "source_data_sha256": hashlib.sha256(DATA_PATH.read_bytes()).hexdigest(),
        "source_rows": [int(row["source_row"]) for row in rows],
        "parameter_order": list(PARAMETER_ORDER),
        "fixed_parameters": {
            "V0_cm3_mol_atoms": 6.414,
            "K0_prime": 4.17,
            "theta0_k": 417.0,
            "Tr_k": 300.0,
        },
        "normalization": "Debye energy per physical atom and volume per mole of atoms; equivalent to n=2 and Z=1 for a B2 FeSi cell.",
        "objective": "unweighted_pressure_residuals",
        "uncertainty": "Conditional local 1-sigma errors: covariance=(J.T J)^-1*SSE/(114-3). Observations are treated as independent with equal residual variance. These errors omit fixed-input uncertainty, run correlations, calibrant systematics and model discrepancy. Sensitivity cases are not a probability distribution.",
        "solver": {
            "method": "scipy.optimize.least_squares",
            "bounds": [[10.0, 0.01, -5.0], [500.0, 5.0, 10.0]],
            "x_scale": "jac",
            "ftol": 1e-11,
            "xtol": 1e-11,
            "gtol": 1e-11,
        },
        "models": {},
    }
    for model in ("bm3", "vinet"):
        published = parameters("fesi_b2", model)

        def fit(mask, start, *, weighted=False, fixed_override=None):
            coeff = published.copy()
            for i, value in (fixed_override or {}).items():
                coeff[i] = value

            def predict(x):
                trial = coeff.copy()
                trial[[1, 4, 5]] = x
                return pressure("fesi_b2", v, t, model, trial)

            def residual(x):
                r = predict(x) - p
                return (r / sigma if weighted else r)[mask]

            solver = least_squares(
                residual,
                start,
                bounds=output["solver"]["bounds"],
                x_scale="jac",
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
            if not solver.success:
                raise RuntimeError(solver.message)
            residuals = predict(solver.x) - p
            covariance = (
                np.linalg.pinv(solver.jac.T @ solver.jac)
                * np.sum(solver.fun**2)
                / (int(mask.sum()) - 3)
            )
            covariance = (covariance + covariance.T) / 2
            errors = np.sqrt(np.diag(covariance))
            result = {
                "parameters": dict(zip(PARAMETER_ORDER, solver.x.tolist())),
                "standard_errors": dict(zip(PARAMETER_ORDER, errors.tolist())),
                "covariance": covariance.tolist(),
                "correlation": (covariance / np.outer(errors, errors)).tolist(),
                "observations": int(mask.sum()),
                "degrees_of_freedom": int(mask.sum()) - 3,
                "pressure_rmse_gpa": float(np.sqrt(np.mean(residuals[mask] ** 2))),
                "residual_range_gpa": [
                    float(residuals[mask].min()),
                    float(residuals[mask].max()),
                ],
                "solver_success": bool(solver.success),
                "jacobian_rank": int(np.linalg.matrix_rank(solver.jac)),
                "active_bounds": solver.active_mask.tolist(),
            }
            return result, residuals, solver.x

        full, residuals, optimum = fit(all_rows, published[[1, 4, 5]])
        full["source_checkpoints"] = [
            {
                "source_row": int(rows[i]["source_row"]),
                "pressure_gpa": float(p[i]),
                "temperature_k": float(t[i]),
                "cell_volume_a3": float(rows[i]["b2_volume_a3"]),
                "refit_pressure_gpa": float(p[i] + residuals[i]),
            }
            for i in (0, len(rows) // 2, len(rows) - 1)
        ]
        multi = [
            fit(all_rows, [published[1], gamma, q])[0]
            for gamma, q in ((1.2, 0.3), (1.2, 3.0), (2.5, 1.7), (4.0, 0.3), (4.0, 3.0))
        ]
        weighted = fit(all_rows, optimum, weighted=True)[0]
        weighted["qualification"] = (
            "Pressure-sigma-only sensitivity; ignores predictor error and P-T covariance, and is not the production objective."
        )
        fixed_sensitivity = {}
        for label, override in (
            ("V0_6_435_cm3_mol_atoms", {0: 6.435}),
            ("theta0_300_k", {3: 300.0}),
            ("theta0_600_k", {3: 600.0}),
            ("K0_prime_3_67", {2: 3.67}),
            ("K0_prime_4_67", {2: 4.67}),
        ):
            fixed_sensitivity[label] = fit(all_rows, optimum, fixed_override=override)[
                0
            ]
        blocks = []
        # Preserve source order for tied pressures across NumPy sort backends.
        for indices in np.array_split(np.argsort(p, kind="stable"), 3):
            mask = all_rows.copy()
            mask[indices] = False
            block, r, _ = fit(mask, optimum)
            block["held_out_source_rows"] = [
                int(rows[i]["source_row"]) for i in indices
            ]
            block["held_out_pressure_range_gpa"] = [
                float(p[indices].min()),
                float(p[indices].max()),
            ]
            block["held_out_pressure_rmse_gpa"] = float(
                np.sqrt(np.mean(r[indices] ** 2))
            )
            blocks.append(block)
        output["models"][model] = {
            "fit": full,
            "multistart": multi,
            "pressure_sigma_only_sensitivity": weighted,
            "fixed_parameter_sensitivity": fixed_sensitivity,
            "fixed_sensitivity_qualification": "V0=6.435 is the nearby Sata comparator in the source table. Theta0=300/600 K and K0-prime=4.17+/-0.5 are analyst-selected stress tests, not source errors or priors.",
            "pressure_block_holdouts": blocks,
            "holdout_qualification": "Three equal-count pressure blocks; omitted blocks test prediction across the observed pressure range. Run identifiers and inter-observation correlations are unavailable, so this is not validation on independent experiments.",
        }
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(reproduce(), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
