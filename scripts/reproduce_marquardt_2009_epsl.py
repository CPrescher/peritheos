"""Validate EPSL HS/LS BM3 fits; never replace published coefficients."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
RECORD = "mg090fe010o_marquardt_2009b_hs_bm3"
DATASET = "mg090fe010o_marquardt_2009_table2_pv"
LS_RECORDS = tuple(
    f"mg090fe010o_marquardt_2009b_ls_bm3_s1_{i:02d}" for i in range(1, 13)
)


def bm3_pressure(volume, v0, k0, kp):
    """Independent Eulerian finite-strain BM3 expression, GPa and cell A^3."""
    f = 0.5 * ((v0 / np.asarray(volume)) ** (2 / 3) - 1)
    return 3 * k0 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (kp - 4) * f)


def reproduce():
    path = ROOT / "peritheos/data/datasets/mg090fe010o-marquardt-2009-table2-pv.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    hs = [row for row in rows if float(row["pressure_gpa"]) < 45]
    p = np.array([float(row["pressure_gpa"]) for row in hs])
    v = np.array([float(row["volume_a3_conventional_cell"]) for row in hs])
    sigma = np.array([float(row["volume_sigma_a3"]) for row in hs])

    def volumes(k0, kp):
        return np.array(
            [
                brentq(lambda x: bm3_pressure(x, 75.62, k0, kp) - pressure, 35, 76)
                for pressure in p
            ]
        )

    # Diagnostic volume-residual fit at fixed measured V0; no invented Brillouin prior.
    fit = least_squares(
        lambda x: (volumes(*x) - v) / sigma,
        [158.2, 3.98],
        bounds=([50, 0], [400, 10]),
    )
    published_v = volumes(158.2, 3.98)
    return {
        "record_identifier": RECORD,
        "purpose": "validation_only_pv_refit_without_Brillouin_constraint",
        "selection": "Table 2 P < 45 GPa; V0 anchored at 75.62 A^3",
        "observations": len(hs),
        "source_orders": [int(row["source_order"]) for row in hs],
        "published_pressure_rmse_gpa": float(
            np.sqrt(np.mean((bm3_pressure(v, 75.62, 158.2, 3.98) - p) ** 2))
        ),
        "published_volume_rmse_a3": float(np.sqrt(np.mean((published_v - v) ** 2))),
        "published_max_abs_volume_residual_sigma": float(
            np.max(np.abs((published_v - v) / sigma))
        ),
        "validation_parameters": {
            "V0": 75.62,
            "K0": float(fit.x[0]),
            "K0_prime": float(fit.x[1]),
        },
        "solver_success": bool(fit.success),
        "qualification": "P-V-only validation cannot reconstruct the source's unspecified ambient Brillouin constraint weight. No mixed-spin/LS rows or Science SOM rows are fit.",
    }


def reproduce_supplement():
    """Check Table S1's LS models against observations and derived KT outputs."""
    data = ROOT / "peritheos/data/datasets"
    with (data / "mg090fe010o-marquardt-2009-table2-pv.csv").open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if float(r["pressure_gpa"]) > 63]
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    v = np.array([float(r["volume_a3_conventional_cell"]) for r in rows])
    sigma = np.array([float(r["volume_sigma_a3"]) for r in rows])
    with (data / "mg090fe010o-marquardt-2009-table-s1-ls-models.csv").open(
        newline=""
    ) as stream:
        models = [r for r in csv.DictReader(stream) if r["source_row"].isdigit()]
    results = []
    for row in models:
        v0 = float(row["v0_ls_a3_conventional_cell"])
        k0 = float(row["k0_gpa"])
        kp = float(row["k0_prime"])

        def volume(pressure, modulus):
            return brentq(lambda x: bm3_pressure(x, v0, modulus, kp) - pressure, 35, v0)

        def residual(modulus):
            return np.array([volume(pi, modulus) for pi in p]) - v

        fit = least_squares(lambda x: residual(x[0]) / sigma, [k0])
        checkpoints = []
        for pressure, column in [
            (69, "kt_69_gpa"),
            (70.1, "kt_70_1_gpa"),
            (75.8, "kt_75_8_gpa"),
            (81.2, "kt_81_2_gpa"),
        ]:
            x = volume(pressure, k0)
            f = 0.5 * ((v0 / x) ** (2 / 3) - 1)
            a = 1.5 * (kp - 4)
            kt = k0 * (1 + 2 * f) ** 2.5 * (1 + (7 + 2 * a) * f + 9 * a * f * f)
            checkpoints.append(
                {
                    "pressure_gpa": pressure,
                    "calculated_kt_gpa": kt,
                    "printed_kt_gpa": float(row[column]),
                }
            )
        results.append(
            {
                "source_row": int(row["source_row"]),
                "published_k0_gpa": k0,
                "validation_k0_gpa": float(fit.x[0]),
                "solver_success": bool(fit.success),
                "published_volume_rmse_a3": float(np.sqrt(np.mean(residual(k0) ** 2))),
                "ls_kt_checkpoints": checkpoints,
            }
        )
    return {
        "purpose": "validation_only_fixed_V0_and_K0_prime_LS_refits",
        "source_orders": [int(r["source_order"]) for r in rows],
        "models": results,
        "ambient_table_s2_ks_gpa": 165,
        "ambient_table_s2_kt_conversion_gpa": 165 / (1 + 31.2e-6 * 1.524 * 300),
        "qualification": "Only six LS XRD rows constrain these fits. S1 KT columns are derived checks, not observations. Mixed-spin columns and min/max/avg summary rows are excluded from independent EOS checks. No statistical errors are published for the S1 coefficients, so uncertainty parity is unestablished.",
    }


if __name__ == "__main__":
    print(
        json.dumps(
            {"high_spin": reproduce(), "supplement": reproduce_supplement()}, indent=2
        )
    )
