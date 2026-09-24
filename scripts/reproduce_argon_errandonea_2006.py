"""Independent BM3 equation check and partial-figure diagnostics, not a refit claim."""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/data/argon-errandonea-2006-reproduction.json"
RECORD = "argon_fcc_errandonea_2006_bm3"


def bm3(volume, v0=143.0, k0=6.5, kp=5.1):
    """Eulerian strain expansion; pressure in GPa, volume in any shared basis."""
    strain = ((v0 / np.asarray(volume)) ** (2 / 3) - 1) / 2
    return 3 * k0 * strain * (1 + 2 * strain) ** 2.5 * (1 + 1.5 * (kp - 4) * strain)


def reproduce():
    from peritheos import Material, get_material_document

    path = (
        ROOT / "peritheos/data/datasets/argon-fcc-errandonea-2006-figure5-digitized.csv"
    )
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    volume = np.array([float(r["volume_a3"]) for r in rows])
    pressure = np.array([float(r["pressure_gpa"]) for r in rows])
    model = Material.from_eosmat(get_material_document("argon_fcc")).get_eos_record(
        RECORD
    )
    calculated = bm3(volume)
    refit = least_squares(
        lambda pars: bm3(volume, *pars) - pressure,
        [143, 6.5, 5.1],
        bounds=([80, 0.01, 1], [300, 50, 15]),
        x_scale="jac",
        max_nfev=5000,
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    grid = np.linspace(48, 143, 51)

    def rmse(residual):
        return float(np.sqrt(np.mean(residual**2)))

    return {
        "record_id": RECORD,
        "equation_status": "verified",
        "reproduction_status": "not_reproduced",
        "reason": "Partial figure observations; full raw data, fit weights and covariance unavailable.",
        "observation_count": len(rows),
        "published_parameters": {"V0": 143.0, "K0": 6.5, "K0_prime": 5.1},
        "native_max_difference_gpa": float(
            np.max(np.abs(model.pressure(grid, 300) - bm3(grid)))
        ),
        "basis_invariance_max_difference_gpa": float(
            np.max(np.abs(bm3(grid) - bm3(grid / 4, v0=143 / 4)))
        ),
        "published_partial_data_rmse_gpa": rmse(calculated - pressure),
        "diagnostic_equal_pressure_weight_refit": {
            "parameters": refit.x.tolist(),
            "rmse_gpa": rmse(bm3(volume, *refit.x) - pressure),
            "solver_success": bool(refit.success),
            "active_bound_mask": refit.active_mask.tolist(),
            "note": "Eight digitized points only; V0 reaches the analyst's 300 A^3 bound. Poorly constrained, not a reproduction or replacement of source coefficients.",
        },
        "observations": [
            {
                "pressure_gpa": float(p),
                "volume_a3": float(v),
                "calculated_pressure_gpa": float(c),
                "residual_gpa": float(c - p),
            }
            for p, v, c in zip(pressure, volume, calculated)
        ],
    }


if __name__ == "__main__":
    REPORT.write_text(json.dumps(reproduce(), indent=2) + "\n")
    print(REPORT)
