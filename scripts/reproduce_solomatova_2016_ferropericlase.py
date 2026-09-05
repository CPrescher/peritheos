"""Check the Solomatova et al. (2016) spin-crossover reference branches."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
TABLE1 = (
    ROOT / "peritheos/data/datasets/mg0490fe0483ti0027o-solomatova-2016-table1-pv.csv"
)
TABLE3 = (
    ROOT
    / "peritheos/data/datasets/mg0490fe0483ti0027o-solomatova-2016-table3-crossover-grid.csv"
)

TABLE_7 = {
    "mg090_marq_hs": (75.55, 159.0, 3.96),
    "mg090_marq_ls": (74.59, 159.0, 4.00),
    "mg083_lin_hs": (75.94, 160.0, 4.04),
    "mg083_lin_ls": (72.29, 190.0, 4.00),
    "mg075_mao_hs": (76.34, 160.0, 4.28),
    "mg075_mao_ls": (73.74, 174.0, 4.00),
    "mg065_chen_hs": (77.10, 162.0, 3.99),
    "mg065_chen_ls": (73.77, 171.0, 4.00),
    "mg061_fei_hs": (77.49, 161.0, 4.25),
    "mg061_fei_ls": (74.83, 162.0, 4.00),
    "mg061_zhuravlev_hs": (77.41, 160.0, 4.07),
    "mg061_zhuravlev_ls": (73.56, 170.0, 4.00),
    "fp48_hs": (77.29, 160.0, 4.12),
    "fp48_ls": (73.64, 173.0, 4.00),
    "mg040_lin_hs": (77.90, 159.0, 3.82),
    "mg040_lin_ls": (73.83, 169.0, 4.00),
}


def bm3_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def _load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def reproduce() -> dict[str, object]:
    table1 = _load(TABLE1)
    table3 = _load(TABLE3)
    raw_pre = [row for row in table1 if float(row["pressure_gpa"]) <= 44.4]
    grid_post = [row for row in table3 if float(row["pressure_gpa"]) >= 88.0]
    hs = TABLE_7["fp48_hs"]
    ls = TABLE_7["fp48_ls"]
    hs_residual = bm3_pressure(
        np.array([float(row["volume_a3_conventional_cell"]) for row in raw_pre]), *hs
    ) - np.array([float(row["pressure_gpa"]) for row in raw_pre])
    ls_residual = bm3_pressure(
        np.array([float(row["volume_a3_conventional_cell"]) for row in grid_post]), *ls
    ) - np.array([float(row["pressure_gpa"]) for row in grid_post])
    ratios = np.array([1.0, 0.9, 0.8, 0.7])
    return {
        "accepted_reference_branches": len(TABLE_7),
        "table1_observations": len(table1),
        "table3_grid_rows": len(table3),
        "fp48_endmember_diagnostics": {
            "high_spin_pre_crossover_rows": len(raw_pre),
            "high_spin_pressure_rmse_gpa": float(np.sqrt(np.mean(hs_residual**2))),
            "low_spin_post_crossover_rows": len(grid_post),
            "low_spin_pressure_rmse_gpa": float(np.sqrt(np.mean(ls_residual**2))),
        },
        "source_equation_checkpoints": {
            key: bm3_pressure(ratios * v0, v0, k0, k0_prime).tolist()
            for key, (v0, k0, k0_prime) in TABLE_7.items()
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
