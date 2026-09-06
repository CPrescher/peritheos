"""Reproduce the Mookherjee et al. (2019) Phase Egg HP BM3 curve."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).parents[1]
MATERIAL = ROOT / "peritheos/data/materials/phase_egg.eosmat"
RECORD = "phase_egg_mookherjee_2019_bm3_hp_1"
DATASET = "phase_egg_mookherjee_2019_supplement_hp_pv"


def bm3_pressure(volume, v0, k0, k0_prime):
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def reproduce():
    document = json.loads(MATERIAL.read_text(encoding="utf-8"))
    record = next(x for x in document["eos_records"] if x["identifier"] == RECORD)
    dataset = next(x for x in document["datasets"] if x["identifier"] == DATASET)
    rows = np.asarray(dataset["rows"], dtype=float)
    volume, observed = rows.T
    parameters = record["eos"]["parameters"]
    calculated = bm3_pressure(
        volume, parameters["V0"], parameters["K0"], parameters["K0_prime"]
    )
    residual = calculated - observed
    fitted = least_squares(
        lambda values: bm3_pressure(volume, *values) - observed,
        [parameters["V0"], parameters["K0"], parameters["K0_prime"]],
    )
    return {
        "observations": len(rows),
        "published_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_max_abs_residual_gpa": float(np.max(np.abs(residual))),
        "diagnostic_unweighted_pv_fit": {
            "V0": float(fitted.x[0]),
            "K0": float(fitted.x[1]),
            "K0_prime": float(fitted.x[2]),
            "pressure_rmse_gpa": float(
                np.sqrt(np.mean((bm3_pressure(volume, *fitted.x) - observed) ** 2))
            ),
        },
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
