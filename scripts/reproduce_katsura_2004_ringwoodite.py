#!/usr/bin/env python3
"""Audit the Katsura et al. (2004) Mg2SiO4 ringwoodite thermal EOS.

Primary source: Journal of Geophysical Research: Solid Earth 109, B12209,
doi:10.1029/2004JB003094. Table 2 is a direct transcription of the official
Wiley full-text table. The publication defines ``n`` as the number of atoms
per formula unit but does not print the numerical value used in its fit.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_rt_eos, fit_thermal_eos

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "docs" / "data" / "ringwoodite-katsura-2004-table2-pvt.csv"

V0_CELL_A3 = 524.8
FORMULA_UNITS_PER_CELL = 8.0
ATOMS_PER_FORMULA = 7.0
PUBLIC_TO_MODEL_SCALE = 6.02214076e23 * 1.0e-25 / FORMULA_UNITS_PER_CELL


def load_data() -> dict[str, np.ndarray]:
    """Load the directly transcribed Table 2 observations."""
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "temperature": np.array([float(row["temperature_k"]) for row in rows]),
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "pressure_uncertainty": np.array(
            [float(row["pressure_uncertainty_gpa"]) for row in rows]
        ),
        "volume_ratio": np.array([float(row["normalized_volume"]) for row in rows]),
        "volume_ratio_uncertainty": np.array(
            [float(row["normalized_volume_uncertainty"]) for row in rows]
        ),
    }


def curve_metrics(model: MieGruneisenDebye, data: dict[str, np.ndarray]) -> dict:
    """Return pressure residuals for a specified MGD model."""
    volume = data["volume_ratio"] * V0_CELL_A3 * PUBLIC_TO_MODEL_SCALE
    residual = model.pressure(volume, data["temperature"]) - data["pressure"]
    heated = data["temperature"] > 304.0
    return {
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "heated_pressure_rmse_gpa": float(np.sqrt(np.mean(residual[heated] ** 2))),
        "heated_mean_pressure_residual_gpa": float(np.mean(residual[heated])),
        "maximum_absolute_pressure_residual_gpa": float(np.max(np.abs(residual))),
    }


def make_model(
    *, n: float, debye_temperature_law: str = "integrated_gruneisen"
) -> MieGruneisenDebye:
    """Construct the candidate mapping with published rounded coefficients."""
    model_v0 = V0_CELL_A3 * PUBLIC_TO_MODEL_SCALE
    return MieGruneisenDebye(
        BM3(model_v0, 182.0, 4.6),
        Tr=300.0,
        theta0=846.0,
        gamma0=1.93,
        q=3.5,
        n=n,
        debye_temperature_law=debye_temperature_law,
    )


def reproduce() -> dict:
    """Return source-curve checks and explicit diagnostic refits."""
    data = load_data()
    model_v0 = V0_CELL_A3 * PUBLIC_TO_MODEL_SCALE
    volume = data["volume_ratio"] * model_v0
    room_temperature = data["temperature"] <= 304.0

    bm3_refit = fit_rt_eos(
        BM3,
        volume[room_temperature],
        data["pressure"][room_temperature],
        initial={"K0_prime": 4.6},
        fixed={"V0": model_v0, "K0": 182.0},
    )
    thermal_refit = fit_thermal_eos(
        MieGruneisenDebye,
        BM3(model_v0, 182.0, 4.6),
        volume,
        data["temperature"],
        data["pressure"],
        initial={"theta0": 846.0, "gamma0": 1.93, "q": 3.5},
        fixed={"Tr": 300.0, "n": ATOMS_PER_FORMULA},
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={
            "theta0": (1.0, 3000.0),
            "gamma0": (0.1, 5.0),
            "q": (-10.0, 20.0),
        },
    )

    # This deliberately nonphysical normalization is a diagnostic only. It is
    # not a proposal for an executable record: Mg2SiO4 has seven atoms per
    # formula unit, while n=5 is neither printed nor chemically defensible.
    hidden_normalization_diagnostic = make_model(n=5.0)

    return {
        "primary_data": {
            "rows": int(len(data["pressure"])),
            "pressure_range_gpa": [
                float(np.min(data["pressure"])),
                float(np.max(data["pressure"])),
            ],
            "temperature_range_k": [
                float(np.min(data["temperature"])),
                float(np.max(data["temperature"])),
            ],
            "volume_ratio_range": [
                float(np.min(data["volume_ratio"])),
                float(np.max(data["volume_ratio"])),
            ],
            "room_temperature_rows": int(np.count_nonzero(room_temperature)),
        },
        "published_mgd_with_chemical_n_7": curve_metrics(
            make_model(n=ATOMS_PER_FORMULA), data
        ),
        "published_mgd_with_n_7_variable_exponent_theta": curve_metrics(
            make_model(
                n=ATOMS_PER_FORMULA,
                debye_temperature_law="variable_exponent",
            ),
            data,
        ),
        "unweighted_300_304_k_bm3_refit": {
            "K0_prime": bm3_refit.parameters["K0_prime"],
            "K0_prime_standard_error": bm3_refit.standard_errors["K0_prime"],
            "pressure_rmse_gpa": float(np.sqrt(np.mean(bm3_refit.residuals**2))),
        },
        "unweighted_mgd_refit_with_chemical_n_7": {
            "theta0": thermal_refit.parameters["theta0"],
            "gamma0": thermal_refit.parameters["gamma0"],
            "q": thermal_refit.parameters["q"],
            "theta0_standard_error": thermal_refit.standard_errors["theta0"],
            "gamma0_standard_error": thermal_refit.standard_errors["gamma0"],
            "q_standard_error": thermal_refit.standard_errors["q"],
            "pressure_rmse_gpa": float(np.sqrt(np.mean(thermal_refit.residuals**2))),
        },
        "n_5_hidden_normalization_diagnostic_only": curve_metrics(
            hidden_normalization_diagnostic, data
        ),
        "decision": (
            "blocked: the printed definition requires n=7 for Mg2SiO4, but the "
            "published coefficients miss Table 2 severely under both Peritheos "
            "Debye-temperature laws; n=5 nearly restores the curve but is neither "
            "reported nor chemically valid"
        ),
    }


def main() -> None:
    """Print deterministic JSON diagnostics."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
