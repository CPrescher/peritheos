"""Reproduce the Mookherjee et al. (2015) 3.65 A-phase EOS audit.

The official MSA supplement contains the complete experimental P-V table but
does not state the exact source row mask, regression objective, covariance, or
weights.  We therefore report (1) residuals of the printed BM3/BM4 curves to
all observations and (2) a transparent ODR diagnostic on compression rows,
using the printed P and V uncertainties.  It validates coefficient identity;
it is not represented as an exact reconstruction of the private fit setup.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy import odr

ROOT = Path(__file__).parents[1]
DATA = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / ("mgsioh6-365a-mookherjee-2015-supplement-table1-pv.csv")
)

PUBLISHED = {
    "experimental_bm3": {
        "order": 3,
        "V0": 194.52,
        "K0": 83.0,
        "K0_prime": 4.9,
    },
    "experimental_bm4_sensitivity": {
        "order": 4,
        "V0": 194.52,
        "K0": 77.0,
        "K0_prime": 7.9,
        "K0_double_prime": -0.7,
    },
    "gga_bm4_model_crystal": {
        "order": 4,
        "V0": 202.02,
        "K0": 80.0,
        "K0_prime": 3.4,
        "K0_double_prime": -0.05,
    },
}


def finite_strain_pressure(
    volume: np.ndarray | float,
    v0: float,
    k0: float,
    k0_prime: float,
    k0_double_prime: float | None = None,
) -> np.ndarray:
    """Evaluate primary-source Equations 1-3 independently of Peritheos."""
    volume = np.asarray(volume, dtype=float)
    strain = 0.5 * ((volume / v0) ** (-2.0 / 3.0) - 1.0)
    if k0_double_prime is None:
        quadratic = 0.0
    else:
        quadratic = (
            k0 * k0_double_prime + (k0_prime - 4.0) * (k0_prime - 3.0) + 35.0 / 9.0
        ) * strain**2
    normalized_pressure = k0 + 1.5 * k0 * ((k0_prime - 4.0) * strain + quadratic)
    return 3.0 * strain * (1.0 + 2.0 * strain) ** 2.5 * normalized_pressure


def _load() -> dict[str, np.ndarray]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "pressure_sigma": np.array(
            [float(row["pressure_sigma_gpa"] or "nan") for row in rows]
        ),
        "volume": np.array([float(row["volume_a3_conventional_cell"]) for row in rows]),
        "volume_sigma": np.array([float(row["volume_sigma_a3"]) for row in rows]),
        "compression": np.array([row["path"] == "compression" for row in rows]),
    }


def _curve_diagnostics(data: dict[str, np.ndarray], key: str) -> dict[str, float]:
    parameters = PUBLISHED[key]
    calculated = finite_strain_pressure(
        data["volume"],
        parameters["V0"],
        parameters["K0"],
        parameters["K0_prime"],
        parameters.get("K0_double_prime"),
    )
    residual = calculated - data["pressure"]
    return {
        "observations": int(residual.size),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "pressure_maximum_absolute_residual_gpa": float(np.max(np.abs(residual))),
    }


def _compression_odr(data: dict[str, np.ndarray], order: int) -> dict[str, object]:
    mask = data["compression"] & np.isfinite(data["pressure_sigma"])
    volume = data["volume"][mask]
    pressure = data["pressure"][mask]
    volume_sigma = data["volume_sigma"][mask]
    pressure_sigma = data["pressure_sigma"][mask]

    def model(beta: np.ndarray, x: np.ndarray) -> np.ndarray:
        if order == 3:
            return finite_strain_pressure(x, 194.52, beta[0], beta[1])
        return finite_strain_pressure(x, 194.52, beta[0], beta[1], beta[2])

    initial = [83.0, 4.9] if order == 3 else [77.0, 7.9, -0.7]
    fit = odr.ODR(
        odr.RealData(volume, pressure, sx=volume_sigma, sy=pressure_sigma),
        odr.Model(model),
        beta0=initial,
        maxit=1000,
    ).run()
    parameter_names = (
        ["K0", "K0_prime"] if order == 3 else ["K0", "K0_prime", "K0_double_prime"]
    )
    return {
        "observations": int(volume.size),
        "fixed_parameters": {"V0": 194.52},
        "parameters": {
            name: float(value) for name, value in zip(parameter_names, fit.beta)
        },
        "standard_errors": {
            name: float(value) for name, value in zip(parameter_names, fit.sd_beta)
        },
        "reduced_weighted_residual_variance": float(fit.res_var),
    }


def reproduce() -> dict[str, object]:
    """Return deterministic curve, refit, and convention diagnostics."""
    data = _load()
    theoretical = PUBLISHED["gga_bm4_model_crystal"]
    checkpoints = np.array([140.0, 160.0, 180.0, 200.0, 220.0])
    return {
        "official_supplement": {
            "observations": int(data["pressure"].size),
            "compression_observations_with_uncertainties": int(
                np.sum(data["compression"] & np.isfinite(data["pressure_sigma"]))
            ),
            "experiments": ["I", "II", "III"],
        },
        "published_curve_residuals_all_rows": {
            key: _curve_diagnostics(data, key)
            for key in ("experimental_bm3", "experimental_bm4_sensitivity")
        },
        "compression_odr_diagnostic": {
            "experimental_bm3": _compression_odr(data, 3),
            "experimental_bm4_sensitivity": _compression_odr(data, 4),
        },
        "gga_bm4_source_equation_checkpoints": {
            "volumes_a3": checkpoints.tolist(),
            "pressures_gpa": finite_strain_pressure(
                checkpoints,
                theoretical["V0"],
                theoretical["K0"],
                theoretical["K0_prime"],
                theoretical["K0_double_prime"],
            ).tolist(),
        },
    }


def main() -> None:
    """Print the reproduction as stable JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
