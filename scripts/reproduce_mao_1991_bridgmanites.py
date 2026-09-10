"""Reproduce the Mao et al. (1991) bridgmanite fit and precursor result."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
FINAL_DATASET = ROOT / "peritheos/data/datasets/mao-1991-tables2-3-bridgmanite-pv.csv"
PRECURSOR_DATASET = (
    ROOT / "peritheos/data/datasets/mao-1989-carnegie-table8-bridgmanite-pv.csv"
)
PUBLISHED = {
    "MgSiO3": {"V0": 162.49, "V0_error": 0.07},
    "Mg0.9Fe0.1SiO3": {"V0": 162.79, "V0_error": 0.08},
    "Mg0.8Fe0.2SiO3": {"V0": 163.53, "V0_error": 0.10},
}
ANDERSON_1989_GOLD = {"V0": 67.79, "K0": 166.65, "K0_prime": 5.4823}


def bm3_pressure(volume, v0, k0=261.0, kp=4.0):
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def murnaghan_pressure(v_over_v0, k0, kp):
    ratio = np.asarray(v_over_v0, dtype=float)
    return (k0 / kp) * (ratio ** (-kp) - 1.0)


def _scaled_standard_errors(result, observations):
    degrees_of_freedom = observations - result.x.size
    covariance = np.linalg.inv(result.jac.T @ result.jac)
    covariance *= np.sum(result.fun**2) / degrees_of_freedom
    return np.sqrt(np.diag(covariance))


def _fixed_k0_fit(pressure, ratio, family):
    if family == "birch":
        unit_pressure = bm3_pressure(ratio, 1.0, k0=1.0, kp=4.0)
    elif family == "murnaghan":
        unit_pressure = murnaghan_pressure(ratio, 1.0, 4.0)
    else:  # pragma: no cover - internal programming error
        raise ValueError(f"Unsupported family: {family}")
    k0 = float(np.dot(unit_pressure, pressure) / np.dot(unit_pressure, unit_pressure))
    residuals = k0 * unit_pressure - pressure
    degrees_of_freedom = pressure.size - 1
    k0_error = float(
        np.sqrt(
            np.sum(residuals**2)
            / degrees_of_freedom
            / np.dot(unit_pressure, unit_pressure)
        )
    )
    return {
        "K0_gpa": k0,
        "K0_error_gpa": k0_error,
        "K0_prime": 4.0,
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "objective": "unweighted_pressure_residuals",
    }


def _load_precursor():
    with PRECURSOR_DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    pressure = np.asarray([float(row["pressure_gpa"]) for row in rows])
    ratio = np.asarray([float(row["v_over_v0"]) for row in rows])
    return rows, pressure, ratio


def _load_final():
    with FINAL_DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    selected = [row for row in rows if row["fit_included"] == "1"]
    pressure = np.asarray([float(row["pressure_gpa"]) for row in selected])
    volume = np.asarray([float(row["volume_a3"]) for row in selected])
    ratio = volume / PUBLISHED["Mg0.9Fe0.1SiO3"]["V0"]
    return rows, selected, pressure, ratio


def _anderson_gold_diagnostic(rows):
    pressure = np.asarray([float(row["pressure_gpa"]) for row in rows])
    volume = np.asarray([float(row["gold_lattice_a"]) ** 3 for row in rows])
    bundled_pressure = bm3_pressure(
        volume,
        ANDERSON_1989_GOLD["V0"],
        ANDERSON_1989_GOLD["K0"],
        ANDERSON_1989_GOLD["K0_prime"],
    )
    bundled_residuals = bundled_pressure - pressure
    effective_v0 = least_squares(
        lambda parameter: (
            bm3_pressure(
                volume,
                parameter[0],
                ANDERSON_1989_GOLD["K0"],
                ANDERSON_1989_GOLD["K0_prime"],
            )
            - pressure
        ),
        x0=np.asarray([ANDERSON_1989_GOLD["V0"]]),
    )
    effective_residuals = (
        bm3_pressure(
            volume,
            effective_v0.x[0],
            ANDERSON_1989_GOLD["K0"],
            ANDERSON_1989_GOLD["K0_prime"],
        )
        - pressure
    )
    return {
        "reference_eos_record": "gold_anderson_1989_bm3_1",
        "observations": len(rows),
        "bundled_parameters": ANDERSON_1989_GOLD,
        "bundled_pressure_rmse_gpa": float(np.sqrt(np.mean(bundled_residuals**2))),
        "bundled_pressure_residual_range_gpa": [
            float(np.min(bundled_residuals)),
            float(np.max(bundled_residuals)),
        ],
        "effective_V0_a3": float(effective_v0.x[0]),
        "effective_V0_pressure_rmse_gpa": float(
            np.sqrt(np.mean(effective_residuals**2))
        ),
        "interpretation": (
            "The calibration is executable and all Au lattice parameters are "
            "available. The remaining small systematic offset is consistent with "
            "the Anderson record's V0 being derived from a rounded ambient density."
        ),
    }


def reproduce() -> dict[str, object]:
    final_rows, selected, final_pressure, final_ratio = _load_final()
    final_fit = _fixed_k0_fit(final_pressure, final_ratio, "birch")
    published_final_residuals = (
        bm3_pressure(final_ratio, 1.0, k0=261.0, kp=4.0) - final_pressure
    )

    rows, pressure, ratio = _load_precursor()
    free_murnaghan = least_squares(
        lambda parameters: murnaghan_pressure(ratio, *parameters) - pressure,
        x0=np.asarray([275.0, 3.7]),
    )
    free_errors = _scaled_standard_errors(free_murnaghan, pressure.size)
    final_residuals = bm3_pressure(ratio, 1.0, k0=261.0, kp=4.0) - pressure

    fractions = np.array([1.0, 0.95, 0.90])
    curves = {
        formula: dict(
            zip(
                map(str, fractions),
                map(float, bm3_pressure(fractions * values["V0"], values["V0"])),
            )
        )
        for formula, values in PUBLISHED.items()
    }
    counts = {
        formula: sum(row["formula"] == formula for row in rows) for formula in PUBLISHED
    }
    return {
        "shared_K0_gpa": 261.0,
        "fixed_K0_prime": 4.0,
        "published_final": PUBLISHED,
        "curves_gpa": curves,
        "final": {
            "dataset": FINAL_DATASET.relative_to(ROOT).as_posix(),
            "observations": len(final_rows),
            "fit_observations": len(selected),
            "fit_composition": "Mg0.9Fe0.1SiO3",
            "selection": "all nine 300 K Fe10 rows in final-article Table 2",
            "reported_weights": None,
            "reproduction_objective": "unweighted_pressure_residuals",
            "fixed_kp_birch": final_fit,
            "published_curve": {
                "K0_gpa": 261.0,
                "pressure_rmse_gpa": float(
                    np.sqrt(np.mean(published_final_residuals**2))
                ),
                "maximum_absolute_residual_gpa": float(
                    np.max(np.abs(published_final_residuals))
                ),
            },
            "composition_counts": {
                formula: sum(row["formula"] == formula for row in final_rows)
                for formula in PUBLISHED
            },
            "anderson_1989_gold_recalculation": _anderson_gold_diagnostic(final_rows),
        },
        "precursor": {
            "dataset": PRECURSOR_DATASET.relative_to(ROOT).as_posix(),
            "observations": len(rows),
            "composition_counts": counts,
            "selection": "all 12 Table 8 room-temperature rows",
            "reported_weights": None,
            "fixed_kp_birch": _fixed_k0_fit(pressure, ratio, "birch"),
            "fixed_kp_murnaghan": _fixed_k0_fit(pressure, ratio, "murnaghan"),
            "free_murnaghan": {
                "K0_gpa": float(free_murnaghan.x[0]),
                "K0_error_gpa": float(free_errors[0]),
                "K0_prime": float(free_murnaghan.x[1]),
                "K0_prime_error": float(free_errors[1]),
                "pressure_rmse_gpa": float(np.sqrt(np.mean(free_murnaghan.fun**2))),
                "objective": "unweighted_pressure_residuals",
            },
            "final_261_curve_on_precursor_rows": {
                "pressure_rmse_gpa": float(np.sqrt(np.mean(final_residuals**2))),
                "maximum_absolute_residual_gpa": float(np.max(np.abs(final_residuals))),
            },
            "is_final_publication_fit_input": False,
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
