"""Audit and refit the Solomatova et al. (2016) crossover branches.

The Table 7 BM3 coefficients are endmembers of a coupled MINUTI spin model,
not independent fits to pressure windows. This reproduction therefore fits
recoverable compression trajectories with the same three Fe2+ spin levels and
Table 2 priors. It never treats Table 7 or derived grids as observations.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets"
TABLE1 = DATA / "mg0490fe0483ti0027o-solomatova-2016-table1-pv.csv"
TABLE3 = DATA / "mg0490fe0483ti0027o-solomatova-2016-table3-crossover-grid.csv"
MAO = DATA / "mg075fe025o-mao-2011-figure1-300k-digitized.csv"
MARQUARDT = DATA / "mg090fe010o-marquardt-2009-table2-pv.csv"
CHEN = DATA / "mg065fe035o-chen-2012-table1-pv.csv"
FEI = DATA / "mg061fe039o-fei-2007-table-s2-pv.csv"
ZHURAVLEV = DATA / "mg061fe039o-zhuravlev-2010-table5-pv.csv"
LIN = DATA / "mg083fe017o-lin-2005-figure2-300k-digitized.csv"

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

TRANSITIONS = {
    "mg090_marq": (52.1, 2.3),
    "mg083_lin": (49.2, 8.2),
    "mg075_mao": (62.6, 15.9),
    "mg065_chen": (63.7, 13.3),
    "mg061_fei": (56.6, 17.2),
    "mg061_zhuravlev": (73.5, 10.3),
    "fp48": (68.8, 18.4),
    "mg040_lin": (78.6, 25.1),
}

RECORD_SOURCES = {
    "mg090_marq": {
        "doi": "10.1016/j.epsl.2009.08.017",
        "status": "similar",
        "dataset_identifiers": ["mg090fe010o_marquardt_2009_table2_pv"],
    },
    "mg083_lin": {
        "doi": "10.1038/nature03825",
        "status": "not_refittable",
        "dataset_identifiers": ["mg083fe017o_lin_2005_figure2_300k_digitized"],
        "reason": "The primary figure supplies 43 V/V0,HS observations but not the absolute V0,HS normalization, so all six absolute Table 7 coefficients are not independently identifiable.",
    },
    "mg075_mao": {
        "doi": "10.1029/2011GL049915",
        "status": "parity_not_achieved",
        "dataset_identifiers": ["mg075fe025o_mao_2011_figure1_300k_digitized"],
    },
    "mg065_chen": {
        "doi": "10.1029/2012JB009162",
        "status": "similar",
        "dataset_identifiers": ["mg065fe035o_chen_2012_table1_pv"],
    },
    "mg061_fei": {
        "doi": "10.1029/2007GL030712",
        "status": "similar",
        "dataset_identifiers": ["mg061fe039o_fei_2007_table_s2_pv"],
    },
    "mg061_zhuravlev": {
        "doi": "10.1007/s00269-009-0347-6",
        "status": "similar",
        "dataset_identifiers": ["mg061fe039o_zhuravlev_2010_table5_pv"],
    },
    "fp48": {
        "doi": "10.2138/am-2016-5510",
        "status": "similar",
        "dataset_identifiers": ["mg0490fe0483ti0027o_solomatova_2016_table1_pv"],
    },
    "mg040_lin": {
        "doi": "10.1038/nature03825",
        "status": "not_refittable",
        "reason": "Lin et al. report the 84-102 GPa spectroscopic interval and an approximately 1.6% volume drop at 95 GPa, but no P-V series for this composition; those constraints are not observations.",
    },
}

_KB_EV_K = 8.617333262145e-5
_GPA_A3_PER_EV = 160.2176634
_TEMPERATURE_K = 300.0
_DEGENERACIES = np.array([15.0, 18.0, 1.0])
_ALPHA = np.array([0.0, 1.0, 2.0])
_UNPAIRED_FRACTION = np.array([1.0, 0.5, 0.0])


def bm3_pressure(
    volume: np.ndarray | float, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def bm3_energy_ev(
    volume: np.ndarray | float, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    volume = np.asarray(volume, dtype=float)
    strain = 0.5 * ((v0 / volume) ** (2.0 / 3.0) - 1.0)
    energy_gpa_a3 = 4.5 * v0 * k0 * strain**2 * (1.0 + (k0_prime - 4.0) * strain)
    return energy_gpa_a3 / _GPA_A3_PER_EV


def spin_populations(omega_ev: np.ndarray | float) -> np.ndarray:
    omega = np.asarray(omega_ev, dtype=float)
    logits = np.log(_DEGENERACIES) - np.expand_dims(omega, -1) * _ALPHA / (
        _KB_EV_K * _TEMPERATURE_K
    )
    logits -= np.max(logits, axis=-1, keepdims=True)
    weights = np.exp(logits)
    return weights / np.sum(weights, axis=-1, keepdims=True)


def mixed_pressure(
    volume: np.ndarray | float, parameters: np.ndarray, fe_fraction: float
) -> np.ndarray:
    v0_hs, k0_hs, kp_hs, v0_ls, k0_ls, omega0 = parameters
    volume = np.asarray(volume, dtype=float)
    elastic_difference = bm3_energy_ev(volume, v0_ls, k0_ls, 4.0) - bm3_energy_ev(
        volume, v0_hs, k0_hs, kp_hs
    )
    # Four cation sites occur in the conventional B1 cell; alpha_LS is two.
    omega = omega0 + elastic_difference / (8.0 * fe_fraction)
    low_spin_equivalent = spin_populations(omega) @ _ALPHA / 2.0
    pressure_hs = bm3_pressure(volume, v0_hs, k0_hs, kp_hs)
    pressure_ls = bm3_pressure(volume, v0_ls, k0_ls, 4.0)
    return pressure_hs + low_spin_equivalent * (pressure_ls - pressure_hs)


def _load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _initial_parameters(key: str, fe_fraction: float) -> np.ndarray:
    hs = TABLE_7[f"{key}_hs"]
    ls = TABLE_7[f"{key}_ls"]
    midpoint = TRANSITIONS[key][0]
    volume_midpoint = brentq(
        lambda volume: (
            0.5 * (bm3_pressure(volume, *hs) + bm3_pressure(volume, *ls)) - midpoint
        ),
        30.0,
        hs[0],
    )
    omega_at_half = brentq(
        lambda omega: spin_populations(omega) @ _UNPAIRED_FRACTION - 0.5,
        -1.0,
        1.0,
    )
    elastic_difference = bm3_energy_ev(volume_midpoint, *ls) - bm3_energy_ev(
        volume_midpoint, *hs
    )
    omega0 = omega_at_half - elastic_difference / (8.0 * fe_fraction)
    return np.array([hs[0], hs[1], hs[2], ls[0], ls[1], omega0])


def _transition_pressures(parameters: np.ndarray, fe_fraction: float) -> list[float]:
    pressures = []
    for unpaired_fraction in (0.8, 0.5, 0.2):

        def objective(volume: float) -> float:
            difference = bm3_energy_ev(
                volume, parameters[3], parameters[4], 4.0
            ) - bm3_energy_ev(volume, parameters[0], parameters[1], parameters[2])
            omega = parameters[5] + difference / (8.0 * fe_fraction)
            return float(
                spin_populations(omega) @ _UNPAIRED_FRACTION - unpaired_fraction
            )

        grid = np.linspace(35.0, 83.0, 1000)
        values = [objective(volume) for volume in grid]
        roots = [
            brentq(objective, left, right)
            for left, right, f_left, f_right in zip(
                grid[:-1], grid[1:], values[:-1], values[1:]
            )
            if f_left * f_right < 0.0
        ]
        if not roots:
            raise RuntimeError("spin-population target is not bracketed")
        pressures.append(float(mixed_pressure(roots[-1], parameters, fe_fraction)))
    return pressures


def _fit_absolute(
    key: str,
    path: Path,
    fe_fraction: float,
    volume_column: str,
    selection_column: str | None = None,
    fix_v0_hs: bool = False,
    cubic_lattice_column: bool = False,
    priors: tuple[float, float, float, float, float, float] = (
        160.0,
        5.0,
        170.0,
        20.0,
        4.0,
        0.5,
    ),
) -> dict[str, Any]:
    rows = _load(path)
    if selection_column is not None:
        rows = [
            row
            for row in rows
            if row[selection_column]
            == ("compression" if selection_column == "compression_path" else "1")
        ]
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row[volume_column]) for row in rows])
    if cubic_lattice_column:
        volume = volume**3
    initial = _initial_parameters(key, fe_fraction)

    if fix_v0_hs:
        fixed_v0 = initial[0]
        free_initial = initial[1:]

        def expand(free: np.ndarray) -> np.ndarray:
            return np.r_[fixed_v0, free]

        lower = [80.0, 1.0, 60.0, 80.0, -3.0]
        upper = [300.0, 8.0, 90.0, 400.0, 3.0]
        x_scale = [30.0, 1.0, 3.0, 30.0, 0.2]
    else:
        free_initial = initial

        def expand(free: np.ndarray) -> np.ndarray:
            return free

        lower = [65.0, 80.0, 1.0, 60.0, 80.0, -3.0]
        upper = [90.0, 300.0, 8.0, 90.0, 400.0, 3.0]
        x_scale = [3.0, 30.0, 1.0, 3.0, 30.0, 0.2]

    def residuals(free: np.ndarray) -> np.ndarray:
        parameters = expand(free)
        calculated_pressure = mixed_pressure(volume, parameters, fe_fraction)
        step = 0.002
        dp_dv = (
            mixed_pressure(volume + step, parameters, fe_fraction)
            - mixed_pressure(volume - step, parameters, fe_fraction)
        ) / (2.0 * step)
        approximate_volume_residual = (calculated_pressure - pressure) / np.abs(dp_dv)
        # Table 2 priors and prior windows. The original MINUTI input files and
        # row weights are not published, so observations receive equal weight.
        prior_residuals = np.array(
            [
                (parameters[1] - priors[0]) / priors[1],
                (parameters[4] - priors[2]) / priors[3],
                (parameters[2] - priors[4]) / priors[5],
            ]
        )
        return np.r_[approximate_volume_residual, prior_residuals]

    fit = least_squares(
        residuals,
        free_initial,
        bounds=(lower, upper),
        x_scale=x_scale,
        max_nfev=500,
    )
    parameters = expand(fit.x)
    data_residuals = residuals(fit.x)[: len(rows)]
    transition_pressures = _transition_pressures(parameters, fe_fraction)
    return {
        "observations": len(rows),
        "selection": (
            f"{selection_column}=1" if selection_column else "all observation rows"
        ),
        "objective": "equal-weight approximate volume residuals plus Table 2 priors",
        "solver_success": bool(fit.success),
        "rmse_volume_a3": float(np.sqrt(np.mean(data_residuals**2))),
        "parameters": {
            "V0_HS": float(parameters[0]),
            "K0_HS": float(parameters[1]),
            "K0_prime_HS": float(parameters[2]),
            "V0_LS": float(parameters[3]),
            "K0_LS": float(parameters[4]),
            "K0_prime_LS": 4.0,
        },
        "transition_pressure_20_percent_ls_gpa": transition_pressures[0],
        "transition_pressure_50_percent_ls_gpa": transition_pressures[1],
        "transition_pressure_80_percent_ls_gpa": transition_pressures[2],
        "transition_width_20_80_gpa": transition_pressures[2] - transition_pressures[0],
    }


def _fit_lin_normalized() -> dict[str, Any]:
    rows = _load(LIN)
    # The public plot omits absolute V0,HS. Fixing it to Table 7 makes this a
    # shape diagnostic, not an independent absolute-volume refit.
    anchor_v0 = TABLE_7["mg083_lin_hs"][0]
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = anchor_v0 * np.array([float(row["volume_ratio_v_v0_hs"]) for row in rows])
    initial = _initial_parameters("mg083_lin", 0.17)

    def expand(free: np.ndarray) -> np.ndarray:
        return np.r_[anchor_v0, free[0], free[1], free[2], free[3], free[4]]

    def residuals(free: np.ndarray) -> np.ndarray:
        parameters = expand(free)
        calculated_pressure = mixed_pressure(volume, parameters, 0.17)
        step = 0.002
        dp_dv = (
            mixed_pressure(volume + step, parameters, 0.17)
            - mixed_pressure(volume - step, parameters, 0.17)
        ) / (2.0 * step)
        return np.r_[
            (calculated_pressure - pressure) / np.abs(dp_dv),
            (parameters[1] - 160.0) / 5.0,
            (parameters[4] - 170.0) / 20.0,
            (parameters[2] - 4.0) / 0.5,
        ]

    fit = least_squares(
        residuals,
        initial[1:],
        bounds=(
            [80.0, 1.0, 60.0, 80.0, -3.0],
            [300.0, 8.0, 90.0, 400.0, 3.0],
        ),
        x_scale=[30.0, 1.0, 3.0, 30.0, 0.2],
        max_nfev=500,
    )
    parameters = expand(fit.x)
    transitions = _transition_pressures(parameters, 0.17)
    return {
        "observations": len(rows),
        "status": "normalization_anchored_diagnostic_only",
        "anchor": "V0_HS fixed to Table 7; not an independent absolute-volume refit",
        "solver_success": bool(fit.success),
        "rmse_volume_a3": float(np.sqrt(np.mean(residuals(fit.x)[: len(rows)] ** 2))),
        "parameters": {
            "V0_HS_fixed": float(parameters[0]),
            "K0_HS": float(parameters[1]),
            "K0_prime_HS": float(parameters[2]),
            "V0_LS": float(parameters[3]),
            "K0_LS": float(parameters[4]),
            "K0_prime_LS": 4.0,
        },
        "transition_pressure_50_percent_ls_gpa": transitions[1],
        "transition_width_20_80_gpa": transitions[2] - transitions[0],
    }


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

    refits = {
        "mg090_marq": _fit_absolute(
            "mg090_marq",
            MARQUARDT,
            0.10,
            "volume_a3_conventional_cell",
            "used_in_solomatova_2016_fit",
            priors=(161.4, 3.0, 170.0, 20.0, 4.0, 0.5),
        ),
        "mg075_mao": _fit_absolute(
            "mg075_mao", MAO, 0.25, "cell_volume_a3", fix_v0_hs=True
        ),
        "mg065_chen": _fit_absolute(
            "mg065_chen",
            CHEN,
            0.35,
            "volume_a3_conventional_cell",
            "used_in_solomatova_2016_fit",
        ),
        "mg061_zhuravlev": _fit_absolute(
            "mg061_zhuravlev",
            ZHURAVLEV,
            0.39,
            "volume_a3_conventional_cell",
            "used_in_solomatova_2016_fit",
        ),
        "mg061_fei": _fit_absolute(
            "mg061_fei",
            FEI,
            0.39,
            "sample_lattice_a_angstrom",
            "compression_path",
            cubic_lattice_column=True,
        ),
        "fp48": _fit_absolute(
            "fp48", TABLE1, 0.483, "volume_a3_conventional_cell", fix_v0_hs=True
        ),
    }
    for key, result in refits.items():
        result["ledger_status"] = RECORD_SOURCES[key]["status"]

    return {
        "accepted_reference_branches": len(TABLE_7),
        "table1_observations": len(table1),
        "table3_grid_rows": len(table3),
        "source_audit": RECORD_SOURCES,
        "coupled_refits": refits,
        "lin_normalized_shape_diagnostic": _fit_lin_normalized(),
        "fp48_endmember_diagnostics": {
            "high_spin_pre_crossover_rows": len(raw_pre),
            "high_spin_pressure_rmse_gpa": float(np.sqrt(np.mean(hs_residual**2))),
            "low_spin_post_crossover_rows": len(grid_post),
            "low_spin_pressure_rmse_gpa": float(np.sqrt(np.mean(ls_residual**2))),
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
