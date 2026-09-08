#!/usr/bin/env python3
"""Audit the accepted tranche-B EOS records and held coesite-V candidate."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from peritheos.eos.rt import BM3 as BM3Model

ROOT = Path(__file__).resolve().parents[1]
PHASE_H_DATA = (
    ROOT / "peritheos/data/datasets/phase-h-tsuchiya-mookherjee-2015-table1-pv.csv"
)
COESITE_IV_DATA = (
    ROOT / "peritheos/data/datasets/coesite-iv-bykova-2018-table10-calc-pv.csv"
)
COESITE_V_DATA = (
    ROOT / "peritheos/data/datasets/coesite-v-bykova-2018-table10-calc-pv.csv"
)
COESITE_I_II_DATA = (
    ROOT / "peritheos/data/datasets/coesite-i-ii-cernok-2014-table1-pv.csv"
)
COESITE_II_III_DATA = (
    ROOT / "peritheos/data/datasets/coesite-ii-iii-bykova-2018-table2-pv.csv"
)

BM3 = {
    "stishovite_akber_knutson_2002_vibc_300k_bm3": (46.53, 313.2, 4.4),
    "ca_perovskite_akber_knutson_2002_vibc_pnma_300k_bm3": (45.90, 228.0, 4.3),
    "ca_perovskite_akber_knutson_2002_vibc_pm3m_300k_bm3": (45.92, 229.0, 4.3),
    "ca_perovskite_akber_knutson_2002_vib_pm3m_300k_bm3": (44.75, 305.0, 3.5),
    "phase_h_tsuchiya_mookherjee_2015_gga_static_bm3": (58.9, 147.5, 4.9),
    "ca_perovskite_sherman_1993_basis_b_static_bm3": (44.96, 301.0, 4.0),
    "stishovite_schoelmerich_2020_shock_300k_bm3": (46.5, 307.0, 4.66),
    "bridgmanite_liu_2011_gga_static_bm3": (162.88, 241.0, 4.0),
    "ca_perovskite_liu_2007_lda_static_bm3": (45.46, 240.0, 4.15),
    "coesite_i_iii_bykova_2018_300k_bm3": (547.20, 103.0, 3.02),
    "coesite_iv_bykova_2018_am05_static_bm3_refit": (438.6702, 168.52, 3.34),
}

COESITE_V_HELD_CANDIDATE = (427.3921208, 185.26, 3.10)

BM2 = {
    "ca_perovskite_wang_weidner_1994_bm2": (45.83, 280.0),
    "fe088sio3_bridgmanite_ismailova_2016_300k_bm2": (178.98, 190.0),
    "klb1_mg_perovskite_ricolleau_2009_bm2_alphakt": (164.0, 245.0),
    "klb1_ca_perovskite_ricolleau_2009_bm2_alphakt": (45.60, 244.0),
    "klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt": (76.44, 158.0),
    "klb1_ferropericlase_ricolleau_2009_low_spin_300k_bm2": (74.04, 170.0),
}

# V0, K0, alpha0, alpha1, dK/dT for the source's AlphaKT form.  The
# isothermal reference temperature is 300 K in every case.
THERMAL_BM2 = {
    "klb1_mg_perovskite_ricolleau_2009_bm2_alphakt": (
        164.0,
        245.0,
        3.19e-5,
        8.8e-9,
        -0.036,
    ),
    "klb1_ca_perovskite_ricolleau_2009_bm2_alphakt": (
        45.60,
        244.0,
        3.06e-5,
        8.7e-9,
        -0.035,
    ),
    "klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt": (
        76.44,
        158.0,
        2.20e-5,
        3.61e-8,
        -0.034,
    ),
}

VINET = {
    "wollastonite_sagatova_2021_gga_300k_vinet": (410.48, 91.5, 5.5),
    "pseudowollastonite_sagatova_2021_gga_300k_vinet": (811.65, 74.3, 4.7),
    "breyite_sagatova_2021_gga_300k_vinet": (380.069, 70.38, 4.3),
    "casio2o5_titanite_sagatova_2021_gga_300k_vinet": (329.324, 155.2, 4.2),
    "larnite_sagatova_2021_gga_300k_vinet": (353.053, 100.96, 4.53),
    "ca_perovskite_tetragonal_sagatova_2021_lda_300k_vinet": (178.808, 244.0, 4.2),
    "ca_perovskite_tetragonal_sagatova_2021_gga_300k_vinet": (187.312, 212.0, 4.2),
}


def bm3_pressure(volume, v0: float, k0: float, kp: float):
    eta = (v0 / np.asarray(volume)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def bm2_pressure(volume, v0: float, k0: float):
    return bm3_pressure(volume, v0, k0, 4.0)


def vinet_pressure(volume, v0: float, k0: float, kp: float):
    x = (np.asarray(volume) / v0) ** (1.0 / 3.0)
    return 3.0 * k0 * (1.0 - x) / (x * x) * np.exp(1.5 * (kp - 1.0) * (1.0 - x))


def thermal_bm2_pressure(
    volume,
    temperature: float,
    v0: float,
    k0: float,
    alpha0: float,
    alpha1: float,
    dk_dt: float,
):
    tr = 300.0
    vt = v0 * np.exp(
        alpha0 * (temperature - tr) + 0.5 * alpha1 * (temperature**2 - tr**2)
    )
    kt = k0 + dk_dt * (temperature - tr)
    return bm2_pressure(volume, vt, kt)


def reproduce() -> dict[str, object]:
    with PHASE_H_DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    phase_h_fit = bm3_pressure(
        volume, *BM3["phase_h_tsuchiya_mookherjee_2015_gga_static_bm3"]
    )
    residual = phase_h_fit - pressure
    coesite_checks = {}
    for name, path, parameters in (
        (
            "coesite_iv",
            COESITE_IV_DATA,
            BM3["coesite_iv_bykova_2018_am05_static_bm3_refit"],
        ),
        ("coesite_v_held_candidate", COESITE_V_DATA, COESITE_V_HELD_CANDIDATE),
    ):
        with path.open(newline="", encoding="utf-8") as stream:
            source_rows = list(csv.DictReader(stream))
        source_pressure = np.array([float(row["pressure_gpa"]) for row in source_rows])
        source_volume = np.array(
            [float(row["volume_a3_conventional_cell"]) for row in source_rows]
        )
        source_residual = bm3_pressure(source_volume, *parameters) - source_pressure
        coesite_checks[name] = {
            "observations": len(source_rows),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(source_residual**2))),
            "max_abs_pressure_residual_gpa": float(np.max(np.abs(source_residual))),
        }

    combined_rows = []
    for path in (COESITE_I_II_DATA, COESITE_II_III_DATA):
        with path.open(newline="", encoding="utf-8") as stream:
            combined_rows.extend(csv.DictReader(stream))
    combined_pressure = np.array([float(row["pressure_gpa"]) for row in combined_rows])
    combined_volume = np.array(
        [float(row["volume_a3_z16_equivalent_cell"]) for row in combined_rows]
    )
    combined_published = bm3_pressure(
        combined_volume, *BM3["coesite_i_iii_bykova_2018_300k_bm3"]
    )
    combined_fit = least_squares(
        lambda parameters: (
            bm3_pressure(combined_volume, *parameters) - combined_pressure
        ),
        x0=BM3["coesite_i_iii_bykova_2018_300k_bm3"],
        bounds=((273.6, 10.3, -10.0), (820.8, 1030.0, 20.0)),
        max_nfev=5000,
    )
    combined_refit_residual = (
        bm3_pressure(combined_volume, *combined_fit.x) - combined_pressure
    )

    coesite_v_parameters = COESITE_V_HELD_CANDIDATE
    coesite_v_anchor_volume = 342.5716

    def v0_from_anchor_pressure(anchor_pressure: float) -> float:
        return float(
            brentq(
                lambda v0: float(
                    bm3_pressure(
                        coesite_v_anchor_volume,
                        v0,
                        coesite_v_parameters[1],
                        coesite_v_parameters[2],
                    )
                    - anchor_pressure
                ),
                400.0,
                460.0,
            )
        )

    v0_rounding_interval = np.array(
        [v0_from_anchor_pressure(56.5), v0_from_anchor_pressure(57.5)]
    )
    central_model = BM3Model(
        V0=coesite_v_parameters[0],
        K0=coesite_v_parameters[1],
        K0_prime=coesite_v_parameters[2],
    )
    pressure_grid = np.linspace(26.0, 64.0, 39)
    volume_grid = np.asarray(central_model.volume(pressure_grid), dtype=float)
    rounding_pressure_shifts = np.vstack(
        [
            bm3_pressure(
                volume_grid,
                v0,
                coesite_v_parameters[1],
                coesite_v_parameters[2],
            )
            - pressure_grid
            for v0 in v0_rounding_interval
        ]
    )
    return {
        "accepted_record_count": len(BM3) + len(BM2) + len(VINET),
        "bm3_zero_pressure_max_abs_gpa": max(
            abs(float(bm3_pressure(v[0], *v))) for v in BM3.values()
        ),
        "bm2_zero_pressure_max_abs_gpa": max(
            abs(float(bm2_pressure(v[0], *v))) for v in BM2.values()
        ),
        "vinet_zero_pressure_max_abs_gpa": max(
            abs(float(vinet_pressure(v[0], *v))) for v in VINET.values()
        ),
        "thermal_bm2_reference_zero_pressure_max_abs_gpa": max(
            abs(float(thermal_bm2_pressure(v[0], 300.0, *v)))
            for v in THERMAL_BM2.values()
        ),
        "compressed_curve_checkpoints_gpa": {
            **{
                key: float(bm3_pressure(0.9 * value[0], *value))
                for key, value in BM3.items()
            },
            **{
                key: float(bm2_pressure(0.9 * value[0], *value))
                for key, value in BM2.items()
            },
            **{
                key: float(vinet_pressure(0.9 * value[0], *value))
                for key, value in VINET.items()
            },
        },
        "thermal_curve_checkpoints_2000k_gpa": {
            key: float(thermal_bm2_pressure(0.9 * value[0], 2000.0, *value))
            for key, value in THERMAL_BM2.items()
        },
        "phase_h_table1": {
            "observations": len(rows),
            "published_bm3_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "published_bm3_max_abs_pressure_residual_gpa": float(
                np.max(np.abs(residual))
            ),
        },
        "bykova_2018_table10": coesite_checks,
        "bykova_2018_combined_coesite_i_ii_iii": {
            "observations": len(combined_rows),
            "phase_counts": {
                phase: sum(row["phase"] == phase for row in combined_rows)
                for phase in ("coesite-I", "coesite-II", "coesite-III")
            },
            "pressure_range_gpa": [
                float(np.min(combined_pressure)),
                float(np.max(combined_pressure)),
            ],
            "published_bm3_pressure_rmse_gpa": float(
                np.sqrt(np.mean((combined_published - combined_pressure) ** 2))
            ),
            "partial_unweighted_refit": {
                "V0": float(combined_fit.x[0]),
                "K0": float(combined_fit.x[1]),
                "K0_prime": float(combined_fit.x[2]),
                "pressure_rmse_gpa": float(
                    np.sqrt(np.mean(combined_refit_residual**2))
                ),
                "max_abs_pressure_residual_gpa": float(
                    np.max(np.abs(combined_refit_residual))
                ),
            },
        },
        "bykova_2018_coesite_v_rounding_sensitivity": {
            "anchor_pressure_interval_gpa": [56.5, 57.5],
            "v0_interval_a3": [float(value) for value in v0_rounding_interval],
            "maximum_abs_pressure_shift_26_64_gpa": float(
                np.max(np.abs(rounding_pressure_shifts))
            ),
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
