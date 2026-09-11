#!/usr/bin/env python3
"""Normalize Ricolleau et al. (2009) Table S1 and audit its KLB-1 EOS fits.

The source table is retained as a lossless, line-ending-normalized TSV
transcription. This script expands parenthetical uncertainties without changing
the printed values, then applies only fit selections that can be stated
explicitly. In particular, the three 300 K observations within the
ferropericlase volume collapse are retained but are not assigned to either
limiting spin branch.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import warnings
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW_TABLE = ROOT / "peritheos/data/datasets/klb1-ricolleau-2009-table-s1.tsv"
NORMALIZED_TABLE = ROOT / "peritheos/data/datasets/klb1-ricolleau-2009-table-s1-pvt.csv"
DEFAULT_JSON = ROOT / "docs/data/ricolleau-2009-klb1-eos-refit.json"

PAIR = re.compile(r"(?<![\w.])(\d+(?:\.\d+)?)\s*\((\d+(?:\.\d+)?)\)")

COLUMNS = (
    "source_line",
    "temperature_k",
    "pressure_gpa",
    "pressure_sigma_gpa",
    "mg_perovskite_volume_a3",
    "mg_perovskite_volume_sigma_a3",
    "ferropericlase_volume_a3",
    "ferropericlase_volume_sigma_a3",
    "ca_perovskite_volume_a3",
    "ca_perovskite_volume_sigma_a3",
    "pressure_medium_volume_a3",
    "pressure_medium_volume_sigma_a3",
    "pressure_medium_volume_2_a3",
    "pressure_medium_volume_2_sigma_a3",
    "gold_volume_a3",
    "gold_volume_sigma_a3",
)


def _sigma(value: str, parenthetical: str) -> str:
    """Expand a crystallographic parenthetical uncertainty losslessly."""
    if "." in parenthetical:
        return parenthetical
    decimal_places = len(value.partition(".")[2])
    return f"{int(parenthetical) * 10.0**-decimal_places:.{decimal_places}f}"


def normalized_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for line_number, line in enumerate(
        RAW_TABLE.read_text(encoding="utf-8").splitlines(), start=1
    ):
        pairs = PAIR.findall(line)
        if len(pairs) == 1 and "and" in line and rows:
            value, uncertainty = pairs[0]
            rows[-1]["pressure_medium_volume_2_a3"] = value
            rows[-1]["pressure_medium_volume_2_sigma_a3"] = _sigma(value, uncertainty)
            continue
        if len(pairs) != 6:
            continue
        first_pair = PAIR.search(line)
        assert first_pair is not None
        temperatures = re.findall(r"\b\d{3,4}\b", line[: first_pair.start()])
        if not temperatures:
            raise ValueError(f"no temperature on source line {line_number}")
        values: list[str] = []
        for value, uncertainty in pairs:
            values.extend((value, _sigma(value, uncertainty)))
        row = dict.fromkeys(COLUMNS, "")
        row.update(
            dict(
                zip(
                    (
                        "pressure_gpa",
                        "pressure_sigma_gpa",
                        "mg_perovskite_volume_a3",
                        "mg_perovskite_volume_sigma_a3",
                        "ferropericlase_volume_a3",
                        "ferropericlase_volume_sigma_a3",
                        "ca_perovskite_volume_a3",
                        "ca_perovskite_volume_sigma_a3",
                        "pressure_medium_volume_a3",
                        "pressure_medium_volume_sigma_a3",
                        "gold_volume_a3",
                        "gold_volume_sigma_a3",
                    ),
                    values,
                )
            )
        )
        row["source_line"] = str(line_number)
        row["temperature_k"] = temperatures[-1]
        rows.append(row)
    if len(rows) != 153:
        raise ValueError(f"expected 153 P-V-T rows, found {len(rows)}")
    return rows


def write_normalized(path: Path = NORMALIZED_TABLE) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized_rows())


def _read_numeric() -> dict[str, np.ndarray]:
    with NORMALIZED_TABLE.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in COLUMNS
        if name
        not in {
            "source_line",
            "pressure_medium_volume_2_a3",
            "pressure_medium_volume_2_sigma_a3",
        }
    }


def bm2_pressure(volume: np.ndarray, v0: float, k0: float) -> np.ndarray:
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5)


def thermal_bm2_pressure(
    volume: np.ndarray,
    temperature: np.ndarray,
    v0: float,
    k0: float,
    dk_dt: float,
    alpha0: float,
    alpha1: float,
) -> np.ndarray:
    tr = 300.0
    reference_volume = v0 * np.exp(
        alpha0 * (temperature - tr) + 0.5 * alpha1 * (temperature**2 - tr**2)
    )
    reference_bulk_modulus = k0 + dk_dt * (temperature - tr)
    return bm2_pressure(volume, reference_volume, reference_bulk_modulus)


def _odr_module():
    # scipy.odr is deprecated in SciPy 1.17 but remains the only bundled solver
    # that exposes the exact errors-in-variables covariance used in this audit.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        from scipy import odr

    return odr


def _static_fit(
    pressure: np.ndarray,
    pressure_sigma: np.ndarray,
    volume: np.ndarray,
    volume_sigma: np.ndarray,
    mask: np.ndarray,
    *,
    fixed_k0: float,
    initial_v0: float,
) -> dict[str, Any]:
    odr = _odr_module()
    model = odr.Model(lambda beta, x: bm2_pressure(x, beta[0], fixed_k0))
    data = odr.RealData(
        volume[mask], pressure[mask], sx=volume_sigma[mask], sy=pressure_sigma[mask]
    )
    result = odr.ODR(data, model, beta0=[initial_v0]).run()
    predicted = bm2_pressure(volume[mask], result.beta[0], fixed_k0)
    return {
        "observations": int(np.sum(mask)),
        "pressure_range_gpa": [
            float(np.min(pressure[mask])),
            float(np.max(pressure[mask])),
        ],
        "fixed_parameters": {"K0": fixed_k0, "K0_prime": 4.0},
        "parameters": {"V0": float(result.beta[0])},
        "standard_errors": {"V0": float(math.sqrt(result.cov_beta[0, 0]))},
        "residual_variance": float(result.res_var),
        "vertical_pressure_rmse_gpa": float(
            np.sqrt(np.mean((predicted - pressure[mask]) ** 2))
        ),
        "solver_stop_reason": list(result.stopreason),
    }


def _thermal_fit(
    pressure: np.ndarray,
    pressure_sigma: np.ndarray,
    volume: np.ndarray,
    volume_sigma: np.ndarray,
    temperature: np.ndarray,
    mask: np.ndarray,
    *,
    fixed_v0: float,
    fixed_k0: float,
    fitted_k0: bool,
    initial: list[float],
) -> dict[str, Any]:
    odr = _odr_module()

    if fitted_k0:
        names = ("K0", "dK_dT", "alpha0", "alpha1")

        def evaluate(beta: np.ndarray, x: np.ndarray) -> np.ndarray:
            return thermal_bm2_pressure(
                x[0], x[1], fixed_v0, beta[0], beta[1], beta[2], beta[3]
            )

    else:
        names = ("dK_dT", "alpha0", "alpha1")

        def evaluate(beta: np.ndarray, x: np.ndarray) -> np.ndarray:
            return thermal_bm2_pressure(
                x[0], x[1], fixed_v0, fixed_k0, beta[0], beta[1], beta[2]
            )

    x = np.vstack((volume[mask], temperature[mask]))
    sx = np.vstack((volume_sigma[mask], np.ones(int(np.sum(mask)))))
    data = odr.RealData(x, pressure[mask], sx=sx, sy=pressure_sigma[mask])
    result = odr.ODR(
        data,
        odr.Model(evaluate),
        beta0=initial,
        ifixx=np.vstack(
            (
                np.ones(int(np.sum(mask)), dtype=int),
                np.zeros(int(np.sum(mask)), dtype=int),
            )
        ),
    ).run()
    predicted = evaluate(result.beta, x)
    return {
        "observations": int(np.sum(mask)),
        "pressure_range_gpa": [
            float(np.min(pressure[mask])),
            float(np.max(pressure[mask])),
        ],
        "temperature_range_k": [
            float(np.min(temperature[mask])),
            float(np.max(temperature[mask])),
        ],
        "fixed_parameters": {"V0": fixed_v0, "K0_prime": 4.0}
        | ({} if fitted_k0 else {"K0": fixed_k0}),
        "parameters": dict(zip(names, map(float, result.beta))),
        "standard_errors": {
            name: float(math.sqrt(result.cov_beta[index, index]))
            for index, name in enumerate(names)
        },
        "residual_variance": float(result.res_var),
        "vertical_pressure_rmse_gpa": float(
            np.sqrt(np.mean((predicted - pressure[mask]) ** 2))
        ),
        "solver_stop_reason": list(result.stopreason),
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _gold_pressure_recalculation(
    pressure: np.ndarray, gold_volume: np.ndarray, temperature: np.ndarray
) -> dict[str, Any]:
    from peritheos import Material, get_material_document

    record = Material.from_eosmat(
        get_material_document("gold"),
        record_identifiers=["gold_fei_2007_vinet_2"],
    ).eos_records[0]
    recalculated = np.asarray(record.pressure(gold_volume, temperature), dtype=float)
    residual = recalculated - pressure
    return {
        "reference_eos_record": "gold_fei_2007_vinet_2",
        "observations": int(pressure.size),
        "recalculated_minus_published_mean_gpa": float(np.mean(residual)),
        "recalculated_minus_published_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "recalculated_minus_published_max_abs_gpa": float(np.max(np.abs(residual))),
        "qualification": (
            "All raw inputs needed for re-reduction are present. The executable "
            "published Fei et al. parameterization closely tracks but does not "
            "exactly reproduce the printed pressures; no undocumented correction "
            "or alternate implementation is inferred."
        ),
    }


def reproduce() -> dict[str, Any]:
    data = _read_numeric()
    pressure = data["pressure_gpa"]
    pressure_sigma = data["pressure_sigma_gpa"]
    temperature = data["temperature_k"]
    room_temperature = temperature == 300.0
    heated = temperature > 300.0
    high_spin = room_temperature & (pressure < 50.0)
    crossover = room_temperature & (pressure >= 50.0) & (pressure < 70.0)
    low_spin = room_temperature & (pressure >= 70.0)

    fp_static_high = _static_fit(
        pressure,
        pressure_sigma,
        data["ferropericlase_volume_a3"],
        data["ferropericlase_volume_sigma_a3"],
        high_spin,
        fixed_k0=158.0,
        initial_v0=76.44,
    )
    fp_static_low = _static_fit(
        pressure,
        pressure_sigma,
        data["ferropericlase_volume_a3"],
        data["ferropericlase_volume_sigma_a3"],
        low_spin,
        fixed_k0=170.0,
        initial_v0=74.04,
    )
    mg = _thermal_fit(
        pressure,
        pressure_sigma,
        data["mg_perovskite_volume_a3"],
        data["mg_perovskite_volume_sigma_a3"],
        temperature,
        np.ones(pressure.shape, dtype=bool),
        fixed_v0=164.0,
        fixed_k0=245.0,
        fitted_k0=True,
        initial=[245.0, -0.036, 3.19e-5, 8.8e-9],
    )
    ca = _thermal_fit(
        pressure,
        pressure_sigma,
        data["ca_perovskite_volume_a3"],
        data["ca_perovskite_volume_sigma_a3"],
        temperature,
        np.ones(pressure.shape, dtype=bool),
        fixed_v0=45.60,
        fixed_k0=244.0,
        fitted_k0=True,
        initial=[244.0, -0.035, 3.06e-5, 8.7e-9],
    )
    fp_thermal = _thermal_fit(
        pressure,
        pressure_sigma,
        data["ferropericlase_volume_a3"],
        data["ferropericlase_volume_sigma_a3"],
        temperature,
        heated,
        fixed_v0=76.44,
        fixed_k0=158.0,
        fitted_k0=False,
        initial=[-0.034, 2.20e-5, 3.61e-8],
    )

    return {
        "format": "peritheos.ricolleau-2009-klb1-eos-audit",
        "format_version": 1,
        "generated_with": "scripts/audit_ricolleau_2009_klb1_eos.py",
        "source": {
            "doi": "10.1029/2008GL036759",
            "location": "Official Supporting Information Table S1",
            "download_name": "grl25568-sup-0002-ts01.txt",
            "raw_table": RAW_TABLE.relative_to(ROOT).as_posix(),
            "raw_table_sha256": _sha256(RAW_TABLE),
            "normalized_table": NORMALIZED_TABLE.relative_to(ROOT).as_posix(),
            "normalized_table_sha256": _sha256(NORMALIZED_TABLE),
        },
        "observations": {
            "total": int(pressure.size),
            "room_temperature": int(np.sum(room_temperature)),
            "heated": int(np.sum(heated)),
            "secondary_pressure_medium_volumes": sum(
                bool(row["pressure_medium_volume_2_a3"]) for row in normalized_rows()
            ),
        },
        "protocol": {
            "equation": "BM2 (K0_prime=4) with K(T)=K0+dK_dT*(T-300 K)",
            "reference_volume_law": (
                "V0(T)=V0*exp[alpha0*(T-300 K)+0.5*alpha1*(T^2-300^2 K^2)]"
            ),
            "objective": (
                "orthogonal distance regression using printed 1-sigma pressure "
                "and phase-volume uncertainties; temperature held fixed"
            ),
            "uncertainty": (
                "unscaled parameter covariance from the printed observation "
                "sigmas; no residual-variance rescaling"
            ),
            "qualification": (
                "The paper specifies BM2, fixed coefficients, and alpha(T), but "
                "does not state its regression objective, weights, treatment of "
                "temperature uncertainty, complete row selection, or covariance "
                "scaling. This is a source-constrained reconstruction, not an "
                "assertion of the authors' undocumented numerical implementation."
            ),
        },
        "pressure_recalculation": _gold_pressure_recalculation(
            pressure, data["gold_volume_a3"], temperature
        ),
        "spin_selection": {
            "high_spin_limiting_branch": {
                "rule": "300 K and reported pressure below 50 GPa",
                "observations": int(np.sum(high_spin)),
                "basis": (
                    "The article places the room-temperature volume contraction "
                    "at about 50 GPa; these seven rows lie wholly below it."
                ),
            },
            "unassigned_crossover": {
                "rule": "300 K and reported pressure from 50 to below 70 GPa",
                "observations": int(np.sum(crossover)),
                "pressures_gpa": sorted(map(float, pressure[crossover])),
                "basis": (
                    "These three printed rows occupy the observed volume-collapse "
                    "interval and are retained without inventing a spin assignment."
                ),
            },
            "low_spin_limiting_branch": {
                "rule": "300 K and reported pressure at or above 70 GPa",
                "observations": int(np.sum(low_spin)),
                "basis": (
                    "The next seven 300 K rows form the separated high-pressure "
                    "limiting branch. The 70 GPa boundary labels the observed "
                    "sampling gap; it is not claimed as a published transition cutoff."
                ),
            },
        },
        "record_refits": {
            "klb1_mg_perovskite_ricolleau_2009_bm2_alphakt": {
                "stage": "all 153 Table S1 P-V-T rows",
                "fit": mg,
            },
            "klb1_ca_perovskite_ricolleau_2009_bm2_alphakt": {
                "stage": "all 153 Table S1 P-V-T rows",
                "fit": ca,
            },
            "klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt": {
                "stage": (
                    "V0 from the seven limiting high-spin 300 K rows; thermal "
                    "coefficients from all 136 heated rows with published V0 and K0 fixed"
                ),
                "static_fit": fp_static_high,
                "thermal_fit": fp_thermal,
            },
            "klb1_ferropericlase_ricolleau_2009_low_spin_300k_bm2": {
                "stage": "seven limiting low-spin 300 K rows",
                "fit": fp_static_low,
            },
        },
        "conclusion": {
            "static_spin_branches": (
                "The errors-in-variables reconstruction reproduces both published "
                "V0 values and their 0.02 A^3 uncertainties to printed rounding."
            ),
            "thermal_fits": (
                "The source-constrained reconstruction is uncertainty-compatible "
                "with every published thermal coefficient except ferropericlase "
                "alpha0, which is just outside combined 2-sigma but remains within "
                "the ledger's 30% thermal-coefficient similarity criterion."
            ),
            "irreducible_blocker": (
                "Exact numerical parity for the authors' thermal regressions cannot "
                "be claimed without their row exclusions (if any), weighting or "
                "objective, temperature-error treatment, and covariance convention."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write CSV and JSON")
    parser.add_argument("--output", type=Path, default=DEFAULT_JSON)
    args = parser.parse_args()
    if args.write:
        write_normalized()
        result = reproduce()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
