#!/usr/bin/env python3
"""Conditional fcc-Fe pressure recovery from Funamori/Boehler/Basinski.

BM3, interpolation at 1400 K, and the use of Boehler's mean-expansion
coefficient with a 300 K volume are explicit audit choices, not recovered
Seagle code. No coefficient is adjusted to Seagle or Campbell targets.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import brentq

if __package__:
    from .reconstruct_seagle_2008_pressures import hcp_pressure
else:
    from reconstruct_seagle_2008_pressures import hcp_pressure

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/data/seagle-2008-fcc-reconstruction.json"
OUTPUT = ROOT / "docs/data/seagle-2008-completed-pressure-audit.csv"
CELL_TO_MOLAR = Avogadro * 1e-24 / 4
# Basinski Table 2 brackets 1400 K with these lattice spacings, in kX.
# The historical 1.00202 conversion is documented by NBS Monograph 25-3 p3.
A_1400_KX = 3.6524 + (1400 - 1347) / (1457 - 1347) * (3.6622 - 3.6524)
V0_1400 = (A_1400_KX * 1.00202) ** 3 * CELL_TO_MOLAR


def fcc_pressure(cell_volume, temperature, delta=6.5):
    """GPa from a four-atom fcc cell in A3; Boehler mean-expansion convention."""
    v, t = np.broadcast_arrays(
        np.asarray(cell_volume, float), np.asarray(temperature, float)
    )

    def scalar(cell, temp):
        measured = cell * CELL_TO_MOLAR

        def thermal_volume(v300, target):
            alpha_mean = 7.70e-5 * (v300 / 6.835) ** delta
            return v300 * (1 + alpha_mean * (target - 300))

        if temp < 300:
            raise ValueError("FCC reconstruction is only supported at T >= 300 K")
        v300 = brentq(
            lambda w: thermal_volume(w, temp) - measured,
            measured * 0.25,
            measured,
            xtol=1e-13,
        )
        ratio = thermal_volume(v300, 1400) / V0_1400
        return (
            1.5
            * 120
            * (ratio ** (-7 / 3) - ratio ** (-5 / 3))
            * (1 + 0.75 * (5 - 4) * (ratio ** (-2 / 3) - 1))
        )

    return np.array([scalar(vi, ti) for vi, ti in zip(v.flat, t.flat)]).reshape(v.shape)


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def calibration_checks():
    checks = []
    for row in read_csv(ROOT / "docs/data/seagle-2008-fcc-pressure-checkpoints.csv"):
        v, t, p, sp = [
            float(row[k])
            for k in [
                "cell_volume_a3",
                "temperature_k",
                "published_pressure_gpa",
                "published_pressure_uncertainty_gpa",
            ]
        ]
        predicted = float(fcc_pressure(v, t))
        checks.append(
            {
                **row,
                "reconstructed_pressure_gpa": predicted,
                "difference_gpa": predicted - p,
                "within_reported_uncertainty": abs(predicted - p) <= sp,
            }
        )
    funamori = []
    for row in read_csv(ROOT / "docs/data/funamori-1996-fcc-table1.csv"):
        v, t, p = [
            float(row[k]) for k in ["cell_volume_a3", "temperature_k", "pressure_gpa"]
        ]
        funamori.append(
            {
                **row,
                "reconstructed_pressure_gpa": float(fcc_pressure(v, t)),
                "difference_gpa": float(fcc_pressure(v, t)) - p,
            }
        )
    residuals = np.array([r["difference_gpa"] for r in checks])
    return {
        "reference_volume_1400_cm3_mol": V0_1400,
        "lattice_parameter_1400_kx": A_1400_KX,
        "seagle_checkpoints": checks,
        "seagle_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "seagle_mean_difference_gpa": float(np.mean(residuals)),
        "funamori_table1": funamori,
        "qualification": "Conditional BM3/mean-expansion reconstruction; all eight Seagle comparisons use printed pressure errors. Funamori reports no row-wise pressure error here, so its differences are not classified as within author uncertainty. Agreement does not identify the exact historical interpolation or thermal reduction.",
    }


def completed_rows():
    rows = read_csv(
        ROOT
        / "peritheos/data/datasets/feo-seagle-2008-supplement-volume-temperature.csv"
    )
    out = []
    for row in rows:
        t = float(row["temperature_k"])
        item = {
            "source_row": row["source_row"],
            "iron_phase": row["iron_phase"],
            "temperature_k": t,
            "reconstructed_pressure_gpa": "",
            "measurement_only_pressure_error_gpa": "",
            "fcc_delta_6_pressure_gpa": "",
            "fcc_delta_7_pressure_gpa": "",
            "status": "missing_iron_volume",
        }
        if row["iron_unit_cell_volume_a3"]:
            v = float(row["iron_unit_cell_volume_a3"])
            model = fcc_pressure if row["iron_phase"] == "fcc" else hcp_pressure
            item["status"] = (
                "conditional_fcc_funamori_boehler_basinski"
                if model is fcc_pressure
                else "hcp_seagle_2006_density_text"
            )
            item["reconstructed_pressure_gpa"] = float(model(v, t))
            sv = row["iron_unit_cell_volume_uncertainty_a3"]
            if sv:
                dpdv = float(model(v + 1e-4, t) - model(v - 1e-4, t)) / 2e-4
                dpdt = float(model(v, t + 0.01) - model(v, t - 0.01)) / 0.02
                item["measurement_only_pressure_error_gpa"] = float(
                    np.hypot(dpdv * float(sv), dpdt * 150)
                )
            if model is fcc_pressure:
                item["fcc_delta_6_pressure_gpa"] = float(model(v, t, delta=6))
                item["fcc_delta_7_pressure_gpa"] = float(model(v, t, delta=7))
        out.append(item)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rows = completed_rows()
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    outputs = {
        OUTPUT: stream.getvalue(),
        REPORT: json.dumps(calibration_checks(), indent=2, sort_keys=True) + "\n",
    }
    for path, text in outputs.items():
        if args.check:
            if not path.exists() or path.read_text() != text:
                raise SystemExit(f"Stale fcc reconstruction: {path.name}")
        else:
            path.write_text(text)
    print("79 reconstructed pressures; two missing iron volumes")


if __name__ == "__main__":
    main()
