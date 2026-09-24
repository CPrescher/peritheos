"""Audit the original argon polynomial; never treat generated curves as data."""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/data/argon-anderson-swenson-1975-reproduction.json"
FACTOR = 4e24 / 6.02214076e23


def rows(name):
    with (ROOT / "peritheos/data/datasets" / name).open() as stream:
        return list(csv.DictReader(stream))


def source_pressure(v, row):
    """Direct source equation in original molar/kbar units, converted to GPa."""
    return (
        sum(
            float(row[f"a{n}_kbar_cm{3 * n}_mol{n}"]) / np.asarray(v) ** n
            for n in (3, 5, 7, 9)
        )
        / 10
    )


def rms(values):
    return float(np.sqrt(np.mean(np.square(values))))


def reproduce():
    table = rows("argon-anderson-swenson-1975-table1.csv")
    raw = rows("argon-anderson-1972-appendix-a.csv")
    material = Material.from_eosmat(get_material_document("argon_fcc"))
    records = [
        r for r in material.eos_records if "anderson_swenson_1975" in r.identifier
    ]
    result = {
        "status": "equation_reproduced",
        "raw_table_cells": len(raw),
        "usable_observations": sum(r["usable"] == "1" for r in raw),
        "independent_refit_status": "not_reproduced",
        "isotherms": [],
    }
    for row, record in zip(table, records):
        v0 = float(row["v0_cm3_mol"])
        vs = np.linspace(
            brentq(lambda v: source_pressure(v, row) - 2, 0.70 * v0, v0), v0, 41
        )
        predicted = record.pressure(vs * FACTOR)
        native_root = record.volume(predicted) / FACTOR
        zero = brentq(lambda v: source_pressure(v, row), 0.99 * v0, 1.01 * v0)
        bulk0 = record.eos.bulk_modulus(zero * FACTOR) * 10
        t = float(row["temperature_k"])
        # Diagnostic at nearest nominal T within 1.1 K. No claim this is an
        # exact common-temperature reduction, and 67.9 K is never interpolated.
        selected = [
            r
            for r in raw
            if r["usable"] == "1" and abs(float(r["temperature_k"]) - t) <= 1.100001
        ]
        raw_residual, corrected_residual, large_holder_residual = [], [], []
        for r in selected:
            pressure = float(r["pressure_kbar"]) / 10
            expected = brentq(
                lambda v: source_pressure(v, row) - pressure, 0.70 * v0, 1.01 * v0
            )
            v = float(r["reported_volume_cm3_mol"])
            correction = (
                {1: 0.253, 2: 0.449}[int(r["sample"])]
                if float(r["holder_diameter_inch"]) == 0.25
                else 0
            )
            raw_residual.append(v - expected)
            corrected_residual.append(v + correction - expected)
            if float(r["holder_diameter_inch"]) != 0.25:
                large_holder_residual.append(v - expected)
        result["isotherms"].append(
            {
                "record_id": record.identifier,
                "temperature_k": t,
                "python_max_difference_gpa": float(
                    max(abs(predicted - source_pressure(vs, row)))
                ),
                "inverse_max_difference_cm3_mol": float(max(abs(native_root - vs))),
                "pressure_at_rounded_v0_gpa": float(source_pressure(v0, row)),
                "zero_pressure_root_cm3_mol": zero,
                "calculated_bulk_at_zero_kbar": bulk0,
                "published_bulk_at_zero_kbar": float(row["k0_kbar"]),
                "published_volume_rms_cm3_mol": float(row["volume_fit_rms_cm3_mol"]),
                "comparison_point_count": len(selected),
                "unadjusted_volume_rmse_cm3_mol": rms(raw_residual),
                "literal_table11_corrected_volume_rmse_cm3_mol": rms(
                    corrected_residual
                ),
                "larger_holder_only_volume_rmse_cm3_mol": rms(large_holder_residual),
            }
        )
    vs = np.linspace(17.37, 22.56, 51)
    reduced = 2.86 * (
        -0.553408 * (22.56 / vs) ** 3
        + 0.6068095 * (22.56 / vs) ** 5
        - 0.0534018 * (22.56 / vs) ** 7
    )
    result["eq9_vs_table1_max_difference_gpa"] = float(
        max(abs(reduced - source_pressure(vs, table[0])))
    )
    result["diagnostic_note"] = (
        "Observed-volume residuals use nearest nominal T within 1.1 K; actual T preserved. "
        "The source used isobaric reductions to common T, not reconstructed here. Table11 corrections applied literally in the diagnostic only; "
        "small-holder discrepancy is unresolved, so no original-fit reproduction or coefficient uncertainties are claimed."
    )
    return result


if __name__ == "__main__":
    REPORT.write_text(json.dumps(reproduce(), indent=2) + "\n")
    print(REPORT)
