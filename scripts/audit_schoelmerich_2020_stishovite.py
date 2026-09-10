#!/usr/bin/env python3
"""Audit Schoelmerich et al. (2020) Table 1 shock states.

The source-owned shock observations are retained as a transcription. The
reported 300 K EOS is deliberately not cataloged because its corrected fit
states and complete thermal-reduction protocol were not published numerically.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "peritheos/data/datasets/stishovite-schoelmerich-2020-table1-shock.csv"
OUTPUT = ROOT / "docs/data/schoelmerich-2020-stishovite-audit.json"

RHO0_G_CM3 = 4.30
MOLAR_MASS_G_MOL = 60.083
FORMULA_UNITS_PER_CELL = 2.0


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _finite(rows: list[dict[str, str]], name: str) -> np.ndarray:
    return np.asarray([float(row[name]) for row in rows if row[name]], dtype=float)


def audit() -> dict[str, Any]:
    table_rows = _rows(TABLE)
    shock_rows = [row for row in table_rows if row["particle_velocity_km_s"]]

    particle = _finite(shock_rows, "particle_velocity_km_s")
    shock = _finite(shock_rows, "shock_velocity_km_s")
    density = _finite(shock_rows, "density_g_cm3")
    pressure = _finite(shock_rows, "hugoniot_pressure_gpa")
    volume = _finite(shock_rows, "volume_a3_conventional_cell")
    reported_energy = _finite(shock_rows, "reported_internal_energy_change")

    slope, intercept = np.polyfit(particle, shock, 1)
    rh_pressure = RHO0_G_CM3 * particle * shock
    density_from_volume = (
        MOLAR_MASS_G_MOL * FORMULA_UNITS_PER_CELL / (0.602214076 * volume)
    )
    # GPa * cm3/g is numerically MJ/kg.
    rh_energy_mj_kg = 0.5 * pressure * (1.0 / RHO0_G_CM3 - 1.0 / density)
    rh_energy_kj_mol = rh_energy_mj_kg * MOLAR_MASS_G_MOL

    return {
        "format": "peritheos.schoelmerich-2020-stishovite-audit",
        "format_version": 1,
        "audit_date": "2026-09-08",
        "source": {
            "doi": "10.1038/s41598-020-66340-y",
            "article_url": "https://www.nature.com/articles/s41598-020-66340-y",
            "supplement_url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-020-66340-y/MediaObjects/41598_2020_66340_MOESM1_ESM.pdf",
            "license": "CC BY 4.0",
            "article_pdf_sha256": "b138b100e1d5fa79eb83bd3ba84fc87517fdd285e7aefc28b4f9fc6f60f8a836",
            "supplement_pdf_sha256": "2b1b43493480e9a92aeaf25bd43738f73f796c971660ff304ac9de4f966800d3",
        },
        "reported_table": {
            "dataset": TABLE.name,
            "sha256": _sha256(TABLE),
            "rows": len(table_rows),
            "finite_shock_rows": len(shock_rows),
            "row_origin_counts": {
                value: sum(row["state_origin"] == value for row in table_rows)
                for value in sorted({row["state_origin"] for row in table_rows})
            },
            "rounded_us_up_ordinary_least_squares": {
                "slope": float(slope),
                "intercept_km_s": float(intercept),
            },
            "supplement_unrounded_us_up_relation": {
                "slope": 1.199576,
                "intercept_km_s": 9.166700,
            },
            "rankine_hugoniot_pressure_check": {
                "rho0_g_cm3": RHO0_G_CM3,
                "calculated_pressure_gpa": rh_pressure.tolist(),
                "reported_minus_calculated_gpa": (pressure - rh_pressure).tolist(),
                "max_abs_difference_gpa": float(np.max(np.abs(pressure - rh_pressure))),
            },
            "density_from_cell_volume_check": {
                "calculated_density_g_cm3": density_from_volume.tolist(),
                "reported_minus_calculated_g_cm3": (
                    density - density_from_volume
                ).tolist(),
                "max_abs_difference_g_cm3": float(
                    np.max(np.abs(density - density_from_volume))
                ),
            },
            "internal_energy_unit_audit": {
                "printed_table_heading": "E-E0 (kJ/mol)",
                "rankine_hugoniot_value_if_mj_kg": rh_energy_mj_kg.tolist(),
                "rankine_hugoniot_value_if_kj_mol": rh_energy_kj_mol.tolist(),
                "reported_value": reported_energy.tolist(),
                "rmse_treating_reported_values_as_mj_kg": float(
                    np.sqrt(np.mean((reported_energy - rh_energy_mj_kg) ** 2))
                ),
                "rmse_treating_reported_values_as_kj_mol": float(
                    np.sqrt(np.mean((reported_energy - rh_energy_kj_mol) ** 2))
                ),
                "finding": (
                    "The printed kJ/mol heading is dimensionally inconsistent. "
                    "The numerical entries reproduce the Rankine-Hugoniot energy "
                    "in MJ/kg (equivalently kJ/g); the CSV preserves the literal "
                    "reported values and labels the source unit defect in metadata."
                ),
            },
        },
        "catalog_decision": {
            "outcome": "direct_refit_unavailable",
            "production_eos_record": False,
            "published_parameters": {
                "V0_a3": 46.5,
                "K0_gpa": 307.0,
                "K0_error_gpa": 4.0,
                "K0_prime": 4.66,
                "K0_prime_error": 0.15,
            },
            "reason": (
                "The five shock-corrected 300 K states are shown only graphically. "
                "The source also omits the exact correction inputs, corrected-state "
                "uncertainties and covariance, EosFit objective and weights, and an "
                "explicit LCLS-233 selection. Figure digitization cannot provide an "
                "independent reproduction of coefficients plotted in the same figure."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = audit()
    rendered = json.dumps(result, indent=2) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"{OUTPUT} is stale; rerun without --check")
        return
    OUTPUT.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
