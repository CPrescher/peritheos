#!/usr/bin/env python3
"""Reproduce and validate the two Dewaele-anchored diamond composites.

This is deliberately a reconstruction, not a coefficient refit.  It checks the
published thermal branches against source observations and verifies the complete
Peritheos composition algebra, including caloric increments.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.constants import Avogadro, electron_volt

from peritheos.materials import (
    DIAMOND_BENEDICT_2014,
    DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED,
    DIAMOND_CORREA_2008,
    DIAMOND_CORREA_2008_DEWAELE_ANCHORED,
)

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "peritheos" / "data" / "datasets"
DEFAULT_OUTPUT = ROOT / "docs" / "data" / "diamond-thermal-composite-reconstruction.json"
ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR = Avogadro * 1.0e-25
EV_PER_ATOM_TO_J_PER_MOL = electron_volt * Avogadro


def _rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [
            {
                key: float(value)
                for key, value in row.items()
                if key not in {"series"} and value != ""
            }
            | ({"series": row["series"]} if "series" in row else {})
            for row in csv.DictReader(stream)
        ]


def _rmse(values: list[float]) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def _pressure_diagnostics(eos: Any, rows: list[dict[str, float]]) -> dict[str, Any]:
    residuals = []
    grouped: dict[float, list[tuple[dict[str, float], float]]] = defaultdict(list)
    for row in rows:
        volume = row["volume_a3_per_atom"] * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
        prediction = float(eos.pressure(volume, row["temperature_k"]))
        residuals.append(prediction - row["pressure_gpa"])
        grouped[row["volume_a3_per_atom"]].append((row, prediction))

    increment_residuals = []
    for states in grouped.values():
        states.sort(key=lambda item: item[0]["temperature_k"])
        reference_row, reference_prediction = states[0]
        for row, prediction in states[1:]:
            observed = row["pressure_gpa"] - reference_row["pressure_gpa"]
            predicted = prediction - reference_prediction
            increment_residuals.append(predicted - observed)

    return {
        "observations": len(rows),
        "isochores": len(grouped),
        "absolute_pressure_rmse_gpa": _rmse(residuals),
        "absolute_pressure_bias_gpa": float(np.mean(residuals)),
        "absolute_pressure_max_abs_residual_gpa": float(np.max(np.abs(residuals))),
        "isochoric_pressure_increment_comparisons": len(increment_residuals),
        "pressure_increment_rmse_gpa": _rmse(increment_residuals),
        "pressure_increment_bias_gpa": float(np.mean(increment_residuals)),
        "pressure_increment_max_abs_residual_gpa": float(
            np.max(np.abs(increment_residuals))
        ),
    }


def _energy_increment_diagnostics(
    eos: Any, rows: list[dict[str, float]]
) -> dict[str, Any]:
    grouped: dict[float, list[tuple[dict[str, float], float]]] = defaultdict(list)
    for row in rows:
        volume = row["volume_a3_per_atom"] * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
        prediction = float(eos.internal_energy(volume, row["temperature_k"]))
        prediction /= EV_PER_ATOM_TO_J_PER_MOL
        grouped[row["volume_a3_per_atom"]].append((row, prediction))

    residuals = []
    for states in grouped.values():
        states.sort(key=lambda item: item[0]["temperature_k"])
        reference_row, reference_prediction = states[0]
        for row, prediction in states[1:]:
            observed = (
                row["internal_energy_ev_per_atom"]
                - reference_row["internal_energy_ev_per_atom"]
            )
            predicted = prediction - reference_prediction
            residuals.append(predicted - observed)
    return {
        "internal_energy_increment_comparisons": len(residuals),
        "internal_energy_increment_rmse_ev_per_atom": _rmse(residuals),
        "internal_energy_increment_bias_ev_per_atom": float(np.mean(residuals)),
        "internal_energy_increment_max_abs_residual_ev_per_atom": float(
            np.max(np.abs(residuals))
        ),
        "qualification": (
            "Fixed-volume increments remove the arbitrary electronic-structure "
            "energy zero; absolute energies from different calculations are not compared."
        ),
    }


def _identity_diagnostics(
    source: Any, anchored: Any, rows: list[dict[str, float]]
) -> dict[str, Any]:
    pressure_errors = []
    energy_increment_errors = []
    reference_errors = []
    for row in rows:
        volume = row["volume_a3_per_atom"] * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
        temperature = row["temperature_k"]
        expected_pressure = (
            anchored.rt_eos.pressure(volume)
            + source.pressure(volume, temperature)
            - source.pressure(volume, anchored.Tr)
        )
        pressure_errors.append(anchored.pressure(volume, temperature) - expected_pressure)
        reference_errors.append(
            anchored.pressure(volume, anchored.Tr) - anchored.rt_eos.pressure(volume)
        )
        source_increment = source.internal_energy(volume, temperature) - source.internal_energy(
            volume, anchored.Tr
        )
        anchored_increment = anchored.internal_energy(
            volume, temperature
        ) - anchored.internal_energy(volume, anchored.Tr)
        energy_increment_errors.append(anchored_increment - source_increment)
    return {
        "states_checked": len(rows),
        "max_abs_composition_pressure_error_gpa": float(
            np.max(np.abs(pressure_errors))
        ),
        "max_abs_reference_isotherm_error_gpa": float(
            np.max(np.abs(reference_errors))
        ),
        "max_abs_thermal_energy_increment_error_j_per_mol": float(
            np.max(np.abs(energy_increment_errors))
        ),
    }


def reproduce() -> dict[str, Any]:
    correa_rows = _rows(
        DATA_ROOT / "diamond-correa-2008-figure8-dft-md-vector-digitized.csv"
    )
    benedict_rows = _rows(
        DATA_ROOT / "diamond-benedict-2014-supplement-solid-dft-md.csv"
    )
    correa_source = DIAMOND_CORREA_2008.eos
    benedict_source = DIAMOND_BENEDICT_2014.eos
    correa_anchored = DIAMOND_CORREA_2008_DEWAELE_ANCHORED.eos
    benedict_anchored = DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED.eos
    return {
        "format": "peritheos.diamond-thermal-composite-reconstruction",
        "format_version": 1,
        "classification": "source-equation reconstruction, not independent refit",
        "coefficient_optimization_performed": False,
        "composition": (
            "P(V,T)=P_Dewaele(V,298 K)+P_theory(V,T)-P_theory(V,298 K)"
        ),
        "correa_2008": {
            "source_dataset": "diamond_correa_2008_figure8_dft_md_vector_digitized",
            "source_model_pressure_validation": _pressure_diagnostics(
                correa_source, correa_rows
            ),
            "complete_composite_identity": _identity_diagnostics(
                correa_source, correa_anchored, correa_rows
            ),
            "qualification": (
                "Figure 8 supplies direct pressure checkpoints but not the upstream "
                "cold-energy grid, phonon DOS/moments, regression weights, or covariance."
            ),
        },
        "benedict_2014": {
            "source_dataset": "diamond_benedict_2014_supplement_solid_dft_md",
            "source_model_pressure_validation": _pressure_diagnostics(
                benedict_source, benedict_rows
            ),
            "source_model_caloric_validation": _energy_increment_diagnostics(
                benedict_source, benedict_rows
            ),
            "complete_composite_identity": _identity_diagnostics(
                benedict_source, benedict_anchored, benedict_rows
            ),
            "qualification": (
                "The exact supplementary solid DFT-MD table is a downstream model "
                "validation grid, not the cold-curve and phonon fitting inputs."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(reproduce(), indent=2) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale generated file: {args.output}")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
