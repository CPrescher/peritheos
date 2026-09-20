#!/usr/bin/env python3
"""Audit-only recovery of the hcp-Fe pressures used by Seagle (2008).

Seagle (2006), equations 1-4, Table 2 and paragraph 18 provide conflicting
reference volumes. Both are evaluated; no coefficient is fitted here.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, R
from scipy.special import roots_legendre

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINTS = ROOT / "docs/data/seagle-hcp-pressure-checkpoints.csv"
OUTPUT = ROOT / "docs/data/seagle-2008-reconstructed-pressures.csv"
NODES, WEIGHTS = roots_legendre(48)
NODES = (NODES + 1) / 2
WEIGHTS = WEIGHTS / 2


def hcp_pressure(cell_volume, temperature, reference="density_text"):
    """GPa from a two-atom hcp cell (A3), with explicit source corrections.

    beta=9.10e-8 kJ/g/K2 is the dimensionally consistent interpretation of
    the misprinted electronic coefficient. The printed equation multiplies
    electronic as well as vibrational energy by gamma(V).
    """
    if reference not in {"density_text", "table_2"}:
        raise ValueError("Unknown Seagle reference-volume convention")
    v0 = 55.845 / 8.30 if reference == "density_text" else 6.687
    v = np.asarray(cell_volume, dtype=float) * Avogadro * 1e-24 / 2
    t = np.asarray(temperature, dtype=float)
    ratio = v / v0
    gamma = 2.4 * ratio**1.2
    theta = 380 * np.exp((2.4 - gamma) / 1.2)

    def energy(temp):
        z = (theta / temp)[..., None] * NODES
        return 9 * R * temp * np.sum(WEIGHTS * NODES**2 * z / np.expm1(z), axis=-1)

    cold = (
        1.5
        * 164.8
        * (ratio ** (-7 / 3) - ratio ** (-5 / 3))
        * (1 + 0.75 * (5.33 - 4) * (ratio ** (-2 / 3) - 1))
    )
    vibrational = gamma * (energy(t) - energy(300.0)) / v / 1000
    electronic = gamma * (55.845 / v) * 9.10e-8 / 2 * ratio**1.34 * (t**2 - 300**2)
    return cold + vibrational + electronic


def calibration_checks():
    """Compare against independent printed outputs, not Campbell fit targets."""
    with CHECKPOINTS.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    out = {}
    for source in sorted({row["source"] for row in rows}):
        selected = [r for r in rows if r["source"] == source]
        v, t, p, sp = np.array(
            [
                [float(r[k]) for r in selected]
                for k in [
                    "iron_cell_volume_a3",
                    "temperature_k",
                    "published_pressure_gpa",
                    "published_pressure_uncertainty_gpa",
                ]
            ]
        )
        out[source] = {"observations": len(selected)}
        for reference in ["density_text", "table_2"]:
            residual = hcp_pressure(v, t, reference) - p
            out[source][reference] = {
                "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
                "max_absolute_difference_gpa": float(np.max(abs(residual))),
                "within_reported_pressure_uncertainty": int(
                    np.sum(abs(residual) <= sp)
                ),
            }
    return out


def reconstructed_rows(assignments):
    """Return all 81 audit rows; unavailable individual pressures stay empty."""
    source = (
        ROOT
        / "peritheos/data/datasets/feo-seagle-2008-supplement-volume-temperature.csv"
    )
    with source.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    nominal = {r["source_row"]: r["nominal_pressure_gpa"] for r in assignments}
    result = []
    for row in rows:
        index = int(row["source_row"])
        t = float(row["temperature_k"])
        out = {
            "source_row": index,
            "iron_phase": row["iron_phase"],
            "temperature_k": t,
            "nominal_series_pressure_gpa": nominal[index],
            "reconstructed_pressure_gpa": "",
            "table_2_reference_pressure_gpa": "",
            "propagated_measurement_uncertainty_gpa": "",
            "status": "fcc_calibration_incomplete"
            if row["iron_phase"] == "fcc"
            else "missing_iron_volume",
        }
        if row["iron_phase"] == "hcp" and row["iron_unit_cell_volume_a3"]:
            v = float(row["iron_unit_cell_volume_a3"])
            out["status"] = "hcp_seagle_2006_density_text"
            out["reconstructed_pressure_gpa"] = float(hcp_pressure(v, t))
            out["table_2_reference_pressure_gpa"] = float(hcp_pressure(v, t, "table_2"))
            if row["iron_unit_cell_volume_uncertainty_a3"]:
                sv = float(row["iron_unit_cell_volume_uncertainty_a3"])
                dpdv = (hcp_pressure(v + 1e-4, t) - hcp_pressure(v - 1e-4, t)) / 2e-4
                dpdt = (hcp_pressure(v, t + 0.01) - hcp_pressure(v, t - 0.01)) / 0.02
                # Independent input errors are an audit assumption; excludes EOS systematics.
                out["propagated_measurement_uncertainty_gpa"] = float(
                    np.hypot(dpdv * sv, dpdt * 150.0)
                )
        result.append(out)
    return result


def encode_rows(rows):
    """Deterministic CSV for the separate, derived audit artifact."""
    import io

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def series_checks(rows):
    """Compare reconstructed point pressures with the seven approximate isobars."""
    result = []
    for nominal in sorted({r["nominal_series_pressure_gpa"] for r in rows}):
        values = [
            r["reconstructed_pressure_gpa"]
            for r in rows
            if r["nominal_series_pressure_gpa"] == nominal
            and r["reconstructed_pressure_gpa"] != ""
        ]
        result.append(
            {
                "nominal_pressure_gpa": nominal,
                "hcp_rows": len(values),
                "reconstructed_mean_gpa": float(np.mean(values)) if values else None,
                "reconstructed_range_gpa": [min(values), max(values)]
                if values
                else None,
            }
        )
    return result
