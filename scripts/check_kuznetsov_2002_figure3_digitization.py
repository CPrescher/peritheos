#!/usr/bin/env python3
"""Check the hcp-Pb marker digitization from Kuznetsov et al. (2002), Fig. 3."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"
CSV_PATH = DATA / "lead-hcp-kuznetsov-2002-figure3-digitized.csv"

# Tick centers measured in the native 1713 x 1312 pixel Figure 3 image extracted
# from the source PDF. Pixel y increases downwards.
X_TICKS = np.array([441.5, 637.5, 834.0, 1028.5, 1224.5, 1419.0, 1613.5])
X_VALUES = np.array([0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
Y_TICKS = np.array([967.5, 855.0, 742.5, 630.5, 518.5, 406.5, 294.0, 181.5, 69.5])
Y_VALUES = np.array([0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0])
V0_FCC_A3_PER_ATOM = 30.307


def _linear_calibration(pixels: np.ndarray, values: np.ndarray) -> tuple[float, float]:
    slope, intercept = np.polyfit(pixels, values, 1)
    return float(slope), float(intercept)


def check() -> dict[str, object]:
    x_slope, x_intercept = _linear_calibration(X_TICKS, X_VALUES)
    y_slope, y_intercept = _linear_calibration(Y_TICKS, Y_VALUES)
    with CSV_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    maximum_rounding_difference = 0.0
    for row in rows:
        ratio = x_slope * float(row["pixel_x"]) + x_intercept
        pressure = y_slope * float(row["pixel_y"]) + y_intercept
        volume = ratio * V0_FCC_A3_PER_ATOM
        maximum_rounding_difference = max(
            maximum_rounding_difference,
            abs(ratio - float(row["v_over_v0_fcc"])),
            abs(pressure - float(row["reduced_pressure_gpa"])),
            abs(volume - float(row["atomic_volume_a3_per_atom"])),
        )

    return {
        "dataset": CSV_PATH.name,
        "sha256": hashlib.sha256(CSV_PATH.read_bytes()).hexdigest(),
        "rows": len(rows),
        "open_symbols": sum(row["symbol_fill"] == "open" for row in rows),
        "solid_symbols": sum(row["symbol_fill"] == "solid" for row in rows),
        "pressure_range_gpa": [
            min(float(row["reduced_pressure_gpa"]) for row in rows),
            max(float(row["reduced_pressure_gpa"]) for row in rows),
        ],
        "axis_calibration": {
            "v_over_v0_fcc": [x_slope, x_intercept],
            "reduced_pressure_gpa": [y_slope, y_intercept],
        },
        "maximum_stored_rounding_difference": maximum_rounding_difference,
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
