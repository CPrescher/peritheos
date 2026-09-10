#!/usr/bin/env python3
"""Reproduce Funamori et al.'s (1998) two MgAl2O4 BM2 constructions.

Each source EOS is constrained by a recovered ambient cell and one in-situ
high-pressure cell.  This is not a regression: the source fixes the separately
measured V0, assumes K0'=4, and obtains K0 as the normalized pressure of the
single compressed state.  The script reproduces that calculation from the
printed lattice axes and also propagates the ruby/Pt pressure bracket.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"
CASES = {
    "cafe2o4_type": {
        "path": DATA / "mgal2o4-cafe2o4-funamori-1998-text-pv.csv",
        "V0": 240.3,
        "K0": 211.0,
    },
    "cati2o4_type": {
        "path": DATA / "mgal2o4-cati2o4-funamori-1998-text-pv.csv",
        "V0": 240.3,
        "K0": 206.0,
    },
}


def bm2_pressure(volume: np.ndarray | float, v0: float, k0: float) -> np.ndarray:
    """Return second-order Birch-Murnaghan pressure in GPa."""
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5)


def normalized_pressure(volume: float, pressure: float, v0: float) -> float:
    """Return Funamori et al.'s second-order Birch normalized pressure."""
    strain = 0.5 * ((v0 / volume) ** (2.0 / 3.0) - 1.0)
    return pressure / (3.0 * strain * (1.0 + 2.0 * strain) ** 2.5)


def reproduce() -> dict[str, dict[str, float | int]]:
    """Return deterministic endpoint-construction diagnostics."""
    result = {}
    for phase, case in CASES.items():
        with case["path"].open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        high = rows[1]
        lattice_volume = (
            float(high["a_angstrom"])
            * float(high["b_angstrom"])
            * float(high["c_angstrom"])
        )
        pressure = float(high["pressure_gpa"])
        ruby_pressure = float(high["pressure_ruby_after_heating_gpa"])
        pt_pressure = float(high["pressure_pt_after_heating_gpa"])
        endpoint_k0 = normalized_pressure(lattice_volume, pressure, case["V0"])
        ruby_k0 = normalized_pressure(lattice_volume, ruby_pressure, case["V0"])
        pt_k0 = normalized_pressure(lattice_volume, pt_pressure, case["V0"])
        published_pressure = float(bm2_pressure(lattice_volume, case["V0"], case["K0"]))
        result[phase] = {
            "lattice_product_volume_a3": lattice_volume,
            "source_pressure_gpa": pressure,
            "published_curve_pressure_gpa": published_pressure,
            "published_curve_pressure_residual_gpa": published_pressure - pressure,
            "fixed_v0_endpoint_k0_gpa": endpoint_k0,
            "ruby_endpoint_k0_gpa": ruby_k0,
            "pt_endpoint_k0_gpa": pt_k0,
            "pressure_bracket_half_range_k0_gpa": abs(ruby_k0 - pt_k0) / 2.0,
            "informative_finite_pressure_observations": 1,
        }
    return result


def main() -> None:
    """Print the reproduction as stable JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
