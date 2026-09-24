"""Audit the Ross supporting study; this does not reproduce the 1986 EOS.

The available full text is the distinct 1985 UCRL-93030 conference precursor.
Only its clearly resolved Figure 1 DAC symbols and Eq. (2) pair potential are
represented. A pair potential is not a thermodynamic pressure-volume EOS.
Run ``python -m scripts.reproduce_argon_ross`` to regenerate the audit report.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/argon-ross-1985-precursor-figure1-subset.csv"
REPORT = ROOT / "docs/data/argon-ross-1986-reproduction.json"
AVOGADRO = 6.02214076e23


def precursor_exp6_energy_kelvin(radius_angstrom, alpha=13.2):
    """UCRL-93030 Eq. (2), phi/kB, on the discussed r >= 2 A interval.

    This is neither a pressure evaluator nor a reproduction of a Monte Carlo
    isotherm. The short-distance exp-6 catastrophe is outside this domain.
    Alpha 13.0 is the paper's comparison potential, not a statistical error.
    """
    if not math.isfinite(radius_angstrom) or radius_angstrom < 2.0:
        raise ValueError("The audited precursor potential requires r >= 2 A")
    if alpha not in (13.0, 13.2):
        raise ValueError("Only the two printed precursor alpha values are audited")
    ratio = radius_angstrom / 3.85
    return (
        122.0
        / (alpha - 6.0)
        * (6.0 * math.exp(alpha * (1.0 - ratio)) - alpha / ratio**6)
    )


def pixel_to_pv(x, y):
    """Affine Figure 1 crop calibration; returns GPa and cm3/mol of atoms."""
    # 240 dpi pdftoppm rendering, cropped to (580, 170, 1280, 1140).
    volume = 8.0 + (x - 60.0) * 11.0 / (620.0 - 60.0)
    pressure = (920.0 - y) * 80.0 / (920.0 - 125.0)
    return pressure, volume


def molar_to_cell(volume_cm3_mol):
    """cm3/mol of Ar atoms -> A3/four-atom conventional fcc cell."""
    return volume_cm3_mol * 4e24 / AVOGADRO


def audit():
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    errors = []
    for row in rows:
        pressure, molar = pixel_to_pv(float(row["pixel_x"]), float(row["pixel_y"]))
        errors.extend(
            [
                abs(pressure - float(row["pressure_gpa"])),
                abs(molar - float(row["molar_volume_cm3_mol"])),
                abs(molar_to_cell(molar) - float(row["volume_a3"])),
            ]
        )
    return {
        "study_identifier": "argon_ross_1986_supporting_study",
        "requested_doi": "10.1063/1.451346",
        "reproduction_status": "not_reproduced",
        "matching_journal_pdf_obtained": False,
        "source_report": "UCRL-93030 (1985), a distinct conference precursor",
        "dataset_identifier": "argon_ross_1985_precursor_figure1_subset",
        "dataset_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
        "digitized_static_points": len(rows),
        "complete_dataset": False,
        "temperature_k": 293.0,
        "pressure_range_gpa": [
            min(float(r["pressure_gpa"]) for r in rows),
            max(float(r["pressure_gpa"]) for r in rows),
        ],
        "conversion_max_absolute_roundoff": max(errors),
        "coordinate_reading_bound_pixels": 10,
        "pressure_reading_bound_gpa": 10 * 80 / 795,
        "molar_volume_reading_bound_cm3_mol": 10 * 11 / 560,
        "reading_bounds_are_experimental_uncertainties": False,
        "precursor_pair_potential_checks": [
            {
                "alpha": alpha,
                "minimum_radius_angstrom": 3.85,
                "energy_at_minimum_kelvin": precursor_exp6_energy_kelvin(3.85, alpha),
                "energy_at_2_angstrom_kelvin": precursor_exp6_energy_kelvin(2.0, alpha),
            }
            for alpha in (13.2, 13.0)
        ],
        "eos_fit_performed": False,
        "pressure_residual_rmse_gpa": None,
        "reason": "Journal full text unavailable; precursor Monte Carlo EOS not reproduced. No surrogate EOS fitted.",
    }


if __name__ == "__main__":
    REPORT.write_text(json.dumps(audit(), indent=2) + "\n")
    print(REPORT)
