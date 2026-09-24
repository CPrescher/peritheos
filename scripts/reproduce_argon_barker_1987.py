"""Audit Barker's transcribed table and digitized figure; not an EOS solver.

Run from the repository root with ``python -m scripts.reproduce_argon_barker_1987``.
Interpolation only compares marks in a published figure. It is not a fit or an
independent Monte Carlo reproduction and must not be registered as an EOS.
"""

import csv
import json
from bisect import bisect_right
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/argon-barker-1987-reproduction.json"
AVOGADRO = 6.02214076e23
GAS_CONSTANT = 8.31446261815324


def molar_to_atomic_a3(volume_cm3_mol):
    return volume_cm3_mol * 1e24 / AVOGADRO


def pixel_to_source(x, y):
    """Linear axes in figure5-detail.png, cropped from PDF page 4 at 300 dpi."""
    return 8 + (x - 150) * 14 / 808, (918 - y) * 450 / 804


def load_rows(name):
    with (DATA / name).open() as stream:
        return list(csv.DictReader(stream))


def audit():
    table = load_rows("argon-barker-1987-table1-liquid.csv")
    experiment = next(r for r in table if r["determination_method"] == "experimental")
    experimental_z = float(experiment["pv_over_nkt"])
    liquid = []
    for row in table:
        z = float(row["pv_over_nkt"])
        t = float(row["temperature_k"])
        vm = float(row["molar_volume_cm3_mol"])
        liquid.append(
            {
                "label": row["label"],
                "determination_method": row["determination_method"],
                "pressure_gpa": z * GAS_CONSTANT * t / (vm * 1e-6) / 1e9,
                "molar_internal_energy_j_mol": float(row["u_over_nkt"])
                * GAS_CONSTANT
                * t,
                "relative_pressure_difference_percent": 100 * (z / experimental_z - 1),
            }
        )
    points = load_rows("argon-barker-1987-figure5-selected.csv")
    curve = sorted(
        (float(r["molar_volume_cm3_mol"]), float(r["pressure_gpa"]))
        for r in points
        if r["series"] == "bfw_at_theoretical_curve"
    )
    volumes = [p[0] for p in curve]
    comparisons = []
    excluded = 0
    for row in points:
        if row["series"] != "ross_1986_experimental_comparison":
            continue
        v = float(row["molar_volume_cm3_mol"])
        p = float(row["pressure_gpa"])
        if not volumes[0] <= v <= volumes[-1]:
            excluded += 1
            continue
        index = min(max(bisect_right(volumes, v), 1), len(curve) - 1)
        v1, p1 = curve[index - 1]
        v2, p2 = curve[index]
        interpolated = p1 + (p2 - p1) * (v - v1) / (v2 - v1)
        comparisons.append(
            {
                "point_id": row["point_id"],
                "molar_volume_cm3_mol": v,
                "digitized_experiment_pressure_gpa": p,
                "interpolated_digitized_theory_pressure_gpa": interpolated,
                "theory_minus_experiment_gpa": interpolated - p,
            }
        )
    return {
        "study_identifier": "argon_barker_1987",
        "status": "not_reproduced",
        "diagnostic_kind": "published_table_and_digitized_figure_comparison",
        "independent_model_reproduction": False,
        "fit_performed": False,
        "liquid_table": liquid,
        "figure5": {
            "theory_sample_count": len(curve),
            "selected_experimental_point_count": sum(
                r["series"] == "ross_1986_experimental_comparison" for r in points
            ),
            "experimental_points_outside_sampled_curve": excluded,
            "comparisons": comparisons,
            "note": "Selected visually separable circles only, not a complete dataset. "
            "Linear interpolation of digitized published theory is a graphical "
            "diagnostic, not independent EOS evaluation. Pixel bounds are not "
            "experimental error bars or Monte Carlo sampling uncertainties.",
        },
    }


if __name__ == "__main__":
    REPORT.write_text(json.dumps(audit(), indent=2) + "\n")
    print(REPORT)
