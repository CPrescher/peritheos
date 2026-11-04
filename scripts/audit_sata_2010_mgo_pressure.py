"""Check Sata's MgO pressure reduction against the original Speziale source.

Repository-only diagnostic; no observations or catalog coefficients are changed.
Run with --write to archive the report and derived-coordinate table, or --check
to verify them. Original source tables and the catalog are never rewritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from scripts.audit_sata_2010_fe3c import TOKENS, observations, pressure, statistics
from scripts.reproduce_speziale_2001_mgo import pressure as thermal_pressure

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/data/sata-2010-mgo-pressure-check.json"
POINTS = ROOT / "docs/data/sata-2010-fe3c-recalculated-points.csv"
ANCHORS = np.array([74.71, 160.2, 3.99])
HALF_DIGITS = np.array([0.005, 0.05, 0.005])


def bm3(volume, parameters=ANCHORS):
    v0, k0, kp = parameters
    return pressure(volume, [k0, kp], v0)


def vinet(volume, kp=4.12):
    """Speziale section 5.3 Vinet control, with the printed V0 adopted."""
    x = (np.asarray(volume) / 74.71) ** (1 / 3)
    return 3 * 160.2 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def logarithmic(volume):
    """Section 5.3 logarithmic control, with the printed V0 adopted."""
    f = np.log(74.71 / np.asarray(volume)) / 3
    return 3 * 160.2 * f * np.exp(3 * f) * (1 + 1.5 * (4.21 - 2) * f)


def recalculated_sample_points():
    """Pair the same Fe3C volume with original and independently reduced MgO P.

    Derived pressure coordinates only: no Fe3C refit, replacement observation,
    inferred total error, Li recalibration or fitted MgO coefficients.
    """
    rows = observations()[1]
    recalculated = bm3(np.array([r["mgo_volume_a3"] for r in rows]))
    source_prediction = pressure(np.array([r["volume_a3"] for r in rows]))
    return [
        {
            "source_row": r["source_row"],
            "run": r["run"],
            "measurement_id": r["measurement_id"],
            "temperature_k": r["temperature_k"],
            "fe3c_volume_a3": r["volume_a3"],
            "fe3c_volume_reported_width_a3": r["error_volume_a3"],
            "mgo_volume_a3": r["mgo_volume_a3"],
            "mgo_volume_reported_width_a3": r["error_mgo_volume_a3"],
            "published_pressure_gpa": r["pressure_gpa"],
            "published_pressure_reported_width_gpa": r["error_pressure_gpa"],
            "recalculated_pressure_gpa": float(p),
            "recalculated_minus_published_pressure_gpa": float(p - r["pressure_gpa"]),
            "published_fe3c_curve_pressure_gpa": float(pred),
            "published_pressure_residual_gpa": float(pred - r["pressure_gpa"]),
            "recalculated_pressure_residual_gpa": float(pred - p),
        }
        for r, p, pred in zip(rows, recalculated, source_prediction)
    ]


def points_csv(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def rounding_minimax(volume, published):
    """One common coefficient set; allow every P and V its own half digit.

    The last variable is the maximum remaining interval gap in GPa. Positive
    means that even this liberal rounding box leaves an uncovered row. This is
    a bounded numerical diagnostic, not a coefficient-confidence calculation.
    """

    def interval_constraints(z):
        parameters = ANCHORS + HALF_DIGITS * z[:3]
        low = bm3(volume + 0.0005, parameters)
        high = bm3(volume - 0.0005, parameters)
        return np.r_[
            z[3] - (low - published - 0.005), z[3] - (published - 0.005 - high)
        ]

    fits = []
    for start in ([0, 0, 0, 0.4], [-1, -1, -1, 1], [1, 1, 1, 1]):
        fit = minimize(
            lambda z: z[3],
            start,
            method="SLSQP",
            # The pressure and interval expressions are analytic. Complex-step
            # derivatives avoid cancellation at the tight rounding-box optimum
            # on older SciPy builds, without changing the objective or bounds.
            jac="cs",
            bounds=[(-1, 1)] * 3 + [(-1, 2)],
            constraints={"type": "ineq", "fun": interval_constraints},
            options={"ftol": 1e-12, "maxiter": 1000},
        )
        if not fit.success or min(interval_constraints(fit.x)) < -1e-8:
            raise RuntimeError(f"Rounding minimax failed: {fit.message}")
        fits.append(fit)
    best = min(fits, key=lambda fit: fit.fun)
    return {
        "coefficient_box_V0_K0_Kprime": np.column_stack(
            [ANCHORS - HALF_DIGITS, ANCHORS + HALF_DIGITS]
        ).tolist(),
        "volume_half_digit_a3": 0.0005,
        "pressure_half_digit_gpa": 0.005,
        "parameters_at_best_bound": (ANCHORS + HALF_DIGITS * best.x[:3]).tolist(),
        "remaining_max_interval_gap_gpa": float(best.fun),
        "multistart_gap_spread_gpa": float(np.ptp([fit.fun for fit in fits])),
        "qualification": "Nearest rounding only; independent closed coefficient bounds are deliberately liberal, including K0 although source K0 is fixed. No probability, confidence interval or recovered author coefficients.",
    }


def reproduce():
    rows = observations()[1]
    v = np.array([r["mgo_volume_a3"] for r in rows])
    a = np.array([r["mgo_a_angstrom"] for r in rows])
    p = np.array([r["pressure_gpa"] for r in rows])
    literal = bm3(v)
    eta = (74.71 / v) ** (1 / 3)
    base = 1.5 * 160.2 * (eta**7 - eta**5)
    sensitivity = base * 0.75 * (eta**2 - 1)
    implied_kp = 4 + float(sensitivity @ (p - base) / (sensitivity @ sensitivity))
    diagnostic = bm3(v, [74.71, 160.2, implied_kp])
    data_rounding = (
        np.maximum(abs(bm3(v + 0.0005) - literal), abs(bm3(v - 0.0005) - literal))
        + 0.005
    )
    corners_low = bm3(v + 0.0005, ANCHORS - HALF_DIGITS)
    corners_high = bm3(v - 0.0005, ANCHORS + HALF_DIGITS)
    cases = {
        "literal_BM3": literal,
        "BM3_from_lattice_a_cubed": bm3(a**3),
        "BM2_Kprime_4_diagnostic": bm3(v, [74.71, 160.2, 4]),
        "source_unconstrained_BM3": bm3(v, [74.53, 170, 3.59]),
        "source_Vinet_Kprime_4_12_adopted_V0": vinet(v),
        "Vinet_same_BM3_coefficients_diagnostic": vinet(v, 3.99),
        "source_logarithmic_Kprime_4_21_adopted_V0": logarithmic(v),
        "Kprime_only_diagnostic": diagnostic,
    }
    thermal_difference = max(
        float(np.max(abs(thermal_pressure(v, 300, **kwargs) - literal)))
        for kwargs in ({}, {"q1": 0}, {"gamma_law": "local_power"})
    )
    return {
        "scope": "Original Speziale coefficients and Sata MgO pressure-reduction rounding/formulation controls; no catalog or source-data changes.",
        "source": {
            "doi": "10.1029/2000JB900318",
            "url": "https://duffy.princeton.edu/sites/g/files/toruqf616/files/speziale_et_al-2001-jgrse.pdf",
            "pdf_sha256": "c0431bdddf8f4bc8269a7e0d1b79d30d92a38b4a778d9edddad9cb4c45b1f94b",
            "locations": "Section 4.2, p.518; Table 2, p.520; section 5.3, pp.520-521; equations 1-2, p.516.",
            "preferred_BM3": {
                "V0_a3": 74.71,
                "V0_error_a3": 0.01,
                "K0_fixed_gpa": 160.2,
                "adopted_K0_error_gpa": 0.2,
                "Kprime": 3.99,
                "Kprime_error": 0.01,
            },
            "alternative_qualification": "Vinet and logarithmic coefficients are source-published alternatives; their V0 is not restated in section 5.3, so controls explicitly adopt 74.71 A3. They are not claimed as Sata's specified BM3 reduction.",
        },
        "source_table_sha256": hashlib.sha256(TOKENS.read_bytes()).hexdigest(),
        "cases": {
            name: {"predicted_gpa": pred.tolist(), **statistics(pred - p)}
            for name, pred in cases.items()
        },
        "source_rows": [
            {
                "run": r["run"],
                "measurement_id": r["measurement_id"],
                "volume_a3": float(vi),
                "published_pressure_gpa": float(pi),
                "literal_pressure_gpa": float(li),
                "data_rounding_bound_gpa": float(di),
                "coefficient_and_volume_rounding_low_gpa": float(lo),
                "coefficient_and_volume_rounding_high_gpa": float(hi),
            }
            for r, vi, pi, li, di, lo, hi in zip(
                rows, v, p, literal, data_rounding, corners_low, corners_high
            )
        ],
        "rounding": {
            "fixed_coefficient_data_rounding_max_gpa": float(max(data_rounding)),
            "all_targets_in_separate_row_coefficient_boxes": bool(
                np.all((p + 0.005 >= corners_low) & (p - 0.005 <= corners_high))
            ),
            "common_coefficient_set": rounding_minimax(v, p),
            "interpretation": "Coefficient display precision can generate the mismatch's magnitude pointwise, but a common nearest-rounding box does not replay all nine rows at input/output display precision. Statistical coefficient widths are distinct from rounding bounds.",
        },
        "thermal_controls_at_300K_max_difference_gpa": thermal_difference,
        "recalculated_fe3c_coordinates": {
            "kind": "recalculated_pressure_coordinates",
            "scope": "Nine Sata Fe3C points only, using measured paired MgO volumes, fixed published Speziale BM3 coefficients and unchanged Fe3C volumes at 300 K. Li points are not recalculated.",
            "qualification": "Shown alongside original observations. Recalculated pressures are derived quantities, not independently measured pressures or established corrections to the publication. No total recalculated-pressure uncertainties are inferred; source widths remain identified with their original quantities. No Fe3C refit, coefficient replacement or change to catalog recommendation/parity.",
            "csv": str(POINTS.relative_to(ROOT)),
            "rows": recalculated_sample_points(),
            "published_curve_on_original_nine": statistics(
                pressure(np.array([r["volume_a3"] for r in rows])) - p
            ),
            "published_curve_on_recalculated_nine": statistics(
                pressure(np.array([r["volume_a3"] for r in rows])) - literal
            ),
        },
        "Kprime_only_diagnostic": {
            "value": implied_kp,
            "inside_source_error_width": bool(abs(implied_kp - 3.99) <= 0.01),
            "inside_nearest_rounding_half_digit": bool(abs(implied_kp - 3.99) <= 0.005),
            "all_rows_overlap_P_V_display_rounding": bool(
                np.all(
                    (bm3(v + 0.0005, [74.71, 160.2, implied_kp]) <= p + 0.005)
                    & (bm3(v - 0.0005, [74.71, 160.2, implied_kp]) >= p - 0.005)
                )
            ),
            "qualification": "Fit to Sata reduction targets with V0/K0 fixed; demonstrates sensitivity and consistency, not historical coefficients, a replacement MgO EOS or a pressure recalibration.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    if args.write:
        REPORT.write_text(json.dumps(report, indent=2) + "\n")
        POINTS.write_text(points_csv(report["recalculated_fe3c_coordinates"]["rows"]))
    elif args.check:
        archived = json.loads(REPORT.read_text())

        # Optimizer last digits may vary between BLAS/SciPy builds.
        def close(a, b):
            if isinstance(a, dict):
                return a.keys() == b.keys() and all(close(a[k], b[k]) for k in a)
            if isinstance(a, list):
                return len(a) == len(b) and all(close(x, y) for x, y in zip(a, b))
            if isinstance(a, float):
                return bool(np.isclose(a, b, atol=1e-8, rtol=1e-8))
            return a == b

        if not close(report, archived):
            raise SystemExit("Sata MgO pressure-check report needs regeneration")
        if POINTS.read_text() != points_csv(
            report["recalculated_fe3c_coordinates"]["rows"]
        ):
            raise SystemExit("Sata recalculated-coordinate CSV needs regeneration")
        print("Sata MgO source, formulation and rounding checks passed")
    else:
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
