"""Independent Walker KCl fit and paired-NaCl input audit (no guessed scale)."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/walker-2002-reproduction.json"
BIRCH_THERMAL_GPA_K = 0.00286  # Birch (1986), original publisher abstract
BIRCH_REFERENCE_K = 298.15  # 25 Celsius; distinct from Walker's KCl reference


def load_rows(filename):
    path = ROOT / "peritheos/data/datasets" / filename
    return list(csv.DictReader(path.read_text().splitlines())), hashlib.sha256(
        path.read_bytes()
    ).hexdigest()


def values(rows, name):
    return np.array([float(row[name]) for row in rows])


def bm3(volume, v0, k0, kp, literal_printed_signs=False):
    f = ((v0 / volume) ** (2 / 3) - 1) / 2
    sign = -1 if literal_printed_signs else 1
    return 3 * k0 * f * (1 + sign * 2 * f) ** 2.5 * (1 + sign * 1.5 * f * (kp - 4))


def fit_metrics(prediction, observed):
    residual = prediction - observed
    return {
        "sum_squared_pressure_residuals_gpa2": float(residual @ residual),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_absolute_residual_gpa": float(np.max(np.abs(residual))),
    }


def paired_inputs(rows, file_column):
    """Only known thermal increments; a full independent replay is blocked."""
    result = []
    for index, row in enumerate(rows):
        temperature = float(row["temperature_celsius"]) + 273.15
        printed_pressure = (
            float(row["pressure_kbar"]) * 0.1 if row["pressure_kbar"] else None
        )
        result.append(
            {
                "source_row_index": index,
                "calibrant_file": row[file_column],
                "row_kind": row.get("row_kind", "sample_observation"),
                "temperature_k": temperature,
                "nacl_lattice_a_angstrom": float(row["nacl_lattice_a_angstrom"]),
                "nacl_conventional_cell_volume_a3": float(
                    row["nacl_lattice_a_angstrom"]
                )
                ** 3,
                "reported_pressure_gpa": printed_pressure,
                "birch_abstract_thermal_increment_gpa": BIRCH_THERMAL_GPA_K
                * (temperature - BIRCH_REFERENCE_K),
                "outside_birch_abstract_temperature_range": (
                    temperature < 298.15 or temperature > 773.15
                ),
                "independently_replayed_pressure_gpa": None,
            }
        )
    return result


def reproduce():
    b1_all, b1_hash = load_rows("kcl-walker-2002-table1-pvt.csv")
    b1 = [r for r in b1_all if r["included_in_fit"] == "1"]
    b2, b2_hash = load_rows("kcl-walker-2002-table2-pvt.csv")
    v = values(b1, "b1_kcl_cell_volume_a3")
    p = values(b1, "pressure_kbar") * 0.1
    dt = values(b1, "temperature_celsius") + 273.15 - 296.15
    shape = bm3(v, 249.53, 1, 5)
    design = np.column_stack((shape, dt))
    fitted = np.linalg.lstsq(design, p, rcond=None)[0]
    nonlinear = least_squares(
        lambda z: bm3(v, 249.53, z[0], 5) + z[0] * z[1] * dt - p,
        [17.7, 0.00011],
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
    )
    literal_shape = bm3(v, 249.53, 1, 5, literal_printed_signs=True)
    literal_design = np.column_stack((literal_shape, dt))
    literal = np.linalg.lstsq(literal_design, p, rcond=None)[0]
    published = shape * 17.7 + 0.00195 * dt
    # Figure 1 prints 37.50 cm3/mol, inconsistent with Tables 1/3's 249.53 A3.
    # This is a separate source-backed sensitivity case, not a catalog change.
    figure_v0 = 37.50 * 4e24 / 6.02214076e23
    figure_shape = bm3(v, figure_v0, 1, 5)
    figure_design = np.column_stack((figure_shape, dt))
    figure_fit = np.linalg.lstsq(figure_design, p, rcond=None)[0]
    b1_result = {
        "observations": len(b1),
        "source_row_indices": [
            i for i, r in enumerate(b1_all) if r["included_in_fit"] == "1"
        ],
        "objective": "unweighted sum of squared pressure residuals; Walker page 808",
        "fixed_parameters": {"V0": 249.53, "K0_prime": 5, "Tr": 296.15},
        "joint_fitted_parameters": {
            "K0": float(fitted[0]),
            "alpha_KT": float(fitted[1]),
            "alpha0": float(fitted[1] / fitted[0]),
        },
        "published": fit_metrics(published, p),
        "joint_refit": fit_metrics(design @ fitted, p),
        "normal_equation_max_absolute_gradient": float(
            np.max(np.abs(design.T @ (design @ fitted - p)))
        ),
        "direct_K0_alpha0_solver": {
            "success": bool(nonlinear.success),
            "K0": float(nonlinear.x[0]),
            "alpha_KT": float(np.prod(nonlinear.x)),
        },
        "literal_printed_BE1_diagnostic": {
            "K0": float(literal[0]),
            "alpha_KT": float(literal[1]),
            **fit_metrics(literal_design @ literal, p),
            "qualification": "Inconsistent signs with positive strain; not a catalog model.",
        },
        "fixed_published_K0_diagnostic": {
            "K0": 17.7,
            "alpha_KT": float(dt @ (p - shape * 17.7) / (dt @ dt)),
            "qualification": "Conditional diagnostic; source does not specify this B1 staging.",
        },
        "figure1_reference_volume_diagnostic": {
            "source_molar_V0_cm3_per_mol": 37.50,
            "conventional_cell_V0_a3": figure_v0,
            "K0": float(figure_fit[0]),
            "alpha_KT": float(figure_fit[1]),
            "alpha0": float(figure_fit[1] / figure_fit[0]),
            **fit_metrics(figure_design @ figure_fit, p),
            "qualification": "Figure 1's rounded 37.50 cm3/mol differs from "
            "Tables 1/3's 249.53 A3 (37.57 cm3/mol). This volume choice "
            "recovers K0 to printed rounding and beta within its printed "
            "error width, but the author workbook is unavailable; no exact "
            "solver-input or uncertainty parity is established. Table-based "
            "catalog reference volume remains unchanged.",
        },
        "status": "similar_not_exact",
        "qualification": "Joint fit follows the stated objective and fixed B1 inputs. "
        "The B2-only staging footnote does not prescribe B1 staging. "
        "Excel workbook, solver settings, unrounded observations and source covariance "
        "are unavailable. Published coefficients and error widths remain unchanged.",
    }
    v = values(b2, "b2_kcl_cell_volume_a3")
    p = values(b2, "pressure_kbar") * 0.1
    dt = values(b2, "temperature_celsius") + 273.15 - 296.15
    f = ((53.53 / v) ** (2 / 3) - 1) / 2
    shape = 3 * f * (1 + 2 * f) ** 2.5
    b2_fits = {}
    for label, selection in (
        ("eight_23_24_celsius_reference_rows", dt <= 1),
        ("seven_exact_23_celsius_rows_sensitivity", dt == 0),
    ):
        # Algebraic equivalent of fitting K0 and K0' at fixed V0.
        design = np.column_stack(
            (shape[selection], 1.5 * f[selection] * shape[selection])
        )
        coefficients = np.linalg.lstsq(design, p[selection], rcond=None)[0]
        k0, kp = coefficients[0], 4 + coefficients[1] / coefficients[0]
        cold = bm3(v, 53.53, k0, kp)
        slope = dt @ (p - cold) / (dt @ dt)
        b2_fits[label] = {
            "cold_source_row_indices": np.flatnonzero(selection).tolist(),
            "K0": float(k0),
            "K0_prime": float(kp),
            "alpha_KT": float(slope),
            **fit_metrics(cold + slope * dt, p),
        }
    return {
        "audit_date": "2026-10-08",
        "source_resources": {
            "table1_sha256": b1_hash,
            "table2_sha256": b2_hash,
            "walker_final_pdf_sha256": "01585fbd045709f3c08516aeedb30dcb553aec1d3073542c5450a72a6c3ff4b6",
        },
        "b1_fit": b1_result,
        "b2_staged_fit": {
            "fixed_V0_a3": 53.53,
            "Tr_k": 296.15,
            "qualification": "Table 3 B2 footnote: cold elastic fit then thermal fit. "
            "Eight-row diagnostic treats the 24 Celsius row as a reference-isotherm "
            "point in the first stage; its actual 297.15 K is retained in the thermal "
            "stage and paired-NaCl input ledger. Seven-row result is sensitivity only.",
            "fits": b2_fits,
        },
        "nacl_pressure_replay": {
            "ancestry": "Walker page 806: Birch (1986) NaCl-B1 BE2 thermal EOS, both phases",
            "status": "blocked_missing_verified_reference_parameters",
            "reference_temperature_k": BIRCH_REFERENCE_K,
            "thermal_coefficient_gpa_per_k": BIRCH_THERMAL_GPA_K,
            "thermal_evidence": "Birch original publisher abstract; DOI 10.1029/JB091iB05p04949",
            "qualification": "Known thermal increment is only one term of the scale. "
            "No independent full pressure is calculated without the adjusted 25 Celsius "
            "BE2 coefficients and verified volume/reference normalization. Thermal "
            "increments are neither observations nor a ready dataset reduction. "
            "Bracketed ambient pressures in Table 1 are imposed zeroes. "
            "600 Celsius lies outside the abstract's 25-500 Celsius domain.",
            "table1_rows": paired_inputs(b1_all, "spectrum"),
            "table2_rows": paired_inputs(b2, "nacl_file"),
            "table2_calibrant_reference_anchors": [
                {
                    "nacl_file": "r34439",
                    "temperature_k": 296.15,
                    "nacl_lattice_a_angstrom": 5.6414,
                    "nacl_lattice_a_esd_angstrom": 0.0014,
                },
                {
                    "nacl_file": "r35101",
                    "temperature_k": 297.15,
                    "nacl_lattice_a_angstrom": 5.6468,
                    "nacl_lattice_a_esd_angstrom": 0.0004,
                },
            ],
        },
    }


def check_report(actual, expected, path="report"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for key in actual:
            check_report(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for index, (a, e) in enumerate(zip(actual, expected)):
            check_report(a, e, f"{path}[{index}]")
    elif isinstance(actual, float):
        assert np.isclose(actual, expected, rtol=1e-7, atol=1e-10), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    if args.check:
        check_report(result, json.loads(OUTPUT.read_text()))
    else:
        OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "b1": result["b1_fit"]["joint_fitted_parameters"],
                "nacl_pressure_replay": result["nacl_pressure_replay"]["status"],
            }
        )
    )


if __name__ == "__main__":
    main()
