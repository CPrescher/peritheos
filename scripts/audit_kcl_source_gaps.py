"""Primary-method and coordinate-precision diagnostics for Tateno/Campbell KCl.

This report does not replace the published EOS or recover unspecified author
weights, covariance, or pressure reductions. Article PDFs are not redistributed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts.reproduce_kcl_variants import (
    A3_PER_MOLAR,
    bm3,
    column,
    metrics,
    mgd,
    reproduce,
    rows,
    vinet,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/kcl-source-gaps-audit.json"
SOURCE = ROOT / "peritheos/data/datasets/kcl_variant_sources"


def campbell_diagnostics():
    data = rows("kcl-campbell-1991-table1-compression.csv")
    p = column(data, "pressure_gpa")
    ratio = column(data, "normalized_volume")
    sigma_p = column(data, "pressure_uncertainty_gpa")
    sigma_ratio = column(data, "normalized_volume_uncertainty")
    # Campbell page 496: Jeanloz effective strain uses ambient *B1* V01.
    g = (ratio ** (-2 / 3) - 1) / 2
    divisor = 3 * (1 + 2 * g) ** 2.5
    stress = p / divisor
    sigma_g_pressure = sigma_p / divisor
    sigma_g_volume = (5 / 3) * stress * sigma_ratio / ratio
    design = np.column_stack((np.ones(len(data)), g))
    fits = {}
    for name, sigma in (
        ("equal_normalized_stress_weights", np.ones(len(data))),
        ("pressure_only_normalized_stress_weights", sigma_g_pressure),
        (
            "pressure_and_volume_response_weights",
            np.hypot(sigma_g_pressure, sigma_g_volume),
        ),
    ):
        intercept, slope = np.linalg.lstsq(
            design / sigma[:, None], stress / sigma, rcond=None
        )[0]
        # Algebraic BM2 mapping: G=b+m*g; m=K02*r0^(7/3).
        r0 = (1 - 2 * intercept / slope) ** (-1.5)
        k0 = slope / r0 ** (7 / 3)
        fits[name] = {
            "V02_over_V01": float(r0),
            "K02_gpa": float(k0),
            "intercept_gpa": float(intercept),
            "slope_gpa": float(slope),
            "pressure_residuals": metrics(bm3(ratio, r0, k0, 4), p),
            "within_printed_parameter_error_widths": bool(
                abs(r0 - 0.8483) <= 0.0057 and abs(k0 - 28.7) <= 0.6
            ),
        }
    # Analytic scale conversion at the *reported mean* pressure, not inversion
    # of any individual measured ruby line or recovery of five-point averages.
    ruby_ratio = (1 + 5 * p / 1904) ** (1 / 5)
    mao86 = 1904 / 7.665 * (ruby_ratio**7.665 - 1)
    return {
        "source_calibration": "ruby_mao_1978",
        "source_location": "Experimental procedure, pages 495-496; reference 5, page 499",
        "observations": len(data),
        "source_volume_basis": {
            "B1_molar_volume_cm3_mol": 37.521,
            "B1_lattice_a_angstrom": 6.2931,
            "qualification": "Printed source normalizations; distinct from the executable Campbell-Dewaele composite V01=62.36 A3/formula unit.",
        },
        "published_BM2": {"V02_over_V01": 0.8483, "K02_gpa": 28.7, "K02_prime": 4},
        "published_pressure_residuals": metrics(bm3(ratio, 0.8483, 28.7, 4), p),
        "fits": fits,
        "rows": [
            {
                "run": row["run"],
                "reported_mean_mao1978_gpa": float(pp),
                "effective_strain_g": float(gg),
                "normalized_stress_G_gpa": float(ss),
                "diagnostic_mean_coordinate_mao1986_gpa": float(converted),
            }
            for row, pp, gg, ss, converted in zip(data, p, g, stress, mao86)
        ],
        "exact_source_regression_reproduced": False,
        "qualification": "The source specifies weighted linear G-versus-g regression but not numerical weights, predictor-response covariance or parameter-error confidence. The two explicit response-weight choices recover both coefficients within printed widths. They are conditional diagnostics; volume error also affects g and correlates it with G, which these response-only fits do not model. Source pressure errors are spatial standard deviations of five readings, not standard errors of their mean. Conversion of the printed mean does not recover the mean of five nonlinearly converted readings. Raw rows and the published composite EOS are unchanged.",
    }


def tateno_precision_diagnostics():
    data = rows("kcl-tateno-2019-official-table-s1.csv")
    rounded = rows("kcl-tateno-2019-table-s1-pvt.csv")
    v = column(data, "kcl_unit_cell_volume_a3") / A3_PER_MOLAR
    t = column(data, "temperature_k")
    vp = column(data, "platinum_unit_cell_volume_a3")
    rounded_vp = column(rounded, "platinum_unit_cell_volume_a3")

    def holmes(volumes):
        return vinet(volumes, 60.4000884, 798.31 / 3, 1 + 7.2119 / 1.5) + (
            0.0069426 * (t - 300)
        )

    pressure = holmes(vp)
    rounded_pressure = holmes(rounded_vp)
    fits = {}
    for kind, initial in (
        ("mgd", [17.4, 5.77, 1.8, 0.7]),
        ("linear", [17.7, 5.73, 0.0033]),
    ):
        names = ["K0", "K0_prime"] + (
            ["gamma0", "q"] if kind == "mgd" else ["alpha_KT"]
        )

        def model(z):
            if kind == "mgd":
                return mgd(v, t, [54.5 / A3_PER_MOLAR, *z])
            return vinet(v, 54.5 / A3_PER_MOLAR, z[0], z[1]) + z[2] * (t - 300)

        fit = least_squares(
            lambda z: model(z) - rounded_pressure,
            initial,
            bounds=([1, 2, 0.1, 0.01], [50, 10, 8, 4])
            if kind == "mgd"
            else (-np.inf, np.inf),
            x_scale="jac",
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        identifier = "kcl_b2_tateno_2019_holmes_vinet_" + (
            "mgd" if kind == "mgd" else "linear_thermal"
        )
        fits[kind] = {
            "full_workbook_Pt_volume_fit": reproduce()[identifier]["fitted_parameters"],
            "rounded_Pt_volume_only_fit": dict(zip(names, fit.x.tolist())),
            "rounded_Pt_volume_only_residuals": metrics(model(fit.x), rounded_pressure),
            "solver_success": bool(fit.success),
        }
    # Inspect the original OpenXML package without evaluating any invented EOS.
    with zipfile.ZipFile(SOURCE / "6779TableS1 revised.xlsx") as archive:
        ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        package = {
            "sheets": [
                s.attrib["name"] for s in workbook.findall("s:sheets/s:sheet", ns)
            ],
            "formula_cells": len(sheet.findall(".//s:f", ns)),
            "external_links": [n for n in archive.namelist() if "externalLink" in n],
        }
    return {
        "workbook_package": package,
        "observations": len(data),
        "coordinate_precision_sensitivity": {
            "maximum_absolute_pressure_shift_gpa": float(
                np.max(abs(rounded_pressure - pressure))
            ),
            "fits": fits,
            "qualification": "Only Pt volume is replaced by its corrected two-decimal transcription; source T and KCl volumes remain at workbook precision. Equal pressure weights in both cases. This is sensitivity to stored precision, not recovered author rounding or extra experimental precision.",
        },
        "outside_holmes_thermal_approximation_excel_rows": [
            int(row["source_excel_row"]) for row, temp in zip(data, t) if temp >= 2000
        ],
        "exact_source_regression_reproduced": False,
        "qualification": "Final article specifies simultaneous elastic/thermal fitting with fixed V0 and theta0, but names only the Holmes scale. Equation 12 versus full thermodynamic Pt reduction, any correction above 2000 K, author input precision, numerical weights and covariance remain unspecified. The official supplemental PDF is Figure S1 only. No Holmes pressure column, formulas, code or fit protocol are deposited in the inspected archive. Separate author Sokolova coordinates and conditional Holmes coordinates remain in the existing reproduction report.",
    }


def audit():
    return {
        "audit_date": "2026-10-08",
        "resource_sha256": {
            name: hashlib.sha256(
                (ROOT / "peritheos/data/datasets" / name).read_bytes()
            ).hexdigest()
            for name in (
                "kcl-campbell-1991-table1-compression.csv",
                "kcl-tateno-2019-official-table-s1.csv",
                "kcl-tateno-2019-table-s1-pvt.csv",
            )
        },
        "campbell_1991": campbell_diagnostics(),
        "tateno_2019": tateno_precision_diagnostics(),
    }


def assert_report_equal(actual, expected):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for key in actual:
            assert_report_equal(actual[key], expected[key])
    elif isinstance(actual, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            assert_report_equal(a, b)
    elif isinstance(actual, float):
        assert np.isclose(actual, expected, rtol=1e-6, atol=1e-9)
    else:
        assert actual == expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = audit()
    if args.check:
        assert_report_equal(result, json.loads(OUTPUT.read_text()))
    else:
        OUTPUT.write_text(
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    print("Campbell calibration resolved; exact source regressions remain conditional.")


if __name__ == "__main__":
    main()
