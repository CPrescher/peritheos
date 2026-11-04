"""Replay the author-supplied Bykova I/II/III export, without catalog edits.

--extract reads the original XLSX with openpyxl (use the bundled workspace
Python) and DAT, preserving values, cached formulas, cell addresses and tokens.
The normal scientific replay needs only the dedicated JSON transcription,
NumPy and SciPy. Neither mode writes an XLSX or modifies a source file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "output/bykova-2018-author-files-2026-10-25"
INPUT = ROOT / "docs/data/bykova-2018-author-fit-inputs.json"
RESULT = ROOT / "docs/data/bykova-2018-author-fit-results.json"
PLOT = ROOT / "docs/literature-reproductions/assets/bykova-2018-author-fit.svg"
DIGESTS = {
    "coesite_final.xlsx": "a451c3cc7728f66d33b8a14bac27c53bca8a20f38e964f1975db6c6474fc76ab",
    "coeII-plot.dat": "dbe88e69780a1dff041dadd3a9547a1e95bc2bf4bf8d5fcf10c4da76af8ca560",
}
PUBLISHED = [34.20, 103.0, 3.02]
AUTHOR_CURVE = [34.2038, 102.75, 3.023]
RECORD_ID = "coesite_i_iii_bykova_2018_300k_bm3"
DATASET_ID = "coesite_i_iii_bykova_2018_author_fit_export"
CSV_PATH = (
    ROOT / "peritheos/data/datasets/coesite-i-iii-bykova-2018-author-fit-export.csv"
)
REPORT = "docs/literature-reproductions/bykova-2018-author-fit.md"
FINDING = (
    "Central-coefficient parity is achieved; the author's fit values are "
    "reproduced at reported precision."
)
QUALIFICATION = (
    FINDING + " Published central coefficients are reproduced on the "
    "recovered 22-row EoSFit export: 21 measured entries plus the derived "
    "Cernok (2014) fitted ambient anchor, not 22 independent measurements. "
    "All three BM3 coefficients vary, but the anchor's tiny sigma_P effectively "
    "pins V0 in the author-communicated pressure-only 1/sigma_P^2 objective; "
    "sigma_V is preserved but not propagated. Current-study V/Z is rounded "
    "to 0.1 A3 per formula unit. Original EoSFit settings, sigma_V treatment "
    "and covariance/uncertainty parity remain unverified. The older 24-row "
    "unweighted failure is a separate historical diagnostic, not the source "
    "fit-input reconstruction. Published coefficients are unchanged; this "
    "finding does not establish physical accuracy or recommendation."
)
CURRENT_LIT_ROWS = [2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16]
CURRENT_SUM_ROWS = [2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17]


def verify_archives() -> None:
    for name, digest in DIGESTS.items():
        path = ARCHIVE / name
        if path.exists():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name


def extract() -> dict[str, Any]:
    """Read originals only; use cached formula results as archived evidence."""
    import openpyxl

    verify_archives()
    raw = (ARCHIVE / "coeII-plot.dat").read_text(encoding="utf-8", errors="replace")
    points, curve, headers = [], [], []
    mode = "points"
    for line_number, line in enumerate(raw.splitlines(), 1):
        if line.startswith("EoS curve:"):
            mode = "curve"
        tokens = line.split()
        try:
            values = [float(x) for x in tokens]
        except ValueError:
            if line.strip():
                headers.append({"line": line_number, "text": line})
            continue
        if (mode == "points" and len(values) == 4) or (
            mode == "curve" and len(values) == 2
        ):
            item = {"line": line_number, "raw": line, "tokens": tokens}
            (points if mode == "points" else curve).append(item)
    assert len(points) == 22 and len(curve) == 100
    formulas = openpyxl.load_workbook(ARCHIVE / "coesite_final.xlsx", data_only=False)
    cached = openpyxl.load_workbook(ARCHIVE / "coesite_final.xlsx", data_only=True)
    cells: dict[str, dict[str, Any]] = {}

    def retain(sheet: str, address: str) -> dict[str, Any]:
        key = f"{sheet}!{address}"
        f, c = formulas[sheet][address], cached[sheet][address]
        cell = {
            "source": key,
            "formula": f.value if f.data_type == "f" else None,
            "value": c.value,
            "number_format": f.number_format,
        }
        cells[key] = cell
        # Keep direct upstream references, but not unrelated workbook contents.
        if f.data_type == "f" and f.value.startswith("="):
            import re

            for target_sheet, target_cell in re.findall(
                r"([A-Za-z0-9_]+)!\$?([A-Z]+\$?\d+)", f.value
            ):
                target = target_cell.replace("$", "")
                if f"{target_sheet}!{target}" not in cells:
                    retain(target_sheet, target)
        return cell

    for i, point in enumerate(points):
        point["point"] = i + 1
        if i < 14:
            row, sum_row = CURRENT_LIT_ROWS[i], CURRENT_SUM_ROWS[i]
            addresses = [f"{col}{row}" for col in ("V", "W", "X", "Y")]
            point["phase"] = "coesite-II" if i < 4 else "coesite-III"
            point["origin"] = "measured_current_study"
            point["sum_row"] = sum_row
            for col in ("B", "C", "P", "Q", "R", "S", "T", "U", "V", "W"):
                retain("SUM", f"{col}{sum_row}")
        else:
            row = i - 12
            addresses = [f"{col}{row}" for col in ("AA", "AB", "AC", "AD")]
            point["phase"] = "coesite-I"
            point["origin"] = (
                "derived_literature_ambient_eos_output"
                if i == 14
                else "measured_literature"
            )
            for col in ("N", "O", "P", "Q", "R", "S", "T"):
                retain("LIT", f"{col}{row}")
        point["workbook_cells_p_sigma_p_v_sigma_v"] = [
            retain("LIT", address)["source"] for address in addresses
        ]
    # Retain the three excluded records, not just a count or a rejection reason.
    for row in (7, 10, 11):
        cols = ("V", "W", "X", "Y") if row == 7 else ("AA", "AB", "AC", "AD")
        for col in cols:
            retain("LIT", f"{col}{row}")
    for col in ("B", "C", "Q", "R", "S", "T", "U", "V", "W"):
        retain("SUM", f"{col}8")
    return {
        "source_doi": "10.1038/s41467-018-07265-z",
        "received_from_user": "2026-10-25",
        "source_sha256": DIGESTS,
        "dat_headers": headers,
        "volume_basis": "angstrom^3 per SiO2 formula unit; pressure and sigma_P in GPa",
        "point_columns": [
            "pressure_gpa",
            "volume_a3_per_formula_unit",
            "sigma_p_gpa",
            "sigma_v_a3_per_formula_unit",
        ],
        "points": points,
        "calculated_curve": curve,
        "workbook_sheet_names": formulas.sheetnames,
        "workbook_cells": cells,
        "workbook_reading": "Original formulas and cached values retained; workbook not recalculated, refreshed, saved or edited.",
        "author_communicated_saved_curve": AUTHOR_CURVE,
        "author_communicated_replay": [34.2040, 102.73, 3.024],
        "author_communicated_objective": "sum[(P_model(V) - P_export)^2 / sigma_P_export^2], all three BM3 parameters free",
    }


def pressure(volume: Any, coefficients: Any) -> Any:
    import numpy as np

    v0, k0, kp = coefficients
    f = ((v0 / np.asarray(volume)) ** (2 / 3) - 1) / 2
    return 3 * k0 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (kp - 4) * f)


def fit_points(points: list[dict[str, Any]], weighted: bool = True) -> dict[str, Any]:
    import numpy as np
    from scipy.optimize import least_squares

    a = np.array([[float(t) for t in p["tokens"]] for p in points])
    p, v, sp = a[:, 0], a[:, 1], a[:, 2]
    fit = least_squares(
        lambda x: (pressure(v, x) - p) / (sp if weighted else 1),
        AUTHOR_CURVE,
        x_scale="jac",
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-11,
        max_nfev=3000,
    )
    assert fit.success
    r = pressure(v, fit.x) - p
    return {
        "parameters_v0_per_formula_unit_k0_gpa_kprime": fit.x.tolist(),
        "observations": len(points),
        "measured_rows": sum(
            p["origin"] != "derived_literature_ambient_eos_output" for p in points
        ),
        "pressure_rms_gpa": float(np.sqrt(np.mean(r**2))),
        "maximum_absolute_pressure_residual_gpa": float(np.max(np.abs(r))),
        "sum_squared_objective": float(2 * fit.cost),
        "residuals_gpa": r.tolist(),
        "objective": "pressure residuals divided by exported sigma_P"
        if weighted
        else "unweighted pressure residuals",
        "volume_errors_used_in_objective": False,
        "free_parameters": ["V0", "K0", "K0_prime"],
    }


def replay(source: dict[str, Any]) -> dict[str, Any]:
    import numpy as np
    from scipy.optimize import least_squares

    verify_archives()
    points = source["points"]
    cells = source["workbook_cells"]
    assert len(points) == 22
    for p in points:
        assert len(p["tokens"]) == 4
        exported = [float(t) for t in p["tokens"]]
        original = [cells[k]["value"] for k in p["workbook_cells_p_sigma_p_v_sigma_v"]]
        # Export coordinate precision plus legacy EoSFit sigma formatting.
        np.testing.assert_allclose(exported[0], original[0], rtol=0, atol=0.00051)
        np.testing.assert_allclose(exported[2], original[1], rtol=0, atol=0.0000051)
        decimals = 1 if p["origin"] == "measured_current_study" else 3
        rounded = float(
            Decimal(str(original[2])).quantize(
                Decimal(10) ** -decimals, rounding=ROUND_HALF_UP
            )
        )
        assert abs(exported[1] - rounded) < 1e-9
        if p["origin"] == "measured_current_study":
            rounded_sigma = float(
                Decimal(str(original[3])).quantize(
                    Decimal("0.01"), rounding=ROUND_HALF_UP
                )
            )
            assert abs(exported[3] - rounded_sigma) < 2e-6
            row = p["sum_row"]
            np.testing.assert_allclose(
                original,
                [cells[f"SUM!{c}{row}"]["value"] for c in ("B", "C", "S", "T")],
                rtol=0,
                atol=1e-12,
            )
    fit = fit_points(points)
    # Preserve demonstrated one-decimal current-study V/Z, while checking
    # underlying workbook precision before the shortened DAT plot export.
    precision_points = []
    for point in points:
        p, sp, v, sv = [
            cells[k]["value"] for k in point["workbook_cells_p_sigma_p_v_sigma_v"]
        ]
        if point["origin"] == "measured_current_study":
            v = float(Decimal(str(v)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
        precision_points.append(dict(point, tokens=[str(x) for x in (p, v, sp, sv)]))
    precision_fit = fit_points(precision_points)
    measured = [
        p for p in points if p["origin"] != "derived_literature_ambient_eos_output"
    ]
    measured_control = fit_points(measured)
    unweighted_control = fit_points(points, weighted=False)
    curve = source["calculated_curve"]
    a = np.array([[float(t) for t in p["tokens"]] for p in curve])
    decode = least_squares(
        lambda x: pressure(a[:, 1], x) - a[:, 0],
        AUTHOR_CURVE,
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-11,
        max_nfev=3000,
    )
    assert decode.success
    curve_checks = []
    fitted_coefficients = fit["parameters_v0_per_formula_unit_k0_gpa_kprime"]
    for row in curve:
        pp, vv = [float(t) for t in row["tokens"]]
        dp, dv = [
            float(Decimal(10) ** Decimal(Decimal(t).as_tuple().exponent) / 2)
            for t in row["tokens"]
        ]
        checks = {}
        for name, coefficients in (
            ("author_communicated_curve", AUTHOR_CURVE),
            ("weighted_replay", fitted_coefficients),
            (
                "workbook_precision_reconstruction",
                precision_fit["parameters_v0_per_formula_unit_k0_gpa_kprime"],
            ),
        ):
            lower = float(pressure(vv + dv, coefficients) - (pp + dp))
            upper = float(pressure(vv - dv, coefficients) - (pp - dp))
            checks[name] = {
                "pressure_residual_gpa": float(pressure(vv, coefficients) - pp),
                "coordinate_rounding_residual_interval_gpa": [lower, upper],
                "crosses_export_rounding_box": lower <= 0 <= upper,
            }
        curve_checks.append({"line": row["line"], **checks})
    from peritheos.eos.rt.bm import BM3

    native = BM3(
        V0=fitted_coefficients[0],
        K0=fitted_coefficients[1],
        K0_prime=fitted_coefficients[2],
    )
    volumes = np.array([float(p["tokens"][1]) for p in points])
    np.testing.assert_allclose(
        native.pressure(volumes),
        pressure(volumes, fitted_coefficients),
        rtol=1e-12,
        atol=1e-10,
    )
    half_printed_units = np.array([0.005, 0.5, 0.005])
    agreement = abs(np.array(fitted_coefficients) - PUBLISHED) <= half_printed_units
    assert all(agreement)
    source_residuals = pressure(volumes, AUTHOR_CURVE) - np.array(
        [float(p["tokens"][0]) for p in points]
    )
    return {
        "audit_date": "2026-10-25",
        "source_sha256": DIGESTS,
        "parameter_agreement": {
            "status": "parity",
            "scope": "parameter_agreement",
            "finding": FINDING,
            "central_coefficient_parity_achieved": bool(all(agreement)),
            "author_fit_values_reproduced_at_reported_precision": bool(
                np.all(
                    abs(
                        np.array(
                            precision_fit[
                                "parameters_v0_per_formula_unit_k0_gpa_kprime"
                            ]
                        )
                        - AUTHOR_CURVE
                    )
                    <= [0.00005, 0.005, 0.0005]
                )
            ),
            "basis": "All three replayed central coefficients round to the publication's printed central values; not a confidence/covariance comparison.",
            "published": PUBLISHED,
            "half_printed_units": half_printed_units.tolist(),
            "each_parameter_within_printing_interval": agreement.tolist(),
        },
        "author_workflow": {
            "status": "not_assessed",
            "scope": "author_workflow",
            "qualification": FINDING
            + " The recovered author input selection and communicated pressure weighting establish central-value reproduction. Original EoSFit settings and covariance were not supplied; uncertainty reproduction and every software-configuration detail remain unverified.",
        },
        "export_phase_counts": {
            phase: sum(p["phase"] == phase for p in points)
            for phase in ("coesite-I", "coesite-II", "coesite-III")
        },
        "author_protocol_pressure_weighted_replay": fit,
        "workbook_precision_reconstruction": {
            **precision_fit,
            "definition": "Underlying workbook P and sigma_P and literature V/Z; current-study V/Z rounded to one decimal as demonstrated by its number format and DAT export. This is a file-supported reconstruction of precision before plot export, not a recovered original EoSFit settings file.",
            "input_p_v_sigma_p_sigma_v": [
                [float(t) for t in p["tokens"]] for p in precision_points
            ],
        },
        "sensitivity_controls_not_replacement_eos": {
            "same_22_unweighted": unweighted_control,
            "same_pressure_weights_without_derived_anchor_21": measured_control,
        },
        "derived_anchor": {
            "export_point": 15,
            "dat_line": points[14]["line"],
            "export_p_v_sigma_p_sigma_v": [float(t) for t in points[14]["tokens"]],
            "workbook_origin": "LIT!N2:T2 (2014 fitted ambient values), copied into LIT!AA2:AD2",
            "qualification": "Earlier fitted EOS output, not an independent measured ambient observation. Tiny sigma_P effectively fixes V0 in the pressure-only replay; sigma_V is not propagated in that objective.",
        },
        "saved_calculated_curve": {
            "rows": len(curve),
            "decoded_bm3_coefficients_not_a_measurement_refit": decode.x.tolist(),
            "decoding_pressure_rms_gpa": float(np.sqrt(np.mean(decode.fun**2))),
            "decoding_maximum_absolute_pressure_residual_gpa": float(
                np.max(abs(decode.fun))
            ),
            "checks": curve_checks,
            "author_communicated_curve_compatible_rows": sum(
                c["author_communicated_curve"]["crosses_export_rounding_box"]
                for c in curve_checks
            ),
            "weighted_replay_compatible_rows": sum(
                c["weighted_replay"]["crosses_export_rounding_box"]
                for c in curve_checks
            ),
            "workbook_precision_reconstruction_compatible_rows": sum(
                c["workbook_precision_reconstruction"]["crosses_export_rounding_box"]
                for c in curve_checks
            ),
            "author_communicated_coefficients_on_22_export_rows": {
                "pressure_rms_gpa": float(np.sqrt(np.mean(source_residuals**2))),
                "maximum_absolute_pressure_residual_gpa": float(
                    np.max(abs(source_residuals))
                ),
                "residuals_gpa": source_residuals.tolist(),
            },
        },
        "qualifications": [
            "21 measured rows plus one derived ambient anchor, not 22 independent measurements.",
            "Workbook cached formulas were inspected without refreshing or editing the workbook.",
            "Current-study V/Z rounded to 0.1 A3; sigma_V rounded to 0.01 A3 in export.",
            "Published coefficients and original 24-row diagnostic retained; catalog integration is separate from this read-only scientific replay.",
            "Exclusion identities established by supplied files; individual scientific rejection reasons remain only partly recalled by the author.",
            "Repeated measurements and the prior-EOS ambient anchor are not established as statistically independent observations.",
        ],
    }


def author_rows() -> list[list[Any]]:
    """Literal DAT numbers plus explicit row-kind and source-address columns."""
    source = json.loads(INPUT.read_text())
    rows = []
    for point in source["points"]:
        p, v, sp, sv = map(float, point["tokens"])
        refs = point["workbook_cells_p_sigma_p_v_sigma_v"]
        rows.append(
            [
                point["point"],
                point["phase"],
                point["origin"],
                p,
                sp,
                v,
                sv,
                v * 16,
                sv * 16,
                point["line"],
                refs[0] + ":" + refs[-1].split("!")[1],
            ]
        )
    return rows


def bundled_csv_text() -> str:
    import csv
    import io

    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(
        [
            "export_point",
            "phase",
            "observation_origin",
            "pressure_gpa",
            "pressure_sigma_gpa",
            "volume_a3_per_formula_unit",
            "volume_sigma_a3_per_formula_unit",
            "volume_a3_z16_equivalent_cell",
            "volume_sigma_a3_z16_equivalent_cell",
            "dat_line",
            "workbook_range",
        ]
    )
    writer.writerows(author_rows())
    return stream.getvalue()


def bundled_dataset() -> dict[str, Any]:
    """Checksummed author inputs; measured and derived rows stay explicit."""
    columns = [
        {
            "name": "export_point",
            "quantity": "row_identifier",
            "unit": "dimensionless",
            "role": "flag",
        },
        {
            "name": "phase",
            "quantity": "phase_identifier",
            "unit": "text",
            "role": "flag",
        },
        {
            "name": "observation_origin",
            "quantity": "source_kind",
            "unit": "text",
            "role": "flag",
        },
    ]
    for name, quantity, unit, role, of in (
        ("pressure_gpa", "pressure", "GPa", "value", None),
        ("pressure_sigma_gpa", "pressure", "GPa", "standard_deviation", "pressure_gpa"),
        (
            "volume_a3_per_formula_unit",
            "formula_unit_volume",
            "angstrom^3/formula_unit",
            "value",
            None,
        ),
        (
            "volume_sigma_a3_per_formula_unit",
            "formula_unit_volume",
            "angstrom^3/formula_unit",
            "standard_deviation",
            "volume_a3_per_formula_unit",
        ),
        (
            "volume_a3_z16_equivalent_cell",
            "conventional_unit_cell_volume",
            "angstrom^3/conventional_unit_cell",
            "value",
            None,
        ),
        (
            "volume_sigma_a3_z16_equivalent_cell",
            "conventional_unit_cell_volume",
            "angstrom^3/conventional_unit_cell",
            "standard_deviation",
            "volume_a3_z16_equivalent_cell",
        ),
    ):
        col = {"name": name, "quantity": quantity, "unit": unit, "role": role}
        if of:
            col["of"] = of
        columns.append(col)
    columns += [
        {
            "name": "dat_line",
            "quantity": "source_line",
            "unit": "dimensionless",
            "role": "flag",
        },
        {
            "name": "workbook_range",
            "quantity": "source_cells",
            "unit": "text",
            "role": "flag",
        },
    ]
    return {
        "identifier": DATASET_ID,
        "kind": "pressure_volume",
        "description": "Original author-provided EoSFit Group 1 inputs: 21 measured entries and one explicitly derived literature ambient EOS anchor; not 22 independent measurements.",
        "reference": {
            "authors": ["Bykova", "et al."],
            "year": 2018,
            "source": "Nature Communications",
            "doi": "10.1038/s41467-018-07265-z",
        },
        "source_location": "Author-provided coeII-plot.dat, Group 1 lines 9-30; matched to coesite_final.xlsx LIT/SUM cells; received 2026-10-25.",
        "columns": columns,
        "resource": {
            "path": "datasets/" + CSV_PATH.name,
            "sha256": hashlib.sha256(bundled_csv_text().encode()).hexdigest(),
            "media_type": "text/csv",
        },
        "used_by_eos_records": [RECORD_ID],
        "uncertainty": {
            "reported_as": "EoSFit Sigma X and Sigma Y widths; original confidence convention and covariance not supplied.",
            "covariance": "not_reported",
            "fit_use": "Only exported sigma_P enters the communicated replay; sigma_V remains preserved, including the derived anchor's volume width.",
        },
        "provenance": {
            "type": "author_provided_fit_export",
            "received_on": "2026-10-25",
            "original_export": {
                "path": "datasets/coesite-i-iii-bykova-2018-author-eosfit-plot.dat",
                "sha256": DIGESTS["coeII-plot.dat"],
            },
            "original_workbook": {
                "name": "coesite_final.xlsx",
                "sha256": DIGESTS["coesite_final.xlsx"],
                "archived_repository_path": "output/bykova-2018-author-files-2026-10-25/coesite_final.xlsx",
            },
            "method": "Verbatim numeric DAT inputs retained; Z=16 equivalents are exact multiplications. Original DAT tokens and its separately labelled 100 calculated curve rows remain in the checksummed raw export. Workbook formulas, cached values and source addresses are retained in the dedicated audit transcription.",
            "derived_ambient_anchor": {
                "export_point": 15,
                "dat_line": 23,
                "origin": "Cernok (2014) Table 1 / thesis Table 5.1 fitted ambient output, copied through workbook LIT!N2:T2 to AA2:AD2.",
                "independently_measured": False,
            },
            "excluded_source_rows": [
                "coesite-III: workbook LIT!V7:Y7, P=35.9246 GPa, V/Z=26.15042083333333 A3",
                "literature coesite-II: 27.52 and 30.94 GPa, LIT!AA10:AD11",
            ],
            "audit_report": REPORT,
        },
        "notes": "Point 15 is derived EOS output and must not be presented as an independent measured ambient observation. Current-study V/Z was rounded to 0.1 A3 per formula unit; Sigma Y to 0.01 A3. Plot-export shortening and legacy floating-point tails are preserved. The 100 calculated curve rows in the raw DAT are not fit observations. Earlier 24 tabulated observations remain separately available, without being substituted for this author selection.",
    }


def ledger_outcome(record: dict[str, Any]) -> dict[str, Any]:
    """Reproduce the author's central fit values, with no reconstructed errors."""
    import numpy as np

    parameters = record.get("eos", {}).get("parameters")
    if parameters is None:
        parameters = {
            k: record["parameters"][k]["value"] for k in ("V0", "K0", "K0_prime")
        }
    assert parameters == {"V0": 547.2, "K0": 103.0, "K0_prime": 3.02}
    source = json.loads(INPUT.read_text())
    result = replay(source)
    fit = result["author_protocol_pressure_weighted_replay"]
    p = np.array([float(r["tokens"][0]) for r in source["points"]])
    v = np.array([float(r["tokens"][1]) for r in source["points"]])
    src = pressure(v, PUBLISHED) - p
    comparison = []
    for name, published, error, refit, tolerance in zip(
        ("V0", "K0", "K0_prime"),
        (547.2, 103, 3.02),
        (0.16, 2, 0.15),
        np.array(fit["parameters_v0_per_formula_unit_k0_gpa_kprime"]) * [16, 1, 1],
        (0.08, 0.5, 0.005),
    ):
        comparison.append(
            {
                "parameter": name,
                "published": published,
                "published_error": error,
                "refit": float(refit),
                "refit_error": None,
                "difference": float(refit - published),
                "relative_difference": float((refit - published) / published),
                "similar": True,
                "within_combined_2sigma": None,
                "published_precision_tolerance": tolerance,
                "within_published_precision": bool(abs(refit - published) <= tolerance),
            }
        )
    return {
        "status": "parity",
        "parity_basis": "published_decimal_precision",
        "comparison_scope": FINDING
        + " The recovered 22 author-export inputs include a derived ambient anchor. Uncertainty/covariance parity and every original EoSFit configuration detail remain unverified.",
        "central_coefficient_parity_achieved": result["parameter_agreement"][
            "central_coefficient_parity_achieved"
        ],
        "author_fit_values_reproduced_at_reported_precision": result[
            "parameter_agreement"
        ]["author_fit_values_reproduced_at_reported_precision"],
        "dataset_identifiers": [DATASET_ID],
        "observations": 22,
        "fit_kind": "isothermal_pv",
        "objective": "weighted_pressure_residuals",
        "selection": "All 22 original EoSFit Group 1 inputs: 21 measured entries plus point 15, a derived literature ambient EOS anchor.",
        "columns": {
            "pressure": "pressure_gpa",
            "volume": "volume_a3_z16_equivalent_cell",
            "temperature": None,
        },
        "weighting": "Sum[(P_BM3 - P_export)^2 / sigma_P_export^2]; sigma_V preserved but not used. Derived anchor sigma_P=0.0001 GPa effectively pins V0.",
        "absolute_sigma": True,
        "free_parameters": ["V0", "K0", "K0_prime"],
        "fixed_parameters": [],
        "observed_pressure_range_gpa": [float(min(p)), float(max(p))],
        "observed_volume_range": [float(min(v) * 16), float(max(v) * 16)],
        "parameters": comparison,
        "rmse_gpa": fit["pressure_rms_gpa"],
        "published_rmse_gpa": float(np.sqrt(np.mean(src**2))),
        "degrees_of_freedom": 19,
        "solver_success": True,
        "parameter_covariance": None,
        "covariance_convention": "Reconstruction errors/covariance and original confidence convention are not established. No combined-two-sigma or uncertainty-parity claim; absolute_sigma describes the declared pressure weights only.",
        "qualification": QUALIFICATION,
        "reproduction_report": REPORT,
        "original_author_workflow_parity": False,
        "uncertainty_parity_achieved": False,
        "measured_entries": 21,
        "derived_anchor_entries": 1,
        "author_input_provenance": {
            "source_sha256": DIGESTS,
            "derived_anchor": result["derived_anchor"],
            "excluded_points": bundled_dataset()["provenance"]["excluded_source_rows"],
        },
        "saved_source_curve_check": {
            k: result["saved_calculated_curve"][k]
            for k in (
                "rows",
                "author_communicated_curve_compatible_rows",
                "workbook_precision_reconstruction_compatible_rows",
            )
        },
        "historical_diagnostic": {
            "report": "docs/literature-reproductions/bykova-2018-other-eos-audit.md",
            "evidence": "docs/data/bykova-2018-other-eos-audit.json",
            "observations": 24,
            "status": "parity_not_achieved",
            "objective": "unweighted_pressure_residuals",
        },
    }


def plot(source: dict[str, Any], result: dict[str, Any]) -> None:
    import matplotlib
    import numpy as np
    from scipy.optimize import brentq

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, rx) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )
    values = np.array([[float(t) for t in p["tokens"]] for p in source["points"]])
    for phase, marker, color in (
        ("coesite-I", "o", "#377d63"),
        ("coesite-II", "s", "#755095"),
        ("coesite-III", "^", "#c97524"),
    ):
        mask = np.array(
            [
                p["phase"] == phase
                and p["origin"] != "derived_literature_ambient_eos_output"
                for p in source["points"]
            ]
        )
        ax.errorbar(
            values[mask, 0],
            values[mask, 1],
            xerr=values[mask, 2],
            yerr=values[mask, 3],
            fmt=marker,
            ms=4,
            color=color,
            label=phase + " measured",
        )
        rr = np.array(
            result["author_protocol_pressure_weighted_replay"]["residuals_gpa"]
        )[mask]
        rx.plot(values[mask, 0], rr, marker, ms=4, color=color)
    ax.errorbar(
        values[14, 0],
        values[14, 1],
        yerr=values[14, 3],
        fmt="*",
        ms=12,
        color="#286390",
        label="Derived ambient anchor (earlier EOS)",
    )
    curve = np.array(
        [[float(t) for t in p["tokens"]] for p in source["calculated_curve"]]
    )
    ax.plot(
        curve[:, 0],
        curve[:, 1],
        color="#555555",
        lw=3,
        alpha=0.5,
        label="Saved EoSFit curve (calculated)",
    )
    coefficients = result["author_protocol_pressure_weighted_replay"][
        "parameters_v0_per_formula_unit_k0_gpa_kprime"
    ]
    ps = np.linspace(0, 37.5, 180)
    vv = [brentq(lambda v: pressure(v, coefficients) - p, 25, 35) for p in ps]
    ax.plot(ps, vv, color="#222222", lw=1.2, label="22-row pressure-weighted replay")
    ax.set_ylabel("Volume per SiO2 (Å³)")
    source_residuals = result["saved_calculated_curve"][
        "author_communicated_coefficients_on_22_export_rows"
    ]["residuals_gpa"]
    measured_mask = np.array(
        [
            p["origin"] != "derived_literature_ambient_eos_output"
            for p in source["points"]
        ]
    )
    rx.plot(
        values[measured_mask, 0],
        np.array(source_residuals)[measured_mask],
        "x",
        ms=6,
        color="#777777",
        alpha=0.7,
        label="Author curve coefficients",
    )
    rx.plot([], [], "o", ms=4, color="#222222", label="Replay (phase colors)")
    rx.legend(fontsize=8)
    rx.axhline(0, color="#777777", lw=0.8)
    rx.set_ylabel("P(model) − P(export) (GPa)")
    rx.set_xlabel("Pressure (GPa)")
    ax.set_title("Bykova (2018): recovered 22-row fit export")
    ax.legend(fontsize=8)
    for panel in (ax, rx):
        panel.grid(alpha=0.2)
    fig.text(
        0.5,
        0.015,
        "21 measured points + one derived ambient anchor. Replay uses sigma_P only; displayed sigma_V is not used in fitting.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOT)
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    fig.savefig(ARCHIVE / "author-fit-check.png", dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.extract:
        source = extract()
        INPUT.write_text(json.dumps(source, indent=2, ensure_ascii=False) + "\n")
        print("Original-file transcription written; source files unchanged.")
        return
    source = json.loads(INPUT.read_text())
    result = replay(source)
    if args.write:
        RESULT.write_text(json.dumps(result, indent=2) + "\n")
    if args.check:
        from scripts.audit_bykova_2018_other_eos import compare

        compare(result, json.loads(RESULT.read_text()))
        print(
            "Author export replay and saved evidence agree; published central-value parity confirmed."
        )
    if args.plot:
        plot(source, result)
    if not any((args.write, args.check, args.plot)):
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
