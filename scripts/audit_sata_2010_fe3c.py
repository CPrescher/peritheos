"""Audit original Sata (2010) Fe3C observations and the 300 K BM3 curve.

Repository maintenance only. Published pressures and coefficients are never
replaced by the conditional fits or pressure-reduction sensitivities below.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts.reproduce_li_2002_fe3c import rows as li_rows

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
TOKENS = DATA / "cementite-sata-2010-table1-tokens.csv"
NORMALIZED = DATA / "cementite-sata-2010-table1.csv"
SELECTED = DATA / "cementite-sata-2010-selected-pv.csv"
SOURCE = DATA / "cementite-sata-2010-source.json"
PEAKS = DATA / "cementite-sata-2010-table4-peaks.csv"
REPORT = ROOT / "docs/data/sata-2010-fe3c-audit.json"
RECORD_ID = "cementite_nonmagnetic_sata_2010_bm3"
DATASET_ID = "cementite_sata_2010_table1"
SELECTED_ID = "cementite_sata_2010_selected_pv"
VR = 9.341 * 16  # Source atomic volume; 16 atoms per conventional Z=4 Fe3C cell.
PARAMETERS = np.array([290.0, 3.76])
SOURCE_ERRORS = [13.0, 0.18]
USAGE_RECOMMENDATION = "qualified_use"
RECOMMENDATION_REASON = (
    "Retain as a nondefault published EOS at 300 K and 27.2–186.6 GPa with "
    "explicit uncertainty/data-scatter qualifications. Published coefficients "
    "are reproduced conditionally within their reported error widths, but the "
    "source confidence convention and covariance are unspecified. Observation "
    "scatter exceeds some tabulated coordinate-error estimates. "
    "The unchanged curve has pressure RMS 1.106733 GPa and maximum "
    "absolute residual 2.591677 GPa on the joint 12-row selection. Six of nine "
    "Sata observations exceed their first-order combined P/V error widths; the "
    "largest ratio is 9.02 under an uncorrelated-coordinate assumption. These "
    "ratios are descriptive, not formal sigma significances. Conditional "
    "parameter parity and primary-source verification do not establish the "
    "complete experimental uncertainty budget or original-workflow reproduction. "
    "The recalculated-pressure refit remains a separate diagnostic."
)
QUALIFICATION = (
    "Conditional parameter agreement for the published Sata (2010) 300 K BM3 "
    "curve, using all nine Table 1 Fe3C states plus the three Li (2002) measured "
    "states above 25 GPa explicitly named in paragraph 15. The catalog retains "
    "V0=9.341 A^3/atom (149.456 A^3/cell), K0=290(13) GPa and K0-prime=3.76(18). "
    "The common comparison fixes the printed V0 and refines only K0/K0-prime; "
    "the separate literal mBM3 replay refines Pr/Kr/Kr-prime at the selected Vr. "
    "Original weights, confidence and covariance are unspecified. Diagnostic "
    "errors are RSS/dof-scaled and omit fixed-input and calibration uncertainty. "
    "Source MgO and Li ruby/neon pressure bases remain distinct; no pressure "
    "column is recalibrated. The printed MgO coefficients miss Sata pressures "
    "by up to 0.353 GPa, so the exact pressure reducer is unverified. This is "
    "neither author-workflow reproduction nor independent pressure-scale "
    "validation. Litasov's later 298 K Vinet recalibration is not this EOS. "
    + RECOMMENDATION_REASON
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_token(token):
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(?:\((\d+)\))?", token)
    if match is None:
        raise ValueError(f"Invalid original number: {token!r}")
    value, error = match.groups()
    return float(value), None if error is None else int(error) * 10 ** -len(
        value.partition(".")[2]
    )


def observations():
    with TOKENS.open(newline="", encoding="utf-8") as stream:
        raw = list(csv.DictReader(stream))
    normalized = []
    for index, row in enumerate(raw, 1):
        item = {
            "source_row": index,
            "run": int(row["run"][-1]),
            "measurement_id": int(row["measurement_id"]),
        }
        for field, name in [
            ("a_raw", "a_angstrom"),
            ("b_raw", "b_angstrom"),
            ("c_raw", "c_angstrom"),
            ("volume_raw", "volume_a3"),
            ("mgo_a_raw", "mgo_a_angstrom"),
            ("mgo_volume_raw", "mgo_volume_a3"),
            ("pressure_raw", "pressure_gpa"),
            ("vc_vfe_raw", "vc_vfe"),
        ]:
            item[name], item["error_" + name] = parse_token(row[field])
        item["temperature_k"] = 300.0
        normalized.append(item)
    return raw, normalized


def selected_rows():
    """Retain the exact source coordinates; Li pressures are not MgO pressures."""
    selected = []
    for row in observations()[1]:
        selected.append(
            {
                "source": 0,
                "source_row": row["source_row"],
                "run": row["run"],
                "measurement_id": row["measurement_id"],
                "pressure_gpa": row["pressure_gpa"],
                "error_pressure_gpa": row["error_pressure_gpa"],
                "volume_a3": row["volume_a3"],
                "error_volume_a3": row["error_volume_a3"],
                "temperature_k": 300.0,
            }
        )
    for row in li_rows()[1]:
        if row["kind"] == "measured" and row["pressure_gpa"] > 25:
            selected.append(
                {
                    "source": 1,
                    "source_row": int(row["row_id"].split("-")[-1]),
                    "run": None,
                    "measurement_id": None,
                    "pressure_gpa": row["pressure_gpa"],
                    "error_pressure_gpa": row["pressure_gpa_error_width"],
                    "volume_a3": row["volume_a3"],
                    "error_volume_a3": row["volume_a3_error_width"],
                    "temperature_k": 300.0,
                }
            )
    return selected


def csv_text(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def pressure(volume, parameters=PARAMETERS, v0=VR):
    eta = (v0 / np.asarray(volume)) ** (1 / 3)
    k0, kp = parameters
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def modified_pressure(volume, parameters, vr=VR):
    pr, kr, kp = parameters
    f = ((vr / np.asarray(volume)) ** (2 / 3) - 1) / 2
    return (1 + 2 * f) ** 2.5 * (
        pr + (3 * kr - 5 * pr) * f + (4.5 * kr * (kp - 4) + 17.5 * pr) * f**2
    )


def modified_bulk(volume, parameters, vr=VR):
    pr, kr, kp = parameters
    f = ((vr / np.asarray(volume)) ** (2 / 3) - 1) / 2
    s = 1 + 2 * f
    b, c = 3 * kr - 5 * pr, 4.5 * kr * (kp - 4) + 17.5 * pr
    h, dh = pr + b * f + c * f**2, b + 2 * c * f
    return s**2.5 * (5 * h + s * dh) / 3


def data_agreement(rows):
    """Describe measured-coordinate residuals without inferring confidence/covariance."""
    v, p, sp, sv = [
        np.array([r[name] for r in rows])
        for name in (
            "volume_a3",
            "pressure_gpa",
            "error_pressure_gpa",
            "error_volume_a3",
        )
    ]
    residual = pressure(v) - p
    volume_width = modified_bulk(v, [0, *PARAMETERS]) / v * sv
    combined_width = np.hypot(sp, volume_width)
    ratio = np.abs(residual) / combined_width
    own = np.array([r["source"] == 0 for r in rows])
    return {
        "status": "poor",
        "coordinates": "Published pressure evaluated at each original measured volume; P(model)-P(observed), GPa.",
        "combined_width_expression": "sqrt(sP^2 + (K(V)*sV/V)^2)",
        "qualification": "First-order propagation assuming uncorrelated P/V coordinate errors. Source confidence, coordinate/calibration covariance and fitted-parameter covariance are unavailable. No formal sigma significance, reduced chi-square or predictive confidence band is inferred; coefficients are held at their printed central values.",
        "rows": [
            {
                "source": row["source"],
                "source_row": row["source_row"],
                "pressure_gpa": float(pi),
                "pressure_residual_gpa": float(di),
                "reported_pressure_width_gpa": float(si),
                "volume_width_as_pressure_gpa": float(vi),
                "combined_coordinate_width_gpa": float(ci),
                "absolute_residual_over_pressure_width": float(abs(di) / si),
                "absolute_residual_over_combined_width": float(ri),
            }
            for row, pi, di, si, vi, ci, ri in zip(
                rows, p, residual, sp, volume_width, combined_width, ratio
            )
        ],
        "all_12": {
            **statistics(residual),
            "outside_combined_width": int(np.sum(ratio > 1)),
            "max_absolute_residual_over_combined_width": float(max(ratio)),
        },
        "sata_9": {
            **statistics(residual[own]),
            "outside_pressure_width": int(np.sum(np.abs(residual[own]) > sp[own])),
            "outside_combined_width": int(np.sum(ratio[own] > 1)),
        },
        "usage_recommendation": USAGE_RECOMMENDATION,
        "recommendation_reason": RECOMMENDATION_REASON,
    }


def statistics(residual):
    return {
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_gpa": float(np.max(np.abs(residual))),
        "residual_gpa": residual.tolist(),
    }


def linear_solution(design, observed):
    beta, _, rank, _ = np.linalg.lstsq(design, observed, rcond=None)
    if rank != design.shape[1]:
        raise ValueError("Rank-deficient source replay")
    residual = design @ beta - observed
    dof = len(observed) - design.shape[1]
    covariance = np.linalg.inv(design.T @ design) * float(residual @ residual) / dof
    return beta, covariance, dof


def bm3_fit():
    rows = selected_rows()
    v, p = [np.array([r[name] for r in rows]) for name in ("volume_a3", "pressure_gpa")]
    f = ((VR / v) ** (2 / 3) - 1) / 2
    design = (1 + 2 * f)[:, None] ** 2.5 * np.column_stack(
        [3 * f - 18 * f**2, 4.5 * f**2]
    )
    beta, covariance, dof = linear_solution(design, p)
    jacobian = np.array([[1, 0], [-beta[1] / beta[0] ** 2, 1 / beta[0]]])
    physical_covariance = jacobian @ covariance @ jacobian.T
    values = [float(beta[0]), float(beta[1] / beta[0])]
    return {
        "parameters": values,
        "standard_errors": np.sqrt(np.diag(physical_covariance)).tolist(),
        "covariance": physical_covariance.tolist(),
        "dof": dof,
        **statistics(pressure(v, values) - p),
    }


def modified_fit(rows, mode="equal_pressure"):
    v, p, sp, sv = [
        np.array([r[name] for r in rows])
        for name in (
            "volume_a3",
            "pressure_gpa",
            "error_pressure_gpa",
            "error_volume_a3",
        )
    ]
    f = ((VR / v) ** (2 / 3) - 1) / 2
    design = (1 + 2 * f)[:, None] ** 2.5 * np.column_stack(
        [1 - 5 * f + 17.5 * f**2, 3 * f - 18 * f**2, 4.5 * f**2]
    )
    if mode in ("equal_pressure", "pressure_width"):
        scale = sp if mode == "pressure_width" else np.ones_like(sp)
        beta, covariance, dof = linear_solution(design / scale[:, None], p / scale)
        jacobian = np.array(
            [[1, 0, 0], [0, 1, 0], [0, -beta[2] / beta[1] ** 2, 1 / beta[1]]]
        )
        values = np.array([beta[0], beta[1], beta[2] / beta[1]])
        covariance = jacobian @ covariance @ jacobian.T
    else:

        def loss(values):
            width = np.sqrt(sp**2 + (modified_bulk(v, values) * sv / v) ** 2)
            return (modified_pressure(v, values) - p) / width

        result = least_squares(loss, [0, 290, 3.76], xtol=1e-12, ftol=1e-12, gtol=1e-11)
        values = result.x
        dof = len(v) - 3
        covariance = (
            np.linalg.inv(result.jac.T @ result.jac)
            * float(result.fun @ result.fun)
            / dof
        )
    return {
        "parameters_Pr_Kr_Kprime": values.tolist(),
        "standard_errors": np.sqrt(np.diag(covariance)).tolist(),
        "covariance": covariance.tolist(),
        "dof": dof,
        **statistics(modified_pressure(v, values) - p),
    }


def reproduce():
    rows = selected_rows()
    own = observations()[1]
    v, p = [np.array([r[name] for r in rows]) for name in ("volume_a3", "pressure_gpa")]
    marker_v = np.array([r["mgo_volume_a3"] for r in own])
    marker_a = np.array([r["mgo_a_angstrom"] for r in own])
    replay = pressure(marker_v, [160.2, 3.99], 74.71)
    published_marker_pressure = p[:9]
    eta = (74.71 / marker_v) ** (1 / 3)
    base, slope = 1.5 * 160.2 * (eta**7 - eta**5), 0.75 * (eta**2 - 1)
    implied_kp = 4 + float(
        np.dot(base * slope, published_marker_pressure - base)
        / np.dot(base * slope, base * slope)
    )
    source_residual = pressure(v) - p
    with PEAKS.open(newline="", encoding="utf-8") as stream:
        peak_rows = list(csv.DictReader(stream))
    peak_checks = []
    for row in peak_rows:
        hkl = [int(row[k]) for k in ("h", "k", "l")]
        observed_d, error_d = parse_token(row["d_observed_raw"])
        calculated_d = 1 / np.sqrt(
            sum(
                (index / axis) ** 2
                for index, axis in zip(hkl, (3.94656, 4.53255, 5.95085))
            )
        )
        peak_checks.append(
            {
                "hkl_source_setting": hkl,
                "d_observed_angstrom": observed_d,
                "error_d_angstrom": error_d,
                "d_calculated_source_angstrom": float(row["d_calculated_raw"]),
                "delta_d_source_angstrom": float(row["delta_d_raw"]),
                "d_from_source_cell_angstrom": float(calculated_d),
            }
        )
    return {
        "record_identifier": RECORD_ID,
        "temperature_ref_k": 300,
        "published_parameters": {
            "V0_cell_a3": VR,
            "V0_atom_a3": 9.341,
            "K0_gpa": 290,
            "K0_prime": 3.76,
        },
        "published_reference_fit": {
            "Vr_cell_a3": VR,
            "Pr_gpa": 0.0,
            "Pr_error_gpa": 1.6,
            "Kr_gpa": 290,
            "Kr_error_gpa": 13,
            "Kprime_r": 3.76,
            "Kprime_r_error": 0.18,
        },
        "source_rows": own,
        "table4_diffraction_checks": peak_checks,
        "selection": rows,
        "source_replay": statistics(source_residual),
        "source_own_nine_replay": statistics(source_residual[:9]),
        "data_agreement": data_agreement(rows),
        "fixed_V0_bm3_diagnostic": bm3_fit(),
        "modified_bm3_replays": {
            mode: modified_fit(rows, mode)
            for mode in ("equal_pressure", "pressure_width", "effective_variance")
        },
        "source_only_modified_bm3": modified_fit(rows[:9]),
        "mgo_pressure_replay": {
            "parameters": {"V0": 74.71, "K0": 160.2, "K0_prime": 3.99},
            "published_pressure_gpa": published_marker_pressure.tolist(),
            "replayed_pressure_gpa": replay.tolist(),
            **statistics(replay - published_marker_pressure),
            "printed_marker_volume_minus_a_cubed_a3": (marker_v - marker_a**3).tolist(),
            "diagnostic_implied_Kprime": implied_kp,
            "qualification": "Implied Kprime is a one-coefficient pressure-reducer diagnostic, not recovered author input, a replacement MgO EOS or a recalibration of observations. Source and calculated pressures remain separate. No exact author-pressure replay is claimed.",
        },
        "later_litasov_comparison_not_integrated": {
            "model": "Vinet",
            "temperature_k": 298,
            "V0_cell_a3": 148.5,
            "K0_gpa": 316,
            "K0_prime": 3.5,
            "source": "Litasov (2013), Table 2; later recalibrated comparison, not original Sata coefficients",
        },
        "qualifications": QUALIFICATION,
        "source_files": {
            p.relative_to(ROOT).as_posix(): digest(p)
            for p in (
                TOKENS,
                PEAKS,
                ROOT / "docs/data/li-2002-fe3c/table1-source-tokens.csv",
            )
        },
    }


def check_manuscript(path):
    text = subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True)
    section = text[text.index("Fe3C #1") :]
    section = section[: section.index("Fe0.95O")]
    actual, run = [], None
    for line in section.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        if tokens[:1] == ["Fe3C"]:
            run = "Fe3C " + tokens[1]
            tokens = tokens[2:]
        if len(tokens) == 9 and tokens[0].isdigit():
            actual.append(
                dict(
                    zip(
                        [
                            "run",
                            "measurement_id",
                            "a_raw",
                            "b_raw",
                            "c_raw",
                            "volume_raw",
                            "mgo_a_raw",
                            "mgo_volume_raw",
                            "pressure_raw",
                            "vc_vfe_raw",
                        ],
                        [run, *tokens],
                    )
                )
            )
    expected = observations()[0]
    if actual != expected:
        raise ValueError(
            "Manuscript Fe3C tokens differ from the final-HTML factual transcription"
        )
    return {
        "rows": len(actual),
        "numerical_tokens": 8 * len(actual),
        "sha256": digest(path),
        "version": "2009 manuscript; corroborates Fe3C tokens only, not the final publication of record",
    }


def ledger_outcome(record):
    if record["identifier"] != RECORD_ID:
        raise ValueError("Wrong source record")
    fit = bm3_fit()
    comparisons = []
    for name, source, width, value, error in zip(
        ("K0", "K0_prime"),
        PARAMETERS,
        SOURCE_ERRORS,
        fit["parameters"],
        fit["standard_errors"],
    ):
        comparisons.append(
            {
                "parameter": name,
                "published": float(source),
                "published_error": width,
                "refit": value,
                "refit_error": error,
                "difference": value - float(source),
                "relative_difference": abs(value - source) / source,
                "within_reported_error": bool(abs(value - source) <= width),
                "within_combined_2sigma": None,
                "similar": bool(abs(value - source) <= width),
            }
        )
    if not all(p["within_reported_error"] and p["similar"] for p in comparisons):
        raise ValueError("Conditional source-width parameter agreement failed")
    return {
        "status": "parity",
        "parity_basis": "within_reported_parameter_errors",
        "conditional": True,
        "qualification": QUALIFICATION,
        "reason": QUALIFICATION,
        "parameters": comparisons,
        "dataset_identifiers": [SELECTED_ID],
        "observations": 12,
        "fit_kind": "isothermal_pv",
        "objective": "equal pressure residuals at measured volumes",
        "selection": "All nine Sata Table 1 Fe3C rows and all three measured Li (2002) Table 1 rows with P>25 GPa; original pressures retained.",
        "weighting": "Equal pressure weights; no source weights recovered",
        "absolute_sigma": False,
        "free_parameters": ["K0", "K0_prime"],
        "fixed_parameters": ["V0"],
        "fixed_values": {"V0": VR},
        "parameter_covariance": fit["covariance"],
        "degrees_of_freedom": fit["dof"],
        "solver_success": True,
        "covariance_convention": "RSS/(12-2) scaled inverse normal matrix, conditional on fixed printed V0; excludes calibration and predictor uncertainty. No source confidence or covariance inferred.",
        "rmse_gpa": fit["rmse_gpa"],
        "published_rmse_gpa": reproduce()["source_replay"]["rmse_gpa"],
        "coefficient_comparison_status": "parity",
        "data_agreement_status": "poor",
        "published_max_abs_pressure_residual_gpa": reproduce()["source_replay"][
            "max_abs_gpa"
        ],
        "pressure_rmse_coordinates": "Original measured volumes for published and diagnostic curves; joint 12-row selection.",
        "usage_recommendation": USAGE_RECOMMENDATION,
        "recommendation_reason": RECOMMENDATION_REASON,
        "observed_pressure_range_gpa": [27.2, 186.6],
        "observed_temperature_range_k": [300, 300],
        "author_workflow_reproduced": False,
        "reproduction_report": REPORT.relative_to(ROOT).as_posix(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write only audit data, not catalog records",
    )
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--manuscript", type=Path)
    args = parser.parse_args()
    result = reproduce()
    if args.manuscript:
        print(json.dumps(check_manuscript(args.manuscript), sort_keys=True))
    outputs = {
        NORMALIZED: csv_text(observations()[1]),
        SELECTED: csv_text(selected_rows()),
        REPORT: json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)
        + "\n",
    }
    if args.check:
        for path, payload in outputs.items():
            if not path.is_file() or path.read_text() != payload:
                raise SystemExit(f"Stale Sata Fe3C audit artifact: {path}")
        print("Sata (2010) Fe3C source/selection/fit audit checks passed")
    elif args.write:
        for path, payload in outputs.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload, encoding="utf-8")
    else:
        print(outputs[REPORT], end="")


if __name__ == "__main__":
    main()
