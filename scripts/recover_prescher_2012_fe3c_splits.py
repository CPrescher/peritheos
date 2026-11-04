"""Repository-only tests of exact printed Prescher (2012) split-fit recovery.

No published observation, parameter, or earlier diagnostic is overwritten.
Rounding feasibility is an optimistic sensitivity check, not recovered data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.recover_prescher_2012_fe3c import ROOT, TABLE, VECTORS, bm3, fit

STRAIN = ROOT / "peritheos/data/research/cementite/prescher-2012-strain-vectors.json"
REPORT = ROOT / "docs/data/prescher-2012-fe3c-split-recovery.json"
TARGETS = {"low": [145, 8.5, 3, 0.7], "high": [172, 5.1, 1, 0.1]}
HALF_PRINTED_UNIT = np.array([0.5, 0.05, 0.5, 0.05])


def matches_printed(values, target):
    """Nearest rounding to the explicitly printed digits, not error-bar overlap.

    Lower boundary included, upper excluded (positive half-up convention).
    These tests cannot identify the source's unprinted digits or rounding rule.
    """
    values, target = np.asarray(values), np.asarray(target)
    return (
        (values >= target - HALF_PRINTED_UNIT) & (values < target + HALF_PRINTED_UNIT)
    ).tolist()


def extract_strain(path):
    import pdfplumber

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != "a908c16fdd1754a447d3959a6bf4d2cb1b46b21dbf08888981c6ab400ab77c12":
        raise ValueError("Inspect the geometry before using a different source PDF")
    f_ticks = [196.304, 221.974, 247.641, 273.310]
    F_ticks = [556.180, 542.844, 529.598, 516.353, 503.108]
    ax = np.polyfit(f_ticks, [0, 0.02, 0.04, 0.06], 1)
    ay = np.polyfit(F_ticks, [150, 160, 170, 180, 190], 1)
    rows = []
    with pdfplumber.open(path) as doc:
        page = doc.pages[1]
        for rect in page.rects:
            if not (
                195 < rect["x0"] < 283
                and 496 < rect["top"] < 559
                and rect["fill"]
                and 1.4 < rect["width"] < 1.5
                and 1.4 < rect["height"] < 1.5
            ):
                continue
            x = (rect["x0"] + rect["x1"]) / 2
            y = (rect["top"] + rect["bottom"]) / 2
            ycaps = [
                line["top"]
                for line in page.lines
                if 1.2 < line["width"] < 1.4
                and line["height"] < 0.002
                and abs((line["x0"] + line["x1"]) / 2 - x) < 0.005
            ]
            xcaps = [
                line["x0"]
                for line in page.lines
                if 1.2 < line["height"] < 1.4
                and line["width"] < 0.002
                and abs((line["top"] + line["bottom"]) / 2 - y) < 0.005
                and abs(line["x0"] - x) < 1
            ]
            assert len(ycaps) == len(xcaps) == 2
            xs, ys = sorted(xcaps), sorted(ycaps)
            f, F = float(np.polyval(ax, x)), float(np.polyval(ay, y))
            sF = float((ys[1] - ys[0]) / 2 * abs(ay[0]))
            rows.append(
                {
                    "pdf_x": x,
                    "pdf_top": y,
                    "f": f,
                    "F_gpa": F,
                    "f_width": float((xs[1] - xs[0]) / 2 * ax[0]),
                    "F_width_gpa": sF,
                    "pdf_xcaps": xs,
                    "pdf_ycaps": ys,
                    "pressure_from_f_F_gpa": float(3 * f * (1 + 2 * f) ** 2.5 * F),
                    "pressure_width_if_F_bar_is_pressure_only_gpa": float(
                        3 * f * (1 + 2 * f) ** 2.5 * sF
                    ),
                }
            )
        rows.sort(key=lambda row: row["f"])
        assert len(rows) == 19
        # Ascending strain order is matched to Table S1, not to the 22 main markers.
        for i, row in enumerate(rows, 1):
            row["matching_table_s1_row"] = i
        guide = [
            curve
            for curve in page.curves
            if 195 < curve["x0"] < 285
            and 490 < curve["top"] < 560
            and curve["width"] > 60
            and curve["stroke"]
        ]
        assert len(guide) == 1 and len(guide[0]["pts"]) == 3
    return {
        "format": "peritheos.prescher-2012-strain-inset-vectors",
        "doi": "10.1103/PhysRevB.85.140402",
        "source_sha256": digest,
        "page": 2,
        "figure": "1 upper inset",
        "kind": "digitized_graphic_coordinates_and_error_bars_not_original_inputs",
        "f_ticks_pdf_x": f_ticks,
        "f_tick_values": [0, 0.02, 0.04, 0.06],
        "F_ticks_pdf_top": F_ticks,
        "F_tick_values_gpa": [150, 160, 170, 180, 190],
        "f_axis_affine": ax.tolist(),
        "F_axis_affine": ay.tolist(),
        "rows": rows,
        "gray_guide_polyline_pdf": guide[0]["pts"],
        "qualification": "Plotted widths have no established Gaussian/confidence meaning. Inverting the displayed F width assumes only P contributes and is a diagnostic of the graphic, not recovery of fit sigmas. f and F share coordinate dependence and are not generally independent. The gray polyline is not identified as a stored fit output.",
    }


def effective_variance_fit(p, v, sp, sv, v0):
    """Related-author 2015 workflow lead; not evidence of the 2012 protocol."""
    from scipy.optimize import least_squares

    def residual(b):
        h = 1e-20
        derivative = bm3(v + 1j * h, v0, *b).imag / h
        return (bm3(v, v0, *b) - p) / np.sqrt(sp**2 + (derivative * sv) ** 2)

    opt = least_squares(
        residual, [145, 8.5], x_scale="jac", ftol=1e-12, xtol=1e-12, gtol=1e-12
    )
    inverse = np.linalg.inv(opt.jac.T @ opt.jac)
    sse = float(opt.fun @ opt.fun)
    return {
        "parameters": {"V0": v0, "K0": float(opt.x[0]), "K0_prime": float(opt.x[1])},
        "residual_scaled_errors": np.sqrt(
            np.diag(inverse) * sse / (len(p) - 2)
        ).tolist(),
        "unscaled_inverse_normal_matrix": inverse.tolist(),
        "degrees_of_freedom": len(p) - 2,
        "sum_squared_objective_residuals": sse,
    }


def compare_fit(result, target, use_absolute):
    from scipy.stats import t

    central = [result["parameters"][key] for key in ("K0", "K0_prime")]
    scaled = np.array(result["residual_scaled_errors"][-2:])
    unscaled = np.sqrt(np.diag(result["unscaled_inverse_normal_matrix"]))[-2:]
    errors = {"residual_scaled_1sigma": scaled}
    if use_absolute:
        errors["input_widths_as_absolute_1sigma"] = unscaled
    comparisons = {}
    for label, widths in errors.items():
        for convention, factor in (
            ("", 1),
            ("_2sigma", 2),
            ("_95pct_t_half_width", t.ppf(0.975, result["degrees_of_freedom"])),
        ):
            values = central + (widths * factor).tolist()
            matches = matches_printed(values, target)
            comparisons[label + convention] = {
                "K0_K0prime_errorK0_errorK0prime": values,
                "matches_printed_digits": matches,
                "all_four_match": all(matches),
            }
    return comparisons


def strain_line_fits(strain, pressure):
    from scipy.optimize import curve_fit, least_squares

    f = np.array([row["f"] for row in strain["rows"]])
    F = np.array([row["F_gpa"] for row in strain["rows"]])
    sf = np.array([row["f_width"] for row in strain["rows"]])
    sF = np.array([row["F_width_gpa"] for row in strain["rows"]])
    results = []
    for limit in (22, 24):
        for side in ("low", "high"):
            mask = pressure < limit if side == "low" else pressure >= limit
            n = int(mask.sum())
            fs, ys = f[mask], F[mask]
            for mode in ("equal_F", "F_widths", "independent_f_F_widths"):
                if mode == "independent_f_F_widths":

                    def residual(b):
                        return np.r_[
                            (b[0] * (1 + 1.5 * (b[1] - 4) * b[2:]) - ys) / sF[mask],
                            (b[2:] - fs) / sf[mask],
                        ]

                    opt = least_squares(residual, np.r_[[145, 8.5], fs], x_scale="jac")
                    b = opt.x[:2]
                    covariance = np.linalg.inv(opt.jac.T @ opt.jac)[:2, :2]
                    ss = float(opt.fun @ opt.fun)
                else:
                    sigma = sF[mask] if mode == "F_widths" else np.ones(n)
                    b, covariance = curve_fit(
                        lambda f, k, kp: k * (1 + 1.5 * (kp - 4) * f),
                        fs,
                        ys,
                        p0=[145, 8.5],
                        sigma=sigma,
                        absolute_sigma=True,
                    )
                    ss = float(
                        sum(((b[0] * (1 + 1.5 * (b[1] - 4) * fs) - ys) / sigma) ** 2)
                    )
                result = {
                    "parameters": {"K0": float(b[0]), "K0_prime": float(b[1])},
                    "unscaled_inverse_normal_matrix": covariance.tolist(),
                    "residual_scaled_errors": np.sqrt(
                        np.diag(covariance) * ss / (n - 2)
                    ).tolist(),
                    "degrees_of_freedom": n - 2,
                }
                results.append(
                    {
                        "selection": side + "_" + str(limit),
                        "source_rows": (np.flatnonzero(mask) + 1).tolist(),
                        "objective": mode,
                        **result,
                        "comparisons": compare_fit(
                            result, TARGETS[side], mode != "equal_F"
                        ),
                    }
                )
    return results


def rounding_feasibility(p, v, v0, target):
    """Optimistic rectangle-distance minimization; never generates observations.

    Table P is printed to 0.1 GPa and V to 0.01 A3. Coordinates can move
    independently to any point in their rounding boxes, deliberately favoring
    recovery. Conditional covariance uses the nominal-coordinate Jacobian.
    Thus the output is a numerical error-floor estimate, not an original fit.
    """
    from scipy.optimize import least_squares

    def residual(b):
        low = bm3(v + 0.005, v0, *b) - (p + 0.05)
        high = bm3(v - 0.005, v0, *b) - (p - 0.05)
        return np.where(low > 0, low, np.where(high < 0, high, 0.0))

    center = np.asarray(target[:2])
    opt = least_squares(
        residual,
        center,
        bounds=(center - [0.5, 0.05], center + [0.5, 0.05]),
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
    )
    h = 1e-20
    jacobian = np.column_stack(
        [
            bm3(v, v0, opt.x[0] + 1j * h, opt.x[1]).imag / h,
            bm3(v, v0, opt.x[0], opt.x[1] + 1j * h).imag / h,
        ]
    )
    ss = float(opt.fun @ opt.fun)
    errors = np.sqrt(np.diag(np.linalg.inv(jacobian.T @ jacobian)) * ss / (len(p) - 2))
    return {
        "V0_fixed": v0,
        "central_coefficients_at_numerical_minimum": opt.x.tolist(),
        "minimum_rectangle_distance_sum_squares_gpa2": ss,
        "optimistic_residual_scaled_error_floor_estimates": errors.tolist(),
        "nominal_J_coordinate_approximation": "Nominal V used in information matrix; not a rigorous global covariance bound. No unprinted observations are inferred or saved.",
        "rounding_box_half_widths": {"pressure_gpa": 0.05, "volume_a3": 0.005},
    }


def legacy_pressure_replay(p, v, v0, pressure_width=0.1):
    """Reimplement the SciPy 0.10.1 wrapper's covariance semantics.

    Modern MINPACK is used; this is not execution of an archived SciPy binary.
    Source: scipy/optimize/minpack.py at tag v0.10.1, curve_fit body.
    """
    from scipy.optimize import leastsq

    def residual(b):
        return (bm3(v, v0, *b) - p) / pressure_width

    params, inverse, info, message, status = leastsq(
        residual, [145, 8.5], full_output=True
    )
    if status not in (1, 2, 3, 4):
        raise RuntimeError(message)
    ss = float(sum(info["fvec"] ** 2))
    covariance = inverse * ss / (len(p) - 2)
    return {
        "parameters_K0_K0prime": params.tolist(),
        "errors_from_legacy_curve_fit": np.sqrt(np.diag(covariance)).tolist(),
        "errors_from_unscaled_leastsq_matrix": np.sqrt(np.diag(inverse)).tolist(),
        "reduced_chi_square": ss / (len(p) - 2),
        "source_url": "https://github.com/scipy/scipy/blob/v0.10.1/scipy/optimize/minpack.py",
        "qualification": "SciPy 0.10.1 curve_fit always multiplies covariance by SSE/(N-k); it has no absolute_sigma argument. Replays those semantics using current MINPACK, not the original author executable.",
    }


def parameter_assessments(cases):
    """Assess central parameters separately from uncertainty/workflow recovery."""
    selected = {
        "low": ("table_19", "low_24", "constant_pressure_0.2"),
        "high": (
            "archive_own_shelx_19_with_table_pressures",
            "high_22",
            "constant_pressure_0.1",
        ),
    }
    qualifications = {
        "low": "Conditional on graph-derived fixed V0, inclusion of the 22.6 GPa row, and absolute 0.2 GPa pressure errors. The source prose says below 22 GPa; the figure indicates a change near 24 GPa. K0 uncertainty is not reproduced.",
        "high": "Conditional on graph-derived fixed V0, archived own/SHELX volumes paired with Table S1 pressures, and absolute 0.1 GPa pressure errors. The working p16 cell differs from the final table; source inputs and covariance convention are unconfirmed.",
    }
    assessments = {}
    for side, signature in selected.items():
        index, case = next(
            (i, row)
            for i, row in enumerate(cases)
            if (row["input"], row["selection"], row["objective"]) == signature
            and row["V0_basis"] == "drawn_full_curve_V0"
        )
        comparison = case["comparisons"]["input_widths_as_absolute_1sigma"]
        values = comparison["K0_K0prime_errorK0_errorK0prime"]
        matches = comparison["matches_printed_digits"]
        assessments[side] = {
            "status": "parity" if all(matches[:2]) else "parity_not_achieved",
            "scope": "parameter_agreement",
            "basis": "K0_and_K0_prime_at_published_printed_precision",
            "conditional": True,
            "compared_parameters": ["K0", "K0_prime"],
            "case_index": index,
            "case_signature": dict(zip(("input", "selection", "objective"), signature)),
            "qualifications": [qualifications[side]],
            "uncertainty_comparisons": [
                {
                    "parameter": key,
                    "published": TARGETS[side][j + 2],
                    "reconstructed": values[j + 2],
                    "matches_printed_precision": matches[j + 2],
                    "finding": "matched_at_printed_precision"
                    if matches[j + 2]
                    else "not_reproduced",
                }
                for j, key in enumerate(("K0", "K0_prime"))
            ],
            "author_workflow_reproduced": False,
        }
    return assessments


def replay():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    previous = json.loads(
        (ROOT / "docs/data/prescher-2012-fe3c-recovery.json").read_text()
    )
    graph = json.loads(VECTORS.read_text())
    strain = json.loads(STRAIN.read_text())
    refs = {
        "free_V0": None,
        "drawn_full_curve_V0": previous["source_drawn_curve_reconstruction"][
            "parameters"
        ]["V0"],
        "full_22_marker_pressure_fit_V0": previous["subsets"]["full_range"]["fits"][
            "figure_equal_pressure_free_V0"
        ]["parameters"]["V0"],
        "table_plus_ambient_pressure_fit_V0": previous["row_set_diagnostics"][
            "table_plus_figure_ambient"
        ]["parameters"]["V0"],
    }
    matched = [
        row for row in graph["markers"] if row["matching_table_s1_row"] is not None
    ]
    origins = json.loads(
        (
            ROOT / "peritheos/data/research/cementite/prescher-2012-archive.json"
        ).read_text()
    )["origin_projects"][0]["finite_observation_rows"]
    own_cifs = json.loads(
        (
            ROOT
            / "peritheos/data/research/cementite/prescher-2012-own-shelx-cells.json"
        ).read_text()
    )["rows"]
    inputs = {
        "table_19": (table["pressure_gpa"], table["volume_a3"]),
        "figure_matched_19": (
            np.array([r["pressure_gpa"] for r in matched]),
            np.array([r["volume_a3"] for r in matched]),
        ),
        "figure_all_22": (
            np.array([r["pressure_gpa"] for r in graph["markers"]]),
            np.array([r["volume_a3"] for r in graph["markers"]]),
        ),
        "archive_origin_10": (
            np.array([r["pressure_gpa"] for r in origins]),
            np.array([r["volume_a3"] for r in origins]),
        ),
        "archive_own_shelx_19_with_table_pressures": (
            table["pressure_gpa"],
            np.array([row["volume_a3"] for row in own_cifs]),
        ),
    }
    results = []
    for input_name, (pressure, volume) in inputs.items():
        for limit in (22, 24):
            for side in ("low", "high"):
                mask = pressure < limit if side == "low" else pressure >= limit
                if mask.sum() < 4:
                    continue
                p, v = pressure[mask], volume[mask]
                for ref_label, v0 in refs.items():
                    modes = [
                        "equal_pressure",
                        "equal_volume",
                        "constant_pressure_0.2",
                        "constant_pressure_0.1",
                    ]
                    if input_name == "table_19":
                        modes += [
                            "pressure_widths",
                            "volume_widths",
                            "coordinate_widths",
                        ]
                        if v0 is not None:
                            modes.append("effective_variance")
                    for mode in modes:
                        supplied = mode not in ("equal_pressure", "equal_volume")
                        sp = (
                            table["pressure_error_gpa"][mask]
                            if input_name == "table_19"
                            else None
                        )
                        sv = (
                            table["volume_error_a3"][mask]
                            if input_name == "table_19"
                            else None
                        )
                        actual_mode = mode
                        if mode.startswith("constant_pressure_"):
                            actual_mode = "pressure_widths"
                            sp = np.full(len(p), float(mode.rsplit("_", 1)[1]))
                        try:
                            if mode == "effective_variance":
                                result = effective_variance_fit(p, v, sp, sv, v0)
                            else:
                                result = fit(p, v, actual_mode, v0, sp, sv)
                        except (ValueError, RuntimeError) as exc:
                            results.append(
                                {
                                    "input": input_name,
                                    "selection": side + "_" + str(limit),
                                    "V0_basis": ref_label,
                                    "objective": mode,
                                    "solver_failure": str(exc),
                                    "comparisons": {},
                                    "qualification": "Failed trial/domain evaluation is not proof that no optimum exists; no scientific conclusion is assigned to this case.",
                                }
                            )
                            continue
                        results.append(
                            {
                                "input": input_name,
                                "selection": side + "_" + str(limit),
                                "row_indices_1based_in_named_input": (
                                    np.flatnonzero(mask) + 1
                                ).tolist(),
                                "V0_basis": ref_label,
                                "objective": mode,
                                "parameters": result["parameters"],
                                "degrees_of_freedom": result["degrees_of_freedom"],
                                "sum_squared_objective_residuals": result[
                                    "sum_squared_objective_residuals"
                                ],
                                "comparisons": compare_fit(
                                    result, TARGETS[side], supplied
                                ),
                            }
                        )
    exact = [
        {"case_index": i, "convention": name}
        for i, result in enumerate(results)
        for name, comparison in result["comparisons"].items()
        if comparison["all_four_match"]
    ]
    v0 = refs["drawn_full_curve_V0"]
    rounding = {}
    for side, mask in (
        ("low24", table["pressure_gpa"] < 24),
        ("high22", table["pressure_gpa"] >= 22),
    ):
        rounding[side] = rounding_feasibility(
            table["pressure_gpa"][mask],
            table["volume_a3"][mask],
            v0,
            TARGETS[side[:3] if side.startswith("low") else "high"],
        )
    line_results = strain_line_fits(strain, table["pressure_gpa"])
    return {
        "format": "peritheos.prescher-2012-split-recovery",
        "date": "2026-10-30",
        "targets_K0_K0prime_errorK0_errorK0prime": TARGETS,
        "comparison": "The strict full-result recovery test requires all four quantities to round to their printed digits. Central-parameter parity is assessed separately under parameter_agreement; an uncertainty mismatch does not erase central-value agreement. Source unprinted digits, confidence convention and rounding rule remain unknown.",
        "qualification": "Finite source-motivated protocol diagnostics, not an exhaustive search of arbitrary masks/weights and not recovered author settings. Constant pressure 0.2 comes from the inset bars; 0.1 additionally assumes those bars are two-sigma. Effective variance is a related-author 2015 numerical lead, not a recovered 2012 method. Inset independence assumptions are conditional.",
        "table_sha256": hashlib.sha256(TABLE.read_bytes()).hexdigest(),
        "strain_evidence_sha256": hashlib.sha256(STRAIN.read_bytes()).hexdigest(),
        "reference_volume_bases": refs,
        "cases": results,
        "parameter_assessments": parameter_assessments(results),
        "four_quantity_rounding_matches": exact,
        "strain_line_cases": line_results,
        "strain_line_four_quantity_matches": [
            i
            for i, r in enumerate(line_results)
            if any(c["all_four_match"] for c in r["comparisons"].values())
        ],
        "optimistic_rounding_sensitivity": rounding,
        "legacy_scipy_covariance_replays": {
            "low_table_through_22point6": legacy_pressure_replay(
                table["pressure_gpa"][table["pressure_gpa"] < 24],
                table["volume_a3"][table["pressure_gpa"] < 24],
                v0,
            ),
            "high_own_shelx_cells": legacy_pressure_replay(
                table["pressure_gpa"][table["pressure_gpa"] >= 22],
                np.array([row["volume_a3"] for row in own_cifs])[
                    table["pressure_gpa"] >= 22
                ],
                v0,
            ),
        },
        "author_workflow_reproduced": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract-pdf", type=Path)
    args = parser.parse_args()
    if args.extract_pdf:
        STRAIN.write_text(json.dumps(extract_strain(args.extract_pdf), indent=2) + "\n")
        print(STRAIN)
    else:
        report = replay()
        REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(REPORT)
        print("All-four rounding matches:", report["four_quantity_rounding_matches"])
