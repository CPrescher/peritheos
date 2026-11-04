"""Scott (2001) diagnostic replay; no catalog writes or replacement EOS.

Repository maintenance: uv run --with pillow python -m scripts.audit_scott_2001.
Requires Poppler's pdfimages. SciPy calculations do not use Peritheos fitting.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.optimize import brentq, least_squares, linprog

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "output/cementite-eos-audit-2026-10-25/sources/scott-2001.pdf"
CSV = ROOT / "peritheos/data/datasets/cementite-scott-2001-figure3-digitized.csv"
OUT = ROOT / "output/scott-2001-parity-2026-10-30"
REPORT = ROOT / "docs/data/scott-2001-parity-investigation.json"
EXTRACTION = ROOT / "peritheos/data/research/cementite/scott-2001-figure3-audit.json"
V0 = 155.26
SOURCE = np.array([175.4, 5.1])


def pressure(volume, k0, kp, v0=V0):
    """GPa; volume and v0 use the same conventional-cell A^3 basis."""
    z = (v0 / np.asarray(volume)) ** (2 / 3)
    return 1.5 * k0 * (z**3.5 - z**2.5) * (1 + 0.75 * (kp - 4) * (z - 1))


def volume(p, k0, kp, v0=V0):
    return np.array(
        [
            brentq(lambda v: pressure(v, k0, kp, v0) - pp, v0 * 0.5, v0 * 1.2)
            for pp in np.atleast_1d(p)
        ]
    )


def linear_fit(v, p, sigma=None, v0=V0):
    """Exact linear-power BM3 solution for fixed V0 and pressure residuals."""
    z = (v0 / v) ** (2 / 3)
    a = 1.5 * (z**3.5 - z**2.5)
    design = np.column_stack([a, 0.75 * a * (z - 1)])
    if sigma is not None:
        design = design / sigma[:, None]
        p = p / sigma
    k0, b = np.linalg.lstsq(design, p, rcond=None)[0]
    return np.array([k0, 4 + b / k0])


def replay(v, p, *, objective="pressure", sp=None, sv=None, free_v0=False, v0=V0):
    def residual(a):
        k, q = a[:2]
        vv = a[2] if free_v0 else v0
        if objective == "volume":
            return volume(p, k, q, vv) - v
        if objective == "coordinate_widths":
            latent = a[3:] if free_v0 else a[2:]
            return np.concatenate(
                [(pressure(latent, k, q, vv) - p) / sp, (latent - v) / sv]
            )
        return (pressure(v, k, q, vv) - p) / (1 if sp is None else sp)

    initial = [175.4, 5.1] + ([v0] if free_v0 else [])
    npar = len(initial)
    if objective == "coordinate_widths":
        initial += list(v)
    fit = least_squares(
        residual,
        initial,
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-10,
        max_nfev=3000,
    )
    dof = len(fit.fun) - len(fit.x)
    _, singular, right = np.linalg.svd(fit.jac, full_matrices=False)
    keep = singular > np.finfo(float).eps * max(fit.jac.shape) * singular[0]
    scaled = right[keep, :npar] / singular[keep, None]
    covariance = np.einsum("ki,kj->ij", scaled, scaled)
    covariance *= np.sum(fit.fun**2) / dof
    std = np.sqrt(np.diag(covariance))
    parameters = fit.x[:npar]
    ref_v0 = parameters[2] if free_v0 else v0
    curve_p = np.linspace(0, 73.2, 1001)
    curve_delta = volume(curve_p, *parameters[:2], ref_v0) - volume(curve_p, *SOURCE)
    return {
        "parameters": dict(zip(["K0_gpa", "K0_prime", "V0_a3"], parameters.tolist())),
        "fixed_V0_a3": None if free_v0 else v0,
        "errors_residual_scaled": std.tolist(),
        "covariance_residual_scaled": covariance.tolist(),
        "correlation_K0_K0_prime": float(covariance[0, 1] / std[0] / std[1]),
        "dof": dof,
        "objective_sse": float(np.sum(fit.fun**2)),
        "pressure_rms_gpa": float(
            np.sqrt(np.mean((pressure(v, *parameters[:2], ref_v0) - p) ** 2))
        ),
        "source_pressure_rms_gpa": float(
            np.sqrt(np.mean((pressure(v, *SOURCE) - p) ** 2))
        ),
        "volume_rms_a3": float(
            np.sqrt(np.mean((volume(p, *parameters[:2], ref_v0) - v) ** 2))
        ),
        "max_curve_delta_a3_0_to_73p2_gpa": float(np.max(np.abs(curve_delta))),
        "central_pair_within_printed_widths": bool(
            np.all(np.abs(parameters[:2] - SOURCE) <= [3.5, 0.3])
        ),
        "optimizer_success": bool(fit.success),
        "observations": len(p),
        "objective": objective,
        "unscaled_inverse_normal_covariance": (
            (covariance * dof / np.sum(fit.fun**2)).tolist() if sp is not None else None
        ),
        "uncertainty_qualification": "Conditional inverse-normal-matrix errors. No source error confidence or inter-row covariance recovered; graphical widths are not established measurement sigmas.",
    }


def extract_pixels(rows):
    from PIL import Image

    OUT.mkdir(exist_ok=True, parents=True)
    subprocess.run(
        [
            "pdfimages",
            "-f",
            "3",
            "-l",
            "3",
            "-png",
            str(PDF),
            str(OUT / "page3-extracted"),
        ],
        check=True,
        capture_output=True,
    )
    image_path = OUT / "page3-extracted-000.png"
    native = Image.open(image_path)
    bitmap = np.asarray(native.convert("L")) < 128
    distance = distance_transform_edt(bitmap)
    points = []
    for i, row in enumerate(rows):
        x = round(283 + (1 - row["volume_ratio"]) * 2580)
        y = round(837 - row["pressure_gpa"] * 6.9)
        patch = distance[y - 5 : y + 6, x - 5 : x + 6]
        yy, xx = np.where(patch == patch.max())
        points.append(
            {
                "original_row_zero_based": i,
                "x": float(x - 5 + xx.mean()),
                "y": float(y - 5 + yy.mean()),
            }
        )
    # Visually identified major ticks, centres measured from native bitmap.
    # Top tick stems avoid the inset and overlapping data at the lower border.
    x_ticks = []
    for guess in [283, 412, 542, 672, 801, 931, 1061]:
        profile = bitmap[255:278, guess - 10 : guess + 11].sum(axis=0)
        x_ticks.append(float(np.where(profile == profile.max())[0].mean() + guess - 10))
    y_ticks = []
    for guess in [837, 768, 699, 630, 561, 492, 423, 354, 285]:
        profile = bitmap[guess - 7 : guess + 8, 283:299].sum(axis=1)
        y_ticks.append(float(np.where(profile >= 14)[0].mean() + guess - 7))
    ratio_axis = np.polyfit(x_ticks, np.linspace(1, 0.7, 7), 1)
    p_axis = np.polyfit(y_ticks, np.arange(0, 81, 10), 1)
    for point in points:
        point["pressure_gpa"] = float(np.polyval(p_axis, point["y"]))
        point["volume_ratio"] = float(np.polyval(ratio_axis, point["x"]))
    # No new observations are discovered by seed refinement; preserve that limit.
    evidence = {
        "source_pdf_sha256": hashlib.sha256(PDF.read_bytes()).hexdigest(),
        "original_csv_sha256": hashlib.sha256(CSV.read_bytes()).hexdigest(),
        "native_image_size": list(native.size),
        "native_image_bits_per_component": 1 if native.mode == "1" else 8,
        "pdf_page_zero_based": 2,
        "pdf_structure_inspection": "PyMuPDF inspection of the hash-identified source on 2026-10-30 found one one-bit image and zero drawings on each page. Native bitmap extraction for this replay uses Poppler pdfimages.",
        "method": "Distance-transform maxima in +/-5 pixel windows around the original 21 markers; average tied maxima. Multiple major-tick affine calibrations. Refinement of existing seeds, not an exhaustive independent point census or run assignment.",
        "axis_x_tick_centres": x_ticks,
        "axis_x_tick_ratios": np.linspace(1, 0.7, 7).tolist(),
        "axis_y_tick_centres": y_ticks,
        "axis_y_tick_pressures_gpa": list(range(0, 81, 10)),
        "ratio_from_x_polynomial": ratio_axis.tolist(),
        "pressure_from_y_polynomial": p_axis.tolist(),
        "remaining_systematics": "Raster skew, line/marker overlap, shared axis calibration, and uncertain glyph centres. These are digitization diagnostics, not measurement errors. Ambient marker is not forced to (0,1).",
        "points": points,
        "open_point": {
            "x": 468,
            "y": 726,
            "status": "Approximate visually identified quenched marker; diagnostic inclusion only. Caption identifies 16 GPa; fit membership is not explicitly stated.",
        },
    }
    traces = []
    windows = [(823, 839, 302, 320), (563, 575, 648, 653), (500, 517, 698, 701)]
    for x0, x1, y0, y1 in windows:
        for y in range(y0, y1):
            ids = np.where(bitmap[y, x0:x1])[0] + x0
            runs = [
                r for r in np.split(ids, np.where(np.diff(ids) > 1)[0] + 1) if len(r)
            ]
            if len(runs) != 1:
                continue
            run = runs[0]
            if run[0] == x0 or run[-1] == x1 - 1:
                continue
            x = float(run.mean())
            traces.append(
                {
                    "x": x,
                    "y": y,
                    "width_pixels": len(run),
                    "pressure_gpa": float(np.polyval(p_axis, y)),
                    "volume_ratio": float(np.polyval(ratio_axis, x)),
                }
            )
    evidence["drawn_curve_trace"] = {
        "windows_x0_x1_y0_y1": windows,
        "points": traces,
        "qualification": "Visually selected isolated portions of the continuous third-order line, avoiding stems and dashed BM2. Correlated drawing samples, not experimental observations; BM3 and Vinet are indistinguishable in the caption.",
    }
    return evidence


def identifiability(v, p, pixels_per_ratio=2580, pixels_per_gpa=6.9):
    """Local minimax stationarity test, not fitted or recovered observations.

    Minimum common pixel displacement bound that can make the source pair
    stationary for equal pressure residuals, allowing both coordinates to move.
    Include the derivative of J with respect to volume in the score derivative.
    """

    def score(vv, pp):
        z = (V0 / vv) ** (2 / 3)
        a = 1.5 * (z**3.5 - z**2.5)
        jac = np.column_stack(
            [pressure(vv, *SOURCE) / SOURCE[0], 0.75 * SOURCE[0] * a * (z - 1)]
        )
        return jac.T @ (pressure(vv, *SOURCE) - pp)

    base = score(v, p)
    columns = []
    for coordinate, step in [("x", -V0 / pixels_per_ratio), ("y", -1 / pixels_per_gpa)]:
        for i in range(len(p)):
            dv, dp = np.zeros(len(p)), np.zeros(len(p))
            (dv if coordinate == "x" else dp)[i] = step * 1e-3
            columns.append((score(v + dv, p + dp) - score(v - dv, p - dp)) / 2e-3)
    design = np.asarray(columns).T
    n = len(columns)
    # |delta coordinate| <= t in native pixels. All movements independent here.
    inequality = np.vstack(
        [
            np.column_stack([np.eye(n), -np.ones(n)]),
            np.column_stack([-np.eye(n), -np.ones(n)]),
        ]
    )
    result = linprog(
        np.r_[np.zeros(n), 1],
        A_ub=inequality,
        b_ub=np.zeros(2 * n),
        A_eq=np.column_stack([design, np.zeros(2)]),
        b_eq=-base,
        bounds=[(None, None)] * n + [(0, None)],
        method="highs",
    )
    assert result.success
    altered_v = v - result.x[: len(p)] * V0 / pixels_per_ratio
    altered_p = p - result.x[len(p) : 2 * len(p)] / pixels_per_gpa
    return {
        "minimum_common_pixel_bound_linearized": float(result.x[-1]),
        "nonlinear_fit_after_linearized_adversarial_movements": linear_fit(
            altered_v, altered_p
        ).tolist(),
        "qualification": "Local score linearization; adversarial coordinate movements, not experimental uncertainty, a probability, or recovered source observations. The ambient point is allowed to move. Independent movements omit shared calibration errors.",
    }


def run():
    with CSV.open() as stream:
        rows = [
            {k: float(v) if v else np.nan for k, v in row.items()}
            for row in csv.DictReader(stream)
        ]
    p = np.array([r["pressure_gpa"] for r in rows])
    v = np.array([r["volume_ratio"] for r in rows]) * V0
    sp = np.array([r["pressure_error_gpa"] for r in rows])
    sv = np.array([r["volume_ratio_error"] for r in rows]) * V0
    mask = np.isfinite(sp) & np.isfinite(sv)
    evidence = extract_pixels(rows)
    rp = np.array([r["pressure_gpa"] for r in evidence["points"]])
    rv = np.array([r["volume_ratio"] for r in evidence["points"]]) * V0
    cases = {
        "original_equal_pressure_fixed_v0": replay(v, p),
        "original_equal_volume_fixed_v0": replay(v, p, objective="volume"),
        "original_equal_pressure_free_v0": replay(v, p, free_v0=True),
        "original_equal_volume_free_v0": replay(v, p, objective="volume", free_v0=True),
        "resolved_bars_subset_equal_pressure": replay(v[mask], p[mask]),
        "resolved_bars_subset_pressure_widths": replay(v[mask], p[mask], sp=sp[mask]),
        "resolved_bars_subset_coordinate_widths": replay(
            v[mask], p[mask], objective="coordinate_widths", sp=sp[mask], sv=sv[mask]
        ),
        "refined_pixels_equal_pressure_fixed_v0": replay(rv, rp),
        "refined_pixels_equal_volume_fixed_v0": replay(rv, rp, objective="volume"),
        "omit_ambient_equal_pressure": replay(v[1:], p[1:]),
        "original_relative_pressure": replay(v[1:], p[1:], sp=p[1:]),
        "refined_relative_pressure": replay(rv[1:], rp[1:], sp=rp[1:]),
    }
    for vv in [155.12, 155.4]:
        cases[f"fixed_v0_{vv}"] = replay(v, p, v0=vv)
    for label, data_v, data_p in [("original", v, p), ("refined", rv, rp)]:
        op = evidence["open_point"]
        open_p = float(np.polyval(evidence["pressure_from_y_polynomial"], op["y"]))
        open_v = float(np.polyval(evidence["ratio_from_x_polynomial"], op["x"])) * V0
        cases[f"{label}_plus_quenched_equal_pressure"] = replay(
            np.r_[data_v, open_v], np.r_[data_p, open_p]
        )
    # BM2 and constrained-Kp tradeoff supply two extra printed cross-checks.
    kps = np.linspace(2.5, 7, 46)
    tradeoff = {}
    for label, vv, pp, widths in [
        ("original_equal_pressure", v, p, None),
        ("resolved_bars_pressure_widths", v[mask], p[mask], sp[mask]),
        ("original_relative_pressure", v[1:], p[1:], p[1:]),
        ("refined_relative_pressure", rv[1:], rp[1:], rp[1:]),
    ]:
        kk = []
        for kp in kps:
            basis = pressure(vv, 1, kp)
            weights = np.ones(len(pp)) if widths is None else 1 / widths**2
            kk.append(float(np.sum(weights * basis * pp) / np.sum(weights * basis**2)))
        slope, intercept = np.polyfit(kps, kk, 1)
        basis = pressure(vv, 1, 4)
        weights = np.ones(len(pp)) if widths is None else 1 / widths**2
        tradeoff[label] = {
            "bm2_K0_gpa": float(
                np.sum(weights * basis * pp) / np.sum(weights * basis**2)
            ),
            "K0_vs_fixed_Kp_linear_regression": {
                "intercept_gpa": float(intercept),
                "slope_gpa": float(slope),
            },
            "fixed_Kp_grid": kps.tolist(),
            "constrained_K0_gpa": kk,
        }
    curve_points = evidence["drawn_curve_trace"]["points"]
    cp = np.array([r["pressure_gpa"] for r in curve_points])
    cv = np.array([r["volume_ratio"] for r in curve_points]) * V0
    curve_check = {
        "samples": len(cp),
        "drawn_pressure_range_gpa": [float(cp.min()), float(cp.max())],
        "qualification": evidence["drawn_curve_trace"]["qualification"],
        "source_comparison": {},
        "conditional_curve_fit": replay(cv, cp, objective="volume"),
    }
    for label, pair in [
        ("published", SOURCE),
        ("original_equal_pressure", linear_fit(v, p)),
    ]:
        delta = (volume(cp, *pair) - cv) / V0 / evidence["ratio_from_x_polynomial"][0]
        curve_check["source_comparison"][label] = {
            "rms_x_pixels": float(np.sqrt(np.mean(delta**2))),
            "max_abs_x_pixels": float(np.max(np.abs(delta))),
        }
    f = 0.5 * ((V0 / v[1:]) ** (2 / 3) - 1)
    normalized = p[1:] / (3 * f * (1 + 2 * f) ** 2.5)
    slope, intercept = np.polyfit(f, normalized, 1)
    # Seeded uniform perturbations are sensitivity envelopes, not confidence intervals.
    rng = np.random.default_rng(2001)
    perturbations = {}
    for label, dx, dy in [
        ("one_native_pixel", V0 / 2580, 1 / 6.9),
        ("three_native_pixels", 3 * V0 / 2580, 3 / 6.9),
        ("archived_digitization_bounds", V0 * 0.0035, 0.8),
    ]:
        draws = np.array(
            [
                linear_fit(
                    v + rng.uniform(-dx, dx, len(v)), p + rng.uniform(-dy, dy, len(p))
                )
                for _ in range(2000)
            ]
        )
        perturbations[label] = {
            "coordinate_halfwidths_volume_a3_pressure_gpa": [dx, dy],
            "samples": len(draws),
            "parameter_min": draws.min(axis=0).tolist(),
            "parameter_max": draws.max(axis=0).tolist(),
            "central_95_percent_sensitivity_envelope": np.percentile(
                draws, [2.5, 97.5], axis=0
            ).tolist(),
            "qualification": "Independent uniform coordinate perturbations with fixed axes; illustrative sensitivity, not experimental confidence or a posterior probability.",
        }
    # Peritheos is imported only after all independent fits have been calculated.
    from peritheos import get_material_document
    from peritheos.eos.rt import BM3
    from peritheos.fitting import fit_rt_eos
    from scripts.validate_primary_eos_refits import _series

    document = get_material_document("cementite")
    record = next(
        r
        for r in document["eos_records"]
        if r["identifier"] == "cementite_scott_2001_bm3_1"
    )
    dataset = next(
        d
        for d in document["datasets"]
        if d["identifier"] == "cementite_scott_2001_figure3_digitized"
    )
    series = _series(document, record, dataset)
    fits = [
        fit_rt_eos(
            BM3,
            v,
            p,
            initial={"K0": 175.4, "K0_prime": 5.1},
            fixed={"V0": V0},
            absolute_sigma=flag,
        )
        for flag in [True, False]
    ]
    model_grid = np.linspace(V0 * 0.78, V0, 101)
    peritheos_check = {
        "generic_reader_pressure_sigma": None
        if series.pressure_sigma is None
        else series.pressure_sigma.tolist(),
        "generic_reader_volume_sigma": None
        if series.volume_sigma is None
        else series.volume_sigma.tolist(),
        "max_pressure_formula_difference_gpa": float(
            np.max(
                np.abs(
                    BM3(V0=V0, K0=175.4, K0_prime=5.1).pressure(model_grid)
                    - pressure(model_grid, *SOURCE)
                )
            )
        ),
        "linear_power_solution": linear_fit(v, p).tolist(),
        "absolute_sigma_true_parameters": fits[0].parameters,
        "absolute_sigma_true_covariance": np.asarray(fits[0].covariance).tolist(),
        "absolute_sigma_false_covariance": np.asarray(fits[1].covariance).tolist(),
    }
    result = {
        "runtime": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
        },
        "scope": "Diagnostic reproduction audit only; no new EOS or catalog assessment change.",
        "source_coefficients": {
            "V0_a3": V0,
            "K0_gpa": 175.4,
            "K0_prime_body": 5.1,
            "K0_prime_abstract": 5.2,
            "printed_errors": [0.14, 3.5, 0.3],
        },
        "source_fixed_v0_status": "Ambient volume measured independently; fixed/free constraint not explicitly stated in inspected final article.",
        "original_input_sha256": evidence["original_csv_sha256"],
        "source_pdf_sha256": evidence["source_pdf_sha256"],
        "resolved_bars_original_rows_zero_based": np.flatnonzero(mask).tolist(),
        "cases": cases,
        "peritheos_implementation_check": peritheos_check,
        "additional_published_fit_checks": tradeoff,
        "drawn_curve_equation_check": curve_check,
        "normalized_stress_equal_weight_diagnostic": {
            "K0_gpa": float(intercept),
            "K0_prime": float(4 + slope / (1.5 * intercept)),
            "qualification": "Equal F weights, not the authors' recovered regression.",
        },
        "visible_pressure_bar_relative_widths": {
            "minimum": float(np.min(sp[mask] / p[mask])),
            "maximum": float(np.max(sp[mask] / p[mask])),
            "median": float(np.median(sp[mask] / p[mask])),
            "qualification": "Motivates testing relative-pressure loss, without establishing that original errors were proportional to pressure.",
        },
        "published_tradeoff_gpa": {
            "intercept": 226.7,
            "slope": -9.97,
            "bm2_K0": 186.5,
            "grid_in_source": "Kp from 2.5 to 7; original spacing unknown",
        },
        "digitization_sensitivity": perturbations,
        "graphical_stationarity_test": identifiability(v, p),
    }
    EXTRACTION.write_text(json.dumps(evidence, indent=2) + "\n")
    REPORT.write_text(json.dumps(result, indent=2) + "\n")
    return result, evidence


def main():
    result, evidence = run()
    from matplotlib import pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), layout="constrained")
    p = np.linspace(0, 73.2, 301)
    for label, pair in [
        ("Published", SOURCE),
        (
            "Equal pressure, 21 points",
            list(
                result["cases"]["original_equal_pressure_fixed_v0"][
                    "parameters"
                ].values()
            ),
        ),
        (
            "Pressure bars, 13 points (diagnostic)",
            list(
                result["cases"]["resolved_bars_subset_pressure_widths"][
                    "parameters"
                ].values()
            ),
        ),
        (
            "Relative pressure, refined (diagnostic)",
            list(result["cases"]["refined_relative_pressure"]["parameters"].values()),
        ),
    ]:
        ax[0].plot(p, volume(p, *pair), label=label)
    raw = np.genfromtxt(CSV, names=True, delimiter=",")
    ax[0].scatter(
        raw["pressure_gpa"],
        raw["volume_ratio"] * V0,
        s=14,
        c="k",
        label="Original digitization",
        zorder=5,
    )
    ax[0].set(xlabel="Pressure (GPa)", ylabel="Conventional-cell volume (Å³)")
    ax[0].legend(fontsize=8)
    names = [
        "original_equal_pressure_fixed_v0",
        "original_equal_volume_fixed_v0",
        "resolved_bars_subset_pressure_widths",
        "resolved_bars_subset_coordinate_widths",
        "refined_pixels_equal_pressure_fixed_v0",
        "refined_relative_pressure",
    ]
    labels = [
        "21 points, pressure",
        "21 points, volume",
        "13 points, P widths",
        "13 points, P/V widths",
        "Refined pixels, pressure",
        "20 refined, relative pressure",
    ]
    for name, label in zip(names, labels):
        c = result["cases"][name]
        ax[1].plot(
            c["parameters"]["K0_gpa"], c["parameters"]["K0_prime"], "o", label=label
        )
    ax[1].errorbar(
        *SOURCE, xerr=3.5, yerr=0.3, fmt="k*", ms=12, label="Published ± printed widths"
    )
    ax[1].set(xlabel="K₀ (GPa)", ylabel="K₀′")
    ax[1].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.2)
    fig.suptitle("Scott (2001): conditional protocol sensitivity; no replacement EOS")
    fig.savefig(OUT / "protocol-comparison.png", dpi=180)
    plt.close(fig)
    from PIL import Image

    fig, ax = plt.subplots(figsize=(9, 7), layout="constrained")
    ax.imshow(Image.open(OUT / "page3-extracted-000.png"), cmap="gray")
    a, b = evidence["ratio_from_x_polynomial"]
    c, d = evidence["pressure_from_y_polynomial"]
    for label, pair, color in [
        ("Published BM3", SOURCE, "#0072b2"),
        (
            "Original equal-pressure replay",
            linear_fit(raw["volume_ratio"] * V0, raw["pressure_gpa"]),
            "#d55e00",
        ),
    ]:
        pp = np.linspace(0, 77.6, 400)
        xx = (volume(pp, *pair) / V0 - b) / a
        yy = (pp - d) / c
        ax.plot(xx, yy, color=color, lw=1.2, label=label)
    trace = evidence["drawn_curve_trace"]["points"]
    ax.scatter(
        [t["x"] for t in trace],
        [t["y"] for t in trace],
        s=12,
        c="#cc79a7",
        label="Separate drawing trace",
        zorder=5,
    )
    ax.set(
        xlim=(265, 1075),
        ylim=(850, 235),
        xlabel="Native image x (pixels)",
        ylabel="Native image y (pixels)",
        title="Scott Figure 3: source-image audit (native geometry)",
    )
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.82), fontsize=8)
    fig.savefig(OUT / "source-curve-overlay.png", dpi=180)
    plt.close(fig)
    print(REPORT.relative_to(ROOT))
    for name, c in result["cases"].items():
        print(name, c["parameters"])


if __name__ == "__main__":
    main()
