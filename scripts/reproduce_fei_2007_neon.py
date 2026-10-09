"""Ne-only Fei (2007) equation, Figure 5, and qualified fit reconstruction.

Run ``python -m scripts.reproduce_fei_2007_neon``. Published coefficients are
never replaced by diagnostics. Raster curve samples are validation outputs,
never measured fitting inputs. No Pt fitting or global regeneration is done.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/fei-2007-neon-reproduction.json"
SOURCE = DATA / "neon-fei-2007-source.json"
RECORDS = {"vinet": "neon_fcc_fei_2007_vinet_2", "bm3": "neon_fcc_fei_2007_bm3_1"}
# V0 [cell A3], K0 [GPa], K0', gamma0, q, theta0 [K].
PARAMETERS = {
    "vinet": (88.967, 1.16, 8.23, 2.05, 0.6, 75.1),
    "bm3": (88.967, 1.43, 8.02, 2.05, 0.6, 75.1),
}
NAMES = ("V0", "K0", "K0_prime", "gamma0", "q", "theta0")


def source_pressure(
    volume, temperature, branch="vinet", coefficients=None, law="printed"
):
    """Independent Eqs. 1-3 using SI molar volume and adaptive quadrature."""
    v0, k0, kp, gamma0, q, theta0 = (
        PARAMETERS[branch] if coefficients is None else coefficients
    )
    volume, temperature = np.broadcast_arrays(
        np.asarray(volume, float), np.asarray(temperature, float)
    )
    ratio = volume / v0
    gamma = gamma0 * ratio**q
    exponent = (
        -gamma * np.log(ratio)
        if law == "printed"
        else (
            -gamma0 * np.log(ratio)
            if abs(q) < 1e-10
            else -gamma0 * np.expm1(q * np.log(ratio)) / q
        )
    )
    theta = theta0 * np.exp(exponent)
    x = ratio ** (1 / 3)
    if branch == "vinet":
        cold = 3 * k0 * (1 - x) * np.exp(1.5 * (kp - 1) * (1 - x)) / x**2
    else:
        z = 1 / x
        cold = 1.5 * k0 * (z**7 - z**5) * (1 + 0.75 * (kp - 4) * (z**2 - 1))

    def energy(t, th):
        y = th / t
        integral = quad(lambda a: a**3 / np.expm1(a), 0, y, epsabs=1e-11, epsrel=1e-11)[
            0
        ]
        return 9 * R * t * integral / y**3

    increment = np.array(
        [energy(t, th) - energy(300, th) for t, th in zip(temperature.flat, theta.flat)]
    ).reshape(volume.shape)
    # Four atoms per conventional fcc cell, one atom per formula unit.
    molar_volume = volume * N_A * 1e-30 / 4
    result = cold + gamma * increment / molar_volume / 1e9
    return float(result) if result.ndim == 0 else result


def rows(name):
    with (DATA / name).open() as f:
        return list(csv.DictReader(f))


def observations():
    finger = rows("neon-finger-1981-table1.csv")
    hemley = rows("neon-hemley-1989-table1-compression.csv")
    cold = rows("neon-fei-2007-figure5-digitized.csv")
    hot = rows("neon-fei-2007-figure5-1000k-digitized.csv")
    combined = finger + hemley + cold + hot
    volume = np.array([float(r["volume_a3_conventional_cell"]) for r in combined])
    pressure = np.array(
        [float(r.get("pressure_gpa_dewaele_2004", r["pressure_gpa"])) for r in combined]
    )
    # Fei says all room-temperature data; exact treatment of Finger's 293 K
    # is undocumented. This diagnostic groups it at 300 K, with sensitivity below.
    temperature = np.r_[np.full(48, 300.0), np.full(6, 1000.0)]
    sp = np.r_[
        [float(r["pressure_uncertainty_kbar"]) / 10 for r in finger],
        [float(r["pressure_uncertainty_gpa_dewaele_2004"]) for r in hemley],
    ]
    sv = np.r_[
        [float(r["volume_uncertainty_a3_conventional_cell"]) for r in finger + hemley]
    ]
    return volume, temperature, pressure, sp, sv


def diagnostic_fit(
    volume,
    temperature,
    target,
    branch,
    free=(1, 2),
    law="printed",
    coefficients=None,
    sigma=None,
):
    base = np.array(PARAMETERS[branch] if coefficients is None else coefficients)
    indices = list(free)
    weights = np.ones_like(target) if sigma is None else sigma

    def unpack(x):
        c = base.copy()
        c[indices] = x
        return c

    fit = least_squares(
        lambda x: (
            (source_pressure(volume, temperature, branch, unpack(x), law) - target)
            / weights
        ),
        base[indices],
        bounds=(np.full(len(indices), -5.0), np.full(len(indices), 20.0)),
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=2000,
    )
    residual = source_pressure(volume, temperature, branch, unpack(fit.x), law) - target
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * (
        fit.fun @ fit.fun / (len(target) - len(indices))
    )
    errors = np.sqrt(np.diag(covariance))
    denominator = np.outer(errors, errors)
    correlation = np.divide(
        covariance, denominator, out=np.zeros_like(covariance), where=denominator != 0
    )
    return {
        "observations": len(target),
        "law": law,
        "free_parameters": [NAMES[i] for i in indices],
        "fixed_parameters": {
            NAMES[i]: float(base[i]) for i in range(6) if i not in indices
        },
        "parameters": dict(zip([NAMES[i] for i in indices], fit.x.tolist())),
        "conditional_residual_scaled_standard_errors": dict(
            zip([NAMES[i] for i in indices], errors.tolist())
        ),
        "conditional_parameter_correlation": correlation.tolist(),
        "published_rmse_gpa": float(
            np.sqrt(
                np.mean(
                    (source_pressure(volume, temperature, branch, base, law) - target)
                    ** 2
                )
            )
        ),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "solver_success": bool(fit.success),
        "uncertainty_note": "Conditional Jacobian/residual-scaled diagnostic errors only. No source covariance, common pressure-scale uncertainty, temperature uncertainty, or digitization axis correlations propagated.",
    }


def curve_check(source):
    """Sample clean curve strokes independently of any predicted coordinates."""
    metadata = source["figure5_hot_digitization"]
    image_path = ROOT / metadata["image_path"]
    assert (
        hashlib.sha256(image_path.read_bytes()).hexdigest() == metadata["image_sha256"]
    )
    pixels = np.asarray(Image.open(image_path).convert("RGB"), dtype=int)
    pc = np.polyfit(*np.array(metadata["axis_calibration_pixels"]["pressure"]).T, 1)
    vc = np.polyfit(*np.array(metadata["axis_calibration_pixels"]["volume"]).T, 1)
    # Visually inspected clean scanlines, clear of plotted markers/text. These
    # positions are selected from the raster, never by agreement with a model.
    scanlines = [550, 610, 650, 740, 910, 1080]
    result = []
    for y in scanlines:
        line = pixels[y]
        masks = {
            300: (line[:, 2] > 150) & (line[:, 0] < 90) & (line[:, 1] < 110),
            1000: (line[:, 0] > 150) & (line[:, 1] < 100) & (line[:, 2] < 100),
            2000: (line[:, 0] > 150) & (line[:, 2] > 150) & (line[:, 1] < 100),
        }
        for t, mask in masks.items():
            xx = np.flatnonzero(
                mask & (np.arange(len(mask)) > 350) & (np.arange(len(mask)) < 1650)
            )
            groups = [
                g for g in np.split(xx, np.where(np.diff(xx) > 2)[0] + 1) if len(g)
            ]
            assert len(groups) == 1, (y, t, groups)
            lo, hi = int(groups[0][0]), int(groups[0][-1])
            x = (lo + hi) / 2
            v = float(np.polyval(vc, y))
            p = float(np.polyval(pc, x))
            predicted = {
                b: {
                    law: float(source_pressure(v, t, b, law=law))
                    for law in ("printed", "integrated")
                }
                for b in PARAMETERS
            }
            sensitivity = (
                (hi - lo) / 2 * abs(pc[0])
                + abs(pc[0])
                + max(
                    abs(
                        source_pressure(v + sign * abs(vc[0]), t)
                        - predicted["vinet"]["printed"]
                    )
                    for sign in (-1, 1)
                )
            )
            result.append(
                {
                    "temperature_k": t,
                    "x_pixel": x,
                    "y_pixel": y,
                    "stroke_x_span": [lo, hi],
                    "volume_a3": v,
                    "figure_pressure_gpa": p,
                    "calculated_pressure_gpa": predicted,
                    "stroke_plus_one_pixel_axis_sensitivity_gpa": float(sensitivity),
                }
            )
    summaries = {}
    for b in PARAMETERS:
        summaries[b] = {}
        for law in ("printed", "integrated"):
            summaries[b][law] = {
                str(t): {
                    "samples": 6,
                    "rmse_gpa": float(
                        np.sqrt(
                            np.mean(
                                [
                                    (
                                        r["calculated_pressure_gpa"][b][law]
                                        - r["figure_pressure_gpa"]
                                    )
                                    ** 2
                                    for r in result
                                    if r["temperature_k"] == t
                                ]
                            )
                        )
                    ),
                    "max_abs_difference_gpa": float(
                        max(
                            abs(
                                r["calculated_pressure_gpa"][b][law]
                                - r["figure_pressure_gpa"]
                            )
                            for r in result
                            if r["temperature_k"] == t
                        )
                    ),
                }
                for t in (300, 1000, 2000)
            }
    separation = []
    for y in scanlines:
        selected = {r["temperature_k"]: r for r in result if r["y_pixel"] == y}
        for t in (1000, 2000):
            delta = (
                selected[t]["figure_pressure_gpa"]
                - selected[300]["figure_pressure_gpa"]
            )
            predicted = source_pressure(selected[t]["volume_a3"], t) - source_pressure(
                selected[t]["volume_a3"], 300
            )
            separation.append(
                {
                    "temperature_k": t,
                    "volume_a3": selected[t]["volume_a3"],
                    "figure_thermal_increment_gpa": float(delta),
                    "calculated_thermal_increment_gpa": float(predicted),
                    "difference_gpa": float(predicted - delta),
                }
            )
    return {
        "source_image": metadata["image_path"],
        "image_sha256": metadata["image_sha256"],
        "axis_calibration_pixels": metadata["axis_calibration_pixels"],
        "scanlines": scanlines,
        "points": result,
        "summary": summaries,
        "thermal_separation": separation,
        "qualification": "Calculated curve validation only, excluded from every fit. Stroke half-width plus illustrative one pixel per calibrated coordinate is a deterministic readout sensitivity, not a confidence interval. Figure 5 does not label its cold branch; Vinet and BM3 are tested separately. Raster curves cannot securely discriminate the two theta laws.",
    }


def reproduce():
    source = json.loads(SOURCE.read_text())
    v, t, p, sp, sv = observations()
    material = Material.from_eosmat(get_material_document("neon_fcc"))
    results = {}
    for branch in PARAMETERS:
        cold = diagnostic_fit(v[:48], t[:48], p[:48], branch)
        # Weighted comparison is restricted to the 35 table rows with source
        # measurement errors. Digitization readout bounds are not source errors.
        source_table_unweighted = diagnostic_fit(v[:35], t[:35], p[:35], branch)
        source_pressure_weighted = diagnostic_fit(
            v[:35], t[:35], p[:35], branch, sigma=sp
        )
        derivative = (
            source_pressure(v[:35] + 1e-4, 300, branch)
            - source_pressure(v[:35] - 1e-4, 300, branch)
        ) / 2e-4
        effective = diagnostic_fit(
            v[:35], t[:35], p[:35], branch, sigma=np.hypot(sp, derivative * sv)
        )
        observed_temp = t[:35].copy()
        observed_temp[:14] = 293
        observed_temperature_sensitivity = diagnostic_fit(
            v[:35], observed_temp, p[:35], branch, sigma=np.hypot(sp, derivative * sv)
        )
        hot = diagnostic_fit(v[48:], t[48:], p[48:], branch, free=(4,))
        hot_integrated = diagnostic_fit(
            v[48:], t[48:], p[48:], branch, free=(4,), law="integrated"
        )
        staged_base = np.array(PARAMETERS[branch])
        staged_base[1:3] = [cold["parameters"][key] for key in ("K0", "K0_prime")]
        staged = diagnostic_fit(
            v[48:], t[48:], p[48:], branch, free=(4,), coefficients=staged_base
        )
        joint = diagnostic_fit(v, t, p, branch, free=(1, 2, 4))
        comparisons = []
        for name, index, source_error in (
            ("K0", 1, 0.14),
            ("K0_prime", 2, 0.31),
            ("q", 4, 0.3),
        ):
            published = PARAMETERS[branch][index]
            fitted = joint["parameters"][name]
            fit_error = joint["conditional_residual_scaled_standard_errors"][name]
            difference = abs(fitted - published)
            relative = difference / abs(published)
            # Match the primary-refit ledger's K0/K0-prime/q similarity limits.
            similar = (
                relative <= 0.15
                if name == "K0"
                else difference <= 1.0 or relative <= 0.20
                if name == "K0_prime"
                else relative <= 0.25
            )
            comparisons.append(
                {
                    "parameter": name,
                    "published": published,
                    "published_error": source_error,
                    "refit": fitted,
                    "refit_error": fit_error,
                    "within_published_error_width": difference <= source_error,
                    "within_combined_2sigma": bool(
                        difference <= 2 * np.hypot(source_error, fit_error)
                    ),
                    "similar": similar,
                }
            )
        fit_status = (
            "parity"
            if all(c["within_combined_2sigma"] and c["similar"] for c in comparisons)
            else "similar"
            if all(c["within_combined_2sigma"] or c["similar"] for c in comparisons)
            else "parity_not_achieved"
        )
        # Uniform readout displacements are correlated deterministic sensitivity
        # cases, not six independent source error bars.
        shifted_q = []
        cal = source["figure5_hot_digitization"]["axis_calibration_pixels"]
        dp = (
            3
            * (cal["pressure"][1][1] - cal["pressure"][0][1])
            / (cal["pressure"][1][0] - cal["pressure"][0][0])
        )
        dv = 3 * abs(
            (cal["volume"][1][1] - cal["volume"][0][1])
            / (cal["volume"][1][0] - cal["volume"][0][0])
        )
        for pv in (-1, 1):
            for vv in (-1, 1):
                shifted_q.append(
                    diagnostic_fit(
                        v[48:] + vv * dv, t[48:], p[48:] + pv * dp, branch, free=(4,)
                    )["parameters"]["q"]
                )
        grid_v = np.array([24, 31, 39, 54, 88.967])[:, None]
        grid_t = np.array([293, 300, 1000, 2000])[None, :]
        expected = source_pressure(grid_v, grid_t, branch)
        actual = material.get_eos_record(RECORDS[branch]).pressure(grid_v, grid_t)
        results[branch] = {
            "record_identifier": RECORDS[branch],
            "status": fit_status,
            "parameter_comparisons": comparisons,
            "published_coefficients": dict(zip(NAMES, PARAMETERS[branch])),
            "cold_48_unweighted": cold,
            "source_tables_35_unweighted": source_table_unweighted,
            "source_tables_35_pressure_weighted": source_pressure_weighted,
            "source_tables_35_effective_error_weighted": effective,
            "source_tables_35_at_original_temperature_sensitivity": observed_temperature_sensitivity,
            "hot_6_q_only": hot,
            "hot_6_integrated_law_q_only_sensitivity": hot_integrated,
            "staged_48_cold_then_6_hot": staged,
            "joint_54_unweighted": joint,
            "hot_q_3px_correlated_readout_sensitivity": [
                float(min(shifted_q)),
                float(max(shifted_q)),
            ],
            "independent_equation_max_difference_gpa": float(
                np.max(abs(actual - expected))
            ),
            "equation_grid": {
                "volume_a3": grid_v[:, 0].tolist(),
                "temperature_k": grid_t[0].tolist(),
                "pressure_gpa": expected.tolist(),
            },
            "cold_individual_parameter_error_widths": {
                key: abs(cold["parameters"][key] - PARAMETERS[branch][i]) / error
                for key, i, error in (("K0", 1, 0.14), ("K0_prime", 2, 0.31))
            },
        }
    finger = rows("neon-finger-1981-table1.csv")
    hemley = rows("neon-hemley-1989-table1-compression.csv")
    mao_pressure = np.array([float(r["pressure_gpa"]) for r in hemley])
    wavelength_ratio = (1 + 7.665 * mao_pressure / 1904) ** (1 / 7.665)
    recalculated = 1904 / 9.5 * (wavelength_ratio**9.5 - 1)
    return {
        "doi": source["doi"],
        "input_sha256": {
            name: meta["sha256"] for name, meta in source["files"].items()
        },
        "status": {
            "fit_reproduction": "parity"
            if all(r["status"] == "parity" for r in results.values())
            else "see_branch_results",
            "reproduction_outcome": "reproduced"
            if all(r["status"] in {"parity", "similar"} for r in results.values())
            else "see_branch_results",
            "source_equations": "independently_checked",
            "curve_reproduction": "see_independent_raster_checks",
            "data_availability": "all_35_prior_table_rows_plus_19_separable_new_figure_markers",
            "exact_original_regression": "not_established",
            "absolute_pressure_accuracy": "not_established",
        },
        "scope": {
            "finger_1981_table_rows": 14,
            "hemley_1989_table_rows": 21,
            "new_fei_300k_digitized_rows": 13,
            "new_fei_1000k_digitized_rows": 6,
            "total_recovered_observations": 54,
            "complete_author_selection": False,
            "paired_new_fei_calibrant_volumes_recovered": False,
            "measured_temperature_range_k": [293, 1000],
            "2000k_curves": "calculated_extrapolation",
        },
        "normalization": {
            "atoms_per_formula_unit": 1,
            "formula_units_per_fcc_cell": 4,
            "cell_a3_to_molar_cm3_mol": N_A * 1e-24 / 4,
            "v0_molar_cm3_mol": 88.967 * N_A * 1e-24 / 4,
            "finger_adopted_v0_molar_cm3_mol": 13.394,
            "finger_v0_modern_conversion_cell_a3": 13.394 / (N_A * 1e-24 / 4),
            "finger_max_printed_vs_lattice_molar_difference_cm3_mol": max(
                abs(
                    float(r["volume_a3_conventional_cell"]) * N_A * 1e-24 / 4
                    - float(r["molar_volume_cm3_mol"])
                )
                for r in finger
            ),
        },
        "hemley_ruby_recalculation_max_rounding_difference_gpa": float(
            max(abs(recalculated - p[14:35]))
        ),
        "results": results,
        "curve_check": curve_check(source),
        "limitations": [
            "Fei original new observations, paired Pt/Au volumes, full source selection, weights, covariance and unrounded inputs not recovered.",
            "Cold unweighted fits can agree with individual published error widths without reproducing the original regression or joint uncertainty region.",
            "Only prior-table source errors permit source-error-weighted comparisons. Hot source errors and temperature errors are absent; no source-error-weighted hot fit is possible.",
            "Effective error weights use slopes of the published model, ignore shared calibration correlations and retain each table's original pressure/volume errors.",
            "Finger 293 K pressures are preserved on the original Piermarini scale; Fei explicitly recalibrates Hemley only. Assigning Finger to 300 K is a documented diagnostic choice.",
            "Curve agreement validates an equation/normalization at plotted conditions; it establishes neither exact original-regression recovery nor absolute pressure accuracy.",
            "Dewaele (2008) note29 discusses a 0 K reinterpretation and contains a self-referential citation typo. Fei Eq.3 and the Figure5 cold curve remain the authority for the retained 300 K subtraction.",
            "Published q error is 0.3; source confidence level and covariance are unspecified. No thermal diagnostic replaces published q=0.6.",
        ],
    }


def main():
    report = reproduce()
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")
    plot_validation(report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "results": {
                    b: {
                        k: r[k]
                        for k in (
                            "cold_48_unweighted",
                            "hot_6_q_only",
                            "joint_54_unweighted",
                            "independent_equation_max_difference_gpa",
                        )
                    }
                    for b, r in report["results"].items()
                },
                "curves": report["curve_check"]["summary"],
            },
            indent=2,
        )
    )


def plot_validation(report):
    """Show source strokes and independently calculated Ne curves/residuals."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curve = report["curve_check"]
    pc = np.polyfit(*np.array(curve["axis_calibration_pixels"]["pressure"]).T, 1)
    vc = np.polyfit(*np.array(curve["axis_calibration_pixels"]["volume"]).T, 1)
    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig = plt.figure(figsize=(16, 9), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.5, 1])
    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(plt.imread(ROOT / curve["source_image"]))
    volumes = np.linspace(24.0, 50.0, 600)
    for branch, style in (("vinet", "--"), ("bm3", ":")):
        for temperature in (300, 1000, 2000):
            pressures = source_pressure(volumes, temperature, branch)
            xx = (pressures - pc[1]) / pc[0]
            yy = (volumes - vc[1]) / vc[0]
            inside = (xx >= 190) & (xx <= 1660)
            ax.plot(
                xx[inside],
                yy[inside],
                style,
                color="black",
                linewidth=1.0,
                label=f"Published {branch.upper()} + MGD"
                if temperature == 300
                else None,
            )
    ax.set_xlim(0, 1683)
    ax.set_ylim(1275, 0)
    ax.axis("off")
    ax.set_title("Figure 5 with independent published-coefficient curves", fontsize=12)
    ax.legend(loc="lower right", framealpha=0.95, fontsize=9)
    residual_ax = fig.add_subplot(gs[0, 1])
    colors = {300: "#2525cc", 1000: "#c52525", 2000: "#c220bd"}
    for temperature, color in colors.items():
        selected = [r for r in curve["points"] if r["temperature_k"] == temperature]
        residual_ax.errorbar(
            [r["figure_pressure_gpa"] for r in selected],
            [
                r["calculated_pressure_gpa"]["vinet"]["printed"]
                - r["figure_pressure_gpa"]
                for r in selected
            ],
            yerr=[r["stroke_plus_one_pixel_axis_sensitivity_gpa"] for r in selected],
            fmt="o",
            color=color,
            capsize=2,
            markersize=4,
            label=f"{temperature} K",
        )
    residual_ax.axhline(0, color=".4", linewidth=0.8)
    residual_ax.set_xlabel("Pressure read from calculated source curve (GPa)")
    residual_ax.set_ylabel("Vinet calculation - source curve (GPa)")
    residual_ax.set_title("Curve checks; deterministic raster readout bounds")
    residual_ax.grid(alpha=0.15)
    residual_ax.legend(fontsize=9)
    hot_ax = fig.add_subplot(gs[1, 1])
    v, t, p, _, _ = observations()
    for branch, color, marker in (("vinet", "#126a84", "o"), ("bm3", "#bf6725", "s")):
        hot_ax.plot(
            p[48:],
            source_pressure(v[48:], t[48:], branch) - p[48:],
            marker + "-",
            color=color,
            label=f"Published {branch.upper()}",
            linewidth=1.0,
        )
    hot_ax.axhline(0, color=".4", linewidth=0.8)
    hot_ax.set_xlabel("Figure 5 hot marker pressure, published Pt scale (GPa)")
    hot_ax.set_ylabel("Calculated - digitized observation (GPa)")
    hot_ax.set_title("Six 1000 K observations; experimental errors unavailable")
    hot_ax.grid(alpha=0.15)
    hot_ax.legend(fontsize=9)
    fig.suptitle(
        "Solid Ne: Fei et al. (2007) equations, source curves and recovered observations",
        fontsize=15,
    )
    fig.supxlabel(
        "Vinet K0=1.16, K0'=8.23; BM3 K0=1.43, K0'=8.02. Common V0=88.967 A^3, gamma0=2.05, q=0.6, theta0=75.1 K; reference 300 K.\n"
        "Printed theta law. Curve samples are validation output, never fitting inputs. Source: DOI 10.1073/pnas.0609013104, Figure 5.",
        fontsize=9,
    )
    fig.savefig(ROOT / "docs/data/fei-2007-neon-figure5-validation.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
