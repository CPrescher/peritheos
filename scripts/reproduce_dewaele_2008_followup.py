"""Scoped source-figure audit; does not alter EOS records or calibrations.

Graph continuation is a diagnostic assumption, never recovered shock data or a
new material EOS. Run with ``python -m scripts.reproduce_dewaele_2008_followup``.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts.reproduce_dewaele_2008_metals import inverse, vinet

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/data/dewaele-2008-followup-source.json"
REPORT = ROOT / "docs/data/dewaele-2008-followup.json"
# Fixed compression shapes used ONLY in the explicitly extrapolative diagnostic.
# (K0 GPa, K0_prime, existing catalog record). Au's 2004 shape is an assumption:
# the Figure 2 Au comparison cites the later 2007 pressure-medium study.
SHAPES = {
    "Al": (74.3, 4.47, "aluminum_dewaele_2004_vinet_1"),
    "Cu": (132.4, 5.32, "copper_dewaele_2004_vinet_1"),
    "Ta": (194, 3.52, "tantalum_dewaele_2004_revised_ruby_vinet"),
    "Au": (167, 6.00, "gold_dewaele_2004_vinet_5"),
    "Ni": (183, 4.99, "nickel_dewaele_2008_table4_us_k0_183_vinet"),
    "Ag": (101, 5.97, "silver_dewaele_2008_table4_us_k0_101_vinet"),
    "Pt": (277, 5.08, "platinum_dewaele_2004_revised_ruby_vinet"),
    "Mo": (261, 4.06, "molybdenum_dewaele_2008_table4_us_k0_261_vinet"),
    "W": (295.2, 4.32, "tungsten_dewaele_2004_vinet_2"),
    "Co": (190, 4.37, "cobalt_hcp_dewaele_2008_table4_us_k0_190_vinet"),
}


def upgraded_pressure(p2004, a=1920.0, b=9.61):
    log_ratio = np.log1p(9.5 * np.asarray(p2004) / 1904.0) / 9.5
    return a / b * np.expm1(b * log_ratio)


def ruby_fit(p, target):
    fit = least_squares(
        lambda q: upgraded_pressure(p, *q) - target,
        [1904.0, 9.5],
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    fitted = upgraded_pressure(p, *fit.x)
    source = upgraded_pressure(p)
    return {
        "A_gpa": float(fit.x[0]),
        "B": float(fit.x[1]),
        "solver_success": bool(fit.success),
        "rms_pressure_residual_gpa": float(np.sqrt(np.mean((fitted - target) ** 2))),
        "published_ruby_rms_pressure_residual_gpa": float(
            np.sqrt(np.mean((source - target) ** 2))
        ),
        "published_vs_fitted_max_pressure_difference_gpa": float(
            np.max(abs(source - fitted))
        ),
        "coefficient_uncertainties": None,
        "author_workflow_reproduced": False,
    }


def fe_comparisons(source):
    outcomes = {}
    for name, c in source["figure4_fe_lda"]["curves"].items():
        rows = np.array(c["pressure_volume_points"])
        rows = rows[(rows[:, 0] >= 17) & (rows[:, 0] <= 197)]
        for state, q, shift in [
            ("static", [9.697, 334, 4.74], 0),
            ("298k", [9.726, 329, 4.75], 1),
        ]:
            p, v = rows.T
            dv = np.array([inverse(x + shift, q) for x in p]) - v
            dp = vinet(v, q) - (p + shift)
            outcomes[f"{name}_{state}"] = {
                "parameters_atomic_a3_gpa_dimensionless": q,
                "source_curve_index": c["curve_index"],
                "graph_pressure_range_gpa": [float(min(p)), float(max(p))],
                "source_thermal_pressure_shift_gpa": shift,
                "checkpoint_count": len(p),
                "volume_difference_range_a3_per_atom": [float(min(dv)), float(max(dv))],
                "max_absolute_volume_difference_a3_per_atom": float(max(abs(dv))),
                "pressure_rms_gpa": float(np.sqrt(np.mean(dp**2))),
                "within_0_04_a3_per_atom_graph_reading_bound": bool(
                    np.all(abs(dv) <= 0.04)
                ),
            }
    return {
        "comparisons": outcomes,
        "finding": "Both Table VI hcp-Fe LDA sets match the thin curve; they fail the thick curve assigned to with-semicore PAW by the source caption. Which source label or coefficient row is wrong remains unresolved.",
        "catalog_changes": False,
        "all_electron_check": source["figure4_fe_lda"]["all_electron_check"],
    }


def continuation(p, shape, adjustment):
    k0, kp = shape[:2]
    x = np.array([inverse(float(z), [1.0, k0, kp]) ** (1 / 3) for z in p])
    return p * np.expm1(adjustment[0] + 1.5 * adjustment[1] * (1 - x))


def input_scale_check():
    """Test both printed coordinate-rounding intervals, excluding He gauges."""
    path = ROOT / "peritheos/data/datasets/dewaele-2008-table2-complete.csv"
    with path.open() as stream:
        rows = [
            r
            for r in csv.DictReader(stream)
            if r["pressure_from_helium_eos"] == "False"
        ]

    def half_digit(token):
        return 0.5 * 10 ** (-len(token.split(".")[1])) if "." in token else 0.5

    def convert(p):
        return 1904 / 9.5 * np.expm1(9.5 / 7.665 * np.log1p(7.665 * p / 1904))

    exceptions = []
    for r in rows:
        p, target = float(r["pressure_mao_token"]), float(r["pressure_revised_token"])
        wp, wt = (
            half_digit(r["pressure_mao_token"]),
            half_digit(r["pressure_revised_token"]),
        )
        if not convert(p - wp) - wt <= target <= convert(p + wp) + wt:
            exceptions.append(
                {
                    "material": r["material"],
                    "run": int(r["run"]),
                    "mao_pressure_token": r["pressure_mao_token"],
                    "revised_pressure_token": r["pressure_revised_token"],
                    "computed_revised_pressure_gpa": float(convert(p)),
                }
            )
    return {
        "non_helium_material_rows": len(rows),
        "within_joint_printed_rounding_bounds": len(rows) - len(exceptions),
        "exceptions": exceptions,
        "qualification": "The two exception material rows share one Ni/Zn pressure pair. Source values are preserved; no correction or cause is inferred.",
    }


def crosscalibration(source):
    curves = {
        m: np.array(c["pressure_difference_points_gpa"])
        for m, c in source["figure2b"]["metals"].items()
    }
    common_max = min(x[-1, 0] for x in curves.values())
    p_common = np.linspace(0, common_max, 201)
    mean_dp = np.mean(
        [np.interp(p_common, xy[:, 0], xy[:, 1]) for xy in curves.values()], axis=0
    )
    common = ruby_fit(p_common, p_common + mean_dp)
    common.update(
        {
            "pressure_range_gpa": [0.0, float(common_max)],
            "objective": "fit the equal-metal mean on 201 uniformly spaced common pressures; linear interpolation; no extrapolation",
            "mean_correction_at_upper_limit_gpa": float(mean_dp[-1]),
        }
    )
    p_each, target_each = [], []
    for xy in curves.values():
        p = np.linspace(0, xy[-1, 0], 201)
        p_each.extend(p)
        target_each.extend(p + np.interp(p, xy[:, 0], xy[:, 1]))
    each = ruby_fit(np.array(p_each), np.array(target_each))
    each["objective"] = (
        "pool 201 uniformly spaced pressures over each metal own plotted range; equal total weight per metal; no extrapolation"
    )
    extensions, models = [], {}
    grid = np.linspace(0, 160, 321)
    for m, xy in curves.items():
        fit = least_squares(
            lambda q: continuation(xy[:, 0], SHAPES[m], q) - xy[:, 1],
            [0.0, 0.1],
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        ext = continuation(grid, SHAPES[m], fit.x)
        extensions.append(ext)
        models[m] = {
            "assumed_dac_shape": {
                "K0_gpa": SHAPES[m][0],
                "K0_prime": SHAPES[m][1],
                "catalog_record": SHAPES[m][2],
            },
            "fitted_graph_adjustment": {
                "log_modulus_ratio": float(fit.x[0]),
                "derivative_difference": float(fit.x[1]),
            },
            "plot_pressure_max_gpa": float(xy[-1, 0]),
            "graph_fit_max_residual_gpa": float(max(abs(fit.fun))),
            "correction_at_160_gpa": float(ext[-1]),
            "extrapolated_to_160_gpa": bool(xy[-1, 0] < 159.9),
        }
    average = np.mean(extensions, axis=0)
    full = ruby_fit(grid, grid + average)
    full.update(
        {
            "objective": "equal-metal mean of graph-derived Vinet-difference continuations, uniform 0–160 GPa grid, 321 pressures",
            "mean_correction_at_160_gpa": float(average[-1]),
            "models": models,
            "qualification": "Diagnostic extrapolation only. Uses assumed compression shapes; Au uses a 2004 shape although Figure 2 cites 2007. Fitted graph adjustments are not published shock coefficients. No uncertainties or precision are assigned outside plotted ranges.",
        }
    )
    return {
        "ten_metals": list(curves),
        "Zn_excluded": True,
        "input_scale_table2_check": input_scale_check(),
        "input_ruby": {"A_gpa": 1904.0, "B": 9.5, "reference": "2004 revised ruby"},
        "published_output_ruby": {
            "A_gpa": 1920.0,
            "B": 9.61,
            "correction_at_160_gpa": float(upgraded_pressure(160) - 160),
        },
        "common_range_mean": common,
        "individual_ranges_pooled": each,
        "full_range_graph_continuation": full,
        "finding": "The published ruby coefficients are not recovered uniquely by these explicit replay choices. Full author-workflow reproduction is unresolved; graph recovery is not independent reconstruction of the reduced-shock inputs.",
        "missing_inputs": [
            "original comparison grid or observation/pair selection",
            "cross-material and within-material weights",
            "handling of metals beyond their DAC ranges",
            "complete original reduced-shock calculations and covariance",
        ],
        "catalog_changes": False,
    }


def reproduce():
    source = json.loads(SOURCE.read_text())
    return {
        "format": "peritheos.dewaele-2008-followup",
        "source_pdf_sha256": source["source_pdf_sha256"],
        "fe_lda": fe_comparisons(source),
        "ruby_crosscalibration": crosscalibration(source),
    }


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    source = json.loads(SOURCE.read_text())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    ax = axes[0]
    for name, c in source["figure4_fe_lda"]["curves"].items():
        xy = np.array(c["pressure_volume_points"])
        ax.plot(
            xy[:, 0],
            xy[:, 1],
            label=name.replace("_", " "),
            lw=1.5 if "thin" in name else 2.8,
        )
    p = np.linspace(0, 197, 150)
    ax.plot(
        p,
        [inverse(x, [9.697, 334, 4.74]) for x in p],
        "--",
        color="black",
        label="Table VI static LDA",
    )
    ax.set(
        xlim=(0, 200),
        xlabel="Pressure (GPa)",
        ylabel="Atomic volume (Å³)",
        title="hcp-Fe: table matches the thin curve",
    )
    ax.legend(fontsize=8)
    ax = axes[1]
    c = report["ruby_crosscalibration"]
    for m, curve in source["figure2b"]["metals"].items():
        xy = np.array(curve["pressure_difference_points_gpa"])
        ax.plot(xy[:, 0], xy[:, 1], lw=0.7, alpha=0.35)
    grid = np.linspace(0, 160, 321)
    ax.plot(
        grid,
        upgraded_pressure(grid) - grid,
        color="black",
        lw=2,
        label="Published ruby: 1920, 9.61",
    )
    common = c["common_range_mean"]
    p = np.linspace(0, common["pressure_range_gpa"][1], 201)
    ax.plot(
        p,
        upgraded_pressure(p, common["A_gpa"], common["B"]) - p,
        lw=2,
        label="Equal-metal common range",
    )
    full = c["full_range_graph_continuation"]
    ax.plot(
        grid,
        upgraded_pressure(grid, full["A_gpa"], full["B"]) - grid,
        "--",
        lw=2,
        label="Extrapolated equal-metal mean",
    )
    ax.set(
        xlabel="2004 ruby pressure (GPa)",
        ylabel="Pressure correction (GPa)",
        title="Ten-metal replay depends on averaging",
    )
    ax.legend(fontsize=8)
    destination = (
        ROOT / "docs/literature-reproductions/assets/dewaele-2008-followup.png"
    )
    destination.parent.mkdir(exist_ok=True)
    fig.savefig(destination, dpi=180)
    plt.close(fig)
    return destination


if __name__ == "__main__":
    result = reproduce()
    REPORT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(REPORT)
    print(plot(result))
