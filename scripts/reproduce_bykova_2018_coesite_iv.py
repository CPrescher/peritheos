"""Replay the conditional coesite-IV refit and check independent figure readings.

Repository-only audit: ``python -m scripts.reproduce_bykova_2018_coesite_iv``.
No EOS coefficient or existing source table is written. Figure markers are
checks only; they never enter the one-parameter fit to Table 10.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
RECORD_ID = "coesite_iv_bykova_2018_am05_static_bm3_refit"
TABLE_ID = "coesite_iv_bykova_2018_table10_calc_pv"
FIGURE_ID = "coesite_iv_bykova_2018_figure6_digitized"
TABLE = ROOT / "peritheos/data/datasets/coesite-iv-bykova-2018-table10-calc-pv.csv"
FIGURE = ROOT / "peritheos/data/datasets/coesite-iv-bykova-2018-figure6-digitized.csv"
DIGITIZATION = ROOT / "docs/data/coesite-iv-bykova-2018-figure6-digitization.json"
RESULT = ROOT / "docs/data/coesite-iv-bykova-2018-audit.json"
REPORT = "docs/literature-reproductions/coesite-iv-bykova-2018.md"
PARAMETERS = {"V0": 438.6702, "K0": 168.52, "K0_prime": 3.34}
QUALIFICATION = (
    "Replays the stored Peritheos one-parameter refit, not a published V0 or "
    "the authors' E-P-V regression. K0=168.52 GPa and K0-prime=3.34 are "
    "adopted from Supplementary Table 11 and held fixed. Three calculated "
    "Table 10 pressures are printed as whole GPa; they are not exact states. "
    "The V0 formal standard error is conditional residual-scatter covariance, "
    "not a source uncertainty or a rounding/systematic error budget. The "
    "unweighted residual sum per degree of freedom is not a statistical "
    "reduced chi-square. Figure 6 supplies separate graphical checks, including "
    "five points beyond the three apparent table counterparts; none is fitted. "
    "See the dedicated coesite-IV audit for reading bounds and limitations."
)


def pressure(volume: Any, v0: float = PARAMETERS["V0"]) -> Any:
    """Independent standard BM3 expression; cell volume in A3, pressure in GPa."""
    eta = (v0 / np.asarray(volume)) ** (1.0 / 3.0)
    return (
        1.5
        * PARAMETERS["K0"]
        * (eta**7 - eta**5)
        * (1.0 + 0.75 * (PARAMETERS["K0_prime"] - 4.0) * (eta**2 - 1.0))
    )


def digitized_rows() -> list[dict[str, float | int]]:
    """Convert archived pixel readings, independently of the EOS curve."""
    source = json.loads(DIGITIZATION.read_text())
    p_ticks = np.array(source["pressure_ticks_gpa_and_x_pixel"])
    v_ticks = np.array(source["volume_ticks_a3_per_formula_unit_and_y_pixel"])
    p_map = np.polyfit(p_ticks[:, 1], p_ticks[:, 0], 1)
    v_map = np.polyfit(v_ticks[:, 1], v_ticks[:, 0], 1)
    bound = source["coordinate_bound_pixels"]
    rows = []
    for marker in source["markers"]:
        v = float(np.polyval(v_map, marker["y_pixel"]))
        rows.append(
            {
                "marker": marker["marker"],
                "figure_x_pixel": marker["x_pixel"],
                "figure_y_pixel": marker["y_pixel"],
                "pressure_gpa": float(np.polyval(p_map, marker["x_pixel"])),
                "volume_a3_per_formula_unit": v,
                "volume_a3_conventional_cell": 16 * v,
                "pressure_digitization_bound_gpa": float(abs(p_map[0]) * bound),
                "volume_digitization_bound_a3_per_formula_unit": float(
                    abs(v_map[0]) * bound
                ),
                "table10_counterpart": int(marker["table10_counterpart"]),
            }
        )
    return rows


def digitized_csv() -> str:
    rows = digitized_rows()
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def check_digitized_csv() -> None:
    """Keep archived pixels exact; compare their calculated coordinates numerically."""
    from scripts.bykova_replay_checks import check_report

    saved = list(csv.DictReader(io.StringIO(FIGURE.read_text())))
    current = digitized_rows()
    if len(saved) != len(current) or not saved:
        raise ValueError("Digitized CSV row count differs")
    identities = {"marker", "figure_x_pixel", "figure_y_pixel", "table10_counterpart"}
    for index, (left, right) in enumerate(zip(saved, current)):
        if left.keys() != right.keys():
            raise ValueError("Digitized CSV columns differ")
        typed = {
            key: int(value) if key in identities else float(value)
            for key, value in left.items()
        }
        check_report(typed, right, f"digitized_rows[{index}]", calculated=True)


def reproduce() -> dict[str, Any]:
    table = np.genfromtxt(TABLE, names=True, delimiter=",")
    p = table["pressure_gpa"]
    v = table["volume_a3_conventional_cell"]
    fit = least_squares(
        lambda x: pressure(v, x[0]) - p,
        [PARAMETERS["V0"]],
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
    )
    dof = len(p) - 1
    scatter_variance = float(np.sum(fit.fun**2) / dof)
    formal_error = float(np.sqrt(scatter_variance / np.sum(fit.jac**2)))
    intervals = [
        [
            brentq(lambda v0: pressure(vi, v0) - (pi + offset), 400, 500)
            for offset in (-0.5, 0.5)
        ]
        for pi, vi in zip(p, v)
    ]
    checks = []
    for row in digitized_rows():
        volume = row["volume_a3_conventional_cell"]
        dp = row["pressure_digitization_bound_gpa"]
        dv = 16 * row["volume_digitization_bound_a3_per_formula_unit"]
        model = float(pressure(volume))
        residual = model - row["pressure_gpa"]
        # Independent rectangular graphical coordinate bounds, not statistical
        # errors: use exact nonlinear pressure extrema, without quadrature.
        residual_bounds = [
            float(pressure(volume + dv) - (row["pressure_gpa"] + dp)),
            float(pressure(volume - dv) - (row["pressure_gpa"] - dp)),
        ]
        checks.append(
            {
                **row,
                "model_pressure_gpa": model,
                "pressure_residual_gpa": residual,
                "residual_interval_gpa": residual_bounds,
                "within_reading_bounds": residual_bounds[0] <= 0 <= residual_bounds[1],
            }
        )
    extra = [row for row in checks if not row["table10_counterpart"]]
    tighter_failures = []
    for row in checks:
        volume = row["volume_a3_conventional_cell"]
        dp = row["pressure_digitization_bound_gpa"] / 2
        dv = 8 * row["volume_digitization_bound_a3_per_formula_unit"]
        low = float(pressure(volume + dv) - (row["pressure_gpa"] + dp))
        high = float(pressure(volume - dv) - (row["pressure_gpa"] - dp))
        if not low <= 0 <= high:
            tighter_failures.append(row["marker"])

    def statistics(rows: list[dict[str, Any]]) -> dict[str, Any]:
        residuals = np.array([row["pressure_residual_gpa"] for row in rows])
        return {
            "points": len(rows),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
            "max_absolute_pressure_residual_gpa": float(np.max(abs(residuals))),
            "all_within_reading_bounds": all(r["within_reading_bounds"] for r in rows),
        }

    return {
        "record_identifier": RECORD_ID,
        "source_pdf_sha256": json.loads(DIGITIZATION.read_text())["source_pdf_sha256"],
        "table_sha256": hashlib.sha256(TABLE.read_bytes()).hexdigest(),
        "stored_coefficients_unchanged": PARAMETERS,
        "fit": {
            "observations": len(p),
            "volume_range_a3": [float(min(v)), float(max(v))],
            "objective": "unweighted_pressure_residual_least_squares",
            "free_parameters": ["V0"],
            "V0_a3": float(fit.x[0]),
            "formal_V0_standard_error_a3": formal_error,
            "residuals_gpa": fit.fun.tolist(),
            "rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
            "stored_curve_rmse_gpa": float(np.sqrt(np.mean((pressure(v) - p) ** 2))),
            "residual_variance_gpa2": scatter_variance,
            "degrees_of_freedom": dof,
            "solver_success": bool(fit.success),
        },
        "rounding_diagnostic": {
            "assumption": "Nearest-integer pressure rounding +/-0.5 GPa, with printed volumes and moduli held fixed; not a reported uncertainty model.",
            "per_row_V0_intervals_a3": intervals,
            "intersection_V0_interval_a3": [
                max(interval[0] for interval in intervals),
                min(interval[1] for interval in intervals),
            ],
            "confidence_interval": False,
        },
        "figure_check": {
            "coefficients_fitted_to_figure": False,
            "all_markers": statistics(checks),
            "additional_markers_only": statistics(extra),
            "points": checks,
            "two_pixel_sensitivity": {
                "markers_outside_bounds": tighter_failures,
                "qualification": "The agreement conclusion depends on the stated analyst reading bound. A tighter +/-2 pixel total bound excludes markers 3 and 4; no sub-resolution or exact-curve agreement is established.",
            },
            "qualification": "Agreement within analyst graphical reading bounds only; not author-fit reproduction or physical validation.",
        },
    }


def ledger_outcome(record: dict[str, Any]) -> dict[str, Any]:
    if record["eos"]["parameters"] != PARAMETERS:
        raise ValueError("Stored coesite-IV coefficients differ from the audited refit")
    result = reproduce()
    fit = result["fit"]
    stored = record["eos"]["parameters"]["V0"]
    difference = fit["V0_a3"] - stored
    reproduced = bool(fit["solver_success"] and abs(difference) <= 0.00005)
    return {
        "status": "parity" if reproduced else "parity_not_achieved",
        "parity_basis": "registered_peritheos_refit_reproduced",
        "comparison_scope": "Stored Peritheos V0 at its registered four-decimal precision; not a published coefficient.",
        "dataset_identifiers": [TABLE_ID],
        "observations": 3,
        "fit_kind": "stored_refit_reproduction",
        "objective": "pressure_residuals",
        "selection": "all three coesite-IV Table 10 rows; Figure 6 excluded",
        "columns": {
            "pressure": "pressure_gpa",
            "volume": "volume_a3_conventional_cell",
            "temperature": None,
        },
        "weighting": "Equal pressure weights; no source standard deviations supplied.",
        "absolute_sigma": False,
        "free_parameters": ["V0"],
        "fixed_parameters": ["K0", "K0_prime"],
        "fixed_values": {"K0": 168.52, "K0_prime": 3.34},
        "observed_pressure_range_gpa": [39.0, 48.0],
        "observed_volume_range": fit["volume_range_a3"],
        "parameters": [
            {
                "parameter": "V0",
                "published": stored,
                "published_error": None,
                "refit": fit["V0_a3"],
                "refit_error": fit["formal_V0_standard_error_a3"],
                "difference": difference,
                "relative_difference": difference / stored,
                "similar": reproduced,
                "within_combined_2sigma": None,
            }
        ],
        "rmse_gpa": fit["rmse_gpa"],
        "published_rmse_gpa": fit["stored_curve_rmse_gpa"],
        "degrees_of_freedom": fit["degrees_of_freedom"],
        "solver_success": fit["solver_success"],
        "covariance_convention": "Residual-scaled conditional one-parameter covariance. Formal V0 error excludes pressure rounding, fixed-modulus uncertainty, DFT systematics and model inadequacy; not a source error or confidence interval.",
        "qualification": QUALIFICATION,
        "reproduction_report": REPORT,
        "figure_check": result["figure_check"],
        "rounding_diagnostic": result["rounding_diagnostic"],
    }


def write_plot(result: dict[str, Any]) -> None:
    """Show the unchanged curve and separately digitized coordinate bounds."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    from scripts.report_plot_themes import dark_svg

    with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": RECORD_ID}):
        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        volume = np.linspace(22.3, 24.15, 200)
        axes[0].plot(
            pressure(16 * volume),
            volume,
            color="#225ea8",
            label="Unchanged Peritheos BM3",
        )
        for counterpart, color, label in [
            (False, "#b35b00", "Additional calculated markers (5)"),
            (True, "#666666", "Apparent Table 10 counterparts (3)"),
        ]:
            rows = [
                r
                for r in result["figure_check"]["points"]
                if bool(r["table10_counterpart"]) == counterpart
            ]
            p = np.array([r["pressure_gpa"] for r in rows])
            v = np.array([r["volume_a3_per_formula_unit"] for r in rows])
            residual = np.array([r["pressure_residual_gpa"] for r in rows])
            bounds = np.array([r["residual_interval_gpa"] for r in rows])
            axes[0].errorbar(
                p,
                v,
                xerr=[r["pressure_digitization_bound_gpa"] for r in rows],
                yerr=[r["volume_digitization_bound_a3_per_formula_unit"] for r in rows],
                fmt="D",
                markersize=5,
                capsize=3,
                color=color,
                label=label,
            )
            axes[1].errorbar(
                p,
                residual,
                yerr=np.array([residual - bounds[:, 0], bounds[:, 1] - residual]),
                fmt="D",
                markersize=5,
                capsize=3,
                color=color,
            )
            for r in rows:
                axes[0].annotate(
                    str(r["marker"]),
                    (r["pressure_gpa"], r["volume_a3_per_formula_unit"]),
                    xytext=(5, 7),
                    textcoords="offset points",
                    fontsize=8,
                )
        axes[0].set_title("Coesite-IV: independent Figure 6 check at 0 K")
        axes[0].set_ylabel("Volume per SiO₂ formula unit (Å³)")
        axes[0].legend(fontsize=9)
        axes[1].axhline(0, color="#225ea8", linewidth=1)
        axes[1].set_ylabel("P(model) − P(figure) (GPa)")
        axes[1].set_xlabel("Digitized pressure (GPa)")
        for ax in axes:
            ax.grid(alpha=0.2)
        fig.text(
            0.5,
            0.02,
            "Bars: ±4-pixel graphical coordinate bounds propagated through BM3.\nNo figure points fitted; bars are not source uncertainties or confidence intervals.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.065, 1, 1))
        path = (
            ROOT
            / "docs/literature-reproductions/assets/coesite-iv-bykova-2018-figure6-check.svg"
        )
        fig.savefig(path, metadata={"Date": None})
        path.with_name(path.stem + "-dark.svg").write_bytes(dark_svg(path.read_bytes()))
        preview = (
            ROOT / "output/coesite-iv-bykova-2018-audit-2026-10-23/figure6-check.png"
        )
        preview.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview, dpi=160)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.check:
        from scripts.bykova_replay_checks import check_report

        check_report(json.loads(RESULT.read_text()), result)
        check_digitized_csv()
    else:
        RESULT.write_text(text)
        FIGURE.write_text(digitized_csv())
    if args.plot:
        write_plot(result)
    print(json.dumps(result["figure_check"]["additional_markers_only"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
