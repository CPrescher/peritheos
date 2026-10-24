"""Check the published joint IV/V EOS without replacing its coefficients.

The reference law is identified against the primary PDF's red vector curve.
The 16 Table 2 measurements and that fitted curve are separate datasets.
Unweighted refits are diagnostic: author weights and constraints are absent.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
TABLE = DATA / "coesite-iv-v-bykova-2018-table2-pv.csv"
VECTOR = DATA / "coesite-iv-v-bykova-2018-figure1-vector.json"
RESULT = ROOT / "docs/data/coesite-iv-v-bykova-2018-experimental.json"
REPORT = "docs/literature-reproductions/coesite-iv-v-bykova-2018-experimental.md"
RECORD_ID = "coesite_iv_v_bykova_2018_experimental_bm2_reference"
TABLE_ID = "coesite_iv_v_bykova_2018_table2_pv"
PARAMETERS = {"V0": 375.04, "K0": 254.0, "P0": 32.7}


def pressure(volume: Any, v0: float = 375.04, k0: float = 254.0) -> Any:
    """Independent expression, cell volume in A3, pressure in GPa."""
    eta = (v0 / np.asarray(volume)) ** (1 / 3)
    return 32.7 + 1.5 * k0 * (eta**7 - eta**5)


def source_curve() -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Sample the archived PDF path, not a newly calculated model curve."""
    source = json.loads(VECTOR.read_text())
    path = source["path"]
    start = np.array(path[1][1])
    control1, control2, end = map(np.array, path[2][1:])
    t = np.linspace(0, 1, 25)[:, None]
    xy = (
        (1 - t) ** 3 * start
        + 3 * (1 - t) ** 2 * t * control1
        + 3 * (1 - t) * t**2 * control2
        + t**3 * end
    )
    xy = np.vstack((path[0][1], xy, path[3][1]))
    xmap = np.polyfit(
        source["pressure_tick_x_pdf_points"], source["pressure_tick_gpa"], 1
    )
    ymap = np.polyfit(
        source["volume_tick_y_pdf_points"], source["volume_tick_a3_per_fu"], 1
    )
    p = np.polyval(xmap, xy[:, 0])
    v = 16 * np.polyval(ymap, xy[:, 1])
    # Source stroke half-width plus maximum axis-map residual in PDF points.
    xerror = max(
        abs(
            np.polyval(
                np.polyfit(
                    source["pressure_tick_gpa"], source["pressure_tick_x_pdf_points"], 1
                ),
                source["pressure_tick_gpa"],
            )
            - source["pressure_tick_x_pdf_points"]
        )
    )
    yerror = max(
        abs(
            np.polyval(
                np.polyfit(
                    source["volume_tick_a3_per_fu"],
                    source["volume_tick_y_pdf_points"],
                    1,
                ),
                source["volume_tick_a3_per_fu"],
            )
            - source["volume_tick_y_pdf_points"]
        )
    )
    bounds = {
        "pressure_bound_gpa": float(
            (source["linewidth_pdf_points"] / 2 + xerror) * abs(xmap[0])
        ),
        "cell_volume_bound_a3": float(
            16 * (source["linewidth_pdf_points"] / 2 + yerror) * abs(ymap[0])
        ),
        "meaning": "Half the source vector stroke width plus maximum tick-map residual, independently in both coordinates. Graphical bounds, not statistical errors. The PDF Bezier is a drawn approximation to the analytic EOS.",
    }
    mask = (p >= 29.9) & (p <= 70.4)
    return p[mask], v[mask], bounds


def reproduce() -> dict[str, Any]:
    with TABLE.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    v = np.array([float(r["volume_a3_conventional_cell"]) for r in rows])
    residuals = pressure(v) - p
    diagnostic = least_squares(
        lambda x: pressure(v, *x) - p, [375.04, 254.0], bounds=([300, 1], [500, 1000])
    )
    curve_p, curve_v, bounds = source_curve()
    curve_residual = pressure(curve_v) - curve_p
    low = (
        pressure(curve_v + bounds["cell_volume_bound_a3"])
        - curve_p
        - bounds["pressure_bound_gpa"]
    )
    high = (
        pressure(curve_v - bounds["cell_volume_bound_a3"])
        - curve_p
        + bounds["pressure_bound_gpa"]
    )
    return {
        "record_identifier": RECORD_ID,
        "published_parameters": PARAMETERS,
        "reference_volume_per_formula_unit_a3": 23.44,
        "reference_pressure_gpa": 32.7,
        "table_sha256": hashlib.sha256(TABLE.read_bytes()).hexdigest(),
        "curve_check": {
            "source_kind": "published_fitted_curve_vector_path",
            "observations": len(curve_p),
            "bounds": bounds,
            "rmse_gpa": float(np.sqrt(np.mean(curve_residual**2))),
            "max_absolute_residual_gpa": float(max(abs(curve_residual))),
            "all_within_graphical_bounds": bool(np.all((low <= 0) & (high >= 0))),
            "checkpoints": [
                {
                    "pressure_gpa": float(p),
                    "volume_a3_conventional_cell": float(v),
                    "residual_gpa": float(r),
                    "residual_interval_gpa": [float(a), float(b)],
                }
                for p, v, r, a, b in zip(curve_p, curve_v, curve_residual, low, high)
            ],
        },
        "observations_check": {
            "observations": len(rows),
            "coesite_iv_rows": sum(r["phase"] == "coesite_iv" for r in rows),
            "coesite_v_rows": sum(r["phase"] == "coesite_v" for r in rows),
            "rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
            "max_absolute_residual_gpa": float(max(abs(residuals))),
            "residuals_gpa": residuals.tolist(),
            "largest_residual_source_row": int(np.argmax(abs(residuals)) + 1),
            "interpretation": "All rows retained, including significant between-run scatter reported in Figure 1. These residuals are data agreement, not PDF curve-reading error or proof of author-regression reproduction. Row 13 is the 44.5 GPa coesite-V observation.",
        },
        "diagnostic_refit": {
            "objective": "Equal-weight pressure residuals on all 16 Table 2 rows, P0=32.7 GPa fixed; V0 and K0 free. This is an analyst diagnostic, not the recovered author protocol.",
            "parameters": {
                "V0": float(diagnostic.x[0]),
                "K0": float(diagnostic.x[1]),
                "P0": 32.7,
            },
            "rmse_gpa": float(np.sqrt(np.mean(diagnostic.fun**2))),
            "adopted_as_replacement": False,
        },
        "author_workflow_reproduced": False,
        "limitations": "Original row selection, objective, weights, fixed/free reference-volume choice and covariance are not specified. Printed errors have no stated confidence convention. No pressure-calibrant paired observations were recovered. Graphical source-equation agreement does not establish physical validation, author-fit parity or a recommended pure-phase EOS.",
    }


def ledger_outcome(record: dict[str, Any]) -> dict[str, Any]:
    if record["eos"]["parameters"] != PARAMETERS:
        raise ValueError("Published experimental coefficients must remain unchanged")
    result = reproduce()
    return {
        "status": "source_reconstruction",
        "dataset_identifiers": [TABLE_ID],
        "observations": 16,
        "fit_kind": "source_curve_validation_only",
        "reconstruction_kind": "finite_pressure_reference_bm2",
        "comparison_scope": "Published Figure 1 curve at graphical resolution, with observation residuals reported separately; not source-author coefficient parity.",
        "fixed_values": {"P0": 32.7},
        "fixed_parameters": ["P0"],
        "published_rmse_gpa": result["observations_check"]["rmse_gpa"],
        "curve_validation": result["curve_check"],
        "diagnostic_refit": result["diagnostic_refit"],
        "qualification": record["scientific_validation"]["note"],
        "reproduction_report": REPORT,
    }


def plot(result: dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    from scripts.report_plot_themes import dark_svg

    with TABLE.open(newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    data = {
        key: np.array(
            [row[key] if key == "phase" else float(row[key]) for row in source_rows]
        )
        for key in [
            "phase",
            "pressure_gpa",
            "pressure_error_gpa",
            "volume_a3_conventional_cell",
            "volume_error_a3_conventional_cell",
        ]
    }
    p, v, _ = source_curve()
    with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": RECORD_ID}):
        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        vv = np.linspace(334, 387, 200)
        axes[0].plot(
            pressure(vv),
            vv / 16,
            color="#d62728",
            label="Published BM2Reference (unchanged)",
        )
        axes[0].plot(p, v / 16, ".", color="black", label="Figure 1 red vector curve")
        diagnostic = result["diagnostic_refit"]["parameters"]
        axes[0].plot(
            pressure(vv, diagnostic["V0"], diagnostic["K0"]),
            vv / 16,
            "--",
            color="#555555",
            label="Equal-weight diagnostic (not adopted)",
        )
        for phase, label, marker, color in [
            ("coesite_iv", "Coesite-IV", "D", "#1f77b4"),
            ("coesite_v", "Coesite-V", "v", "#ff7f0e"),
        ]:
            mask = data["phase"] == phase
            pp, vol = (
                data["pressure_gpa"][mask],
                data["volume_a3_conventional_cell"][mask],
            )
            axes[0].errorbar(
                pp,
                vol / 16,
                xerr=data["pressure_error_gpa"][mask],
                yerr=data["volume_error_a3_conventional_cell"][mask] / 16,
                fmt=marker,
                color=color,
                capsize=2,
                label=label,
            )
            axes[1].plot(pp, pressure(vol) - pp, marker, color=color, label=label)
            axes[1].plot(
                pp,
                pressure(vol, diagnostic["V0"], diagnostic["K0"]) - pp,
                marker,
                markerfacecolor="none",
                markeredgecolor="#555555",
            )
        axes[0].set_ylabel("Volume per SiO₂ formula unit (Å³)")
        axes[0].set_title("Published joint coesite-IV/V experimental EOS")
        axes[0].legend(fontsize=8)
        axes[1].set_ylabel("P(model) − P(table) (GPa)")
        axes[1].set_xlabel("Pressure (GPa)")
        axes[1].axhline(0, color="black", linewidth=0.8)
        for ax in axes:
            ax.grid(alpha=0.2)
        fig.text(
            0.5,
            0.02,
            "Bars: printed coordinate errors; confidence unspecified. Lower: filled = published,\nopen grey = diagnostic residuals at measured volumes. All 16 rows retained.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.065, 1, 1))
        dest = (
            ROOT
            / "docs/literature-reproductions/assets/coesite-iv-v-bykova-2018-experimental.svg"
        )
        fig.savefig(dest, metadata={"Date": None})
        dest.with_name(dest.stem + "-dark.svg").write_bytes(dark_svg(dest.read_bytes()))
        fig.savefig(
            ROOT / "output/coesite-iv-v-experimental-2026-10-24/experimental-check.png",
            dpi=160,
        )
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    if args.check:
        from scripts.bykova_replay_checks import check_report

        check_report(json.loads(RESULT.read_text()), result)
    else:
        RESULT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if args.plot:
        plot(result)
    print(
        json.dumps(
            {
                "curve": {
                    k: v for k, v in result["curve_check"].items() if k != "checkpoints"
                },
                "observations": result["observations_check"]["rmse_gpa"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
