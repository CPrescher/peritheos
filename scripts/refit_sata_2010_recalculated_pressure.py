"""Separate Fe3C diagnostic: recalculated Sata P plus unchanged Li points.

Equal-pressure-weight BM3 at the previous diagnostic fixed V0, with the
original-reference mBM3 as a matched-method check. No catalog/source writes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.audit_sata_2010_fe3c import (
    PARAMETERS,
    ROOT,
    VR,
    csv_text,
    linear_solution,
    modified_fit,
    pressure,
    selected_rows,
    statistics,
)
from scripts.audit_sata_2010_mgo_pressure import recalculated_sample_points
from scripts.report_plot_themes import dark_path, dark_svg

REPORT = ROOT / "docs/data/sata-2010-fe3c-recalculated-pressure-refit.json"
INPUTS = ROOT / "docs/data/sata-2010-fe3c-recalculated-refit-inputs.csv"
RESIDUALS = ROOT / "docs/data/sata-2010-fe3c-recalculated-refit-residuals.csv"
FIGURE = (
    ROOT
    / "docs/literature-reproductions/assets/sata-2010-fe3c-recalculated-pressure-refit.svg"
)


def input_rows():
    """Retain both pressure columns and explicit calibration identities."""
    recalculated = {r["source_row"]: r for r in recalculated_sample_points()}
    rows = []
    for source in selected_rows():
        row = dict(source)
        row["published_pressure_gpa"] = row.pop("pressure_gpa")
        # Keep the width's original meaning; it is not used as a fit weight.
        row["published_pressure_error_width_gpa"] = row.pop("error_pressure_gpa")
        row["fit_pressure_gpa"] = (
            recalculated[row["source_row"]]["recalculated_pressure_gpa"]
            if row["source"] == 0
            else row["published_pressure_gpa"]
        )
        row["pressure_shift_gpa"] = (
            row["fit_pressure_gpa"] - row["published_pressure_gpa"]
        )
        row["pressure_basis"] = (
            "recalculated_Sata_MgO_printed_Speziale_BM3"
            if row["source"] == 0
            else "unchanged_Li_original_ruby_neon"
        )
        rows.append(row)
    return rows


def fit_fixed_v0(volume, observed):
    """Same two-coefficient, equal-pressure objective as the original diagnostic."""
    f = ((VR / np.asarray(volume)) ** (2 / 3) - 1) / 2
    design = (1 + 2 * f)[:, None] ** 2.5 * np.column_stack(
        [3 * f - 18 * f**2, 4.5 * f**2]
    )
    beta, covariance, dof = linear_solution(design, observed)
    transform = np.array([[1, 0], [-beta[1] / beta[0] ** 2, 1 / beta[0]]])
    covariance = transform @ covariance @ transform.T
    parameters = np.array([beta[0], beta[1] / beta[0]])
    return {
        "model": "BM3",
        "fixed_values": {"V0_cell_a3": VR},
        "parameter_order": ["K0_gpa", "K0_prime"],
        "parameters": parameters.tolist(),
        "conditional_standard_errors": np.sqrt(np.diag(covariance)).tolist(),
        "conditional_covariance": covariance.tolist(),
        "degrees_of_freedom": dof,
        **statistics(pressure(volume, parameters) - observed),
    }


def reproduce():
    rows = input_rows()
    volume, original, recalculated = [
        np.array([r[key] for r in rows])
        for key in ("volume_a3", "published_pressure_gpa", "fit_pressure_gpa")
    ]
    original_fit = fit_fixed_v0(volume, original)
    recalculated_fit = fit_fixed_v0(volume, recalculated)
    models = {
        "published_BM3": PARAMETERS,
        "original_pressure_diagnostic_BM3": original_fit["parameters"],
        "recalculated_pressure_diagnostic_BM3": recalculated_fit["parameters"],
    }
    cross = {
        name: {
            "on_original_pressures": statistics(pressure(volume, values) - original),
            "on_recalculated_Sata_unchanged_Li": statistics(
                pressure(volume, values) - recalculated
            ),
        }
        for name, values in models.items()
    }
    original_mbm3 = modified_fit(selected_rows())
    recalculated_mbm3_rows = selected_rows()
    for row, value in zip(recalculated_mbm3_rows, recalculated):
        row["pressure_gpa"] = float(value)
    recalculated_mbm3 = modified_fit(recalculated_mbm3_rows)
    residual_rows = []
    for i, row in enumerate(rows):
        result = {
            k: row[k]
            for k in (
                "source",
                "source_row",
                "run",
                "measurement_id",
                "volume_a3",
                "published_pressure_gpa",
                "fit_pressure_gpa",
            )
        }
        for name, evaluation in cross.items():
            result[name + "_residual_original_gpa"] = evaluation[
                "on_original_pressures"
            ]["residual_gpa"][i]
            result[name + "_residual_recalculated_gpa"] = evaluation[
                "on_recalculated_Sata_unchanged_Li"
            ]["residual_gpa"][i]
        residual_rows.append(result)
    source_files = [
        ROOT / "peritheos/data/datasets/cementite-sata-2010-table1-tokens.csv",
        ROOT / "docs/data/li-2002-fe3c/table1-source-tokens.csv",
        ROOT / "docs/data/sata-2010-fe3c-audit.json",
        ROOT / "docs/data/sata-2010-mgo-pressure-check.json",
    ]
    return {
        "kind": "diagnostic_refit",
        "catalog_integration": False,
        "temperature_k": 300,
        "original_pressure_range_gpa": [float(min(original)), float(max(original))],
        "diagnostic_input_pressure_range_gpa": [
            float(min(recalculated)),
            float(max(recalculated)),
        ],
        "selection": "All nine Sata Fe3C/MgO pairs plus all three cited Li measured rows above 25 GPa. No exclusions or additional observations.",
        "pressure_reduction": {
            "Sata": "Measured MgO cell V, BM3 V0=74.71 A3, K0=160.2 GPa, Kprime=3.99 at 300 K",
            "Li": "Original pressure, volume and source identity retained unchanged; no MgO or Litasov recalibration.",
        },
        "objective": "Minimize sum [P_model(V_i)-P_input_i]^2, equal pressure weights, at original measured Fe3C volumes.",
        "fixed_V0_qualification": "V0=149.456 A3/cell is the same fixed printed anchor as the previous two-coefficient diagnostic, not a newly measured ambient volume or an inferred author constraint.",
        "uncertainty_convention": "Conditional covariance: RSS/(N-p) scaled inverse normal matrix, transformed to the reported coefficients. Coordinate errors, calibration coefficients/covariance and fixed-V0 uncertainty are not propagated. Errors are descriptive residual-scaled estimates, not recovered source errors or predictive confidence bands. Old/new fits share observations; their coefficient differences are paired deterministic shifts, not independent-fit significance tests.",
        "scope_restriction": "Separate requested pressure-recalculation diagnostic only. No original author-workflow or pressure-scale validation, established correction, replacement source EOS, quantitative-use recommendation or new catalog record. Existing published coefficients, source observations, original-data conditional parameter parity and poor-data-agreement/recommendation findings remain unchanged.",
        "inputs": rows,
        "fixed_V0_BM3": {
            "original_pressures": original_fit,
            "recalculated_Sata_unchanged_Li": recalculated_fit,
            "coefficient_shift_recalculated_minus_original": (
                np.array(recalculated_fit["parameters"]) - original_fit["parameters"]
            ).tolist(),
        },
        "reference_mBM3_check": {
            "fixed_Vr_cell_a3": VR,
            "parameter_order": ["Pr_gpa", "Kr_gpa", "Kr_prime"],
            "original_pressures": original_mbm3,
            "recalculated_Sata_unchanged_Li": recalculated_mbm3,
            "coefficient_shift_recalculated_minus_original": (
                np.array(recalculated_mbm3["parameters_Pr_Kr_Kprime"])
                - original_mbm3["parameters_Pr_Kr_Kprime"]
            ).tolist(),
            "qualification": "Matched three-coefficient equal-pressure check at selected Vr, with fitted Pr; do not reinterpret Vr/Kr as zero-pressure V0/K0 when Pr differs from zero.",
        },
        "cross_evaluation": cross,
        "residuals": residual_rows,
        "source_sha256": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source_files
        },
    }


def plot(report, preview=None):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    rows = report["inputs"]
    v, old, new = [
        np.array([r[k] for r in rows])
        for k in ("volume_a3", "published_pressure_gpa", "fit_pressure_gpa")
    ]
    own = np.array([r["source"] == 0 for r in rows])
    vg = np.linspace(min(v), max(v), 250)
    fits = report["fixed_V0_BM3"]
    values = [
        PARAMETERS,
        fits["original_pressures"]["parameters"],
        fits["recalculated_Sata_unchanged_Li"]["parameters"],
    ]
    labels = [
        "Published BM3",
        "Diagnostic fit: original pressures",
        "Diagnostic fit: recalculated Sata + unchanged Li",
    ]
    colors = ["#4e5660", "#245b8a", "#ae4c14"]
    style = {
        "svg.hashsalt": "sata-recalculated-refit-v1",
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
        "font.size": 9,
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(
            2,
            2,
            figsize=(11.8, 7),
            height_ratios=[1.5, 1],
            sharex="col",
            sharey="row",
            layout="constrained",
        )
        fig.suptitle(
            "Fe₃C pressure-recalculation diagnostic: matched BM3 refits", fontsize=13
        )
        for column, (coordinate, key, title) in enumerate(
            [
                (old, "on_original_pressures", "Original Sata + original Li pressures"),
                (
                    new,
                    "on_recalculated_Sata_unchanged_Li",
                    "Recalculated Sata + original Li pressures",
                ),
            ]
        ):
            top, bottom = axes[:, column]
            for pars, label, color, linestyle in zip(
                values, labels, colors, [":", "--", "-"]
            ):
                top.plot(pressure(vg, pars), vg, color=color, ls=linestyle, label=label)
            top.scatter(
                coordinate[own], v[own], s=23, color="#4e5660", label="Nine Sata states"
            )
            top.scatter(
                coordinate[~own],
                v[~own],
                marker="D",
                s=23,
                facecolors="none",
                edgecolors="#4e5660",
                label="Three unchanged Li states",
            )
            for name, label, color, marker in zip(
                report["cross_evaluation"], labels, colors, ["x", "o", "^"]
            ):
                bottom.scatter(
                    coordinate,
                    report["cross_evaluation"][name][key]["residual_gpa"],
                    color=color,
                    marker=marker,
                    s=23,
                    label=label,
                )
            bottom.axhline(0, color="#4e5660", lw=0.8)
            top.set_title(title)
            bottom.set_xlabel("Pressure coordinate used (GPa)")
        axes[0, 0].set_ylabel("Volume (Å³ / conventional cell)")
        axes[1, 0].set_ylabel("P(model) − P(input) (GPa)")
        axes[0, 0].legend(loc="upper right", fontsize=7, frameon=False)
        for ax in axes.flat:
            ax.grid(color="#e7eaee")
            ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False)
        fig.text(
            0.5,
            0.02,
            "All fits: 300 K, equal pressure weights, fixed V₀=149.456 Å³/cell. No new uncertainty bands inferred.\nSeparate diagnostic only; published record and original-data parity are preserved; catalog qualifications are separate.",
            ha="center",
            fontsize=8,
        )
        fig.get_layout_engine().set(rect=(0, 0.085, 1, 0.915))
        fig.savefig(
            FIGURE, metadata={"Date": None, "Description": report["scope_restriction"]}
        )
        FIGURE.with_name(Path(dark_path(FIGURE.name)).name).write_bytes(
            dark_svg(FIGURE.read_bytes())
        )
        if preview:
            preview.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(preview, dpi=160)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--preview", type=Path)
    args = parser.parse_args()
    report = reproduce()
    outputs = {
        REPORT: json.dumps(report, indent=2) + "\n",
        INPUTS: csv_text(report["inputs"]),
        RESIDUALS: csv_text(report["residuals"]),
    }
    if args.write:
        for path, content in outputs.items():
            path.write_text(content)
    elif args.check:
        # Deterministic linear calculations in the same maintained environment.
        for path, content in outputs.items():
            if path.read_text(encoding="utf-8") != content:
                raise SystemExit(
                    f"Recalculated-pressure diagnostic needs regeneration: {path}"
                )
        print("Separate Sata recalculated-pressure refit checks passed")
    elif not args.plot:
        print(outputs[REPORT])
    if args.plot:
        plot(report, args.preview)
        print(FIGURE.relative_to(ROOT))


if __name__ == "__main__":
    main()
