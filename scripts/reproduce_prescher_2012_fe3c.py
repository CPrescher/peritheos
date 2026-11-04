"""Replay the recovered Fe3C table without inventing published V0 values.

Repository maintenance only. The equal-weight baseline and coordinate-width
diagnostic are explicit assumptions, not recovery of the source's objective.
Neither creates a replacement catalog EOS.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from peritheos.eos.rt import BM3
from peritheos.fitting import fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
DATASET_ID = "cementite_prescher_2012_table_s1"
DATA = ROOT / "peritheos/data/datasets/cementite-prescher-2012-table-s1.csv"
REPORT = ROOT / "docs/data/prescher-2012-fe3c-reproduction.json"
SOURCE = {
    "full_range": {"K0": 161.0, "K0_prime": 5.9},
    "below_22_gpa": {"K0": 145.0, "K0_prime": 8.5},
    "above_22_gpa": {"K0": 172.0, "K0_prime": 5.1},
}
SOURCE_ERRORS = {
    "full_range": {"K0": 2.0, "K0_prime": 0.2},
    "below_22_gpa": {"K0": 3.0, "K0_prime": 0.7},
    "above_22_gpa": {"K0": 1.0, "K0_prime": 0.1},
}


def observations():
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    return {name: np.array([float(row[name]) for row in rows]) for name in rows[0]}


def mask_for(name, pressure):
    if name == "full_range":
        return np.ones(pressure.shape, dtype=bool)
    if name == "below_22_gpa":
        return pressure < 22
    if name == "above_22_gpa":
        return pressure >= 22
    raise ValueError(f"Unknown source selection {name}")


def reproduce():
    data = observations()
    results = {}
    for name, source in SOURCE.items():
        mask = mask_for(name, data["pressure_gpa"])
        p = data["pressure_gpa"][mask]
        v = data["volume_a3"][mask]
        fits = {}
        for method in ("equal_pressure_weights", "reported_coordinate_widths"):
            coordinate = method == "reported_coordinate_widths"
            fit = fit_rt_eos(
                BM3,
                v,
                p,
                initial={"V0": 155, **source},
                bounds={"V0": (145, 170), "K0": (50, 400), "K0_prime": (0, 15)},
                pressure_sigma=data["pressure_error_gpa"][mask] if coordinate else None,
                volume_sigma=data["volume_error_a3"][mask] if coordinate else None,
                absolute_sigma=False,
                max_nfev=5000,
            )
            residual = np.asarray(fit.model.pressure(v)) - p
            comparisons = []
            for key, value in source.items():
                reconstructed = float(fit.parameters[key])
                width = SOURCE_ERRORS[name][key]
                comparisons.append(
                    {
                        "parameter": key,
                        "published": value,
                        "published_error": width,
                        "reconstructed": reconstructed,
                        "reconstruction_error": float(fit.standard_errors[key]),
                        "difference": reconstructed - value,
                        "within_printed_source_error_width": abs(reconstructed - value)
                        <= width,
                    }
                )
            fits[method] = {
                "parameters": {k: float(x) for k, x in fit.parameters.items()},
                "parameter_errors": {
                    k: float(x) for k, x in fit.standard_errors.items()
                },
                "source_comparisons": comparisons,
                "K0_and_K0_prime_within_source_widths": all(
                    item["within_printed_source_error_width"] for item in comparisons
                ),
                "observations": int(mask.sum()),
                "free_parameters": list(fit.free_parameters),
                "fixed_parameters": [],
                "source_rows": data["source_row"][mask].astype(int).tolist(),
                "pressure_range_gpa": [float(p.min()), float(p.max())],
                "covariance_order": list(fit.free_parameters),
                "covariance": np.asarray(fit.covariance).tolist(),
                "covariance_convention": "Residual-scaled inverse normal matrix, conditional on the stated objective and row selection; not a recovered source confidence convention.",
                "objective": "errors_in_variables"
                if coordinate
                else "pressure_residuals",
                "weighting": (
                    "Reported last-digit P/V widths treated as numerical coordinate weights; source does not define them as independent standard deviations."
                    if coordinate
                    else "Equal pressure weights; printed coordinate errors not used."
                ),
                "pressure_residuals_gpa_at_observed_volumes": residual.tolist(),
                "pressure_rmse_gpa_at_observed_volumes": float(
                    np.sqrt(np.mean(residual**2))
                ),
                "solver_success": bool(fit.success),
                "solver_message": str(fit.message),
            }
        results[name] = {
            "published_parameters": {"V0": None, **source},
            "published_parameter_errors": {"V0": None, **SOURCE_ERRORS[name]},
            "published_V0_status": "Not stated in inspected main paper or supplied supplement; fixed/shared/refined convention unknown.",
            "published_curve_replay": "unavailable_without_source_V0",
            "author_workflow_reproduced": False,
            "reconstruction_kind": "conditional_diagnostics_not_catalog_refits",
            "fits": fits,
        }
    return {
        "format": "peritheos.prescher-2012-fe3c-reproduction",
        "reference_doi": "10.1103/PhysRevB.85.140402",
        "dataset_identifier": DATASET_ID,
        "observations": len(data["pressure_gpa"]),
        "temperature_k": 293,
        "qualification": "All 19 Table S1 observations are preserved. The ambient point plotted in Figure 1 is not tabulated here; no numerical ambient anchor or published V0 has been synthesized. K0/K0-prime agreement alone does not recover the complete three-parameter EOS or the author's workflow.",
        "source_subsets": results,
    }


def plot(report):
    """Plot only reconstructed curves: the published V0 is unavailable."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = observations()
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    axes[0].errorbar(
        data["pressure_gpa"],
        data["volume_a3"],
        xerr=data["pressure_error_gpa"],
        yerr=data["volume_error_a3"],
        fmt="o",
        color="black",
        ms=3,
        capsize=2,
        label="Table S1 (19 rows)",
    )
    colors = ["#0072B2", "#D55E00", "#009E73"]
    for (name, source), color in zip(report["source_subsets"].items(), colors):
        fit = source["fits"]["reported_coordinate_widths"]
        model = BM3(**fit["parameters"])
        mask = mask_for(name, data["pressure_gpa"])
        pressure = np.linspace(*fit["pressure_range_gpa"], 200)
        label = name.replace("_", " ")
        axes[0].plot(pressure, model.volume(pressure), color=color, label=label)
        axes[1].plot(
            data["pressure_gpa"][mask],
            fit["pressure_residuals_gpa_at_observed_volumes"],
            "o",
            color=color,
            ms=4,
            label=label,
        )
    axes[0].set(xlabel="Pressure (GPa)", ylabel="Volume (Å³ / conventional cell)")
    axes[1].set(xlabel="Pressure (GPa)", ylabel="P(reconstruction) − P(table) (GPa)")
    axes[1].axhline(0, color="0.6", lw=1)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    figure.suptitle("Conditional BM3 reconstructions — published V₀ remains unknown")
    out = (
        ROOT
        / "docs/literature-reproductions/figures/prescher-2012-fe3c-conditional.png"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out, dpi=180)
    plt.close(figure)
    return out


if __name__ == "__main__":
    result = reproduce()
    REPORT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Preserved 19 observations; wrote {REPORT.relative_to(ROOT)}")
