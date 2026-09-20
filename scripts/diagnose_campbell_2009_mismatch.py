#!/usr/bin/env python3
"""Locate conditional Campbell/Seagle fit tension without selecting new EOS values."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

if __package__:
    from .audit_campbell_2009_protocol import dataset
    from .reconstruct_seagle_2008_fcc import (
        CELL_TO_MOLAR,
        V0_1400,
        read_csv,
    )
    from .reproduce_campbell_2009_buffers import SOURCE, pressure, read_rows
else:
    from audit_campbell_2009_protocol import dataset
    from reconstruct_seagle_2008_fcc import (
        CELL_TO_MOLAR,
        V0_1400,
        read_csv,
    )
    from reproduce_campbell_2009_buffers import SOURCE, pressure, read_rows

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/campbell-2009-mismatch-diagnostics.json"
FIGURE = ROOT / "docs/images/campbell-2009-mismatch.png"
PARAMETERS = ["K0_gpa", "gamma0", "q"]
VARIANTS = {
    "baseline": (300, 6.5),
    "delta_6": (300, 6.0),
    "delta_7": (300, 7.0),
    "alpha_at_1400": (1400, 6.5),
    "alpha_at_observed_temperature": ("observed", 6.5),
}


def thermal_variant(cell, temperature, anchor=300, delta=6.5):
    """Mean alpha scaled by V/V0 at a specified common reference temperature.

    These are alternative implementations of an incompletely specified
    historical reduction, not alternative source-published EOS. Keep the
    Boehler mean-expansion law; do not exponentiate an instantaneous alpha.
    """
    cells, temperatures = np.broadcast_arrays(cell, temperature)

    def scalar(c, t):
        v = c * CELL_TO_MOLAR
        tr = t if anchor == "observed" else anchor
        ambient = 6.835 * (1 + 7.70e-5 * (tr - 300))

        def at_reference(a):
            return v * (1 + a * (tr - 300)) / (1 + a * (t - 300))

        alpha = brentq(
            lambda a: a - 7.70e-5 * (at_reference(a) / ambient) ** delta,
            0,
            0.001,
            xtol=1e-15,
        )
        ratio = v * (1 + alpha * 1100) / (1 + alpha * (t - 300)) / V0_1400
        return (
            180
            * (ratio ** (-7 / 3) - ratio ** (-5 / 3))
            * (1 + 0.75 * (ratio ** (-2 / 3) - 1))
        )

    return np.array(
        [scalar(c, t) for c, t in zip(cells.flat, temperatures.flat)]
    ).reshape(cells.shape)


def fit(material, data, fixed_q=None):
    """Same equal-pressure objective and bounds as the existing audit."""
    source = np.array(SOURCE[material])
    indices = [1, 4, 5] if fixed_q is None else [1, 4]
    if fixed_q is not None:
        source[5] = fixed_q

    def candidate(x):
        p = source.copy()
        p[indices] = x
        return p

    def residual(x):
        return pressure(data[0], data[1], candidate(x)) - data[2]

    lower, upper = ([20, 0.01, 0.01], [400, 10, 10])
    solution = least_squares(
        residual,
        source[indices],
        bounds=(lower[: len(indices)], upper[: len(indices)]),
        x_scale="jac",
        ftol=1e-10,
        xtol=1e-10,
        gtol=1e-10,
    )
    p = candidate(solution.x)
    return {
        "observations": data.shape[1],
        "parameters": dict(zip(PARAMETERS, p[[1, 4, 5]].tolist())),
        "sse_gpa2": float(np.sum(solution.fun**2)),
        "rmse_gpa": float(np.sqrt(np.mean(solution.fun**2))),
        "active_bounds": [
            PARAMETERS[i] for i, active in enumerate(solution.active_mask) if active
        ],
        "solver_success": bool(solution.success),
    }


def residual_summary(residual):
    return {
        "observations": len(residual),
        "mean_gpa": float(np.mean(residual)),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "sse_gpa2": float(np.sum(residual**2)),
    }


def metadata(ids):
    current = {
        r["source_row"]: r for r in read_rows("feo-campbell-2009-table-s2-pvt.csv")
    }
    previous = {
        r["source_row"]: r
        for r in read_rows("feo-seagle-2008-supplement-volume-temperature.csv")
    }
    series = {
        r["source_row"]: r["nominal_series_pressure_gpa"]
        for r in read_csv(ROOT / "docs/data/seagle-2008-reconstructed-pressures.csv")
    }
    result = []
    for key in ids:
        study, row = key.split(":")
        r = current[row] if study == "campbell_s2" else previous[row]
        result.append(
            {
                "row_id": key,
                "study": study,
                "iron_phase": r["iron_phase"],
                "series": r["run_id"] if study == "campbell_s2" else series[row],
            }
        )
    return result


def material_diagnostics(material, scale="seagle_2006"):
    data, ids, _ = dataset(material, scale, reconstruct_fcc=True)
    meta = metadata(ids)
    published = pressure(data[0], data[1], SOURCE[material]) - data[2]
    full = fit(material, data)
    p = np.array(SOURCE[material])
    p[[1, 4, 5]] = list(full["parameters"].values())
    refitted = pressure(data[0], data[1], p) - data[2]
    source_masks = {
        "Campbell": np.array([r["study"] == "campbell_s2" for r in meta]),
        "Seagle_fcc": np.array(
            [r["study"] == "seagle" and r["iron_phase"] == "fcc" for r in meta]
        ),
        "Seagle_hcp": np.array(
            [r["study"] == "seagle" and r["iron_phase"] == "hcp" for r in meta]
        ),
    }
    partitions = {
        "study_and_marker": source_masks,
        "temperature": {
            label: (data[1] >= lo) & (data[1] < hi)
            for label, lo, hi in [
                ("below_1500", 0, 1500),
                ("1500_to_2000", 1500, 2000),
                ("2000_to_2500", 2000, 2500),
                ("2500_and_above", 2500, 10000),
            ]
        },
        "seagle_nominal_series": {
            series: np.array(
                [r["study"] == "seagle" and r["series"] == series for r in meta]
            )
            for series in sorted(
                {r["series"] for r in meta if r["study"] == "seagle"}, key=float
            )
        },
    }
    groups = {}
    for partition, masks in partitions.items():
        groups[partition] = {}
        for name, mask in masks.items():
            if not mask.any():
                continue
            groups[partition][name] = {
                "row_ids": [key for key, keep in zip(ids, mask) if keep],
                "published": residual_summary(published[mask]),
                "refitted": residual_summary(refitted[mask]),
                "published_sse_fraction": float(
                    np.sum(published[mask] ** 2) / np.sum(published**2)
                ),
                "without_group": fit(material, data[:, ~mask]),
            }
    profiles = {}
    for name, mask in {"combined": np.ones(len(ids), bool), **source_masks}.items():
        if np.count_nonzero(mask) < 4:
            continue
        subset = data[:, mask]
        profiles[name] = {
            "free_q_fit": fit(material, subset),
            "fixed_q_fits": [
                fit(material, subset, q)
                for q in [0.01, 0.3, 0.6, 1.0, 1.3, 1.6, 2.0, 3.0]
            ],
        }
    return {
        "scale": scale,
        "full_fit": full,
        "groups": groups,
        "q_profiles": profiles,
        "rows": [
            {
                **m,
                "temperature_k": float(t),
                "pressure_gpa": float(p),
                "published_residual_gpa": float(r),
                "refitted_residual_gpa": float(f),
            }
            for m, t, p, r, f in zip(meta, data[1], data[2], published, refitted)
        ],
    }


def thermal_diagnostics():
    checkpoints = read_csv(ROOT / "docs/data/seagle-2008-fcc-pressure-checkpoints.csv")
    seagle = {
        f"seagle:{r['source_row']}": r
        for r in read_rows("feo-seagle-2008-supplement-volume-temperature.csv")
    }
    results = {}
    for name, (anchor, delta) in VARIANTS.items():
        differences, inside = [], 0
        for r in checkpoints:
            difference = float(
                thermal_variant(
                    float(r["cell_volume_a3"]), float(r["temperature_k"]), anchor, delta
                )
            ) - float(r["published_pressure_gpa"])
            differences.append(difference)
            inside += abs(difference) <= float(r["published_pressure_uncertainty_gpa"])
        fits = {}
        for material, scale in [
            ("fe_fcc", "seagle_2006"),
            ("feo", "seagle_2006"),
            ("feo", "dewaele_2006_hypothesis"),
        ]:
            data, ids, _ = dataset(material, scale, reconstruct_fcc=True)
            shifts = []
            for i, key in enumerate(ids):
                if key in seagle and seagle[key]["iron_phase"] == "fcc":
                    r = seagle[key]
                    new_p = float(
                        thermal_variant(
                            float(r["iron_unit_cell_volume_a3"]),
                            data[1, i],
                            anchor,
                            delta,
                        )
                    )
                    shifts.append(new_p - data[2, i])
                    data[2, i] = new_p
            fits[f"{material}_{scale}"] = {
                **fit(material, data),
                "fcc_pressure_shift_mean_gpa": float(np.mean(shifts)),
                "fcc_pressure_shift_range_gpa": [
                    float(min(shifts)),
                    float(max(shifts)),
                ],
            }
        results[name] = {
            "alpha_volume_reference_temperature_k": anchor,
            "delta": delta,
            "checkpoint_differences_gpa": differences,
            "checkpoint_summary": residual_summary(np.array(differences)),
            "checkpoints_within_reported_error": int(inside),
            "fits": fits,
        }
    return results


def reproduce():
    return {
        "format": "peritheos.campbell-2009-mismatch-diagnostics",
        "format_version": 1,
        "qualification": "Diagnostic equal-pressure fits only, with the prior audit bounds. Residual sign is P(Campbell published model) minus input pressure. Group removal changes data coverage and is not evidence that a source group is invalid. Nominal Seagle series assignments are inherited audit inferences from FeO expansion lines, not recovered author run identifiers. Temperature and pressure/source coverage are confounded. Thermal variants are explicit assumptions, not recovered author implementations; no variant is adopted by closeness to Campbell coefficients. No published coefficient or original pressure reconstruction is modified. The Fe/FeO source records are deferred from execution. Fixed-q profiles are sensitivity curves, not confidence intervals or proof of a global optimum.",
        "materials": {m: material_diagnostics(m) for m in ["fe_fcc", "feo"]},
        "feo_dewaele_control": material_diagnostics("feo", "dewaele_2006_hypothesis"),
        "thermal_variants": thermal_diagnostics(),
    }


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.4), layout="constrained")
    colors = {"Campbell": "#2563a6", "Seagle_fcc": "#b66915", "Seagle_hcp": "#9b3278"}
    for row, material in enumerate(["fe_fcc", "feo"]):
        entry = report["materials"][material]
        ax = axes[row, 0]
        for label, group in entry["groups"]["study_and_marker"].items():
            selected = [r for r in entry["rows"] if r["row_id"] in group["row_ids"]]
            ax.scatter(
                [r["temperature_k"] for r in selected],
                [r["published_residual_gpa"] for r in selected],
                s=28,
                alpha=0.85,
                label=label.replace("_", " "),
                color=colors[label],
            )
        ax.axhline(0, color="#555555", lw=0.8)
        ax.set(
            xlabel="Temperature (K)",
            ylabel="Published-model pressure − input (GPa)",
            title=("fcc Fe" if material == "fe_fcc" else "FeO")
            + ": residuals by source",
        )
        ax.legend(fontsize=8)
        ax = axes[row, 1]
        for name, profile in entry["q_profiles"].items():
            fits = profile["fixed_q_fits"]
            ax.plot(
                [r["parameters"]["q"] for r in fits],
                [
                    100 * (r["sse_gpa2"] / profile["free_q_fit"]["sse_gpa2"] - 1)
                    for r in fits
                ],
                marker=".",
                label=name.replace("_", " "),
                color=colors.get(name, "#222222"),
            )
        ax.axvline(
            SOURCE[material][5], color="#555555", ls="--", lw=0.8, label="Published q"
        )
        ax.set(
            xlabel="Fixed q (K₀ and γ₀ refitted)",
            ylabel="SSE increase over free-q fit (%)",
            title="q sensitivity by source",
        )
        ax.legend(fontsize=8)
    fig.suptitle(
        "Campbell–Seagle mismatch diagnostics\nConditional Seagle pressure reconstruction; equal pressure-residual weights",
        fontsize=13,
    )
    FIGURE.parent.mkdir(exist_ok=True)
    fig.savefig(FIGURE, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != text:
            raise SystemExit("Mismatch diagnostic report is stale")
    else:
        OUTPUT.write_text(text)
        plot(report)
    print(
        "Checked source/temperature/series influence, q profiles, and five thermal variants"
    )


if __name__ == "__main__":
    main()
