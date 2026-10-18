"""User-requested independent BM3 fits; no bundled catalog changes.

Run from the repository with its Python environment. All fitting and plot
outputs remain beside this script. Source observations are never modified.
"""

from __future__ import annotations

import csv
import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3
from peritheos.fitting import fit_rt_eos

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
DATA = ROOT / "peritheos/data/datasets/palladium-baty-2024-table-s1-compression.csv"
SHA256 = "513835bba914f6fa1f02646901e0453edda808c6da7bc425fbb71c8cf1c99829"
PUBLISHED = {"V0": 58.88, "K0": 190.0, "K0_prime": 5.3}
BOUNDS = {"V0": (50.0, 70.0), "K0": (10.0, 500.0), "K0_prime": (1.0, 15.0)}


def literal_bm3(v, q):
    eta = q["V0"] / np.asarray(v)
    return (
        1.5
        * q["K0"]
        * (eta ** (7 / 3) - eta ** (5 / 3))
        * (1 + 0.75 * (q["K0_prime"] - 4) * (eta ** (2 / 3) - 1))
    )


def fit(v, p, anchor):
    names = ("K0", "K0_prime")
    return fit_rt_eos(
        BM3,
        volume=v,
        pressure=p,
        initial={name: PUBLISHED[name] for name in names},
        fixed={"V0": anchor},
        bounds={name: BOUNDS[name] for name in names},
        absolute_sigma=False,
        max_nfev=5000,
    )


def metrics(residual):
    return {
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "mean_pressure_residual_gpa": float(np.mean(residual)),
        "max_absolute_pressure_residual_gpa": float(np.max(np.abs(residual))),
    }


def main():
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == SHA256
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 78
    p = np.array([float(row["pressure_gpa"]) for row in rows])
    v = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    finite = p > 0
    assert np.count_nonzero(~finite) == 1
    v0 = float(v[~finite][0])
    fits = {"ambient_anchored": fit(v[finite], p[finite], v0)}
    report = {
        "record_kind": "refit",
        "date": "2026-10-17",
        "source": {
            "doi": "10.1063/5.0179469",
            "table": "Official supplementary Table S1",
            "data_sha256": SHA256,
        },
        "software": {
            name: version(name)
            for name in ("peritheos", "numpy", "scipy", "matplotlib")
        },
        "dataset": "palladium_baty_2024_table_s1_compression",
        "source_rows": 78,
        "model": "BM3",
        "temperature_k": 300.0,
        "volume_basis": "angstrom^3/conventional_fcc_cell_with_4_atoms",
        "observed_pressure_range_gpa": [float(p.min()), float(p.max())],
        "observed_volume_range_a3_cell": [float(v.min()), float(v.max())],
        "source_protocol_reproduced": False,
        "catalog_integrated": False,
        "usage": "Fixed-ambient-V0 fit-quality diagnostic only. Systematic residuals remain unresolved; no single EOS adopted for quantitative scientific use or pressure calibration.",
        "historical_diagnostic_archive": "fit-results-with-free-v0-historical.json",
        "objective": "sum((P_BM3(V_observed) - P_observed)^2), equal weight per fitted row",
        "weights": "No invented row-wise sigmas. Stated global coordinate limits are not treated as statistical weights.",
        "bounds": BOUNDS,
        "published": {
            "parameters": PUBLISHED,
            **metrics(literal_bm3(v, PUBLISHED) - p),
        },
        "refits": {},
    }
    residuals = {"published": literal_bm3(v, PUBLISHED) - p}
    for label, result in fits.items():
        assert result.success, result.message
        anchored = label == "ambient_anchored"
        selected = finite if anchored else np.ones(p.size, dtype=bool)
        q = result.parameters
        residuals[label] = literal_bm3(v, q) - p
        # Independent Python expression and SciPy solver, distinct from native
        # Peritheos fitting and evaluation. The objective/row selection is identical.
        names = result.free_parameters
        independent = least_squares(
            lambda a: (
                literal_bm3(v[selected], {**q, **dict(zip(names, a))}) - p[selected]
            ),
            [PUBLISHED[name] for name in names],
            bounds=(
                [BOUNDS[name][0] for name in names],
                [BOUNDS[name][1] for name in names],
            ),
            x_scale="jac",
            ftol=1e-13,
            xtol=1e-13,
            gtol=1e-13,
        )
        independent_q = {**q, **dict(zip(names, independent.x))}
        assert independent.success
        assert (
            np.max(np.abs(literal_bm3(v, independent_q) - result.model.pressure(v)))
            < 1e-5
        )
        pressure_grid = np.linspace(0, p.max(), 100)
        volumes = result.model.volume(pressure_grid)
        assert np.all(result.model.bulk_modulus(volumes) > 0)
        assert np.max(np.abs(result.model.pressure(volumes) - pressure_grid)) < 1e-7
        # Leave-one-out prediction checks describe same-dataset stability only;
        # unidentified run correlations and calibration systematics remain unknown.
        cv = []
        for idx in np.flatnonzero(selected):
            train = selected.copy()
            train[idx] = False
            trial = fit(v[train], p[train], v0 if anchored else None)
            assert trial.success
            cv.append(float(trial.model.pressure(v[idx]) - p[idx]))
        report["refits"][label] = {
            "identifier": f"palladium_baty_2024_{label}_bm3_refit",
            "parameters": q,
            "free_parameters": names,
            "fixed_parameters": {"V0": v0} if anchored else {},
            "selection": "77 finite-pressure rows fitted; the ambient row fixes V0"
            if anchored
            else "all 78 rows; V0, K0, K0_prime varied",
            "all_78_row_statistics": metrics(residuals[label]),
            "fitted_rows_statistics": metrics(residuals[label][selected]),
            "conditional_standard_errors": {
                name: result.standard_errors[name] for name in names
            },
            "conditional_covariance_parameter_order": names,
            "conditional_covariance": result.covariance.tolist(),
            "conditional_correlation": result.correlation.tolist(),
            "covariance_assumptions": "Local Jacobian; independent equal-variance pressure residuals; SSE/(n_fit-n_free) scaling. For the anchored fit, V0 is conditional/exact in this calculation, not an error-free measurement. These are formal regression errors, not total scientific uncertainties.",
            "degrees_of_freedom": result.degrees_of_freedom,
            "pressure_at_measured_ambient_volume_gpa": float(result.model.pressure(v0)),
            "V0_minus_measured_ambient_volume_a3_cell": q["V0"] - v0,
            "leave_one_out_same_dataset": {
                **metrics(np.array(cv)),
                "rows": len(cv),
                "qualification": "Not independent experimental validation; run identities/calibration systematics unavailable.",
            },
            "independent_scipy_parameters": independent_q,
            "native_scipy_curve_max_difference_gpa": float(
                np.max(np.abs(literal_bm3(v, independent_q) - result.model.pressure(v)))
            ),
        }
    report["ambient_volume_bound_sensitivity"] = []
    for anchor in (v0 - 0.04, v0, v0 + 0.04):
        result = fit(v[finite], p[finite], anchor)
        assert result.success
        report["ambient_volume_bound_sensitivity"].append(
            {
                "V0_fixed_a3_cell": anchor,
                "parameters": result.parameters,
                "fitted_rows_statistics": metrics(result.residuals),
            }
        )
    report["ambient_volume_bound_sensitivity_meaning"] = (
        "Anchor varied by +/- 4*0.01 A^3/cell using the source atomic-volume upper uncertainty limit. Sensitivity bounds, not Gaussian standard deviations or a confidence interval."
    )
    (OUT / "fit-results.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    with (OUT / "residuals.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "source_csv_data_row_1_based",
                "pressure_observed_gpa",
                "atomic_volume_observed_a3",
                *[f"{label}_Pcalc_minus_Pobs_gpa" for label in residuals],
            ]
        )
        for idx in range(p.size):
            writer.writerow(
                [idx + 1, p[idx], v[idx] / 4, *[r[idx] for r in residuals.values()]]
            )

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=(12, 8.5), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1], height_ratios=[1.35, 1])
    full, zoom, resid = (
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, :]),
    )
    config = {
        "published": ("Published experimental BM3", "#3366aa", "-", PUBLISHED),
        "ambient_anchored": (
            "BM3 check: measured V0 fixed",
            "#278352",
            "-.",
            fits["ambient_anchored"].parameters,
        ),
    }
    curve_p = np.linspace(0, p.max(), 500)
    volume_residuals = {
        label: v / 4 - np.asarray(BM3(**q).volume(p)) / 4
        for label, (_, _, _, q) in config.items()
    }
    for label, delta_v in volume_residuals.items():
        if label == "published":
            report["published"]["volume_rmse_a3_atom"] = float(
                np.sqrt(np.mean(delta_v**2))
            )
        else:
            report["refits"][label]["all_78_row_statistics"]["volume_rmse_a3_atom"] = (
                float(np.sqrt(np.mean(delta_v**2)))
            )
    (OUT / "fit-results.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    with (OUT / "volume-residuals.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "pressure_observed_gpa",
                "atomic_volume_observed_a3",
                *[f"{label}_Vobs_minus_Vcalc_a3_atom" for label in volume_residuals],
            ]
        )
        writer.writerows(
            [p[i], v[i] / 4, *[r[i] for r in volume_residuals.values()]]
            for i in range(len(p))
        )
    for ax in (full, zoom):
        ax.scatter(
            p, v / 4, s=18, c="#222222", zorder=5, label="Table S1: 78 observations"
        )
        for name, color, style, q in config.values():
            ax.plot(
                curve_p,
                np.asarray(BM3(**q).volume(curve_p)) / 4,
                color=color,
                ls=style,
                lw=1.9,
                label=name,
            )
        ax.set(xlabel="Pressure (GPa)", ylabel="Atomic volume (Å³/atom)")
        ax.grid(alpha=0.17)
    full.set(
        title="Measured compression and BM3 curves", xlim=(-1, 77), ylim=(11.7, 14.85)
    )
    full.legend(fontsize=8.5, loc="upper right", frameon=False)
    zoom.set(title="Low-pressure detail", xlim=(-0.7, 25), ylim=(13.3, 14.85))
    zoom.scatter([0], [v0 / 4], s=70, facecolors="none", edgecolors="#222222", zorder=6)
    zoom.annotate(
        "Measured ambient V0 retained\n14.72391303 Å³/atom",
        xy=(0, v0 / 4),
        xytext=(12, 14.62),
        fontsize=9,
        ha="center",
        arrowprops={"arrowstyle": "->", "color": "#278352"},
        color="#278352",
    )
    for label, (name, color, style, _) in config.items():
        r = volume_residuals[label]
        resid.scatter(
            p,
            r,
            color=color,
            s=19,
            marker={
                "published": "o",
                "free_reference_volume": "^",
                "ambient_anchored": "s",
            }[label],
            alpha=0.8,
            label=f"{name}: volume RMS {np.sqrt(np.mean(r**2)):.4f} Å³/atom",
        )
    resid.axhline(0, color="#333333", lw=1)
    resid.set(
        xlabel="Observed pressure (GPa)",
        ylabel="Observed − calculated volume (Å³/atom)",
        title="Volume residuals at every observed pressure",
        xlim=(-1, 77),
    )
    resid.grid(alpha=0.17)
    resid.legend(fontsize=8.5, loc="upper right", frameon=False)
    fig.suptitle(
        "Palladium at 300 K — Baty et al. (2024), official Table S1",
        fontsize=15,
        weight="bold",
    )
    fig.supxlabel(
        "Measured V0 fixed in the independent BM3 check. Systematic residuals remain unresolved; no single EOS adopted for quantitative use.",
        fontsize=9,
    )
    fig.savefig(OUT / "fit-comparison.png", dpi=170)
    fig.savefig(OUT / "fit-comparison.svg")
    plt.close(fig)
    print(
        json.dumps(
            {
                label: {
                    "parameters": r["parameters"],
                    "all_rows": r["all_78_row_statistics"],
                    "leave_one_out": r["leave_one_out_same_dataset"],
                }
                for label, r in report["refits"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
