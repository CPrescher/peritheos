"""Cross-check published Pt models against Fei data without refitting.

Independent SI quadrature is the oracle. Original 2004 calibration and the
2007 Au recalibration are separate targets; figure curves are not observations.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from peritheos import get_eos_record
from scripts import audit_matsui_2009_platinum as matsui
from scripts import reproduce_fei_2004_pressure_scales as fei04
from scripts import reproduce_fei_2007_platinum as fei07
from scripts import reproduce_zhu_2025_pressure_standards as zhu

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/data/pt-thermal-eos-fei-crosscheck"
MODELS = {
    "Fei 2004": "platinum_fei_2004_bm3_mgd",
    "Fei 2007": "platinum_fei_2007_vinet_mgd",
    "Matsui 2009": matsui.THERMAL_RECORD,
    "Zhu v3": "platinum_zhu_2025_pvt",
}
COLORS = {
    "Fei 2004": "#808080",
    "Fei 2007": "#272727",
    "Matsui 2009": "#007b98",
    "Zhu v3": "#b44915",
}


def zhu_si(volume, temperature):
    """Derive pressure from the integrated gamma law and excess free energy."""
    v, t = np.broadcast_arrays(
        np.asarray(volume, float), np.asarray(temperature, float)
    )
    p = zhu.THERMAL_RECORDS["platinum_zhu_2025_pvt"]
    r = v / p["V0"]
    gamma = p["gamma0"] * (1 + p["a"] * (r ** p["b"] - 1))
    theta = p["theta0"] * np.exp(
        -p["gamma0"] * ((1 - p["a"]) * np.log(r) + p["a"] * (r ** p["b"] - 1) / p["b"])
    )

    def energy(temp, th):
        y = th / temp
        return (
            9
            * R
            * temp
            / y**3
            * quad(lambda x: x**3 / np.expm1(x), 0, y, epsabs=1e-11)[0]
        )

    increment = np.array(
        [energy(tt, th) - energy(300, th) for tt, th in zip(t.flat, theta.flat)]
    ).reshape(v.shape)
    vm = v * N_A * 1e-30 / 4
    phonon = gamma * increment / vm / 1e9
    # -dF_ex/dV_m, F_ex = -beta0/2 (V_m/V_m0)^m T^2.
    excess = p["beta0"] * p["m"] * r ** p["m"] / (2 * vm) * (t**2 - 300**2) / 1e9
    x = np.cbrt(r)
    cold = 3 * p["K0"] * (1 - x) / x**2 * np.exp(1.5 * (p["K0_prime"] - 1) * (1 - x))
    return cold + phonon + excess


def predict(name, volume, temperature):
    if name == "Fei 2004":
        return fei04.source_pressure(
            np.asarray(volume), np.asarray(temperature), "platinum"
        )
    if name == "Fei 2007":
        return fei07.source_pressure(volume, temperature)
    if name == "Matsui 2009":
        cold, phonon, _, _ = matsui.source_components(volume, temperature)
        return cold + phonon + matsui.electronic_increment(temperature)
    return zhu_si(volume, temperature)


def stats(residual):
    r = np.asarray(residual)
    return {
        "rows": int(r.size),
        "mean_gpa": float(r.mean()),
        "rmse_gpa": float(np.sqrt(np.mean(r**2))),
        "max_abs_gpa": float(np.max(abs(r))),
    }


def main():
    volume, temperature, pressure07, paired, _ = fei07.observations()
    mask_hot = temperature > 300
    mask_cold = np.arange(volume.size) < 36
    pv = np.array([float(r["volume_a3"]) for r in paired])
    pt = np.array([float(r["temperature_k"]) for r in paired])
    p04 = np.array([float(r["pressure_gpa"]) for r in paired])
    reported04 = np.array([float(r["platinum_fei_pressure_gpa"]) for r in paired])
    predictions = {name: predict(name, volume, temperature) for name in MODELS}
    grid_v = np.array([60.38, 0.95 * 60.38, 0.8 * 60.38, 0.7 * 60.38])[:, None]
    grid_t = np.array([300.0, 1000.0, 1473.0, 1873.0, 3000.0])[None, :]
    checks = {}
    for name, identifier in MODELS.items():
        native = get_eos_record(identifier).pressure(grid_v, grid_t)
        oracle = predict(name, grid_v, grid_t)
        difference = float(np.max(abs(native - oracle)))
        assert difference < 1e-8, (name, difference)
        checks[name] = {"native_vs_independent_si_max_gpa": difference}

    # Recheck the integrated theta/gamma identity independently by quadrature.
    identities = {}
    for name in ("Matsui 2009", "Zhu v3"):
        errors = []
        for r in (0.7, 0.8, 0.95):
            if name == "Matsui 2009":

                def gamma(x):
                    return 2.70 * np.exp(x * 1.10)

                theta_ratio = np.exp((2.70 - 2.70 * r**1.10) / 1.10)
            else:

                def gamma(x):
                    return 2.75 * (1 + 0.39 * (np.exp(5.1 * x) - 1))

                theta_ratio = np.exp(
                    -2.75 * ((1 - 0.39) * np.log(r) + 0.39 * (r**5.1 - 1) / 5.1)
                )
            errors.append(abs(np.log(theta_ratio) + quad(gamma, 0, np.log(r))[0]))
        identities[name] = float(max(errors))
        assert max(errors) < 1e-12

    datasets = {}
    for name in MODELS:
        paired_prediction = predict(name, pv, pt)
        datasets[name] = {
            "dewaele_36_cold_rows_revised_ruby": stats(
                predictions[name][mask_cold] - pressure07[mask_cold]
            ),
            "fei_35_hot_rows_2007_au_recalibration": stats(
                predictions[name][mask_hot] - pressure07[mask_hot]
            ),
            "fei_35_hot_rows_original_2004_calibration": stats(
                (paired_prediction - p04)[pt > 300]
            ),
            "fei_35_hot_rows_2004_printed_pt_calculated_outputs": stats(
                (paired_prediction - reported04)[pt > 300]
            ),
        }

    source = json.loads(
        (ROOT / "docs/data/fei-2007-platinum-source-curve-check.json").read_text(
            encoding="utf-8"
        )
    )
    curve_checks = {}
    low = {r["y_pixel"]: r for r in source["points"] if r["temperature_k"] == 1473}
    high = {r["y_pixel"]: r for r in source["points"] if r["temperature_k"] == 1873}
    sv = np.array([low[y]["volume_a3"] for y in low])
    separation = np.array(
        [high[y]["figure_pressure_gpa"] - low[y]["figure_pressure_gpa"] for y in low]
    )
    for name in MODELS:
        curve_checks[name] = {}
        for t in (300, 1473, 1873):
            points = [r for r in source["points"] if r["temperature_k"] == t]
            residual = predict(name, [r["volume_a3"] for r in points], t) - [
                r["figure_pressure_gpa"] for r in points
            ]
            curve_checks[name][str(t)] = stats(residual)
        curve_checks[name]["1873_minus_1473_separation"] = stats(
            predict(name, sv, 1873) - predict(name, sv, 1473) - separation
        )

    examples = []
    for v in (58.0, 55.0, 50.0, 48.0):
        values = {name: float(predict(name, v, 1873)) for name in MODELS}
        increments = {
            name: float(predict(name, v, 1873) - predict(name, v, 300))
            for name in MODELS
        }
        examples.append(
            {
                "volume_a3": v,
                "temperature_k": 1873,
                "pressure_gpa": values,
                "thermal_increment_gpa": increments,
            }
        )

    comparison_grid = {}
    sweep_v = np.linspace(0.69 * 60.38, 60.38, 1000)
    for t in (300, 1000, 2000, 3000):
        pm = predict("Matsui 2009", sweep_v, t)
        valid = (pm >= 0) & (pm <= 250)
        comparison_grid[str(t)] = {
            "max_abs_fei2007_minus_matsui_gpa_where_matsui_0_to_250gpa": float(
                np.max(abs(predict("Fei 2007", sweep_v[valid], t) - pm[valid]))
            )
        }

    # A qualified cold-baseline replacement isolates thermal differences.
    common_cold = predict("Fei 2007", volume, 300)
    common_stats = {
        name: stats(
            common_cold[mask_hot]
            + predictions[name][mask_hot]
            - predict(name, volume[mask_hot], 300)
            - pressure07[mask_hot]
        )
        for name in MODELS
    }
    report = {
        "scope": "Published coefficients only; no refitting or default/status changes. Original 2004 calibration and derived 2007 Au pressures are separate comparison targets. Figure curves are calculated outputs, not measurements. Equal-weight pressure residuals are descriptive, not a statistical ranking of absolute accuracy.",
        "model_identifiers": MODELS,
        "equation_checks": checks,
        "integrated_theta_identity_max_log_error": identities,
        "zhu_source_version": "Released-v3 optimizer/property constants: theta0=240 K, gamma0=2.75, a=0.39, b=5.1, beta0=0.002145 J/mol/K^2, m=0.65. First preprint and standalone calculator differ; no fresh download of original v3 code was available in this check.",
        "data_comparisons": datasets,
        "figure2_curve_comparisons": curve_checks,
        "examples": examples,
        "matsui_published_fei_comparison_claim_grid": comparison_grid,
        "common_fei_2007_cold_baseline_hot_diagnostic": common_stats,
        "common_baseline_qualification": "Replace each model's cold curve with Fei 2007 solely to isolate the thermal increment. These are diagnostic hybrids, not the published models.",
        "input_sha256": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                ROOT / "peritheos/data/datasets/platinum-fei-2004-table2.csv",
                ROOT
                / "peritheos/data/datasets/platinum-dewaele-2004-table1-compression.csv",
                ROOT
                / "peritheos/data/datasets/tsuchiya-kawamura-2002-table1-electronic-pressure.csv",
                ROOT / "docs/data/fei-2007-platinum-source-curve-check.json",
            ]
        },
    }
    OUT.with_suffix(".json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    with OUT.with_suffix(".csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "volume_a3",
                "temperature_k",
                "target_pressure_gpa",
                "target_calibration",
                *MODELS,
            ]
        )
        for index, (v, t, p, values) in enumerate(
            zip(volume, temperature, pressure07, zip(*predictions.values()))
        ):
            writer.writerow(
                [
                    v,
                    t,
                    p,
                    "Dewaele revised ruby"
                    if index < 36
                    else "Fei 2007 Au recalibration",
                    *values,
                ]
            )

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), layout="constrained")
    for name in MODELS:
        color = COLORS[name]
        axes[0, 0].plot(
            pressure07[mask_cold],
            predictions[name][mask_cold] - pressure07[mask_cold],
            ".",
            color=color,
            label=name,
        )
        axes[0, 1].scatter(
            pressure07[mask_hot],
            predictions[name][mask_hot] - pressure07[mask_hot],
            s=16,
            color=color,
            alpha=0.8,
        )
        v = np.linspace(48, 60.38, 170)
        thermal = predict(name, v, 1873) - predict(name, v, 300)
        axes[1, 0].plot(v, thermal, color=color, label=name)
        axes[1, 1].plot(
            v, predict(name, v, 1873) - predict(name, v, 1473), color=color, label=name
        )
    axes[0, 0].set(
        title="36 room-temperature Dewaele observations",
        xlabel="Revised ruby pressure (GPa)",
        ylabel="Model - target pressure (GPa)",
    )
    axes[0, 1].set(
        title="35 hot Fei observations; Au recalibrated to 2007",
        xlabel="Fei 2007 Au pressure (GPa)",
        ylabel="Model - target pressure (GPa)",
    )
    axes[1, 0].set(
        title="Thermal pressure at 1873 K; own cold curve subtracted",
        xlabel="Pt cell volume ($\\AA^3$)",
        ylabel="$P_{1873} - P_{300}$ (GPa)",
    )
    axes[1, 1].scatter(
        sv,
        separation,
        marker="o",
        s=45,
        facecolor="white",
        edgecolor="black",
        label="Fei Figure 2 curve spacing",
        zorder=5,
    )
    axes[1, 1].set(
        title="Fei Figure 2: thermal spacing at high compression",
        xlabel="Pt cell volume ($\\AA^3$)",
        ylabel="$P_{1873} - P_{1473}$ (GPa)",
        xlim=(47.8, 55.3),
    )
    for ax in axes.flat:
        ax.grid(alpha=0.16)
    for ax in axes[0]:
        ax.axhline(0, color="0.4", linewidth=0.8)
    axes[0, 0].legend(fontsize=9)
    axes[1, 1].legend(fontsize=8, loc="lower right")
    fig.suptitle(
        "Pt pressure scales checked against Fei: published coefficients, no refitting",
        fontsize=14,
    )
    fig.supxlabel(
        "Figure samples describe calculated curves, not hot measurements at these volumes. Zhu uses released-v3 optimizer constants.",
        fontsize=9,
    )
    fig.savefig(OUT.with_suffix(".png"), dpi=180)
    plt.close(fig)
    print(
        json.dumps(
            {
                "equations": checks,
                "data": datasets,
                "figure_spacing": {
                    k: v["1873_minus_1473_separation"] for k, v in curve_checks.items()
                },
                "examples": examples,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
