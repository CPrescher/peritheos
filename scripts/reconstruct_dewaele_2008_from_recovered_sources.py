"""Conditional ten-metal ruby reconstruction from recovered numerical sources.

Repository maintenance analysis, not an endorsed calibration or original author
workflow. Writes only the dedicated source-reconstruction output directory.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "output/dewaele-2008-ruby-construction"
OUT = BASE / "source-reconstruction"
METALS = ("Co", "Ni", "Mo", "Ag", "Al", "Cu", "Pt", "Ta", "W", "Au")
DAC_MAX = dict(
    Co=65.2,
    Ni=156.6,
    Mo=97.9,
    Ag=97.9,
    Al=153,
    Cu=153,
    Pt=93.6,
    Ta=93.6,
    W=153,
    Au=93.6,
)
COLORS = {
    "source_stated_range": "#1768ac",
    "source_Co_full_table": "#d17700",
    "literal_published_pairs": "#72767d",
    "common_DAC_support": "#8c4fa3",
}
LABELS = {
    "source_stated_range": "Recovered sources; stated fit ranges",
    "source_Co_full_table": "Recovered sources; Co full table",
    "literal_published_pairs": "Literal published Vinet pairs",
    "common_DAC_support": "Recovered sources; common DAC range",
}


def read_csv(name):
    return list(csv.DictReader((BASE / name).open()))


def save_json(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def save_csv(name, rows):
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def vinet(ratio, k0, kp):
    x = np.cbrt(np.asarray(ratio))
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def invert(function, pressures):
    return np.array(
        [
            brentq(lambda r: float(function(r)) - p, 0.25, 1.0, xtol=1e-14)
            for p in np.atleast_1d(pressures)
        ]
    )


def ruby(log_ratio, a=1920.0, b=9.61):
    return a / b * np.expm1(b * np.asarray(log_ratio))


def ruby_coordinate(pressure, a=1904.0, b=9.5):
    return np.log1p(np.asarray(pressure) * b / a) / b


def fit_vinet(pressure, ratio, k0):
    pressure, ratio = np.asarray(pressure), np.asarray(ratio)
    fit = least_squares(
        lambda q: vinet(ratio, k0, q[0]) - pressure,
        [5.0],
        bounds=(0.0, 15.0),
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
    )
    residual = vinet(ratio, k0, fit.x[0]) - pressure
    return {
        "K0_GPa": float(k0),
        "K0_prime": float(fit.x[0]),
        "samples": len(pressure),
        "sample_pressure_min_GPa": float(pressure.min()),
        "sample_pressure_max_GPa": float(pressure.max()),
        "pressure_RMSE_GPa": float(np.sqrt(np.mean(residual**2))),
        "pressure_max_abs_residual_GPa": float(np.max(np.abs(residual))),
    }


def source_models(n_grid=321, co_k0=184, co_full=False, temperature=293, pt="Holmes"):
    handbook = read_csv("aip-1972-table4d12-co-ni-ag-isotherms.csv")
    hixson = read_csv(
        "hixson-1992-293K-density-isotherms.csv"
        if temperature == 293
        else "hixson-1992-300K-isotherms.csv"
    )
    wang = read_csv("wang-2002-table2-generalized-eos.csv")
    result = {}
    for metal, k0 in [("Co", co_k0), ("Ni", 183), ("Ag", 101)]:
        end = 120.0 if metal == "Co" and co_full else DAC_MAX[metal]
        rows = [
            r
            for r in handbook
            if r["material"] == metal and float(r["P_GPa_derived"]) <= end
        ]
        pressure = np.array([float(r["P_GPa_derived"]) for r in rows])
        ratio = np.array([float(r["V_over_V0_source"]) for r in rows])
        result[metal] = fit_vinet(pressure, ratio, k0)
        result[metal].update(
            source="Kennedy-Keeler 1972 Table4d-12",
            T_K=298.15,
            selection_max_GPa=end,
            table_rows_available=sum(r["material"] == metal for r in handbook),
        )
    for metal, k0, end in [("Mo", 261, 97.9), ("W", 296, 160)]:
        rows = [
            r
            for r in hixson
            if r["material"] == metal and 0 <= float(r["P_GPa"]) <= end
        ]
        result[metal] = fit_vinet(
            [float(r["P_GPa"]) for r in rows],
            [float(r["V_over_V0_derived"]) for r in rows],
            k0,
        )
        result[metal].update(
            source="Hixson-Fritz 1992 tabulated isotherm",
            T_K=temperature,
            selection_max_GPa=end,
        )
    for metal, k0 in [("Cu", 133), ("Ta", 194), ("Au", 167), ("Pt", 277)]:
        row = next(r for r in wang if r["material"] == metal)

        def generalized(r, row=row):
            x = np.cbrt(r)
            y = 1 - x
            exponent = sum(
                float(row[key]) * y**i
                for i, key in enumerate(("eta", "beta", "xi", "delta"), 1)
            )
            return 3 * float(row["B0_GPa"]) * y / x**2 * np.exp(exponent)

        law = generalized
        source = "Wang2002 Eq8/TableII"
        if metal == "Pt" and pt == "Holmes":
            coefficients = json.loads(
                (BASE / "holmes-1989-300K-parameterization.json").read_text()
            )["printed_coefficients"]

            def law(r):
                return vinet(
                    r, coefficients["P_T_GPa"] / 3, 1 + coefficients["eta"] / 1.5
                )

            source = "Holmes1989 Eq11/TableIV exact P_T and eta"
        pressure = np.linspace(0, 160, n_grid)
        ratio = invert(law, pressure)
        assert np.max(np.abs(law(ratio) - pressure)) < 1e-8
        result[metal] = fit_vinet(pressure, ratio, k0)
        result[metal].update(
            source=source,
            T_K=300,
            selection_max_GPa=160,
            sampling="uniform pressure; analyst choice",
        )
    result["Al"] = {
        "K0_GPa": 73.46,
        "K0_prime": 4.52,
        "samples": None,
        "source": "2008 adopted preliminary Dorogokupets2007 pair",
        "fit_replayed": False,
        "reason": "Preliminary full thermal fit inputs unavailable",
    }
    return result


def published_models(co_k0=184):
    dac, rsw, v0 = {}, {}, {}
    for row in read_csv("source-table4-vinet-pairs.csv"):
        name = row["branch"]
        if name.startswith("Co") and name != f"Co{co_k0}":
            continue
        metal = "Co" if name.startswith("Co") else name
        k0 = float(row["K0_gpa_fixed_shared"])
        dac[metal] = (k0, float(row["DAC_K0_prime"]))
        rsw[metal] = (k0, float(row["RSW_K0_prime"]))
        v0[metal] = float(row["V0_atomic_a3_fixed_shared"])
    for row in read_csv("dewaele-2004-table2-vinet-branches.csv"):
        metal = row["material"]
        pair = float(row["K0_GPa"]), float(row["K0_prime"])
        if row["branch"] == "DAC_2004_fixed_K0" and metal != "Au":
            dac[metal] = pair
            v0[metal] = float(row["V0_A3_atom"])
        if row["branch"] == "RSW_Vinet_0_to_160GPa":
            rsw[metal] = pair
    correction = json.loads(
        (BASE / "dewaele-loubeyre-2007-stress-correction.json").read_text()
    )
    gold = next(row for row in correction["coefficients"] if row["material"] == "Au")
    dac["Au"] = gold["K0_GPa"], gold["K0_prime_corrected"]
    v0["Au"] = gold["V0_A3_atom"]
    rsw["Al"] = 73.46, 4.52
    return dac, rsw, v0


def targets(pressure_2004, dac, rsw):
    per_metal = {}
    for metal in METALS:
        ratio = invert(lambda r: vinet(r, *dac[metal]), pressure_2004)
        per_metal[metal] = vinet(ratio, *rsw[metal])
    return np.mean(list(per_metal.values()), axis=0), per_metal


def fit_ruby(p04, target):
    coordinate = ruby_coordinate(p04)
    fit = least_squares(
        lambda q: ruby(coordinate, *q) - target,
        [1904, 9.5],
        bounds=([1000, 1], [3000, 20]),
        x_scale=[1904, 9.5],
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
    )
    residual = ruby(coordinate, *fit.x) - target
    return {
        "A_GPa": float(fit.x[0]),
        "B": float(fit.x[1]),
        "target_mean_RMSE_GPa": float(np.sqrt(np.mean(residual**2))),
        "target_mean_max_abs_residual_GPa": float(np.max(np.abs(residual))),
        "fit_coordinate_2004_range_GPa": [float(min(p04)), float(max(p04))],
        "grid_points": len(p04),
    }


def pair_map(models):
    return {m: (models[m]["K0_GPa"], models[m]["K0_prime"]) for m in METALS}


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    dac, literal, v0 = published_models()
    source = source_models()
    full_co = source_models(co_full=True)
    dense = source_models(n_grid=1601)
    hot = source_models(temperature=300)
    wang_pt = source_models(pt="Wang")
    co190_dac, _, _ = published_models(190)
    co190 = source_models(co_k0=190)
    n_grid = 321
    grid = np.linspace(0, 160, n_grid)
    definitions = {
        "source_stated_range": (dac, pair_map(source), 160),
        "source_Co_full_table": (dac, pair_map(full_co), 160),
        "literal_published_pairs": (dac, literal, 160),
        "common_DAC_support": (dac, pair_map(source), min(DAC_MAX.values())),
        "source_Hixson_300K": (dac, pair_map(hot), 160),
        "source_Pt_Wang": (dac, pair_map(wang_pt), 160),
        "source_Co190": (co190_dac, pair_map(co190), 160),
    }
    fits, curves = {}, {}
    p_published = np.linspace(0, 160, 641)
    coordinate = ruby_coordinate(p_published, 1920, 9.61)
    p04_eval = ruby(coordinate, 1904, 9.5)
    for name, (d, s, end) in definitions.items():
        training = np.linspace(0, end, n_grid)
        target, _ = targets(training, d, s)
        fits[name] = fit_ruby(training, target)
        fit = fits[name]
        curve = ruby(coordinate, fit["A_GPa"], fit["B"])
        delta = curve - p_published
        fits[name].update(
            delta_at_published_160_GPa=float(delta[-1]),
            max_abs_delta_0_to_160_GPa=float(np.max(np.abs(delta))),
            comparison_axis="published2008 pressure at identical wavelength ratio",
        )
        curves[name] = delta
    mean, per_metal = targets(grid, dac, pair_map(source))
    mean_eval, _ = targets(p04_eval, dac, pair_map(source))
    _, literal_metal = targets(grid, dac, literal)
    target_rows = []
    for i, p in enumerate(grid):
        target_rows.append(
            {
                "P2004_GPa": p,
                "ln_lambda_over_lambda0_derived": ruby_coordinate(p),
                "mean_RSW_GPa": mean[i],
                **{f"{m}_RSW_GPa": per_metal[m][i] for m in METALS},
                "all_ten_DAC_measured_domains_overlap": p <= min(DAC_MAX.values()),
            }
        )
    save_csv("ten-metal-targets.csv", target_rows)
    comparison = []
    for i, pressure in enumerate(p_published):
        comparison.append(
            {
                "P_published2008_GPa": pressure,
                "P2004_same_wavelength_GPa": p04_eval[i],
                "mean_RSW_minus_published_GPa": mean_eval[i] - pressure,
                **{f"{name}_minus_published_GPa": curves[name][i] for name in curves},
            }
        )
    save_csv("calibration-comparison.csv", comparison)

    parameter_rows = []
    for metal in METALS:
        parameter_rows.append(
            {
                "metal": metal,
                "V0_A3_atom_shared_assumption": v0[metal],
                "DAC_K0_GPa": dac[metal][0],
                "DAC_K0_prime": dac[metal][1],
                "RSW_literal_K0_GPa": literal[metal][0],
                "RSW_literal_K0_prime": literal[metal][1],
                "RSW_refit_K0_GPa": source[metal]["K0_GPa"],
                "RSW_refit_K0_prime": source[metal]["K0_prime"],
                "RSW_refit_minus_literal_K0_prime": source[metal]["K0_prime"]
                - literal[metal][1],
                "RSW_source": source[metal]["source"],
                "DAC_max_observed_GPa": DAC_MAX[metal],
            }
        )
    save_csv("metal-eos-parameters.csv", parameter_rows)
    save_json(
        "source-isotherm-refits.json",
        {
            "stated_range": source,
            "Co_full_table": full_co["Co"],
            "Co190_stated_range": co190["Co"],
            "Co190_full_table": source_models(co_k0=190, co_full=True)["Co"],
            "Hixson_300K": {m: hot[m] for m in ("Mo", "W")},
            "Pt_Wang": wang_pt["Pt"],
        },
    )

    # Independent published graphical checkpoints are comparison-only. No graph
    # coordinates or output calibration coefficients enter the source refits.
    graph = json.loads(
        (ROOT / "docs/data/dewaele-2008-followup-source.json").read_text()
    )["figure2b"]["metals"]
    graph_comparison = {}
    for metal in METALS:
        points = np.array(graph[metal]["pressure_difference_points_gpa"])
        p, graph_dp = points.T
        ratio = invert(lambda r: vinet(r, *dac[metal]), p)
        predicted = vinet(ratio, *pair_map(source)[metal]) - p
        printed = vinet(ratio, *literal[metal]) - p
        graph_comparison[metal] = {
            "N": len(p),
            "source_refit_RMSE_GPa": float(
                np.sqrt(np.mean((predicted - graph_dp) ** 2))
            ),
            "literal_pairs_RMSE_GPa": float(
                np.sqrt(np.mean((printed - graph_dp) ** 2))
            ),
            "kind": "graphic source-equation comparison; not observations or fitting inputs",
        }
    save_json("figure2-independent-comparison.json", graph_comparison)

    dac_checks = []
    for row in read_csv("pressure-volume-input-inventory.csv"):
        metal = row["metal"]
        if metal == "Au" or row["source_eos_fit_included"] != "True":
            continue
        # Parentheses mark He-gauge pressures in the source, not an uncertainty.
        printed_pressure = float(row["revised_token"].strip("()"))
        residual = (
            float(
                vinet(float(row["volume_token_a3_per_atom"]) / v0[metal], *dac[metal])
            )
            - printed_pressure
        )
        dac_checks.append(
            {
                "metal": metal,
                "source": row["source"],
                "source_row": row["source_material_row"],
                "P2004_printed_GPa": printed_pressure,
                "model_minus_printed_pressure_GPa": residual,
                "He_gauge": row["helium_gauge"],
                "qualification": "EOS residual only; no measured ruby wavelength inferred from He",
            }
        )
    save_csv("DAC-observation-check.csv", dac_checks)

    dense_p = np.linspace(0, 160, 1601)
    dense_target, _ = targets(dense_p, dac, pair_map(dense))
    dense_fit = fit_ruby(dense_p, dense_target)
    coarse_delta = curves["source_stated_range"]
    dense_delta = ruby(coordinate, dense_fit["A_GPa"], dense_fit["B"]) - p_published
    inverse_error = max(
        float(
            np.max(
                np.abs(vinet(invert(lambda r: vinet(r, *dac[m]), grid), *dac[m]) - grid)
            )
        )
        for m in METALS
    )
    identity_error = 0.0
    for m in METALS:
        if dac[m][0] == literal[m][0]:
            ratio = invert(lambda r: vinet(r, *dac[m]), grid)
            x = np.cbrt(ratio)
            identity = grid * np.exp(1.5 * (literal[m][1] - dac[m][1]) * (1 - x))
            identity_error = max(
                identity_error, float(np.max(np.abs(identity - literal_metal[m])))
            )
    assert inverse_error < 1e-8 and identity_error < 1e-8
    figure3 = json.loads((BASE / "figure3-source.json").read_text())
    xy = np.array(figure3["points_mao_pressure_and_pressure_difference_gpa"])
    source_formula_error = (
        ruby(ruby_coordinate(xy[:, 0], b=7.665)) - xy[:, 0] - xy[:, 1]
    )
    verification = {
        "vinet_inverse_max_error_GPa": inverse_error,
        "independent_equal_K0_ratio_identity_max_error_GPa": identity_error,
        "published_2008_formula_vs_Figure3_max_error_GPa": float(
            np.max(np.abs(source_formula_error))
        ),
        "321_vs_1601_grid_max_calibration_change_GPa": float(
            np.max(np.abs(dense_delta - coarse_delta))
        ),
        "dense_grid_fit": dense_fit,
        "source_hashes": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                BASE / f
                for f in (
                    "aip-1972-table4d12-co-ni-ag-isotherms.csv",
                    "hixson-1992-293K-density-isotherms.csv",
                    "hixson-1992-300K-isotherms.csv",
                    "wang-2002-table2-generalized-eos.csv",
                    "dewaele-2004-table2-vinet-branches.csv",
                    "source-table4-vinet-pairs.csv",
                    "holmes-1989-300K-parameterization.json",
                    "dewaele-loubeyre-2007-stress-correction.json",
                    "pressure-volume-input-inventory.csv",
                    "figure3-source.json",
                )
            ]
            + [ROOT / "docs/data/dewaele-2008-followup-source.json"]
        },
    }
    save_json("verification.json", verification)
    summary = {
        "status": "conditional_source_based_reconstruction; author workflow not reproduced",
        "published_ruby": {"A_GPa": 1920, "B": 9.61},
        "fits": fits,
        "protocol": {
            "fixed_V0": "Each reduced isotherm normalized to its own source V0, then matched to DAC V/V0; reference normalization assumption explicit",
            "DAC": "Published ultrasonic-constrained Vinet branches; Au corrected2007; Co184 baseline",
            "RSW": "Recovered room-temperature data refit to fixed-K0 Vinet in stated ranges; Al adopted preliminary pair; Pt Holmes exact analytic curve",
            "source_refit_objective": "Unweighted squared pressure residuals; no missing pointwise errors invented",
            "analytic_sampling": "321 uniform pressure points over0-160GPa for Wang/Holmes; convergence checked at1601",
            "final_averaging": "Arithmetic mean of ten pressures at identical DAC volumes and a common P2004 coordinate; equal total weight per metal",
            "ruby_objective": "Unweighted squared pressure residual of ruby power law relative to that mean on321 uniform P2004 points",
            "coverage": "All-ten numerical extension to160GPa requires continuation of Co DAC above65.2, Mo/Ag above97.9, Pt/Ta/Au above93.6; Ni source table ends120GPa",
            "temperatures": "Handbook298.15K, Hixson293K baseline/300K alternative, Wang/Holmes300K; no unreported temperature conversion",
            "fit_uncertainty": "No statistical interval on A/B; source model errors, shared ancestors and covariance unavailable",
        },
        "limitations": [
            "Exact author objective, grid and weights remain undocumented",
            "Co full-table alternative is a declared range change, not a proven author selection",
            "Pt literal2004 table is preserved, not corrected to agree with graph",
            "Corrected Au observations unavailable; corrected source EOS used",
            "Al preliminary shock reduction not independently replayed",
            "Original shock-to-isotherm reductions not rerun from incomplete thermal/electronic inputs",
            "No new pressure standard or catalogue replacement",
        ],
        "common_DAC_P2004_end_GPa": min(DAC_MAX.values()),
        "common_DAC_end_on_published2008_axis_GPa": float(
            ruby(ruby_coordinate(min(DAC_MAX.values())))
        ),
    }
    save_json("reconstruction.json", summary)
    make_plots(
        p_published,
        p04_eval,
        curves,
        grid,
        per_metal,
        literal_metal,
        graph,
        dac,
        source,
    )
    write_report(
        summary,
        source,
        full_co,
        parameter_rows,
        graph_comparison,
        verification,
        p_published,
        curves,
    )
    print(json.dumps({"fits": fits, "verification": verification}, indent=2))


def make_plots(
    p_published, p04_eval, curves, grid, per_metal, literal_metal, graph, dac, source
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, layout="constrained")
    published_correction = p_published - p04_eval
    for name in COLORS:
        style = "--" if name == "literal_published_pairs" else "-"
        y = curves[name]
        if name == "common_DAC_support":
            end = float(ruby(ruby_coordinate(65.2)))
            mask = p_published <= end
            for axis, values in [(axes[0], y + published_correction), (axes[1], y)]:
                axis.plot(
                    p_published[mask],
                    values[mask],
                    color=COLORS[name],
                    lw=1.7,
                    label=LABELS[name] if axis == axes[0] else None,
                )
                axis.plot(
                    p_published[~mask],
                    values[~mask],
                    color=COLORS[name],
                    lw=1.7,
                    ls=":",
                )
        else:
            axes[0].plot(
                p_published,
                y + published_correction,
                color=COLORS[name],
                ls=style,
                lw=1.8,
                label=LABELS[name],
            )
            axes[1].plot(p_published, y, color=COLORS[name], ls=style, lw=1.8)
    axes[0].plot(
        p_published,
        published_correction,
        color="black",
        lw=2,
        label="Published 2008: A=1920, B=9.61",
    )
    axes[1].axhline(0, color="black", lw=1.4)
    axes[0].set(
        ylabel="Correction relative to 2004 ruby (GPa)",
        title="Dewaele 2008 ruby calibration: reconstruction from numerical sources",
    )
    axes[1].set(
        ylabel="Reconstructed minus published 2008 (GPa)",
        xlabel="Published 2008 pressure at the same ruby wavelength ratio (GPa)",
        xlim=(0, 160),
    )
    for axis in axes:
        axis.grid(alpha=0.2)
    axes[0].legend(loc="upper left", fontsize=8.5)
    axes[1].text(
        0.02,
        0.96,
        "Equal-metal averages; stated choices are conditional. Dotted purple: extrapolation of the common-range fit.",
        transform=axes[1].transAxes,
        va="top",
        fontsize=8.5,
    )
    for suffix in ("png", "svg"):
        fig.savefig(OUT / f"calibration-reconstruction.{suffix}", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(16, 7), layout="constrained")
    for axis, metal in zip(axes.flat, METALS):
        points = np.array(graph[metal]["pressure_difference_points_gpa"])
        axis.scatter(
            *points.T, color="black", s=9, label="Figure2 checkpoints", zorder=3
        )
        axis.plot(
            grid,
            per_metal[metal] - grid,
            color=COLORS["source_stated_range"],
            label="Recovered-source Vinet",
        )
        axis.plot(
            grid,
            literal_metal[metal] - grid,
            color=COLORS["literal_published_pairs"],
            ls="--",
            label="Literal published pair",
        )
        axis.axvspan(DAC_MAX[metal], 160, color="0.93", zorder=0)
        axis.axhline(0, color="0.7", lw=0.6)
        axis.set(
            title=f"{metal}: RSW K₀′={source[metal]['K0_prime']:.4f}",
            xlim=(0, 160),
            xlabel="2004 ruby pressure (GPa)",
            ylabel="P(RSW) − P(2004) (GPa)",
        )
        axis.grid(alpha=0.15)
    axes.flat[0].legend(fontsize=7, loc="upper left")
    fig.suptitle(
        "Per-metal reconstruction; Figure2 is an independent check. Grey: beyond DAC data coverage.",
        fontsize=13,
    )
    for suffix in ("png", "svg"):
        fig.savefig(OUT / f"per-metal-source-comparison.{suffix}", dpi=170)
    plt.close(fig)


def write_report(
    summary, source, full_co, parameters, graph_checks, checks, pressures, curves
):
    lines = [
        "# Dewaele 2008: reconstruction from recovered numerical sources",
        "",
        "## Result and status",
        "",
        "The available source isotherms support a reproducible **conditional reconstruction**. "
        "They do not reproduce the authors' published A=1920 GPa, B=9.61 exactly under "
        "the declared equal-metal averaging and pressure-residual objective. This is not "
        "a replacement pressure standard. The original calibration is not disproven.",
        "",
        "The principal reconstruction uses recovered room-temperature isotherms, "
        "ultrasonic-fixed bulk moduli and the stated fitting ranges. A second source "
        "reconstruction uses the full Co handbook table because it reproduces the "
        "printed Co derivative much more closely. That range change remains a diagnostic "
        "hypothesis. The literal printed Vinet-pair replay is retained separately, including "
        "the inconsistent Pt entry.",
        "",
        "| Construction | A (GPa) | B | Reconstructed minus published at 160 GPa |",
        "| --- | ---: | ---: | ---: |",
        "| Published 2008 | 1920 | 9.61 | 0 |",
    ]
    for name in COLORS:
        fit = summary["fits"][name]
        lines.append(
            f"| {LABELS[name]} | {fit['A_GPa']:.4f} | {fit['B']:.6f} | {fit['delta_at_published_160_GPa']:+.4f} GPa |"
        )
    lines += [
        "",
        "For the common-range fit, values beyond the common DAC domain are "
        "extrapolations. All differences compare identical ruby wavelength ratios; "
        "the horizontal coordinate is the published 2008 pressure, not an independently "
        "known true pressure.",
        "",
        "![Calibration reconstruction](calibration-reconstruction.png)",
        "",
        "## Numerical construction",
        "",
        "1. Preserve the original source tables and published parameters. Fit recovered "
        "shock-reduced P versus V/V0 values to Vinet with K0 fixed to the reported ultrasonic value.",
        "2. Use the published fixed-K0 DAC Vinet branches: 2008 Table IV for Co/Ni/Mo/Ag; "
        "2004 Table II for Al/Cu/Pt/Ta/W; stress-corrected 2007 Au. The main Co branch uses K0=184 GPa.",
        "3. At each P04 in 0–160 GPa, invert each DAC EOS to obtain V/V0, then evaluate "
        "its reduced-shock EOS at that same volume. Take the arithmetic mean of ten pressures. "
        "Zn and Fe are not members of the ten-metal average.",
        "4. Convert P04 to the ruby coordinate using ln(r)=ln(1+9.5*P04/1904)/9.5. "
        "Fit P=A/B*expm1(B*ln(r)) to that mean, with A and B free.",
        "",
        "All pressure-residual fits are unweighted least squares. Each metal receives "
        "equal weight in the final average. Analytic source curves and final ruby fits use "
        "321 uniformly spaced pressure points. These sampling and weighting choices are "
        "explicit analyst assumptions; the papers do not specify them sufficiently for an author-workflow replay.",
        "",
        "Each reduced isotherm is normalized to its own source reference volume and "
        "matched to the DAC V/V0. This adopts a common reference normalization; no "
        "absolute shock-derived ambient atomic volume or thermal correction is invented.",
        "",
        "## Recovered input use and coefficient comparison",
        "",
        "| Metal | Fixed RSW K0 (GPa) | Literal source K0′ | Recovered-source fit K0′ |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in parameters:
        lines.append(
            f"| {row['metal']} | {row['RSW_refit_K0_GPa']:.4g} | {row['RSW_literal_K0_prime']:.4f} | {row['RSW_refit_K0_prime']:.6f} |"
        )
    lines += [
        "",
        "Source selection:",
        "",
        "- **Co/Ni/Ag:** handbook 25 °C (298.15 K) values, with source shock-pressure "
        "limits 65.2/156.6/97.9 GPa respectively. This retains 36/44/42 tabulated "
        "points. Ni's table ends at 120 GPa. All 136 recovered table values remain "
        "preserved; entries outside the selected fit range are not silently included.",
        "- **Mo/W:** the 293 K Hixson–Fritz density tables, fitted to 97.9 and 160 GPa "
        "respectively: 10/17 points. The separate 300 K PVT tables are a declared sensitivity check.",
        "- **Cu/Ta/Au:** Wang Eq.8/Table II generalized universal EOS, sampled over "
        "0–160 GPa then refit to ordinary Vinet. Quartic exponent terms and their signs are retained.",
        "- **Pt:** Holmes Eq.11/Table IV, using exact P_T=798.31 GPa and eta=7.2119, "
        "then refit with K0=277 GPa. Wang's Pt parameterization is a separate source-identity check.",
        "- **Al:** the expressly adopted preliminary Dorogokupets pair 73.46/4.52. "
        "Its original preliminary thermal fit is unavailable, so this is adopted, not independently refitted.",
        "",
        "The recovered raw shock-state tables provide provenance but cannot be substituted "
        "for room-temperature isotherms without their thermal reductions. Those original "
        "reductions were not rerun: intermediate electronic-energy arrays, selected shock "
        "rows and some thermal inputs remain missing. For the same reason, temperatures "
        "293, 298.15 and 300 K are explicitly retained rather than silently equated.",
        "",
        "## Two resolved diagnostics and remaining non-parity",
        "",
        f"**Co range:** the source-stated-range fit gives K0′={source['Co']['K0_prime']:.6f}, "
        f"whereas all 44 handbook Co points through 120 GPa give {full_co['Co']['K0_prime']:.6f}, "
        "which rounds to the printed 4.93. This is evidence for a possible wider Co "
        "fit domain, not proof of the original selection. No continuous cutoff search was performed.",
        "",
        "**Pt source inconsistency:** the 2004 table prints 277/4.66, with ref.23 "
        "pointing to Wang, while Figure3 labels the curve Holmes. The recovered Holmes "
        f"curve gives K0′={source['Pt']['K0_prime']:.6f} when refitted at fixed K0=277 GPa. "
        f"Against the independent 2008 Figure2 Pt checkpoints, the literal pair has "
        f"RMSE {graph_checks['Pt']['literal_pairs_RMSE_GPa']:.3f} GPa, versus "
        f"{graph_checks['Pt']['source_refit_RMSE_GPa']:.3f} GPa for the recovered-source fit. "
        "The printed value is not silently corrected. Close agreement of the literal "
        "ten-metal average with the final ruby curve is therefore insufficient to validate its inputs.",
        "",
        "![Per-metal comparisons](per-metal-source-comparison.png)",
        "",
        "Other small input differences remain, including Mo's recovered-table derivative "
        "versus the printed 4.18. Changing arbitrary weights or row masks to force A/B "
        "parity has not been attempted.",
        "",
        "## Sensitivity and coverage",
        "",
        "| Declared change from principal reconstruction | A (GPa) | B | Difference at 160 GPa |",
        "| --- | ---: | ---: | ---: |",
    ]
    for name, title in [
        ("source_Hixson_300K", "Hixson 300 K tables"),
        ("source_Pt_Wang", "Wang Pt source"),
        ("source_Co190", "Co K0=190 GPa branch"),
    ]:
        fit = summary["fits"][name]
        lines.append(
            f"| {title} | {fit['A_GPa']:.4f} | {fit['B']:.6f} | {fit['delta_at_published_160_GPa']:+.4f} GPa |"
        )
    lines += [
        "",
        "These cases are not a statistical confidence interval or an exhaustive error "
        "budget. No interval for A/B is claimed: source-model errors, shared calibrations "
        "and parameter covariances are not supplied.",
        "",
        "All-ten DAC support overlaps only through P04=65.2 GPa. Full-range comparisons "
        "therefore continue fitted EOS beyond several metals' measured domains, and the "
        "Ni source table is continued above 120 GPa. Corrected Au numerical observations "
        "remain unavailable; its published corrected parameterization is used.",
        "",
        "## Comparison at selected pressures",
        "",
        "| Published 2008 pressure | Principal reconstruction minus published | Co full-table alternative minus published | Literal pairs minus published |",
        "| ---: | ---: | ---: | ---: |",
    ]
    for p in (40, 80, 120, 160):
        i = int(np.argmin(abs(pressures - p)))
        lines.append(
            f"| {p} GPa | {curves['source_stated_range'][i]:+.4f} GPa | {curves['source_Co_full_table'][i]:+.4f} GPa | {curves['literal_published_pairs'][i]:+.4f} GPa |"
        )
    lines += [
        "",
        "## Verification and reproducibility",
        "",
        f"Vinet inversion and an independent equal-K0 pressure-ratio identity agree within "
        f"{max(checks['vinet_inverse_max_error_GPa'], checks['independent_equal_K0_ratio_identity_max_error_GPa']):.2g} GPa. "
        f"Increasing both analytic and ruby-fit grids from321 to1601 changes the final "
        f"0–160 GPa pressure curve by at most {checks['321_vs_1601_grid_max_calibration_change_GPa']:.6f} GPa.",
        "",
        f"The published 1920/9.61 equation agrees with independently recovered Figure3 "
        f"checkpoints to {checks['published_2008_formula_vs_Figure3_max_error_GPa']:.5f} GPa. "
        "That graph is a validation target only and never enters the source-based fit.",
        "",
        "The published DAC EOS are also evaluated against the recovered observation "
        "rows in `DAC-observation-check.csv`; historical uncorrected Au rows are excluded. "
        "He-derived Ni pressures retain their flag and do not acquire invented ruby wavelengths.",
        "",
        "Run from the repository or extracted numerical bundle root:",
        "",
        "```bash",
        "uv run --no-project --with scipy --with matplotlib python scripts/reconstruct_dewaele_2008_from_recovered_sources.py",
        "```",
        "",
        "Inputs are checksummed in `verification.json`. Numerical targets, "
        "parameters, final curves and per-metal graphical comparisons are exported alongside "
        "this report. No catalogue record, source observation, earlier reconstruction "
        "or pressure standard was changed.",
        "",
        "## Primary sources",
        "",
        "- [Dewaele et al. 2008](https://doi.org/10.1103/PhysRevB.78.104102), Sec.IV, TableIV, Figures2–3.",
        "- [Dewaele et al. 2004](https://doi.org/10.1103/PhysRevB.70.094112), TableII, Sec.III.",
        "- Kennedy–Keeler, AIP Handbook, third edition (1972), Sec.4d-2 and Table4d-12, pp.4-96–4-101.",
        "- [Hixson–Fritz 1992](https://doi.org/10.1063/1.351203), TablesIV–VII.",
        "- [Wang et al. 2002](https://doi.org/10.1063/1.1518781), Eq.8/TableII.",
        "- [Holmes et al. 1989](https://doi.org/10.1063/1.344177), Eq.11/TableIV.",
        "- [Dewaele–Loubeyre 2007](https://doi.org/10.1080/08957950701659627), Sec.5.",
        "- [Dorogokupets–Oganov 2007](https://doi.org/10.1103/PhysRevB.75.024115), preliminary Al coefficients, p.024115-4.",
        "",
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines))


if __name__ == "__main__":
    run()
