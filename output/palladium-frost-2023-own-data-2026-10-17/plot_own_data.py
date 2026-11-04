"""Evaluate the published Frost Vinet coefficients without refitting."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from peritheos.eos.rt import Vinet

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SOURCE = (
    ROOT / "peritheos/data/datasets/palladium-frost-2023-tables-i-iii-compression.csv"
)
EXPECTED_SHA256 = "92e4ece39614b1a9f980bff32ae3986083b06bd9673ef0e615f4d27bd39fa75c"


def summary(dv, dp, v):
    return {
        "volume_rmse_a3_atom": float(np.sqrt(np.mean(dv**2))),
        "volume_mean_residual_a3_atom": float(np.mean(dv)),
        "volume_max_absolute_residual_a3_atom": float(np.max(np.abs(dv))),
        "relative_volume_rmse_percent": float(
            100 * np.sqrt(np.mean((dv / (v / 4)) ** 2))
        ),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(dp**2))),
        "pressure_mean_residual_gpa": float(np.mean(dp)),
        "pressure_max_absolute_residual_gpa": float(np.max(np.abs(dp))),
    }


def main():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED_SHA256
    with SOURCE.open() as stream:
        rows = list(csv.DictReader(stream))
    assert Counter(r["loading"] for r in rows) == {"BX18": 26, "BX32": 20, "HP": 47}
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    sp = np.array([float(r["pressure_error_gpa"]) for r in rows])
    v = np.array([float(r["volume_a3_conventional_cell"]) for r in rows])
    sv = np.array([float(r["volume_error_a3_conventional_cell"]) for r in rows])
    a = np.array([float(r["lattice_parameter_a"]) for r in rows])
    sa = np.array([float(r["lattice_parameter_error_a"]) for r in rows])
    np.testing.assert_allclose(v, a**3, rtol=0, atol=5.1e-10)
    np.testing.assert_allclose(sv, 3 * a**2 * sa, rtol=0, atol=5.1e-10)
    parameters = {"V0": 58.678, "K0": 189.3, "K0_prime": 5.473}
    eos = Vinet(**parameters)
    model_v = np.asarray(eos.volume(p))
    dv = (v - model_v) / 4
    dp = np.asarray(eos.pressure(v)) - p
    # First-order propagation of the reported coordinate errors into vertical
    # volume residuals. This does not include EOS-parameter covariance.
    sigma_dv = np.sqrt(sv**2 + (model_v / eos.bulk_modulus(model_v) * sp) ** 2) / 4
    x = (v / parameters["V0"]) ** (1 / 3)
    independent_p = (
        3
        * parameters["K0"]
        * (1 - x)
        / x**2
        * np.exp(1.5 * (parameters["K0_prime"] - 1) * (1 - x))
    )
    np.testing.assert_allclose(independent_p - p, dp, rtol=0, atol=1e-10)
    np.testing.assert_allclose(eos.pressure(model_v), p, rtol=0, atol=1e-8)

    # The supplementary ambient measurement is separate from Tables I–III.
    ambient_v = 3.89088**3 / 4
    ambient_sv = 3 * 3.89088**2 * 0.00002 / 4
    ambient_residual = ambient_v - parameters["V0"] / 4
    report = {
        "source_doi": "10.1063/5.0161038",
        "source_tables": "Supplementary Tables I–III",
        "data_sha256": EXPECTED_SHA256,
        "published_parameters_cell_units": parameters,
        "temperature_k": 300,
        "refit_performed": False,
        "observations": len(rows),
        "observed_pressure_range_gpa": [float(p.min()), float(p.max())],
        "residual_definitions": {
            "volume": "(V_observed - V_published(P_observed))/4, in angstrom^3/atom",
            "pressure": "P_published(V_observed) - P_observed, in GPa",
        },
        "error_bars": "Reported coordinate errors; residual bars combine volume and pressure contributions by independent first-order propagation. Parameter covariance, run correlations, and calibration systematics are not included.",
        "all_93_rows": summary(dv, dp, v),
        "loadings": {},
        "separate_ambient_measurement": {
            "used_in_93_row_statistics": False,
            "a_angstrom": 3.89088,
            "a_error_angstrom": 0.00002,
            "volume_a3_atom": ambient_v,
            "volume_residual_a3_atom": ambient_residual,
        },
    }
    config = {
        "BX18": ("#3575b4", "o", "BX18 · W marker · 26 points"),
        "BX32": ("#da8b22", "s", "BX32 · W marker · 20 points"),
        "HP": ("#26896a", "^", "High-pressure run · Au marker · 47 points"),
    }
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    fig, (ax, residual_ax) = plt.subplots(
        2,
        1,
        figsize=(11.5, 8.6),
        sharex=True,
        height_ratios=[1.65, 1],
        layout="constrained",
    )
    grid = np.linspace(0, 185, 600)
    ax.plot(
        grid,
        eos.volume(grid) / 4,
        color="#202020",
        lw=2,
        label="Published Vinet EOS",
        zorder=2,
    )
    for loading, (color, marker, label) in config.items():
        select = np.array([r["loading"] == loading for r in rows])
        report["loadings"][loading] = {
            "observations": int(select.sum()),
            "pressure_range_gpa": [float(p[select].min()), float(p[select].max())],
            **summary(dv[select], dp[select], v[select]),
        }
        ax.errorbar(
            p[select],
            v[select] / 4,
            xerr=sp[select],
            yerr=sv[select] / 4,
            fmt=marker,
            color=color,
            ms=4,
            elinewidth=0.7,
            capsize=1.5,
            label=label,
            zorder=3,
        )
        residual_ax.errorbar(
            p[select],
            dv[select],
            xerr=sp[select],
            yerr=sigma_dv[select],
            fmt=marker,
            color=color,
            ms=4,
            elinewidth=0.8,
            capsize=1.5,
            zorder=3,
        )
    ax.errorbar(
        [0],
        [ambient_v],
        yerr=[ambient_sv],
        fmt="D",
        mfc="white",
        mec="#666666",
        ecolor="#666666",
        ms=6,
        label="Separate ambient measurement (not in fit table)",
        zorder=4,
    )
    residual_ax.errorbar(
        [0],
        [ambient_residual],
        yerr=[ambient_sv],
        fmt="D",
        mfc="white",
        mec="#666666",
        ecolor="#666666",
        ms=6,
        zorder=4,
    )
    ax.set(
        ylabel="Atomic volume (Å³/atom)",
        title="Published curve against all 93 compression observations",
    )
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    ax.text(
        0.98,
        0.53,
        "V₀ = 58.678 Å³/cell · K₀ = 189.3 GPa · K₀′ = 5.473\nPublished coefficients unchanged; no refit",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
    )
    residual_ax.axhline(0, color="#202020", lw=1)
    residual_ax.set(
        xlabel="Pressure (GPa)",
        ylabel="Observed − calculated volume\n(Å³/atom)",
        title="Volume residuals; error bars propagate reported coordinate errors",
        xlim=(-3, 188),
    )
    residual_ax.text(
        0.98,
        0.95,
        f"Volume RMS: {np.sqrt(np.mean(dv**2)):.4f} Å³/atom\nPressure-equivalent RMS: {np.sqrt(np.mean(dp**2)):.3f} GPa",
        transform=residual_ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )
    for panel in (ax, residual_ax):
        panel.grid(alpha=0.17)
    fig.suptitle(
        "Palladium at 300 K — Frost et al. (2023), published Vinet EOS",
        fontsize=15,
        weight="bold",
    )
    fig.supxlabel(
        "Source: supplementary Tables I–III. Statistics use the 93 compression rows only; the separate ambient value is shown for context.",
        fontsize=9,
    )
    fig.savefig(OUT / "frost-vinet-own-data.png", dpi=180)
    fig.savefig(OUT / "frost-vinet-own-data.svg")
    plt.close(fig)
    (OUT / "comparison.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    with (OUT / "residuals.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "loading",
                "pressure_marker",
                "pressure_gpa",
                "observed_volume_a3_atom",
                "predicted_volume_a3_atom",
                "Vobs_minus_Vcalc_a3_atom",
                "propagated_coordinate_error_a3_atom",
                "Pcalc_minus_Pobs_gpa",
            ]
        )
        writer.writerows(
            [
                r["loading"],
                r["pressure_marker"],
                p[i],
                v[i] / 4,
                model_v[i] / 4,
                dv[i],
                sigma_dv[i],
                dp[i],
            ]
            for i, r in enumerate(rows)
        )
    print(json.dumps(report["all_93_rows"], indent=2))


if __name__ == "__main__":
    main()
