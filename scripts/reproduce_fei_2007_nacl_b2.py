"""Recover Fei (2007) NaCl figure vectors and run conditional EOS diagnostics.

Run ``python -m scripts.reproduce_fei_2007_nacl_b2`` against the archived
vector JSON. Add ``--pdf PATH`` to recover the vectors again (pdfplumber is
required only for extraction). Calculated curves never enter regressions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

from peritheos.materials import NACL_B2_FEI_2007

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
SOURCE = DATA / "nacl-b2-fei-2007-source.json"
OUTPUT = ROOT / "docs/data/fei-2007-nacl-b2-reproduction.json"
PARAMETERS = np.array([41.35, 26.86, 5.25, 1.70, 0.5, 290.0])
NAMES = ("V0", "K0", "K0_prime", "gamma0", "q", "theta0")
PDF_SHA256 = "19109c9685f6862c733ffc14ada22077d61553152b18dab13eb295e8d1da3306"


def pressure(volume, temperature, coefficients=PARAMETERS, law="printed"):
    """Independent SI Vinet/MGD replay, two atoms per NaCl formula unit."""
    v0, k0, kp, gamma0, q, theta0 = coefficients
    volume, temperature = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    logarithm = np.log(volume / v0)
    gamma = gamma0 * np.exp(q * logarithm)
    if law == "printed":
        exponent = -gamma * logarithm
    elif law == "integrated":
        exponent = (
            -gamma0 * logarithm
            if abs(q) < 1e-10
            else -gamma0 * np.expm1(q * logarithm) / q
        )
    else:
        raise ValueError(law)
    theta = theta0 * np.exp(exponent)

    def energy(t, th):
        y = th / t
        integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11, epsrel=1e-11)[
            0
        ]
        return 18 * R * t * integral / y**3

    increment = np.array(
        [
            energy(t, th) - energy(300.0, th)
            for t, th in zip(temperature.flat, theta.flat)
        ]
    ).reshape(volume.shape)
    x = np.exp(logarithm / 3)
    cold = 3 * k0 * (1 - x) * np.exp(1.5 * (kp - 1) * (1 - x)) / x**2
    return cold + gamma * increment / (volume * N_A * 1e-30) / 1e9


def recover(pdf):
    """Extract marker centers and curve vertices from the verified primary PDF."""
    import pdfplumber

    pdf = Path(pdf)
    digest = hashlib.sha256(pdf.read_bytes()).hexdigest()
    if digest != PDF_SHA256:
        raise ValueError(
            "PDF hash differs; inspect layout before reusing vector indices"
        )
    with pdfplumber.open(pdf) as document:
        page = document.pages[2]
        curves = page.curves
        # Major ticks, independently read from the PDF line objects and labels.
        axes = {
            "3": {
                "pressure": [[356.228, 0], [518.843, 120]],
                "volume": [[170.646, 20], [53.223, 40]],
            },
            "4": {
                "pressure": [[361.642, 0], [517.188, 125]],
                "volume": [[659.721, 20], [543.991, 40]],
            },
        }

        def physical(fig, x, y):
            ax = axes[str(fig)]
            (x0, p0), (x1, p1) = ax["pressure"]
            (y0, v0), (y1, v1) = ax["volume"]
            return p0 + (x - x0) * (p1 - p0) / (x1 - x0), v0 + (y - y0) * (v1 - v0) / (
                y1 - y0
            )

        observations = []
        # Figure 3: only the 12 filled circles of this study. The open Sata/Ono
        # symbols are distinct PDF paths and are deliberately not selected.
        # Figure 4: each white-filled cold circle has a duplicate blue outline;
        # select one path per circle. Hot circles remain separate vector paths,
        # including the two touching markers (indices 133 and 134).
        selections = [
            (3, range(84, 96), 300),
            (4, range(102, 126, 2), 300),
            (4, range(126, 138), 1000),
        ]
        for fig, indices, temperature in selections:
            for number, index in enumerate(indices, 1):
                c = curves[index]
                assert len(c["pts"]) == 5 and 2.8 < c["width"] < 3.0
                x, y = (c["x0"] + c["x1"]) / 2, (c["top"] + c["bottom"]) / 2
                p, v = physical(fig, x, y)
                observations.append(
                    {
                        "row_id": f"fig{fig}_{temperature}k_{number:02d}",
                        "figure": fig,
                        "temperature_k": temperature,
                        "pressure_gpa": p,
                        "volume_a3": v,
                        "pdf_curve_index": index,
                        "x_pt": x,
                        "y_pt": y,
                        "marker_box_xyxy_pt": [c["x0"], c["top"], c["x1"], c["bottom"]],
                    }
                )
        calculated_curves = []
        for index, temperature in [(28, 300), (99, 300), (100, 1000), (101, 2000)]:
            fig = 3 if index == 28 else 4
            c = curves[index]
            calculated_curves.append(
                {
                    "figure": fig,
                    "temperature_k": temperature,
                    "pdf_curve_index": index,
                    "linewidth_pt": c["linewidth"],
                    "points": [
                        {
                            "x_pt": x,
                            "y_pt": y,
                            "pressure_gpa": physical(fig, x, y)[0],
                            "volume_a3": physical(fig, x, y)[1],
                        }
                        for x, y in c["pts"]
                    ],
                }
            )
    source = {
        "doi": "10.1073/pnas.0609013104",
        "primary_source": "https://www.pnas.org/doi/10.1073/pnas.0609013104",
        "source_pdf_sha256": digest,
        "source_page": 9184,
        "coordinate_basis": "PDF points from page top left",
        "axis_calibration": axes,
        "extraction": "pdfplumber marker path bounding-box centers; duplicate cold outlines removed. Figure 3 filled study circles only; comparison symbols excluded. Figure 4 touching hot circles are separate vector objects.",
        "precision": "Vector layout precision, not original experimental precision. No original numerical table, paired Pt cells, P/V/T uncertainties, fit weights or covariance recovered. No experimental error bars are plotted. Blank CSV error/calibrant cells mean unavailable, not zero.",
        "duplicate_scope": "Figure 3 and Figure 4 cold series depict the same 12 observations at slightly different graphical coordinates; never concatenate them in a fit.",
        "retrieval_audit": {
            "date": "2026-10-09",
            "publisher_full_text": "Retrieved through web tool; Table 1, Figures 3-4, and Methods verified.",
            "numerical_supplement": "No article-specific numerical supplement recovered. Generic publisher PDF and Supporting Information control is not evidence that a supplement exists. Direct publisher PDF/supplement HTTP requests returned 403; PMC XML and Europe PMC fullTextXML did not supply article XML. No claim of definitive supplement absence.",
            "pdf_copy": "Existing local copy, full 5-page paper; SHA256 matches the independently used primary PDF. Contains model Table 1 but no experimental NaCl/Pt table.",
        },
        "observations": observations,
        "calculated_curves": calculated_curves,
    }
    SOURCE.write_text(json.dumps(source, indent=2) + "\n", encoding="utf-8")
    write_csv(source)
    return source


def write_csv(source):
    fields = [
        "row_id",
        "figure",
        "temperature_k",
        "pressure_gpa",
        "volume_a3",
        "pressure_error_gpa",
        "volume_error_a3",
        "temperature_error_k",
        "platinum_volume_a3",
        "x_pt",
        "y_pt",
        "pdf_curve_index",
    ]
    for fig in (3, 4):
        path = DATA / f"nacl-b2-fei-2007-figure{fig}-digitized.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields)
            writer.writeheader()
            for row in source["observations"]:
                if row["figure"] == fig:
                    writer.writerow({k: row.get(k, "") for k in fields})


def subset(source, figure, temperature):
    rows = [
        r
        for r in source["observations"]
        if r["figure"] == figure and r["temperature_k"] == temperature
    ]
    return tuple(
        np.array([r[k] for r in rows], float)
        for k in ("volume_a3", "temperature_k", "pressure_gpa")
    )


def fit(
    volume,
    temperature,
    target,
    free,
    base=PARAMETERS,
    law="printed",
    residual_direction="pressure",
):
    base = np.array(base, copy=True)

    def unpack(x):
        c = base.copy()
        c[list(free)] = x
        return c

    def residual(x):
        coefficients = unpack(x)
        if residual_direction == "pressure":
            return pressure(volume, temperature, coefficients, law) - target
        predicted = np.array(
            [
                brentq(
                    lambda v: float(pressure(v, t, coefficients, law)) - p,
                    0.3 * coefficients[0],
                    1.2 * coefficients[0],
                )
                for t, p in zip(temperature, target)
            ]
        )
        return predicted - volume

    result = least_squares(
        residual,
        base[list(free)],
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
    )
    covariance = np.linalg.inv(result.jac.T @ result.jac) * (
        result.fun @ result.fun / (len(volume) - len(free))
    )
    return unpack(result.x), {
        "observations": len(volume),
        "law": law,
        "residual_direction": residual_direction,
        "weighting": f"equal {residual_direction} weights; no experimental error weights recovered",
        "free_parameters": [NAMES[i] for i in free],
        "fixed_parameters": {
            NAMES[i]: float(base[i]) for i in range(6) if i not in free
        },
        "parameters": dict(zip([NAMES[i] for i in free], result.x.tolist())),
        "conditional_standard_errors": dict(
            zip([NAMES[i] for i in free], np.sqrt(np.diag(covariance)).tolist())
        ),
        "rmse_gpa" if residual_direction == "pressure" else "rmse_a3": float(
            np.sqrt(np.mean(result.fun**2))
        ),
        "solver_success": bool(result.success),
    }


def reproduce(source=None):
    source = source or json.loads(SOURCE.read_text(encoding="utf-8"))
    cold = subset(source, 3, 300)
    cold4 = subset(source, 4, 300)
    hot = subset(source, 4, 1000)
    c, cold_fit = fit(*cold, free=(1, 2))
    _, cold4_fit = fit(*cold4, free=(1, 2))
    _, hot_fit = fit(*hot, free=(4,))
    retained = np.ones(len(hot[0]), dtype=bool)
    retained[[7, 8]] = False
    _, partial_hot_fit = fit(*(a[retained] for a in hot), free=(4,))
    _, integrated_fit = fit(*hot, free=(4,), law="integrated")
    _, staged_hot_fit = fit(*hot, free=(4,), base=c)
    _, volume_hot_fit = fit(*hot, free=(4,), residual_direction="volume")
    published = {}
    for name, (v, t, p) in [
        ("figure3_300k", cold),
        ("figure4_300k", cold4),
        ("figure4_1000k", hot),
    ]:
        residual = pressure(v, t) - p
        published[name] = {
            "observations": len(v),
            "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "mean_residual_gpa": float(residual.mean()),
            "max_abs_residual_gpa": float(np.max(abs(residual))),
            "pressure_residuals_gpa": residual.tolist(),
        }
    curve_checks = []
    for curve in source["calculated_curves"]:
        v = np.array([p["volume_a3"] for p in curve["points"]])
        p = np.array([p["pressure_gpa"] for p in curve["points"]])
        ax = source["axis_calibration"][str(curve["figure"])]["pressure"]
        pts_per_gpa = (ax[1][0] - ax[0][0]) / (ax[1][1] - ax[0][1])
        for law in ("printed", "integrated"):
            residual = pressure(v, curve["temperature_k"], law=law) - p
            curve_checks.append(
                {
                    "figure": curve["figure"],
                    "temperature_k": curve["temperature_k"],
                    "law": law,
                    "vertices": len(v),
                    "role": "calculated curve, excluded from fits",
                    "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
                    "max_abs_residual_gpa": float(np.max(abs(residual))),
                    "max_horizontal_difference_pt": float(
                        np.max(abs(residual)) * pts_per_gpa
                    ),
                    "printed_linewidth_pt": curve["linewidth_pt"],
                }
            )
    v = np.array([20.5, 23, 27, 35, 41.35])[:, None]
    t = np.array([300, 1000, 2000])[None, :]
    native_difference = float(
        np.max(abs(NACL_B2_FEI_2007.pressure(v, t) - pressure(v, t)))
    )
    # Illustration of coordinate sensitivity, not source-error weighting.
    tick_x = source["axis_calibration"]["4"]["pressure"]
    tick_y = source["axis_calibration"]["4"]["volume"]
    dp = 0.357 * (tick_x[1][1] - tick_x[0][1]) / (tick_x[1][0] - tick_x[0][0])
    dv = 0.357 * abs((tick_y[1][1] - tick_y[0][1]) / (tick_y[1][0] - tick_y[0][0]))
    offset_q = []
    for ps in (-1, 1):
        for vs in (-1, 1):
            _, f = fit(hot[0] + vs * dv, hot[1], hot[2] + ps * dp, free=(4,))
            offset_q.append(f["parameters"]["q"])
    report = {
        "doi": source["doi"],
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "published_parameters": dict(zip(NAMES, PARAMETERS.tolist())),
        "normalization": "V0 and volumes: A^3 per B2 cell, one NaCl FU, two atoms; molar volume = V*N_A*1e-30 m^3/mol(FU); Debye energy uses n=2; reference 300 K.",
        "native_max_difference_gpa": native_difference,
        "published_observation_comparisons": published,
        "cold_fit_figure3": cold_fit,
        "cold_fit_figure4_sensitivity": cold4_fit,
        "hot_q_only_published_cold": hot_fit,
        "hot_excluding_touching_points_sensitivity": partial_hot_fit,
        "hot_q_only_integrated_law_sensitivity": integrated_fit,
        "staged_cold_then_hot": staged_hot_fit,
        "hot_q_only_volume_residual_sensitivity": volume_hot_fit,
        "figure3_vs_figure4_same_cold_rows": {
            "max_pressure_difference_gpa": float(np.max(abs(cold[2] - cold4[2]))),
            "max_volume_difference_a3": float(np.max(abs(cold[0] - cold4[0]))),
            "duplicate_observations_not_pooled": True,
        },
        "curve_checks": curve_checks,
        "coordinate_sensitivity": {
            "bound_pt": 0.357,
            "pressure_bound_gpa": dp,
            "volume_bound_a3": dv,
            "common_offset_q_range": [min(offset_q), max(offset_q)],
            "meaning": "One printed stroke-width coordinate perturbations; not source errors, statistical interval, or confidence on q.",
        },
        "author_fit_reproduced": False,
        "source_error_weighted_fit": {
            "status": "not_possible",
            "reason": "No P/V/T experimental uncertainties recovered; graphical stroke widths are not source-error weights.",
        },
        "missing_for_exact_reproduction": [
            "original unrounded NaCl cell parameters and temperatures",
            "paired measured Pt cells for independent pressure reduction",
            "P/V/T uncertainty columns and calibration covariance",
            "author weights, residual direction and optimization/covariance details",
        ],
        "qualification": "The source describes a least-squares cold EOS and Mie-Gruneisen thermal fit. q is the remaining thermal parameter once source-adopted gamma0/theta0 are fixed. These diagnostics use only study markers, fixed V0, cold K0/K0_prime free followed by hot q free. Conditional errors omit cold-stage and Pt-calibration covariance. Graphical agreement and q proximity do not establish author-fit parity, theta convention uniqueness or absolute pressure accuracy.",
    }
    return report


def plot(source, report):
    import matplotlib.pyplot as plt

    fig, (ax, residual_ax) = plt.subplots(1, 2, figsize=(10, 4.3), layout="constrained")
    grid = np.linspace(20, 41.35, 160)
    for temperature, color in [(300, "#2463aa"), (1000, "#d74632"), (2000, "#b845ab")]:
        ax.plot(
            pressure(grid, temperature),
            grid,
            color=color,
            label=f"Published equation {temperature} K",
        )
    for temperature, color in [(300, "#2463aa"), (1000, "#d74632")]:
        v, t, p = subset(source, 4, temperature)
        ax.scatter(
            p,
            v,
            s=24,
            facecolors="white" if temperature == 300 else color,
            edgecolors=color,
            zorder=4,
        )
        residual_ax.scatter(
            p, pressure(v, t) - p, color=color, s=24, label=f"Figure 4 {temperature} K"
        )
    for curve in source["calculated_curves"]:
        if curve["figure"] == 4:
            ax.scatter(
                [p["pressure_gpa"] for p in curve["points"]],
                [p["volume_a3"] for p in curve["points"]],
                marker="+",
                s=12,
                color="#777777",
                alpha=0.65,
            )
    ax.set(
        xlabel="Pressure (GPa)",
        ylabel="B2 cell volume (A$^3$)",
        xlim=(0, 135),
        ylim=(19, 42),
    )
    residual_ax.axhline(0, color="grey", lw=0.8)
    residual_ax.set(
        xlabel="Figure pressure (GPa)", ylabel="Equation minus marker pressure (GPa)"
    )
    ax.legend(fontsize=8)
    residual_ax.legend(fontsize=8)
    fig.suptitle("Fei 2007 NaCl-B2: published curves and digitized observations")
    fig.savefig(ROOT / "docs/data/fei-2007-nacl-b2-validation.png", dpi=180)
    plt.close(fig)


def ledger_outcome(record):
    """Keep duplicate cold figures separate and retain the staged hot diagnostic."""
    report = reproduce(json.loads(SOURCE.read_text(encoding="utf-8")))
    cold = report["cold_fit_figure3"]
    hot = report["staged_cold_then_hot"]
    parameters = []
    for name, published, error, fit in [
        ("K0", 26.86, 2.9, cold),
        ("K0_prime", 5.25, 0.26, cold),
        ("q", 0.5, 0.3, hot),
    ]:
        value = fit["parameters"][name]
        parameters.append(
            {
                "parameter": name,
                "published": published,
                "published_error": error,
                "refit": value,
                "refit_error": fit["conditional_standard_errors"][name],
                "relative_difference": abs(value - published) / abs(published),
                "similar": bool(abs(value - published) <= error),
                "within_combined_2sigma": None,
            }
        )
    return {
        "status": "parity_not_achieved",
        "fit_kind": "diagnostic_staged_figure_pvt",
        "dataset_identifiers": record["fit_datasets"],
        "observations": 24,
        "free_parameters": ["K0", "K0_prime", "q"],
        "fixed_parameters": ["V0", "gamma0", "theta0", "n", "Tr"],
        "parameters": parameters,
        "solver_success": cold["solver_success"] and hot["solver_success"],
        "reason": "Twelve Figure 3 cold markers set K0/K0-prime, then twelve Figure 4 hot markers set q with the cold fit frozen. The duplicated Figure 4 cold markers and calculated curves are excluded. The staged hot fit does not recover published q; source errors, measured Pt volumes, original weights and covariance are unavailable. Conditional fit errors do not propagate cold-stage uncertainty. Cold and hot residual summaries remain separate.",
        "reproduction": report,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path)
    args = parser.parse_args()
    source = (
        recover(args.pdf)
        if args.pdf
        else json.loads(SOURCE.read_text(encoding="utf-8"))
    )
    report = reproduce(source)
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    plot(source, report)
    print(
        json.dumps(
            {
                "native_max_difference_gpa": report["native_max_difference_gpa"],
                "cold_fit": report["cold_fit_figure3"]["parameters"],
                "hot_fit": report["hot_q_only_published_cold"]["parameters"],
                "staged_hot_fit": report["staged_cold_then_hot"]["parameters"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
