"""Audit the two published Frost Pd EOS against all 93 supplementary rows.

No fitted coefficients from this diagnostic are installed in the catalog.
The source describes uncertainty-aware Marquardt-Levenberg fitting, but does
not define its objective or covariance convention. Exact orthogonal distance
and first-order effective variance are explicit interpretations, not recovered
author code. An unweighted pressure baseline is reported separately.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re

# Retain the archived ODRPACK cross-check without changing its numerical method.
# SciPy 1.17+ warns on import; migration is required before SciPy 1.19 removes it.
import warnings
from pathlib import Path

import numpy as np

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore", message=r"`scipy\.odr` is deprecated.*", category=DeprecationWarning
    )
    from scipy import odr
from scipy.optimize import brentq, least_squares

from peritheos.eos.rt import BM3, Vinet
from peritheos.fitting import fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT / "peritheos/data/datasets/palladium-frost-2023-tables-i-iii-compression.csv"
)
SOURCE_ROWS = ROOT / "docs/data/frost-2023-palladium-source-rows.json"
OUTPUT = ROOT / "docs/data/frost-2023-palladium-audit.json"
DATA_HASH = "92e4ece39614b1a9f980bff32ae3986083b06bd9673ef0e615f4d27bd39fa75c"
NAMES = ("V0", "K0", "K0_prime")
PUBLISHED = {
    "Vinet": [58.678, 189.3, 5.473],
    "BM3": [58.588, 197.5, 4.996],
}
ERRORS = {"Vinet": [0.073, 3.0, 0.063], "BM3": [0.076, 3.3, 0.066]}
MODELS = {"Vinet": Vinet, "BM3": BM3}


def pressure(kind, parameters, volume):
    """Independent literal equations in GPa and conventional-cell A^3."""
    v0, k0, kp = parameters
    x = (np.asarray(volume) / v0) ** (1.0 / 3.0)
    if kind == "Vinet":
        return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))
    eta = 1 / x
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def slope(kind, parameters, volume):
    """Central difference only for the effective-variance diagnostic."""
    return (
        pressure(kind, parameters, volume + 1e-4)
        - pressure(kind, parameters, volume - 1e-4)
    ) / 2e-4


def verify_source(text_path):
    text = text_path.read_text(encoding="utf-8")
    # Rh Table IV is excluded. Retain every printed Pd row and original token.
    text = text.split("B.   Rhodium in Neon")[0]
    rows = [
        dict(zip(("pressure", "a", "delta_a", "loading", "marker"), match))
        for match in re.findall(
            r"^\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(BX18|BX32|HP)\s+(W|Au)\s*$",
            text,
            re.MULTILINE,
        )
    ]
    assert len(rows) == 93
    evidence = {
        "source": "Official supplementary PDF, Section II.A, Tables I-III, pp. 5-7",
        "url": "https://pubs.aip.org/jap/article-supplement/2903099/zip/035901_1_5.0161038.suppl_material/",
        "archive_sha256": "7f03394af010c3eaa6b4346ac6ed0e2273ff86f28ba1cf910471d564873ae015",
        "pdf_sha256": "7e47e13580ecc7faec0a372347217e860aee65c6e32cc279a2c6e66651653f3b",
        "rows": rows,
        "ambient_separate_from_tables": {"a": "3.89088", "delta_a": "0.00002"},
    }
    SOURCE_ROWS.write_text(json.dumps(evidence, indent=2) + "\n")


def observations():
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == DATA_HASH
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    source = json.loads(SOURCE_ROWS.read_text(encoding="utf-8"))["rows"]
    assert len(rows) == len(source) == 93
    for row, original in zip(rows, source):
        for stored, printed in (
            ("pressure_gpa", "pressure"),
            ("lattice_parameter_a", "a"),
            ("lattice_parameter_error_a", "delta_a"),
        ):
            assert float(row[stored]) == float(original[printed])
        assert row["loading"] == original["loading"]
        assert row["pressure_marker"] == original["marker"]

    def values(name):
        return np.array([float(row[name]) for row in rows])

    p = values("pressure_gpa")
    a = values("lattice_parameter_a")
    da = values("lattice_parameter_error_a")
    v, sv = a**3, 3 * a**2 * da
    np.testing.assert_allclose(
        values("volume_a3_conventional_cell"), v, atol=1e-9, rtol=0
    )
    np.testing.assert_allclose(
        values("volume_error_a3_conventional_cell"), sv, atol=1e-9, rtol=0
    )
    np.testing.assert_array_equal(values("pressure_error_gpa"), np.full(93, 0.1))
    return p, v, np.full(93, 0.1), sv


def replay():
    p, v, sp, sv = observations()
    report = {
        "date": "2026-10-17",
        "dataset_sha256": DATA_HASH,
        "observations": len(p),
        "pressure_range_gpa": [float(min(p)), float(max(p))],
        "selection": "All 93 Tables I-III rows; no exclusions; ambient point separate",
        "source_protocol": "Uncertainty-aware nonlinear Marquardt-Levenberg; exact objective, constraints, covariance convention and software not supplied",
        "original_author_workflow_reproduced": False,
        "parameter_error_confidence": "Not stated in the source; two-error-width comparison is conditional",
        "covariance_scope": "Conditional iid coordinate errors; excludes run correlations, calibration systematics and EOS inadequacy",
        "models": {},
    }
    for kind, published in PUBLISHED.items():
        model_type = MODELS[kind]
        native = fit_rt_eos(
            model_type,
            volume=v,
            pressure=p,
            initial=dict(zip(NAMES, published)),
            pressure_sigma=sp,
            volume_sigma=sv,
            absolute_sigma=True,
            max_nfev=5000,
        )
        assert native.success
        fitted = np.array([native.parameters[name] for name in NAMES])
        cross_check = odr.ODR(
            odr.RealData(v, p, sx=sv, sy=sp),
            odr.Model(lambda b, x: pressure(kind, b, x)),
            beta0=published,
            maxit=1000,
            sstol=1e-12,
            partol=1e-12,
        ).run()
        assert cross_check.info in (1, 2, 3)
        np.testing.assert_allclose(cross_check.beta, fitted, rtol=2e-6, atol=2e-5)
        parameter_scale = np.asarray(published)

        def orthogonal_lm_residuals(z):
            return np.r_[
                (pressure(kind, z[:3] * parameter_scale, v + z[3:] * sv) - p) / sp,
                z[3:],
            ]

        literal_lm = least_squares(
            orthogonal_lm_residuals,
            np.r_[np.ones(3), np.zeros(len(p))],
            method="lm",
            max_nfev=30000,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        assert literal_lm.success
        lm_parameters = literal_lm.x[:3] * parameter_scale
        np.testing.assert_allclose(lm_parameters, fitted, rtol=2e-6, atol=2e-5)
        np.testing.assert_allclose(
            pressure(kind, published, v),
            model_type(**dict(zip(NAMES, published))).pressure(v),
            atol=1e-10,
            rtol=0,
        )
        baseline = least_squares(
            lambda b: pressure(kind, b, v) - p,
            published,
            method="lm",
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
        )
        effective = least_squares(
            lambda b: (
                (pressure(kind, b, v) - p)
                / np.sqrt(sp**2 + (slope(kind, b, v) * sv) ** 2)
            ),
            published,
            method="lm",
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        # Evaluate the exact orthogonal objective with source coefficients held
        # fixed, optimizing only latent observation volumes, for a fair comparison.
        fixed = least_squares(
            lambda delta: np.r_[
                (pressure(kind, published, v + delta * sv) - p) / sp, delta
            ],
            np.zeros(len(p)),
            method="lm",
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        published_residual = pressure(kind, published, v) - p
        fitted_residual = pressure(kind, fitted, v) - p
        combined = np.hypot(ERRORS[kind], [native.standard_errors[n] for n in NAMES])
        scaled_errors = np.array([native.standard_errors[n] for n in NAMES]) * np.sqrt(
            native.reduced_chi_square
        )
        scaled_combined = np.hypot(ERRORS[kind], scaled_errors)
        # Source table rounds P to .01 GPa and a to .0001 A. This is a
        # first-order worst-coordinate pressure perturbation bound, not a fit.
        rounding_bound = 0.005 + abs(slope(kind, published, v)) * (
            3 * a_square(v) * 0.00005
        )
        report["models"][kind] = {
            "published": published,
            "published_error_widths": ERRORS[kind],
            "unweighted_pressure_baseline": baseline.x.tolist(),
            "unweighted_pressure_rms_gpa": float(np.sqrt(np.mean(baseline.fun**2))),
            "effective_variance_interpretation": effective.x.tolist(),
            "orthogonal_refit": fitted.tolist(),
            "scipy_odr_independent": cross_check.beta.tolist(),
            "scipy_odr_objective": float(cross_check.sum_square),
            "literal_marquardt_levenberg_same_orthogonal_objective": lm_parameters.tolist(),
            "literal_marquardt_levenberg_objective": float(sum(literal_lm.fun**2)),
            "orthogonal_objective": float(native.chi_square),
            "published_fixed_orthogonal_objective": float(sum(fixed.fun**2)),
            "degrees_of_freedom_refit": native.degrees_of_freedom,
            "reduced_chi_square_refit": float(native.reduced_chi_square),
            "published_pressure_rms_gpa": float(
                np.sqrt(np.mean(published_residual**2))
            ),
            "refit_pressure_rms_gpa": float(np.sqrt(np.mean(fitted_residual**2))),
            "published_pressure_max_abs_gpa": float(max(abs(published_residual))),
            "formal_errors_absolute_coordinate_sigma": [
                native.standard_errors[n] for n in NAMES
            ],
            "conditional_errors_scaled_by_residual_variance": scaled_errors.tolist(),
            "covariance_absolute_coordinate_sigma": native.covariance.tolist(),
            "correlation": native.correlation.tolist(),
            "shift_in_combined_error_widths": (
                (fitted - published) / combined
            ).tolist(),
            "within_combined_two_error_widths": (
                abs(fitted - published) <= 2 * combined
            ).tolist(),
            "shift_in_combined_error_widths_residual_scaled": (
                (fitted - published) / scaled_combined
            ).tolist(),
            "within_combined_two_error_widths_residual_scaled": (
                abs(fitted - published) <= 2 * scaled_combined
            ).tolist(),
            "max_first_order_source_coordinate_rounding_bound_gpa": float(
                max(rounding_bound)
            ),
            "row_diagnostics": [
                {
                    "source_row": i + 1,
                    "pressure_gpa": float(p[i]),
                    "observed_volume_cell_a3": float(v[i]),
                    "volume_error_cell_a3": float(sv[i]),
                    "published_pressure_residual_at_measured_volume_gpa": float(
                        published_residual[i]
                    ),
                    "refit_pressure_residual_at_measured_volume_gpa": float(
                        fitted_residual[i]
                    ),
                    "refit_adjusted_volume_cell_a3": float(native.adjusted_volume[i]),
                    "refit_pressure_residual_at_adjusted_volume_gpa": float(
                        native.residuals[i]
                    ),
                }
                for i in range(len(p))
            ],
        }
    return report


def a_square(volume):
    return volume ** (2 / 3)


def plot_report(report, destination):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p, v, sp, sv = observations()
    with DATA.open() as stream:
        loadings = np.array([r["loading"] for r in csv.DictReader(stream)])
    grid = np.linspace(0, 185, 700)
    fig, axes = plt.subplots(
        2, 2, figsize=(12, 8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
    )

    def volumes(kind, b, pressures):
        return np.array(
            [
                brentq(lambda v: pressure(kind, b, v) - q, 0.5 * b[0], b[0] * 1.01)
                for q in pressures
            ]
        )

    for column, kind in enumerate(PUBLISHED):
        top, bottom = axes[:, column]
        row = report["models"][kind]
        source_v = volumes(kind, row["published"], p)
        fit_v = volumes(kind, row["orthogonal_refit"], p)
        top.plot(
            grid,
            volumes(kind, row["published"], grid) / 4,
            color="black",
            label="Published coefficients",
            lw=1.8,
        )
        top.plot(
            grid,
            volumes(kind, row["orthogonal_refit"], grid) / 4,
            color="#b52666",
            label="Orthogonal replay (diagnostic)",
            lw=1.5,
            ls="--",
        )
        for loading, color in zip(
            ("BX18", "BX32", "HP"), ("#2269ad", "#d77b23", "#299674")
        ):
            mask = loadings == loading
            top.errorbar(
                p[mask],
                v[mask] / 4,
                xerr=sp[mask],
                yerr=sv[mask] / 4,
                fmt="o",
                ms=3.5,
                color=color,
                elinewidth=0.6,
                label=f"{loading} ({sum(mask)} observations)",
            )
            bottom.errorbar(
                p[mask],
                (v[mask] - source_v[mask]) / 4,
                yerr=np.sqrt(
                    sv[mask] ** 2
                    + (sp[mask] / slope(kind, row["published"], source_v[mask])) ** 2
                )
                / 4,
                fmt="o",
                ms=3.5,
                color=color,
                elinewidth=0.6,
            )
            bottom.plot(
                p[mask],
                (v[mask] - fit_v[mask]) / 4,
                linestyle="none",
                marker="x",
                ms=4,
                color=color,
                alpha=0.8,
            )
        top.plot(
            0,
            3.89088**3 / 4,
            marker="D",
            color="#734b19",
            ms=5,
            linestyle="none",
            label="Separate ambient measurement",
        )
        top.set_title(kind)
        top.set_ylabel("Volume (Å³/atom)")
        top.legend(fontsize=8, loc="upper right")
        bottom.axhline(0, color="black", lw=0.8)
        bottom.set_xlabel("Pressure (GPa)")
        bottom.set_ylabel("Observed − curve\nvolume (Å³/atom)")
        bottom.set_ylim(-0.085, 0.085)
        for ax in (top, bottom):
            ax.grid(alpha=0.16)
            ax.set_xlim(-3, 185)
    fig.suptitle(
        "Frost et al. (2023): all 93 Pd observations and uncertainty-aware replay",
        fontsize=14,
    )
    fig.text(
        0.5,
        0.015,
        "Residuals: circles = published; crosses = orthogonal replay. Bars: coordinate errors only; no calibration/run covariance.\nBoth refits are diagnostics; neither is installed as a replacement EOS.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 0.96))
    destination.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination / "frost-pd-audit.png", dpi=180)
    fig.savefig(destination / "frost-pd-audit.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-text", type=Path)
    parser.add_argument("--plot-dir", type=Path)
    args = parser.parse_args()
    if args.source_text:
        verify_source(args.source_text)
    report = replay()
    OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if args.plot_dir:
        plot_report(report, args.plot_dir)
    for kind, row in report["models"].items():
        print(kind, json.dumps(row))


if __name__ == "__main__":
    main()
