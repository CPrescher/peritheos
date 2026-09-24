"""Independent BM3 refit of Chen's digitized Brillouin density markers.

This is not reproduction of the source-author parameters. Graphical error
halfwidths define relative coordinate scales; confidence levels and cross-row
covariances are unknown, so no statistical parameter errors are exported.
"""

import argparse
import json
from functools import lru_cache

import numpy as np
import scipy
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3
from scripts.reproduce_chen_2010_argon import MASS_PER_CELL, ROOT, read_points

RECORD = "argon_fcc_chen_2010_bm3_digitized_refit"
DATASET = "argon_fcc_chen_2010_figure5_brillouin"
REPORT = ROOT / "docs/data/chen-2010-argon-refit.json"
LOWER = np.array([0.5, 0.001, 1.0])
UPPER = np.array([1.99, 50.0, 20.0])


def observations():
    keys = (
        "pressure_gpa",
        "density_g_cm3",
        "pressure_plot_halfwidth_gpa",
        "density_plot_halfwidth_g_cm3",
    )
    return np.array(
        [[float(row[k]) for k in keys] for row in read_points("figure5-brillouin")]
    ).T


def pressure(density, pars):
    """Independent BM3 in density coordinates (rho0, K0, Kprime)."""
    rho0, k0, kp = pars
    f = ((np.asarray(density) / rho0) ** (2 / 3) - 1) / 2
    return 3 * k0 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (kp - 4) * f)


def spacing_weights(p, width=None, offset=0.0):
    if width is None:
        return np.ones(len(p))
    _, inverse, counts = np.unique(
        np.floor((p - offset) / width), return_inverse=True, return_counts=True
    )
    # Each occupied bin receives unit total weight; every source row survives.
    return 1 / counts[inverse]


def fit(data, start=(1.65, 4.0, 5.0), mode="eiv", weights=None):
    p, rho, sp, sr = data
    w = np.sqrt(np.ones(len(p)) if weights is None else weights)
    eiv = mode == "eiv"

    def residual(x):
        if eiv:
            return np.r_[w * (pressure(x[3:], x[:3]) - p) / sp, w * (x[3:] - rho) / sr]
        return w * (pressure(rho, x) - p)

    x0 = np.r_[start, rho] if eiv else np.array(start)
    bounds = (
        np.r_[LOWER, np.full(len(p), 1.0)] if eiv else LOWER,
        np.r_[UPPER, np.full(len(p), 6.0)] if eiv else UPPER,
    )
    result = least_squares(
        residual,
        x0,
        bounds=bounds,
        x_scale="jac",
        ftol=1e-11,
        xtol=1e-11,
        gtol=1e-11,
        max_nfev=800,
    )
    if not result.success or np.any(result.active_mask):
        raise RuntimeError(f"Unconverged or bound-active fit: {result.message}")
    pars = result.x[:3]
    raw = pressure(rho, pars) - p
    adjusted = pressure(result.x[3:], pars) - p if eiv else raw
    return {
        "density_parameters": pars.tolist(),
        "parameters": {
            "V0": float(MASS_PER_CELL / pars[0]),
            "K0": float(pars[1]),
            "K0_prime": float(pars[2]),
        },
        "objective": float(np.sum(result.fun**2)),
        "pressure_rms_original_coordinates_gpa": float(np.sqrt(np.mean(raw**2))),
        "pressure_max_abs_original_coordinates_gpa": float(np.max(np.abs(raw))),
        "pressure_rms_adjusted_coordinates_gpa": float(np.sqrt(np.mean(adjusted**2))),
        "adjusted_density_g_cm3": result.x[3:].tolist() if eiv else None,
        "converged": bool(result.success),
        "bound_active": bool(np.any(result.active_mask)),
        "row_weights": (w**2).tolist(),
    }


def volumes_at(fitted, pressures):
    return np.asarray(BM3(**fitted["parameters"]).volume(pressures))


@lru_cache(maxsize=1)
def reproduce():
    data = observations()
    p, rho, sp, sr = data
    starts = [
        fit(data, start=start) for start in [(1.4, 2, 8), (1.65, 4, 5), (1.85, 8, 4)]
    ]
    primary = min(starts, key=lambda x: x["objective"])
    initial = primary["density_parameters"]
    alternatives = {"equal_pressure": fit(data, initial, mode="pressure")}
    for width, offset in [(0.5, 0), (1.0, 0), (2.0, 0), (1.0, 0.5)]:
        name = f"eiv_bin_{width:g}_gpa_offset_{offset:g}"
        alternatives[name] = fit(
            data, initial, weights=spacing_weights(p, width, offset)
        )
    # Perturb all positions together: this tests coordinate-reading offsets,
    # not an assumed distribution of independent random errors.
    perturbations = {}
    for dp, dr in [(0.05, 0), (-0.05, 0), (0, 0.01), (0, -0.01)]:
        shifted = data.copy()
        shifted[0] += dp
        shifted[1] += dr
        perturbations[f"pressure_{dp:+g}_density_{dr:+g}"] = fit(shifted, initial)
    omitted = {}
    bins = np.floor(p)
    for group in np.unique(bins):
        keep = bins != group
        fitted = fit(data[:, keep], initial)
        withheld = pressure(rho[~keep], fitted["density_parameters"]) - p[~keep]
        omitted[f"omit_{group:g}_to_{group + 1:g}_gpa"] = {
            "parameters": fitted["parameters"],
            "withheld_positions": int(sum(~keep)),
            "withheld_pressure_rms_gpa": float(np.sqrt(np.mean(withheld**2))),
        }
    grid = np.linspace(float(min(p)), float(max(p)), 101)
    base = volumes_at(primary, grid)

    def envelope(fits):
        values = np.array([volumes_at(value, grid) for value in fits.values()])
        return {
            "parameter_min": {
                k: min(f["parameters"][k] for f in fits.values())
                for k in primary["parameters"]
            },
            "parameter_max": {
                k: max(f["parameters"][k] for f in fits.values())
                for k in primary["parameters"]
            },
            "max_volume_difference_from_primary_percent": float(
                100 * np.max(np.abs(values / base - 1))
            ),
        }

    native = BM3(**primary["parameters"])
    native_difference = native.pressure(MASS_PER_CELL / rho) - pressure(rho, initial)
    two_gpa_volume = float(native.volume(2))
    curve = read_points("figure5-curve")
    pc = np.array([float(r["pressure_gpa"]) for r in curve])
    rc = np.array([float(r["density_g_cm3"]) for r in curve])
    within = (pc >= min(p)) & (pc <= max(p))
    return {
        "identifier": RECORD,
        "status": "independent_digitized_data_refit",
        "original_publication_reproduction_status": "not_reproduced",
        "dataset_identifier": DATASET,
        "source_rows": len(p),
        "source_pressure_range_gpa": [float(min(p)), float(max(p))],
        "temperature_k": 290,
        "software": {
            "name": "Peritheos audit using scipy.optimize.least_squares",
            "version": scipy.__version__,
        },
        "primary": primary,
        "multistart_max_parameter_difference": float(
            np.max(np.ptp([f["density_parameters"] for f in starts], axis=0))
        ),
        "weighting_alternatives": alternatives,
        "weighting_sensitivity": envelope(alternatives),
        "common_coordinate_shifts": perturbations,
        "coordinate_shift_sensitivity": envelope(perturbations),
        "leave_one_pressure_bin_out": omitted,
        "deletion_sensitivity": envelope(omitted),
        "native_pressure_max_difference_gpa": float(np.max(np.abs(native_difference))),
        "refit_at_2_gpa": {
            "density_g_cm3": MASS_PER_CELL / two_gpa_volume,
            "K_T_gpa": float(native.bulk_modulus(two_gpa_volume)),
            "K_T_prime": float(native.bulk_modulus_derivative(two_gpa_volume)),
        },
        "published_curve_comparison": {
            "used_in_fit": False,
            "within_data_range_pressure_rms_gpa": float(
                np.sqrt(np.mean((pressure(rc[within], initial) - pc[within]) ** 2))
            ),
        },
        "uncertainty_policy": "No statistical parameter errors or covariance claimed. Graphical halfwidths only define relative error-coordinate scales; cross-row covariance and confidence levels are unknown. Sensitivity envelopes are not confidence intervals.",
        "selection_policy": "All 80 distinct marker positions, each once, ignoring PDF rendering multiplicity. No published line vertices, extrapolated density or reported elastic constants enter the fit.",
    }


def plot(report, target):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p, rho, sp, sr = observations()
    primary = report["primary"]
    grid = np.linspace(min(p), max(p), 250)
    fig, (ax, residual) = plt.subplots(
        2, 1, figsize=(8.5, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.errorbar(
        p,
        rho,
        xerr=sp,
        yerr=sr,
        fmt="o",
        markersize=3,
        color="#52758a",
        ecolor="#b8c8d0",
        elinewidth=0.7,
        label="Brillouin-derived density (digitized)",
    )
    ax.plot(
        grid,
        MASS_PER_CELL / volumes_at(primary, grid),
        color="#b54635",
        linewidth=2,
        label="Independent BM3 refit",
    )
    curve = read_points("figure5-curve")
    ax.plot(
        [float(r["pressure_gpa"]) for r in curve],
        [float(r["density_g_cm3"]) for r in curve],
        "--",
        color="#474747",
        label="Published fit line (comparison only)",
    )
    variants = np.array(
        [
            MASS_PER_CELL / volumes_at(f, grid)
            for f in report["weighting_alternatives"].values()
        ]
    )
    ax.fill_between(
        grid,
        variants.min(axis=0),
        variants.max(axis=0),
        color="#e7a79d",
        alpha=0.5,
        label="Weighting sensitivity; not a confidence interval",
    )
    ax.set_ylabel("Density (g/cm³)")
    ax.set_title(
        "Chen 2010 argon: independent refit of digitized data", loc="left", pad=14
    )
    ax.legend(fontsize=8, loc="lower right")
    residual.axhline(0, color="#777777", linewidth=0.8)
    residual.scatter(
        p, pressure(rho, primary["density_parameters"]) - p, s=12, color="#b54635"
    )
    residual.set_ylabel("ΔP (GPa)")
    residual.set_xlabel("Reported pressure (GPa)")
    residual.text(
        0.02,
        0.92,
        "Residuals at original coordinates",
        transform=residual.transAxes,
        va="top",
        fontsize=8,
    )
    for panel in (ax, residual):
        panel.spines[["top", "right"]].set_visible(False)
        panel.grid(alpha=0.15)
    ax.set_xlim(0, 27)
    fig.tight_layout()
    fig.savefig(target, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    if args.plot:
        plot(report, ROOT / "docs/data/chen-2010-argon-refit.png")
    print(REPORT)
