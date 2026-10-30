"""Independent, repository-only Prescher (2012) source recovery.

Default replay needs NumPy/SciPy only; it never imports Peritheos fitting code.
PDF extraction additionally needs pdfplumber (available in the document runtime).
Digitized markers, drawn model curves, and published table rows stay distinct.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "peritheos/data/datasets/cementite-prescher-2012-table-s1.csv"
VECTORS = ROOT / "peritheos/data/research/cementite/prescher-2012-figure1-vectors.json"
REPORT = ROOT / "docs/data/prescher-2012-fe3c-recovery.json"
SOURCE = {
    "full_range": ([161, 5.9], [2, 0.2]),
    "below_22_gpa": ([145, 8.5], [3, 0.7]),
    "above_22_gpa": ([172, 5.1], [1, 0.1]),
}


def bm3(volume, v0, k0, kp):
    """GPa and conventional-cell Å³; V0 is the zero-pressure reference."""
    eta = (v0 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def volume_at_pressure(pressure, v0, k0, kp):
    from scipy.optimize import brentq

    return np.array(
        [
            brentq(lambda v: bm3(v, v0, k0, kp) - p, 0.65 * v0, 1.03 * v0)
            for p in np.atleast_1d(pressure)
        ]
    )


def selection(name, pressure):
    if name == "below_22_gpa":
        return pressure < 22
    if name == "above_22_gpa":
        return pressure >= 22
    return np.ones(len(pressure), dtype=bool)


def extract_pdf(path):
    """Recover vector geometry with axes calibrated independently of Table S1."""
    import pdfplumber

    if hashlib.sha256(path.read_bytes()).hexdigest() != (
        "a908c16fdd1754a447d3959a6bf4d2cb1b46b21dbf08888981c6ab400ab77c12"
    ):
        # Geometry is version-specific; inspect a different PDF before using it.
        raise ValueError("Unexpected source PDF; inspect its geometry first")
    with pdfplumber.open(path) as pdf:
        page = pdf.pages[1]
        x_ticks = np.array(
            [
                80.660,
                101.670,
                122.679,
                143.780,
                164.789,
                185.799,
                206.809,
                227.818,
                248.920,
                269.929,
                290.938,
            ]
        )
        y_ticks = np.array(
            [
                646.978,
                627.522,
                607.975,
                588.517,
                568.969,
                549.512,
                529.964,
                510.507,
                490.959,
            ]
        )
        # Assert the inspected tick positions still occur in the actual PDF.
        for x in x_ticks[:-1]:
            assert any(
                abs(line["x0"] - x) < 0.002
                and abs(line["bottom"] - 646.978) < 0.002
                and 3 < line["height"] < 3.2
                for line in page.lines
            )
        for y in y_ticks:
            assert any(
                abs(line["top"] - y) < 0.002
                and abs(line["x0"] - 78.559) < 0.002
                and line["width"] > 3
                for line in page.lines
            )
        ax = np.polyfit(x_ticks, np.arange(0, 51, 5), 1)
        ay = np.polyfit(y_ticks, np.arange(120, 161, 5), 1)
        markers = []
        for curve in page.curves:
            if (
                77 < curve["x0"] < 295
                and 490 < curve["top"] < 648
                and curve["fill"]
                and 3.7 < curve["width"] < 4
                and 3.7 < curve["height"] < 4
            ):
                x = (curve["x0"] + curve["x1"]) / 2
                y = (curve["top"] + curve["bottom"]) / 2
                markers.append(
                    {
                        "pdf_x": x,
                        "pdf_top": y,
                        "pressure_gpa": float(np.polyval(ax, x)),
                        "volume_a3": float(np.polyval(ay, y)),
                    }
                )
        markers.sort(key=lambda row: row["pressure_gpa"])
        assert len(markers) == 22
        table = np.genfromtxt(TABLE, delimiter=",", names=True)
        for i, row in enumerate(markers, 1):
            row["marker"] = i
            nearest = int(np.argmin(abs(table["pressure_gpa"] - row["pressure_gpa"])))
            matched = abs(table["pressure_gpa"][nearest] - row["pressure_gpa"]) < 0.12
            row["matching_table_s1_row"] = nearest + 1 if matched else None
        curves = [
            c
            for c in page.curves
            if 77 < c["x0"] < 82
            and c["x1"] > 291
            and c["top"] > 490
            and c["bottom"] < 648
            and c["stroke"]
            and c["width"] > 200
        ]
        assert len(curves) == 1
        return {
            "format": "peritheos.prescher-2012-figure1-vectors",
            "doi": "10.1103/PhysRevB.85.140402",
            "page": 2,
            "figure": 1,
            "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "kind": "digitized_vector_graphics_not_original_numeric_observations",
            "reader": f"pdfplumber {pdfplumber.__version__}",
            "coordinates": "PDF points, x from left and top from page top",
            "pressure_ticks_pdf_x": x_ticks.tolist(),
            "pressure_ticks_gpa": list(range(0, 51, 5)),
            "volume_ticks_pdf_top": y_ticks.tolist(),
            "volume_ticks_a3": list(range(120, 161, 5)),
            "pressure_affine": ax.tolist(),
            "volume_affine": ay.tolist(),
            "precision_limit": "The drawing has approximately 0.09-point coordinate steps (~0.02 GPa, ~0.02 A3). Original numeric inputs and experimental errors cannot be recovered to this precision. Axis and graphics errors are correlated, not independent experimental sigmas.",
            "markers": markers,
            "drawn_fit_curve": [
                {
                    "pdf_x": x,
                    "pdf_top": y,
                    "pressure_gpa": float(np.polyval(ax, x)),
                    "volume_a3": float(np.polyval(ay, y)),
                }
                for x, y in curves[0]["pts"]
            ],
            "curve_qualification": "Calculated source curve, not measured observations. Regression errors on densely drawn vertices are not experimental parameter uncertainties.",
        }


def fit(p, v, mode="equal_pressure", fixed_v0=None, sp=None, sv=None, initial=None):
    """Independent nonlinear LS; widths are conditional numerical weights."""
    from scipy.optimize import curve_fit, least_squares
    from scipy.stats import t

    initial = [155, 160, 6] if initial is None else initial
    free = fixed_v0 is None
    start = initial if free else initial[1:]

    def unpack(b):
        return b[:3] if free else np.r_[fixed_v0, b[:2]]

    npar = 3 if free else 2
    if mode == "coordinate_widths":

        def residual(b):
            return np.r_[(bm3(b[npar:], *unpack(b)) - p) / sp, (b[npar:] - v) / sv]

        result = least_squares(
            residual,
            np.r_[start, v],
            x_scale="jac",
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            max_nfev=5000,
        )
        b = unpack(result.x)
        covariance = np.linalg.inv(result.jac.T @ result.jac)[:npar, :npar]
        ss = float(result.fun @ result.fun)
        success = bool(result.success)
    else:
        inverse = mode in ("equal_volume", "volume_widths")
        fn = volume_at_pressure if inverse else bm3
        x, y = (p, v) if inverse else (v, p)
        sigma = (
            sp if mode == "pressure_widths" else sv if mode == "volume_widths" else None
        )
        params, covariance = curve_fit(
            lambda x, *b: fn(x, *unpack(b)),
            x,
            y,
            p0=start,
            sigma=sigma,
            absolute_sigma=True,
            maxfev=10000,
            ftol=1e-12,
            xtol=1e-12,
        )
        b = unpack(params)
        residuals = fn(x, *b) - y
        ss = float(np.sum((residuals / (1 if sigma is None else sigma)) ** 2))
        success = True
    dof = len(p) - npar
    scaled = covariance * ss / dof
    errors = np.sqrt(np.diag(scaled))
    corr = scaled / np.outer(errors, errors)
    rp = bm3(v, *b) - p
    rv = volume_at_pressure(p, *b) - v
    return {
        "parameters": dict(zip(("V0", "K0", "K0_prime"), map(float, b))),
        "fixed_parameters": {} if free else {"V0": float(fixed_v0)},
        "free_parameters": ["V0", "K0", "K0_prime"] if free else ["K0", "K0_prime"],
        "residual_scaled_errors": errors.tolist(),
        "residual_scaled_covariance": scaled.tolist(),
        "unscaled_inverse_normal_matrix": covariance.tolist(),
        "conditional_95pct_half_widths": (errors * t.ppf(0.975, dof)).tolist(),
        "correlation": corr.tolist(),
        "objective": mode,
        "observations": len(p),
        "degrees_of_freedom": dof,
        "sum_squared_objective_residuals": ss,
        "pressure_rmse_gpa": float(np.sqrt(np.mean(rp**2))),
        "volume_rmse_a3": float(np.sqrt(np.mean(rv**2))),
        "pressure_residuals_gpa": rp.tolist(),
        "volume_residuals_a3": rv.tolist(),
        "success": success,
    }


def polynomial_pressure_fit(p, v):
    """Unconstrained linear replay: BM3 spans V^-3, V^-7/3, V^-5/3.

    This checks the pressure-LS optimum without nonlinear initial values.
    The root near the measured ambient volume identifies the physical branch.
    """
    q = (155 / v) ** (2 / 3)
    design = np.column_stack([q**4.5, q**3.5, q**2.5])
    a, b, c = np.linalg.lstsq(design, p, rcond=None)[0]
    roots = np.roots([a, b, c])
    root = min(
        (r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 0),
        key=lambda r: abs(r - 1),
    )
    v0 = 155 / root**1.5
    first = root**2.5 * (2 * a * root + b)
    second = 5 * root**1.5 * (2 * a * root + b) + 2 * a * root**2.5
    k0 = 2 * root * first / 3
    kp = 2 / 3 + 2 * root * second / (3 * first)
    return [float(v0), float(k0), float(kp)]


def replay():
    import scipy

    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    vectors = json.loads(VECTORS.read_text())
    points = vectors["markers"]
    fp = np.array([r["pressure_gpa"] for r in points])
    fv = np.array([r["volume_a3"] for r in points])
    cp = np.array([r["pressure_gpa"] for r in vectors["drawn_fit_curve"]])
    cv = np.array([r["volume_a3"] for r in vectors["drawn_fit_curve"]])
    curve_fit_result = fit(cp, cv)
    curve_v0 = curve_fit_result["parameters"]["V0"]
    results = {}
    for name, (source, widths) in SOURCE.items():
        fits = {}
        m = selection(name, table["pressure_gpa"])
        p, v = table["pressure_gpa"][m], table["volume_a3"][m]
        for mode in (
            "equal_pressure",
            "equal_volume",
            "pressure_widths",
            "volume_widths",
            "coordinate_widths",
        ):
            for fixed in (False, True):
                key = "table_" + mode + ("_curve_V0_fixed" if fixed else "_free_V0")
                fits[key] = fit(
                    p,
                    v,
                    mode,
                    curve_v0 if fixed else None,
                    table["pressure_error_gpa"][m],
                    table["volume_error_a3"][m],
                )
        m = selection(name, fp)
        for mode in ("equal_pressure", "equal_volume"):
            for fixed in (False, True):
                key = "figure_" + mode + ("_curve_V0_fixed" if fixed else "_free_V0")
                fits[key] = fit(fp[m], fv[m], mode, curve_v0 if fixed else None)
        for value in fits.values():
            pair = np.array([value["parameters"][k] for k in ("K0", "K0_prime")])
            value["difference_from_published_pair"] = (pair - source).tolist()
            value["both_within_printed_widths"] = bool(
                np.all(abs(pair - source) <= widths)
            )
        results[name] = {
            "published_pair": source,
            "published_widths": widths,
            "fits": fits,
            "table_equal_pressure_polynomial_crosscheck": polynomial_pressure_fit(p, v),
            "figure_equal_pressure_polynomial_crosscheck": polynomial_pressure_fit(
                fp[m], fv[m]
            ),
            "table_pressure_start_sensitivity": [
                fit(p, v, initial=initial)["parameters"]
                for initial in ([154, 100, 10], [155, 160, 6], [156, 250, 4])
            ],
        }
    # Isolate row-set effects without changing any published table row.
    unmatched = [r for r in points if r["matching_table_s1_row"] is None]
    ablations = {}
    for label, extra in (
        ("table_plus_figure_ambient", unmatched[:1]),
        ("table_plus_three_figure_only_markers", unmatched),
    ):
        p = np.r_[table["pressure_gpa"], [r["pressure_gpa"] for r in extra]]
        v = np.r_[table["volume_a3"], [r["volume_a3"] for r in extra]]
        ablations[label] = fit(p, v)
    products = table["a_angstrom"] * table["b_angstrom"] * table["c_angstrom"]
    source_curve_data_residuals = {}
    curve_parameters = list(curve_fit_result["parameters"].values())
    for label, p, v in (
        ("table_s1_19_rows", table["pressure_gpa"], table["volume_a3"]),
        ("figure1_22_digitized_markers", fp, fv),
    ):
        residuals = bm3(v, *curve_parameters) - p
        source_curve_data_residuals[label] = {
            "pressure_residuals_gpa": residuals.tolist(),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
            "maximum_absolute_pressure_residual_gpa": float(max(abs(residuals))),
            "qualification": "Recovered drawn model evaluated against observations/markers; distinct from regression residuals on the calculated source-curve vertices.",
        }
    return {
        "format": "peritheos.prescher-2012-independent-recovery",
        "date": "2026-10-30",
        "doi": "10.1103/PhysRevB.85.140402",
        "scipy_version": scipy.__version__,
        "numpy_version": np.__version__,
        "qualification": "Source recovery and conditional diagnostics. No catalog coefficients or observations replaced. Graph-derived V0 and weights are not recovered author constraints. Dense source-curve vertices are calculated, not experimental inputs. Printed source widths have unspecified confidence meaning.",
        "table_sha256": hashlib.sha256(TABLE.read_bytes()).hexdigest(),
        "vectors_sha256": hashlib.sha256(VECTORS.read_bytes()).hexdigest(),
        "table_lattice_product_minus_printed_volume_a3": (
            products - table["volume_a3"]
        ).tolist(),
        "figure_only_markers": unmatched,
        "source_drawn_curve_reconstruction": curve_fit_result,
        "source_drawn_curve_volume_objective": fit(cp, cv, "equal_volume"),
        "source_drawn_curve_against_data": source_curve_data_residuals,
        "row_set_diagnostics": ablations,
        "subsets": results,
        "author_workflow_reproduced": False,
    }


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    p = np.linspace(0, 49, 250)
    drawn = report["source_drawn_curve_reconstruction"]["parameters"]
    axes[0].plot(
        p,
        volume_at_pressure(p, *drawn.values()),
        color="black",
        label="BM3 recovered from published drawn curve",
    )
    axes[0].plot(
        table["pressure_gpa"],
        table["volume_a3"],
        "o",
        ms=4,
        label="Table S1: 19 states",
    )
    extra = report["figure_only_markers"]
    axes[0].plot(
        [r["pressure_gpa"] for r in extra],
        [r["volume_a3"] for r in extra],
        "s",
        mfc="none",
        ms=7,
        label="Figure-only markers (digitized)",
    )
    allfit = report["subsets"]["full_range"]["fits"]
    for key, label, color in (
        ("table_equal_pressure_free_V0", "19 table rows, free V0", "#D55E00"),
        ("figure_equal_pressure_free_V0", "22 figure markers, free V0", "#009E73"),
    ):
        par = allfit[key]["parameters"]
        delta = volume_at_pressure(p, *par.values()) - volume_at_pressure(
            p, *drawn.values()
        )
        axes[1].plot(p, delta, color=color, label=label)
    axes[0].set(xlabel="Pressure (GPa)", ylabel="Volume (Å³ / conventional cell)")
    axes[1].set(
        xlabel="Pressure (GPa)", ylabel="V(replay) − V(drawn source curve) (Å³)"
    )
    axes[1].axhline(0, color="0.6", lw=0.7)
    for ax in axes:
        ax.legend(fontsize=8)
    fig.suptitle("Prescher (2012): Figure 1 and Table S1 contain different row sets")
    dest = (
        ROOT / "docs/literature-reproductions/figures/prescher-2012-fe3c-recovery.png"
    )
    fig.savefig(dest, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract-pdf", type=Path)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.extract_pdf:
        VECTORS.write_text(json.dumps(extract_pdf(args.extract_pdf), indent=2) + "\n")
        print(VECTORS)
    else:
        report = replay()
        REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        if args.plot:
            plot(report)
        print(REPORT)
