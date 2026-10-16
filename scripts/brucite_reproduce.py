"""Recover brucite tables and audit source equations without catalog mutation.

Run ``uv run python scripts/brucite_reproduce.py --write`` to regenerate the
three CSVs and report; ``--check`` verifies the committed results. Optional
``--source-check`` checks the locally downloaded final Fukui PDF against the
transcription. No network access or author-workflow parity is implied.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import re
import subprocess
from decimal import Decimal
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
SOURCE = DATA / "brucite-source.json"
REPORT = ROOT / "docs/data/brucite-reproduction.json"
STAGED = ROOT / "output/material-candidates/brucite/brucite.eosmat"
TABLES = (
    ("fukui_2003", "table1", "brucite-fukui-2003-table1.csv"),
    ("fukui_2003", "table4", "brucite-fukui-2003-table4.csv"),
    ("xia_1998", "table6", "brucite-xia-1998-table6.csv"),
)
COLUMNS = (
    "source_row",
    "temperature_k",
    "temperature_error_k",
    "pressure_gpa",
    "pressure_error_gpa",
    "lattice_a_angstrom",
    "lattice_a_error_angstrom",
    "lattice_c_angstrom",
    "lattice_c_error_angstrom",
    "c_over_a",
    "cell_volume_a3",
    "cell_volume_error_a3",
)


def source():
    return json.loads(SOURCE.read_text(encoding="utf-8"))


def parse_number(token):
    """Preserve decimal precision and unpack last-digit parentheses exactly."""
    match = re.fullmatch(r"(-?\d+(?:\.\d+)?)(?:\((\d+)\))?", token)
    if match is None:
        raise ValueError(token)
    value, digits = match.groups()
    error = ""
    if digits is not None:
        error = format(
            Decimal(digits) * Decimal(10) ** Decimal(value).as_tuple().exponent,
            "f",
        )
    return value, error


def normalized_csv(source_id, table_id):
    table = source()["sources"][source_id][table_id]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    errors = {
        "pressure_gpa": "pressure_error_gpa",
        "lattice_a_angstrom": "lattice_a_error_angstrom",
        "lattice_c_angstrom": "lattice_c_error_angstrom",
        "cell_volume_a3": "cell_volume_error_a3",
    }
    for index, tokens in enumerate(table["rows"], 1):
        if len(tokens) != len(table["columns"]):
            raise ValueError(f"{source_id}/{table_id} row {index}: wrong width")
        row = dict.fromkeys(COLUMNS, "")
        row["source_row"] = str(index)
        for column, token in zip(table["columns"], tokens):
            row[column], error = parse_number(token)
            if column in errors:
                row[errors[column]] = error
        writer.writerow(row)
    return stream.getvalue()


def observations(source_id, table_id):
    return list(csv.DictReader(normalized_csv(source_id, table_id).splitlines()))


def bm3(volume, v0, k0, kp):
    """Direct transcription of Fukui Eq. 2, independent of Peritheos."""
    eta = (v0 / np.asarray(volume)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def comparison_volume(pressure, v0, k0, kp):
    """Invert the brucite comparison curves on their near/compressed branch.

    The modest tensile bracket includes Table 4's -0.06 GPa observation.
    This helper is for these source parameters, not a general EOS solver.
    """
    return np.array(
        [
            brentq(lambda v: bm3(v, v0, k0, kp) - p, 0.5 * v0, 1.05 * v0)
            for p in np.atleast_1d(pressure)
        ]
    )


def bulk_modulus(volume, v0, k0, kp):
    eta = (v0 / np.asarray(volume)) ** (1.0 / 3.0)
    a = 0.75 * (kp - 4)
    return (
        0.5
        * k0
        * (
            (7 * eta**7 - 5 * eta**5) * (1 + a * (eta**2 - 1))
            + (eta**7 - eta**5) * 2 * a * eta**2
        )
    )


def reference_volume(temperature):
    return 40.917 * np.exp(6.69e-5 * (np.asarray(temperature) - 300.0))


def conditional_continuous_pressure(volume, temperature):
    """Audit only: assume Table 6's 300 K values anchor reported slopes.

    The source does not print the separately regressed line intercepts. This
    function therefore must not be represented as an exact author surface.
    """
    delta = np.asarray(temperature) - 300.0
    return bm3(
        volume,
        reference_volume(temperature),
        41.8 - 0.0319 * delta,
        6.6 + 0.00571 * delta,
    )


def metrics(residual):
    residual = np.asarray(residual)
    return {
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "maximum_absolute_pressure_residual_gpa": float(np.max(abs(residual))),
    }


def fit_summary(fit, parameter_order, scale=None):
    """Conditional OLS covariance: residual variance times inverse normal matrix."""
    scale = np.ones(len(fit.x)) if scale is None else np.asarray(scale)
    dof = len(fit.fun) - len(fit.x)
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun @ fit.fun) / dof
    covariance *= scale[:, None] * scale[None, :]
    return {
        "parameters": dict(zip(parameter_order, (fit.x * scale).tolist())),
        "parameter_errors": dict(
            zip(parameter_order, np.sqrt(np.diag(covariance)).tolist())
        ),
        "parameter_covariance": {
            "parameter_order": parameter_order,
            "matrix": covariance.tolist(),
        },
        "degrees_of_freedom": dof,
        "success": bool(fit.success),
        "covariance_scaling": "SSE/(n-p) * inverse(J.T @ J); conditional on fixed V0 and alpha; no calibrant systematic or reference-volume uncertainty included",
        **metrics(fit.fun),
    }


def broad_refits():
    """New Peritheos fits to printed observations, never to plotted curves."""
    rows = observations("fukui_2003", "table4")
    t = np.array([float(r["temperature_k"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    v = np.array([float(r["cell_volume_a3"]) for r in rows])
    hot = t == 873
    fit = least_squares(
        lambda x: bm3(v[hot], reference_volume(873), *x) - p[hot],
        [24, 10],
        bounds=([1, 1], [200, 20]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    result = {"isotherm_873k": fit_summary(fit, ["K0", "K0_prime"])}
    dt = (t - 300) / 1000

    def residual(x):
        return bm3(v, reference_volume(t), x[0] + x[2] * dt, x[1] + x[3] * dt) - p

    fit = least_squares(
        residual,
        [41.8, 6.6, -31.9, 5.71],
        bounds=([1, 1, -60, -8], [200, 20, 60, 20]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    result["continuous_thermal"] = fit_summary(
        fit, ["K0", "K0_prime", "dK_dT", "dKprime_dT"], [1, 1, 0.001, 0.001]
    )
    result["continuous_thermal"]["per_isotherm"] = {
        str(int(temperature)): metrics(fit.fun[t == temperature])
        for temperature in np.unique(t)
    }
    result["objective"] = (
        "Unweighted sum of squared pressure residuals [P_BM3(V_observed,T)-P_observed]^2 in GPa^2; all printed rows; V0(300)=40.917 A^3 and alpha=6.69e-5 K^-1 fixed"
    )
    result["author_workflow_reproduced"] = False
    return result


def figure6c_audit():
    """Manual graphical checkpoints in the original embedded 1292x3452 image.

    Markers test the plot's normalization; solid-line checkpoints test its
    dimensionless equation. Neither set is included in any observation fit.
    """
    audit = source()["sources"]["fukui_2003"]["figure6c_audit"]
    axes = audit["axis_pixels"]
    y0, y025 = axes["ratio_1_y"], axes["ratio_075_y"]
    markers = []
    rows = observations("fukui_2003", "table4")
    for point in audit["filled_marker_centers"]:
        row = rows[point["source_row"] - 1]
        ratio = 1 - 0.25 * (point["y"] - y0) / (y025 - y0)
        markers.append(
            {
                **point,
                "digitized_volume_ratio": ratio,
                "implied_normalization_volume_a3": float(row["cell_volume_a3"]) / ratio,
            }
        )
    curves = []
    for point in audit["solid_curve_checkpoints"]:
        ratio = 1 - 0.25 * (point["y"] - y0) / (y025 - y0)
        pressure = (
            16
            * (point["x"] - axes["pressure_0_x"])
            / (axes["pressure_16_x"] - axes["pressure_0_x"])
        )
        curves.append(
            {
                **point,
                "digitized_volume_ratio": ratio,
                "digitized_pressure_gpa": pressure,
                "printed_BM3_pressure_gpa": float(bm3(ratio, 1, 23.9, 10)),
            }
        )
    values = [m["implied_normalization_volume_a3"] for m in markers]
    return {
        "markers": markers,
        "solid_curve": curves,
        "mean_implied_normalization_volume_a3": float(np.mean(values)),
        "range_implied_normalization_volume_a3": [min(values), max(values)],
        "printed_expansion_reference_volume_a3": float(reference_volume(873)),
        "interpretation": "Solid curve agrees graphically with printed K0=23.9, Kprime=10. Filled observations imply a different V0 normalization, about 42.7 A^3 rather than 42.516 A^3. The mismatch lies within the published table/figure/reference-volume combination; author intent and exact hidden normalization are unrecovered. No corrected source coefficient is inferred; independent fits use printed observations and explicit source V0 only.",
    }


def write_comparison_plot(report):
    """Plot factual table observations and clearly labelled calculated curves."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [
        r for r in observations("fukui_2003", "table4") if r["temperature_k"] == "873"
    ]
    v = np.array([float(r["cell_volume_a3"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    sv = np.array([float(r["cell_volume_error_a3"]) for r in rows])
    sp = np.array([float(r["pressure_error_gpa"]) for r in rows])
    order = np.argsort(v)
    v, p = v[order], p[order]
    sv, sp = sv[order], sp[order]
    fit = report["peritheos_refits"]["isotherm_873k"]["parameters"]
    published = comparison_volume(p, reference_volume(873), 23.9, 10)
    alternate_v0 = comparison_volume(p, 42.70, 23.9, 10)
    refitted = comparison_volume(p, reference_volume(873), fit["K0"], fit["K0_prime"])
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8.2, 7.0),
        sharex=True,
        height_ratios=[2, 1],
        layout="constrained",
    )
    grid = np.linspace(p.min(), p.max(), 200)
    axes[0].errorbar(
        p,
        v,
        xerr=sp,
        yerr=sv,
        fmt="o",
        color="black",
        markersize=4,
        elinewidth=0.9,
        capsize=2,
        label="Table 4 observations + reported errors",
        zorder=3,
    )
    axes[0].plot(
        grid,
        comparison_volume(grid, reference_volume(873), 23.9, 10),
        color="#b64b35",
        label="Published coefficients + printed V0(T)",
    )
    axes[0].plot(
        grid,
        comparison_volume(grid, reference_volume(873), fit["K0"], fit["K0_prime"]),
        color="#176a9b",
        label="Independent Peritheos unweighted refit",
    )
    axes[0].plot(
        grid,
        comparison_volume(grid, 42.70, 23.9, 10),
        "--",
        color="#22804a",
        label="Published K0, K′; assumed V0 = 42.70 Å³",
    )
    axes[0].set(
        ylabel="Cell volume (Å³; Z = 1)",
    )
    fig.suptitle(
        "Brucite at 873 K: testing the reference-volume discrepancy", fontsize=12
    )
    axes[0].legend(fontsize=8)
    axes[1].axhline(0, color="0.55", lw=0.8)
    axes[1].plot(
        p,
        published - v,
        "o-",
        color="#b64b35",
        label=f"Printed V0(T): RMS {np.sqrt(np.mean((published - v) ** 2)):.3f} Å³",
    )
    axes[1].plot(
        p,
        refitted - v,
        "o-",
        color="#176a9b",
        label=f"Peritheos refit: RMS {np.sqrt(np.mean((refitted - v) ** 2)):.3f} Å³",
    )
    axes[1].plot(
        p,
        alternate_v0 - v,
        "s--",
        color="#22804a",
        label=f"V0 = 42.70: RMS {np.sqrt(np.mean((alternate_v0 - v) ** 2)):.3f} Å³",
    )
    axes[1].set(xlabel="Pressure (GPa)", ylabel="Calculated − observed\nvolume (Å³)")
    axes[1].legend(fontsize=8, loc="upper left")
    for ax in axes:
        ax.grid(alpha=0.16)
    fig.savefig(
        ROOT / "docs/data/brucite-873k-comparison.png", dpi=170, bbox_inches="tight"
    )
    plt.close(fig)


def metric_check(rows):
    a = np.array([float(r["lattice_a_angstrom"]) for r in rows])
    c = np.array([float(r["lattice_c_angstrom"]) for r in rows])
    v = np.array([float(r["cell_volume_a3"]) for r in rows])
    difference = np.sqrt(3) / 2 * a**2 * c - v
    return {
        "metric_minus_printed_volume_a3": difference.tolist(),
        "maximum_absolute_difference_a3": float(np.max(abs(difference))),
        "rows_exceeding_printed_volume_error": [
            i
            for i, (d, r) in enumerate(zip(difference, rows), 1)
            if abs(d) > float(r["cell_volume_error_a3"])
        ],
        "interpretation": "Printed a,c,V retained independently; no silent correction.",
    }


def reproduce():
    src = source()["sources"]["fukui_2003"]
    rows = observations("fukui_2003", "table1")
    t = np.array([float(r["temperature_k"]) for r in rows])
    v = np.array([float(r["cell_volume_a3"]) for r in rows])
    sv = np.array([float(r["cell_volume_error_a3"]) for r in rows])
    ambient = {}
    for mode in ("unweighted", "volume_sigma"):
        scale = np.ones_like(sv) if mode == "unweighted" else sv
        fit = least_squares(
            lambda x: (x[0] * np.exp(x[1] * 1e-5 * (t - 300)) - v) / scale,
            [40.917, 6.69],
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
        )
        values = [float(fit.x[0]), float(fit.x[1] * 1e-5)]
        ambient[mode] = {
            "V0_a3": values[0],
            "alpha_per_k": values[1],
            "success": bool(fit.success),
            "delta_over_reported_errors": [
                (values[0] - 40.917) / 0.012,
                (values[1] - 6.69e-5) / 0.18e-5,
            ],
        }

    rows4 = observations("fukui_2003", "table4")
    isotherms = []
    for params in src["parameters"]["isotherms"]:
        temperature = params["temperature_k"]
        group = [r for r in rows4 if float(r["temperature_k"]) == temperature]
        p = np.array([float(r["pressure_gpa"]) for r in group])
        vol = np.array([float(r["cell_volume_a3"]) for r in group])
        sp = np.array([float(r["pressure_error_gpa"]) for r in group])
        sv = np.array([float(r["cell_volume_error_a3"]) for r in group])
        v0 = float(reference_volume(temperature))
        published = [params["K0_gpa"], params["K0_prime"]]
        residual = bm3(vol, v0, *published) - p
        effective_sigma = np.sqrt(
            sp**2 + (bulk_modulus(vol, v0, *published) * sv / vol) ** 2
        )
        fits = {}
        for mode in ("unweighted_pressure", "pressure_sigma", "effective_sigma"):

            def objective(x):
                difference = bm3(vol, v0, *x) - p
                if mode == "pressure_sigma":
                    return difference / sp
                if mode == "effective_sigma":
                    # Diagnostic parameter-dependent variance, not author weights.
                    scale = np.sqrt(sp**2 + (bulk_modulus(vol, v0, *x) * sv / vol) ** 2)
                    return difference / scale
                return difference

            fit = least_squares(
                objective,
                published,
                bounds=([1, 1], [200, 20]),
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
            )
            fits[mode] = {
                "success": bool(fit.success),
                "K0_gpa": float(fit.x[0]),
                "K0_prime": float(fit.x[1]),
                "delta_over_reported_errors": [
                    float((fit.x[0] - published[0]) / params["K0_error_gpa"]),
                    float((fit.x[1] - published[1]) / params["K0_prime_error"]),
                ],
                **metrics(bm3(vol, v0, *fit.x) - p),
            }
        isotherms.append(
            {
                "temperature_k": temperature,
                "observations": len(group),
                "source_rows": [int(r["source_row"]) for r in group],
                "pressure_range_gpa": [float(min(p)), float(max(p))],
                "reference_volume_a3": v0,
                "published_parameters": published,
                "published_curve": {
                    **metrics(residual),
                    "pressure_residuals_gpa": residual.tolist(),
                    "max_absolute_effective_sigma_residual": float(
                        np.max(abs(residual / effective_sigma))
                    ),
                },
                "diagnostic_refits": fits,
                "source_equation_validated": True,
                "author_workflow_reproduced": False,
            }
        )

    # The slope fits test whether rounded slice parameters explain the reported
    # temperature slopes, without assigning an unreported intercept to the EOS.
    slices = src["parameters"]["isotherms"]
    temperatures = np.array([x["temperature_k"] for x in slices]) - 300.0
    trends = {}
    for name, error in (("K0_gpa", "K0_error_gpa"), ("K0_prime", "K0_prime_error")):
        y = np.array([x[name] for x in slices])
        sigma = np.array([x[error] for x in slices])
        for mode in ("unweighted", "reported_sigma"):
            scale = np.ones_like(sigma) if mode == "unweighted" else sigma
            design = np.column_stack((np.ones_like(temperatures), temperatures))
            x = np.linalg.lstsq(design / scale[:, None], y / scale, rcond=None)[0]
            trends[f"{name}_{mode}"] = {
                "intercept_300k": float(x[0]),
                "slope_per_k": float(x[1]),
            }

    csvs = {}
    for owner, table, filename in TABLES:
        payload = normalized_csv(owner, table).encode()
        csvs[filename] = {
            "rows": len(observations(owner, table)),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }
    return {
        "audit_date": "2026-10-16",
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "optimizer": "scipy.optimize.least_squares; explicit tolerances and bounds in script",
        },
        "source_json_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "datasets": csvs,
        "volume_basis": source()["volume_basis"],
        "ambient_published_residuals_a3": (reference_volume(t) - v).tolist(),
        "ambient_diagnostic_refits": ambient,
        "isotherms": isotherms,
        "table6_temperature_trends": trends,
        "figure6c_audit": figure6c_audit(),
        "peritheos_refits": broad_refits(),
        "metric_checks": {
            "fukui_table1": metric_check(rows),
            "fukui_table4": metric_check(rows4),
            "xia_table6": metric_check(observations("xia_1998", "table6")),
        },
        "assessments": {
            "source_equation": "Fukui Eq. 2 is BM3; four source isotherms executable.",
            "coefficients": "Diagnostic fits use explicit assumed weights; see deltas.",
            "author_workflow": "Not reproduced: weights, residual orientation, full precision, covariance and static V0 treatment unreported.",
            "continuous_thermal": "Published-slope reconstruction remains qualified. A separate Peritheos unweighted thermal fit to all 37 printed PVT rows provides a reproducible broad-range alternative with conditional covariance; neither is an author-workflow reproduction.",
            "phase": "Layered lattice-average volume EOS accepted through the published range; H rearrangement is annotated, not treated as a first-order lattice boundary. Sun (2025) primary XRD supports this scope. Ambient H sites are not high-pressure refinements; no phase-equilibrium or dehydration model is implied.",
            "calibration": "Anderson et al. (1989) Au scale identified; paired Au observations absent, so observation-level re-reduction unavailable.",
            "structure": "Direct primary verification resolved using He (2023), Table 4.3, pp. 81-82: fluorine-free hydrogenated Mg(OH)2 neutron refinement. Original Zigan COD archive preserved separately, not claimed directly checked.",
            "refit_records": "Separate Peritheos unweighted 873 K and continuous thermal fits preserve all published coefficients in their original records. Conditional covariance and fixed reference-volume assumptions are explicit.",
        },
        "catalog_scope": source()["catalog_disposition"],
    }


def check_fukui_pdf():
    """Check table tokens in the exact final PDF used for visual verification."""
    src = source()["sources"]["fukui_2003"]
    path = ROOT / "output/material-candidates/brucite/sources/brucite-fukui-2003.pdf"
    if hashlib.sha256(path.read_bytes()).hexdigest() != src["pdf_sha256"]:
        raise ValueError("Fukui source PDF checksum mismatch")
    result = subprocess.run(
        ["pdftotext", "-layout", str(path), "-"],
        check=True,
        capture_output=True,
        text=True,
    )
    # The embedded font maps the printed minus on the first pressure to ')'.
    # This single visually verified substitution does not alter other tokens.
    lines = [
        re.sub(r"\s+", " ", x.replace(")0.06(13)", "-0.06(13)")).strip()
        for x in result.stdout.splitlines()
    ]
    for table in ("table1", "table4"):
        for index, row in enumerate(src[table]["rows"], 1):
            # Text extraction inserts a footnote 'a' after the first datum.
            expected = " ".join(row)
            found = any(
                expected in re.sub(r"(?<=\))a(?= |$)", "", line) for line in lines
            )
            if not found:
                raise ValueError(f"Fukui {table} row {index}: PDF token mismatch")
    return {"pdf_sha256": src["pdf_sha256"], "table1_rows": 11, "table4_rows": 37}


def check_structure_and_phase_pdfs():
    sources = source()["sources"]
    checks = {}
    for name, filename in (
        ("structure", "brucite-thesis-2024.pdf"),
        ("duffy_1995_phase", "brucite-duffy-1995.pdf"),
        ("wu_2026_phase", "brucite-h-disorder-2026.pdf"),
        ("sun_2025_phase", "brucite-sun-2025.pdf"),
    ):
        path = ROOT / "output/material-candidates/brucite/sources" / filename
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != sources[name]["pdf_sha256"]:
            raise ValueError(f"Primary {name} PDF checksum mismatch")
        checks[name] = digest
        if name == "structure":
            result = subprocess.run(
                ["pdftotext", "-f", "89", "-l", "90", "-layout", str(path), "-"],
                check=True,
                capture_output=True,
                text=True,
            )
            for token in (
                "3.14771(3)",
                "4.77163(5)",
                "40.9437(7)",
                "0.2189(2)",
                "0.4108(4)",
                "0.915(2)",
                "0.9528",
            ):
                if token not in result.stdout:
                    raise ValueError(f"He (2023) Table 4.3 token mismatch: {token}")
    return checks


def ledger_outcome(record):
    """Separate original coefficients, independent fits, and reconstruction."""
    if "curved_expansion" in record["identifier"]:
        from scripts.brucite_curved_expansion import ledger_outcome as curved_outcome

        return curved_outcome(record)
    report = reproduce()
    datasets = ["brucite_fukui_2003_table1", "brucite_fukui_2003_table4"]
    common = {
        "dataset_identifiers": datasets,
        "original_author_workflow_parity": False,
        "phase_pure_fit": False,
        "accepted_pressure_gpa": record["validity"]["pressure_gpa"],
        "reproduction_report": "docs/data/brucite-reproduction.json",
    }
    if record.get("record_kind") == "refit":
        key = "continuous_thermal" if "thermal" in record else "isotherm_873k"
        fit = report["peritheos_refits"][key]
        stored = {
            **record["eos"]["parameters"],
            **record.get("thermal", {}).get("parameters", {}),
        }
        comparisons = [
            {
                "parameter": name,
                "published": stored[name],
                "refit": value,
                "published_error": (
                    record.get("parameter_errors", {})
                    | record.get("thermal", {}).get("parameter_errors", {})
                )[name],
                "refit_error": fit["parameter_errors"][name],
                "difference": value - stored[name],
                "relative_difference": abs((value - stored[name]) / stored[name]),
                "within_combined_2sigma": None,
                "similar": bool(np.isclose(value, stored[name], rtol=1e-4, atol=1e-24)),
            }
            for name, value in fit["parameters"].items()
        ]
        rows = observations("fukui_2003", "table4")
        selected = [
            r for r in rows if "thermal" in record or float(r["temperature_k"]) == 873
        ]
        temperatures = np.array([float(r["temperature_k"]) for r in selected])
        volumes = np.array([float(r["cell_volume_a3"]) for r in selected])
        delta = (
            temperatures - 300 if "thermal" in record else np.zeros_like(temperatures)
        )

        def predictions(parameters):
            return bm3(
                volumes,
                reference_volume(temperatures),
                parameters["K0"] + parameters.get("dK_dT", 0) * delta,
                parameters["K0_prime"] + parameters.get("dKprime_dT", 0) * delta,
            )

        # Tiny, correlated slopes are poorly conditioned under finite-difference
        # optimization. Require pressure agreement as well as coefficient agreement.
        curve_difference = float(
            np.max(np.abs(predictions(stored) - predictions(fit["parameters"])))
        )
        if (
            not fit["success"]
            or not all(p["similar"] for p in comparisons)
            or curve_difference > 1e-5
        ):
            raise ValueError("Registered brucite Peritheos refit is stale")
        return {
            **common,
            "status": "parity",
            "parity_basis": "registered_peritheos_refit_reproduced",
            "reproduction_status": "independent_refit_reproduced",
            "original_publication_reproduction_status": "not_reproduced",
            "observations": 37 if "thermal" in record else 8,
            "parameters": comparisons,
            "fit": fit,
            "rmse_gpa": fit["pressure_rmse_gpa"],
            "solver_success": fit["success"],
            "objective": report["peritheos_refits"]["objective"],
            "reason": "New Peritheos unweighted fit to printed Table 4 observations with source V0(T) fixed; reproducible with conditional covariance. No agreement with original author coefficients or workflow is claimed.",
        }
    if "thermal" in record:
        return {
            **common,
            "status": "source_reconstruction",
            "reconstruction_kind": "qualified_300k_anchored_linear_temperature_laws",
            "observations": 48,
            "source_intercepts_reproduced": False,
            "parameters": [],
            "reason": "Literal 300 K anchored reconstruction of reported slopes; separately regressed intercepts and author workflow unrecovered. Explicit-selection record over the common measured pressure interval, 300-873 K. Separate Peritheos refits are available; no phase-pure author-fit claim.",
        }
    temperature = record["validity"]["temperature_k"][0]
    result = next(s for s in report["isotherms"] if s["temperature_k"] == temperature)
    fit = result["diagnostic_refits"]["unweighted_pressure"]
    comparisons = []
    for name, fitted in (("K0", "K0_gpa"), ("K0_prime", "K0_prime")):
        value = record["eos"]["parameters"][name]
        error = record["parameter_errors"][name]
        difference = fit[fitted] - value
        comparisons.append(
            {
                "parameter": name,
                "published": value,
                "published_error": error,
                "refit": fit[fitted],
                "refit_error": None,
                "difference": difference,
                "relative_difference": abs(difference / value),
                "within_reported_error": abs(difference) <= error,
                "within_combined_2sigma": None,
                "similar": abs(difference) <= error,
            }
        )
    agreement = all(p["within_reported_error"] for p in comparisons)
    return {
        **common,
        "status": "parity" if agreement else "parity_not_achieved",
        "parity_basis": "within_reported_parameter_errors",
        "parameters": comparisons,
        "observations": result["observations"],
        "fit_kind": "unweighted_pressure_residuals_fixed_source_V0",
        "published_rmse_gpa": result["published_curve"]["pressure_rmse_gpa"],
        "rmse_gpa": fit["pressure_rmse_gpa"],
        "solver_success": fit["success"],
        "qualification": "Conditional comparison against individual reported errors using every printed same-temperature row and the source expansion reference volume fixed. Unweighted pressure residuals, not author weights. At 873 K the K0 difference exceeds one reported error and the figure implies inconsistent volume normalization. Lattice-average fit includes H rearrangement; no phase-pure or original-workflow claim.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source-check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.source_check:
        print(json.dumps(check_fukui_pdf(), indent=2))
        print(json.dumps(check_structure_and_phase_pdfs(), indent=2))
    report = reproduce()
    if args.plot:
        write_comparison_plot(report)
    products = {
        DATA / filename: normalized_csv(owner, table)
        for owner, table, filename in TABLES
    }
    products[REPORT] = json.dumps(report, indent=2, allow_nan=False) + "\n"
    for path, payload in products.items():
        if args.write:
            path.write_text(payload, encoding="utf-8")
        if args.check:
            if path == REPORT:
                from scripts.brucite_replay_checks import check_report

                check_report(report, json.loads(path.read_text(encoding="utf-8")))
            elif path.read_text(encoding="utf-8") != payload:
                raise SystemExit(f"Outdated brucite artifact: {path}")
    if not (args.check or args.write):
        print(products[REPORT])
    elif args.check:
        print("Brucite tables and reproduction report are current.")


if __name__ == "__main__":
    main()
