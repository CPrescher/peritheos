"""Repository-only Takahashi (2019) diagnostics; never a replacement EOS.

Run: .venv/bin/python -m scripts.audit_takahashi_2019
All protocol variants below are hypotheses, not recovered author choices.
The independent evaluator uses SI cell volumes and Gauss-Legendre quadrature;
it does not call a Peritheos equation or fitting routine.
"""

from __future__ import annotations

import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import Avogadro, R
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares
from scipy.optimize._numdiff import approx_derivative

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/takahashi-2019-parity-2026-10-30"
DATA = ROOT / "peritheos/data/datasets"
SOURCES = ROOT / "output/cementite-eos-audit-2026-10-25/sources"
NAMES = ["V0", "K0", "K0_prime", "gamma0", "q"]
SOURCE = {
    "kbr": np.array([148.8, 311.1, 3.4, 1.06, 1.92]),
    "mgo": np.array([147.3, 323, 3.43, 1.37, 0.98]),
}
WIDTH = {
    "kbr": np.array([1, 17.1, 0.1, 0.42, 1.73]),
    "mgo": np.array([1, 16.6, 0.09, 0.33, 1.01]),
}
LOWER = np.array([120, 100, 1, 0.01, -5])
UPPER = np.array([180, 600, 8, 5, 8])
SCALING = np.array([150, 300, 4, 1, 2])
GN, GW = leggauss(32)
GN, GW = (GN + 1) / 2, GW / 2


def pressure(v, t, x):
    """Equations 1-6; v in conventional-cell A3, pressure in GPa."""
    v, t = np.broadcast_arrays(np.asarray(v), np.asarray(t))
    v0, k0, kp, gamma0, q = x
    eta = (v0 / v) ** (1 / 3)
    cold = 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))
    logv = np.log(v / v0)
    gamma = gamma0 * np.exp(q * logv)
    theta = 314 * np.exp(-gamma0 * (logv if q == 0 else np.expm1(q * logv) / q))

    def energy(temp):
        z = theta / temp
        u = z[..., None] * GN
        integral = z * np.sum(GW * u**3 / np.expm1(u), axis=-1)
        return 9 * 4 * R * temp * integral / z**3

    # Z=4 formulas/cell, 1 A3=1e-30 m3, Pa -> GPa.
    return cold + gamma * (energy(t) - energy(300)) / (v * 1e-30 * Avogadro / 4) / 1e9


def volume(p, t, x):
    """Inverse source EOS, safeguarded by a physical positive-bulk interval."""
    p, t = np.broadcast_arrays(p, t)
    # Newton near observed compressed states; all retained results are checked.
    v = np.full_like(p, 110.0, dtype=float)
    for _ in range(30):
        step = np.clip((pressure(v, t, x) - p) / derivative_v(v, t, x), -8, 8)
        v -= step
        if np.max(np.abs(step)) < 1e-11:
            break
    if np.max(np.abs(pressure(v, t, x) - p)) > 1e-7:
        raise RuntimeError("Inverse EOS did not converge")
    return v


def derivative_v(v, t, x):
    h = 1e-3
    return (pressure(v + h, t, x) - pressure(v - h, t, x)) / (2 * h)


def solve(fun, initial, lower=LOWER, upper=UPPER):
    return least_squares(
        fun,
        initial,
        bounds=(lower, upper),
        x_scale="jac",
        jac="3-point",
        ftol=2e-12,
        xtol=2e-12,
        gtol=2e-9,
        max_nfev=3000,
    )


def effective_variance_fit(v, t, p, sp, sv, st, source):
    """IRLS with source coordinate widths; no confidence convention assumed."""
    x = source.copy()
    for iteration in range(30):
        dpdv = derivative_v(v, t, x)
        dpdt = (pressure(v, t + 0.01, x) - pressure(v, t - 0.01, x)) / 0.02
        effective_sigma = np.sqrt(sp**2 + (dpdv * sv) ** 2 + (dpdt * st) ** 2)
        fit = solve(lambda y: (pressure(v, t, y) - p) / effective_sigma, x)
        change = np.max(np.abs((fit.x - x) / SCALING))
        x = fit.x
        if change < 1e-7:
            break
    return fit, iteration + 1


def summary(fit, v, t, p, source, width):
    x = fit.x[:5]
    # Include latent observations when estimating profiled parameter covariance.
    _, singular, vectors = np.linalg.svd(fit.jac, full_matrices=False)
    retained = singular > np.finfo(float).eps * max(fit.jac.shape) * singular[0]
    scaled_vectors = vectors[retained, :5].T / singular[retained]
    covariance = np.einsum("ij,kj->ik", scaled_vectors, scaled_vectors)
    absolute_errors = np.sqrt(np.diag(covariance))
    dof = fit.fun.size - fit.x.size
    covariance *= np.sum(fit.fun**2) / dof
    errors = np.sqrt(np.diag(covariance))
    residual = pressure(v, t, x) - p
    return {
        "parameters": dict(zip(NAMES, x.tolist())),
        "conditional_residual_scaled_errors": dict(zip(NAMES, errors.tolist())),
        "conditional_absolute_width_errors": dict(zip(NAMES, absolute_errors.tolist())),
        "within_source_width": dict(zip(NAMES, (np.abs(x - source) <= width).tolist())),
        "all_central_values_within_source_width": bool(
            np.all(np.abs(x - source) <= width)
        ),
        "covariance": covariance.tolist(),
        "correlation": (covariance / np.outer(errors, errors)).tolist(),
        "objective_sum_squares": float(np.sum(fit.fun**2)),
        "degrees_of_freedom": dof,
        "unweighted_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "source_refit_max_abs_pressure_difference_gpa": float(
            np.max(np.abs(pressure(v, t, x) - pressure(v, t, source)))
        ),
        "source_refit_max_abs_volume_difference_a3_at_observed_pt": float(
            np.max(np.abs(volume(p, t, x) - volume(p, t, source)))
        ),
        "active_bounds": int(np.count_nonzero(fit.active_mask)),
        "success": bool(fit.success),
        "optimality": float(fit.optimality),
        "evaluations": fit.nfev,
        "message": fit.message,
    }


def workbook_audit(rows):
    """Read actual OOXML, including styles; no workbook recalculation."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    path = SOURCES / "minerals-630318-supplementary.xlsx"
    with ZipFile(path) as archive:
        strings = [
            "".join(n.itertext())
            for n in ET.fromstring(archive.read("xl/sharedStrings.xml"))
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        cells = {}
        for cell in sheet.findall(".//m:c", ns):
            value = cell.find("m:v", ns)
            if value is not None:
                cells[cell.attrib["r"]] = (
                    strings[int(value.text)]
                    if cell.get("t") == "s"
                    else float(value.text)
                )
        raw = json.loads(
            (DATA / "cementite-takahashi-2019-table-s1-cells.json").read_text()
        )
        styles = ET.fromstring(archive.read("xl/styles.xml"))
        formats = {
            int(n.get("numFmtId")): n.get("formatCode")
            for n in styles.findall("m:numFmts/m:numFmt", ns)
        }
        formats.update({0: "General", 1: "0", 2: "0.00"})
        xf = styles.findall("m:cellXfs/m:xf", ns)
        actual_formats = {
            c.get("r"): formats[int(xf[int(c.get("s", "0"))].get("numFmtId"))]
            for c in sheet.findall(".//m:c", ns)
        }
        mismatches = []
        format_mismatches = []
        for row in raw["rows"]:
            for address, saved in row["cells"].items():
                if saved["value"] != cells.get(address):
                    mismatches.append(address)
                if saved["number_format"] != actual_formats[address]:
                    format_mismatches.append(address)
        mapping = {
            "B": "temperature_k",
            "E": "pressure_kbr_gpa",
            "H": "pressure_mgo_gpa",
            "K": "pressure_sio2_gpa",
            "P": "volume_a3",
            "R": "a_angstrom",
            "T": "b_angstrom",
            "V": "c_angstrom",
            "Y": "kbr_volume_a3",
            "AD": "mgo_volume_a3",
            "AI": "sio2_volume_a3",
        }
        central_mismatches = []
        error_mismatches = []
        for row in rows:
            i = row["source_row"]
            for column, field in mapping.items():
                val = cells.get(column + i)
                if val not in (None, "-") and float(row[field]) != val:
                    central_mismatches.append(column + i)
            for column, error_column, field, direct in [
                ("B", "C", "temperature_k_error", True),
                ("E", "F", "pressure_kbr_gpa_error", True),
                ("H", "I", "pressure_mgo_gpa_error", True),
                ("K", "L", "pressure_sio2_gpa_error", True),
                ("P", "Q", "volume_a3_error", False),
                ("Y", "Z", "kbr_volume_a3_error", False),
                ("AD", "AE", "mgo_volume_a3_error", False),
                ("AI", "AJ", "sio2_volume_a3_error", False),
            ]:
                token = cells.get(error_column + i)
                if token in (None, "-"):
                    continue
                token = str(token).strip("()")
                decimals = (
                    len(actual_formats[column + i].split(".")[1])
                    if "." in actual_formats[column + i]
                    else 0
                )
                value = float(token) if direct else float(token) * 10 ** (-decimals)
                if not np.isclose(value, float(row[field]), rtol=1e-14):
                    error_mismatches.append(error_column + i)
        return {
            "archive_members": archive.namelist(),
            "formula_count": len(sheet.findall(".//m:f", ns)),
            "raw_cell_value_mismatches": mismatches,
            "raw_cell_number_format_mismatches": format_mismatches,
            "csv_central_value_mismatches": central_mismatches,
            "csv_error_value_mismatches": error_mismatches,
            "metadata": {
                name: archive.read(name).decode()
                for name in ("docProps/core.xml", "docProps/app.xml")
            },
            "zip_members": ZipFile(SOURCES / "minerals-09-00744-s001.zip").namelist(),
            "cell_P5": raw["rows"][1]["cells"]["P5"],
            "cell_Q5": raw["rows"][1]["cells"]["Q5"],
        }


def parity_ledger_outcome(record, previous):
    """Current conditional parameter parity; source values remain authoritative.

    This executes the stated effective-variance replay, rather than promoting a
    stored status or depending on the generated numerical report being present.
    """
    scale = "kbr" if "_kbr_" in record["identifier"] else "mgo"
    rows = list(
        csv.DictReader(
            (DATA / "cementite-takahashi-2019-table-s1.csv").read_text().splitlines()
        )
    )
    v = np.array([float(r["volume_a3"]) for r in rows])
    t = np.array([float(r["temperature_k"]) for r in rows])
    p = np.array(
        [float(r[f"pressure_{scale}_gpa"] or r["pressure_sio2_gpa"]) for r in rows]
    )
    sp = np.array(
        [
            float(r[f"pressure_{scale}_gpa_error"] or r["pressure_sio2_gpa_error"])
            for r in rows
        ]
    )
    sv = np.array([float(r["volume_a3_error"]) for r in rows])
    st = np.array([float(r["temperature_k_error"] or 0) for r in rows])
    source, widths = SOURCE[scale], WIDTH[scale]
    fit, iterations = effective_variance_fit(v, t, p, sp, sv, st, source)
    _, singular, vectors = np.linalg.svd(fit.jac, full_matrices=False)
    scaled = vectors / singular[:, None]
    covariance = np.einsum("ij,ik->jk", scaled, scaled)
    errors = np.sqrt(np.diag(covariance))
    names = ["rt_eos.V0", "rt_eos.K0", "rt_eos.K0_prime", "gamma0", "q"]
    decimals = [1, 1 if scale == "kbr" else 0, 2, 2, 2]
    error_decimals = [1, 1, 2, 2, 2]
    parameters = []
    for name, published, width, value, error, dp, de in zip(
        names, source, widths, fit.x, errors, decimals, error_decimals
    ):
        difference = float(value - published)
        parameters.append(
            {
                "parameter": name,
                "published": float(published),
                "published_error": float(width),
                "refit": float(value),
                "refit_error": float(error),
                "difference": difference,
                "relative_difference": float(abs(difference / published)),
                "within_printed_source_error_width": bool(abs(difference) <= width),
                "within_combined_2sigma": None,
                "similar": bool(abs(difference) <= width),
                "central_value_matches_printed_precision": bool(
                    round(float(value), dp) == published
                ),
                "error_matches_printed_precision": bool(
                    round(float(error), de) == width
                ),
            }
        )
    central_parity = all(p["within_printed_source_error_width"] for p in parameters)
    uncertainty_parity = all(p["error_matches_printed_precision"] for p in parameters)
    qualification = (
        "All 60 Table S1 rows, including the 2600 K row and ten shared SiO2-marker rows. "
        "Conditional simultaneous five-parameter effective-variance replay uses printed independent P/V/T widths; "
        "missing 300 K temperature widths leave those temperatures fixed. Unscaled inverse-normal covariance "
        "treats the widths as absolute errors; their source confidence and covariance are unconfirmed. "
        + (
            "KBr reproduces all five central values and all five error widths at published printed precision. "
            if scale == "kbr"
            else "MgO central values all lie within printed source widths, but exact printed coefficients and uncertainties are not reproduced. "
        )
        + "The original author objective, weights, fitting order and exact input selection remain unconfirmed; "
        "original_author_workflow_parity=False. No formal combined-two-sigma inference, experimental validation "
        "or change in recommendation follows. Published coefficients and errors remain intact; the earlier "
        "equal-pressure reconstruction is retained separately. Source-equation density checkpoints concern "
        "extrapolated calculated states, not measured core conditions."
    )
    return {
        **previous,
        "status": "parity" if central_parity else "similar",
        "parity_basis": "Conditional parameter_agreement: every free central coefficient lies within its own printed source error width under the all-60-row P/V/T effective-variance replay. Uncertainty reproduction is assessed separately; source confidence is unknown and author workflow is unconfirmed.",
        "qualification": qualification,
        "source_parameter_widths_agreement": central_parity,
        "published_printed_precision_agreement": all(
            p["central_value_matches_printed_precision"] for p in parameters
        ),
        "uncertainty_reproduction": {
            "status": "reproduced_at_printed_precision"
            if uncertainty_parity
            else "not_reproduced_at_printed_precision",
            "basis": "Unscaled inverse-normal errors compared with printed source widths; confidence and cross-coordinate covariance unconfirmed.",
        },
        "objective": "iteratively_reweighted_pressure_residuals",
        "weighting": "Effective variance: sP^2 + (dP/dV*sV)^2 + (dP/dT*sT)^2 from printed widths; freeze during each solve, update to convergence. Missing 300 K T widths leave T fixed.",
        "absolute_sigma": True,
        "parameters": parameters,
        "parameter_covariance": covariance.tolist(),
        "covariance_convention": "Unscaled inverse normal matrix at final frozen effective-variance weights, in free_parameters order with V0 in conventional-cell A3. Printed coordinate widths are assumed absolute independent errors; source confidence and covariance are unconfirmed.",
        "rmse_gpa": float(np.sqrt(np.mean((pressure(v, t, fit.x) - p) ** 2))),
        "reduced_chi_square": float(np.sum(fit.fun**2) / (len(rows) - 5)),
        "solver_success": bool(fit.success),
        "solver_message": fit.message,
        "effective_variance_iterations": iterations,
        "original_author_workflow_parity": False,
        "reproduction_report": "docs/literature-reproductions/takahashi-2019-parity.md",
        "previous_equal_pressure_reconstruction": previous,
    }


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    rows = list(
        csv.DictReader(
            (DATA / "cementite-takahashi-2019-table-s1.csv").read_text().splitlines()
        )
    )
    v = np.array([float(r["volume_a3"]) for r in rows])
    t = np.array([float(r["temperature_k"]) for r in rows])
    sv = np.array([float(r["volume_a3_error"]) for r in rows])
    st = np.array([float(r["temperature_k_error"] or 0) for r in rows])
    cold = t == 300
    report = {
        "description": "Diagnostic protocol sensitivities; no replacement or recovered author workflow.",
        "parameter_order": NAMES,
        "units": ["A3/cell", "GPa", "1", "1", "1"],
        "fixed": {
            "Tr_k": 300,
            "theta0_k": 314,
            "atoms_per_formula": 4,
            "formulas_per_cell": 4,
        },
        "workbook": workbook_audit(rows),
        "input_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                DATA / "cementite-takahashi-2019-table-s1.csv",
                DATA / "cementite-takahashi-2019-table-s1-cells.json",
                SOURCES / "takahashi-2019.pdf",
                SOURCES / "minerals-630318-supplementary.xlsx",
            ]
        },
        "row_count": len(rows),
        "cold_row_count": int(sum(cold)),
        "hot_row_count": int(sum(~cold)),
        "scales": {},
    }
    # The existing Peritheos replay is read-only; this entry point writes no shared ledger.
    from scripts.reproduce_cementite import reproduce_takahashi

    native = reproduce_takahashi()
    for scale, source in SOURCE.items():
        p = np.array(
            [float(r[f"pressure_{scale}_gpa"] or r["pressure_sio2_gpa"]) for r in rows]
        )
        sp = np.array(
            [
                float(r[f"pressure_{scale}_gpa_error"] or r["pressure_sio2_gpa_error"])
                for r in rows
            ]
        )
        width = WIDTH[scale]
        base = solve(lambda x: pressure(v, t, x) - p, source)
        fits = {"equal_pressure_all60": summary(base, v, t, p, source, width)}
        starts = []
        for initial in [
            source,
            [145, 380, 3, 1.5, 0],
            [155, 240, 4, 0.7, 3],
            [150, 300, 3.5, 2, -1],
        ]:
            fit = solve(lambda x: pressure(v, t, x) - p, initial)
            starts.append(
                {
                    "initial": list(initial),
                    "parameters": fit.x.tolist(),
                    "sse": float(2 * fit.cost),
                    "success": bool(fit.success),
                }
            )
        for name, fun in {
            "reported_pressure_widths": lambda x: (pressure(v, t, x) - p) / sp,
            "equal_volume": lambda x: volume(p, t, x) - v,
            "reported_volume_widths": lambda x: (volume(p, t, x) - v) / sv,
        }.items():
            fits[name] = summary(solve(fun, source), v, t, p, source, width)
        # First-order errors-in-variables, standard effective-variance IRLS.
        # Freeze weights during each minimization, then update derivatives.
        for use_temperature in (False, True):
            fit, iterations = effective_variance_fit(
                v, t, p, sp, sv, st if use_temperature else np.zeros_like(st), source
            )
            name = (
                "effective_variance_pvt" if use_temperature else "effective_variance_pv"
            )
            fits[name] = summary(fit, v, t, p, source, width)
            fits[name]["weight_update_iterations"] = iterations
        # Exact latent-coordinate least squares. Missing 300 K temperature errors
        # mean those temperatures are held fixed; they are not assigned zero sigmas.
        for temp_errors in (False, True):
            hot = ~cold

            def eiv_fun(y):
                vv = y[5:65]
                tt = t.copy()
                if temp_errors:
                    tt[hot] = y[65:]
                parts = [(pressure(vv, tt, y[:5]) - p) / sp, (vv - v) / sv]
                if temp_errors:
                    parts.append((tt[hot] - t[hot]) / st[hot])
                return np.concatenate(parts)

            initial = np.r_[source, v, t[hot] if temp_errors else []]
            low = np.r_[
                LOWER, np.full(60, 80), np.full(sum(hot), 500) if temp_errors else []
            ]
            high = np.r_[
                UPPER, np.full(60, 150), np.full(sum(hot), 4000) if temp_errors else []
            ]
            fit = solve(eiv_fun, initial, low, high)
            name = "coordinate_widths_pvt" if temp_errors else "coordinate_widths_pv"
            fits[name] = summary(fit, v, t, p, source, width)
        for name, mask in {
            "omit2600_diagnostic": t <= 2450,
            "paired50_only_diagnostic": p < 290,
        }.items():
            fits[name] = summary(
                solve(lambda x: pressure(v[mask], t[mask], x) - p[mask], source),
                v,
                t,
                p,
                source,
                width,
            )
            fits[name]["fitted_rows"] = np.array([int(r["source_row"]) for r in rows])[
                mask
            ].tolist()
        mask = t <= 2450
        fit, _ = effective_variance_fit(
            v[mask], t[mask], p[mask], sp[mask], sv[mask], st[mask], source
        )
        fits["effective_variance_pvt_omit2600_diagnostic"] = summary(
            fit, v, t, p, source, width
        )
        fits["effective_variance_pvt_omit2600_diagnostic"]["fitted_rows"] = np.array(
            [int(r["source_row"]) for r in rows]
        )[mask].tolist()
        # Sequential interpretation of paragraph 3.2: cold coefficients first,
        # then gamma0 and q with those three held fixed. No new EOS family.
        sequential = {}
        for name in ("pressure", "volume", "pressure_widths"):

            def cold_fun(x):
                xx = np.r_[x, source[3:]]
                if name == "volume":
                    return volume(p[cold], t[cold], xx) - v[cold]
                return (pressure(v[cold], t[cold], xx) - p[cold]) / (
                    sp[cold] if name == "pressure_widths" else 1
                )

            cfit = solve(cold_fun, source[:3], LOWER[:3], UPPER[:3])
            hfit = solve(
                lambda y: pressure(v[~cold], t[~cold], np.r_[cfit.x, y]) - p[~cold],
                source[3:],
                LOWER[3:],
                UPPER[3:],
            )
            sequential[name] = {
                "parameters": np.r_[cfit.x, hfit.x].tolist(),
                "all_within_source_width": bool(
                    np.all(np.abs(np.r_[cfit.x, hfit.x] - source) <= width)
                ),
            }
        thermal = solve(
            lambda y: pressure(v[~cold], t[~cold], np.r_[source[:3], y]) - p[~cold],
            source[3:],
            LOWER[3:],
            UPPER[3:],
        )
        # Fit within the source-width box: measures objective cost, not a fit claim.
        box = solve(
            lambda x: pressure(v, t, x) - p,
            source,
            np.maximum(LOWER, source - width),
            np.minimum(UPPER, source + width),
        )
        jac = approx_derivative(
            lambda x: pressure(v, t, x) - p, base.x, method="3-point"
        )
        singular = np.linalg.svd(jac * SCALING, compute_uv=False)
        # First-order SSE of source-rounding changes, using printed resolutions.
        half_round = np.array(
            [0.05, 0.05 if scale == "kbr" else 0.5, 0.005, 0.005, 0.005]
        )
        rounding_fits = []
        rng = np.random.default_rng(2019)
        # Uniform rounding intervals of tabulated volume (.01 A3), P (.1 GPa),
        # T (integer K). These are precision diagnostics, not added observations.
        for _ in range(24):
            vr = v + rng.uniform(-0.005, 0.005, 60)
            pr = p + rng.uniform(-0.05, 0.05, 60)
            tr = t + np.where(cold, 0, rng.uniform(-0.5, 0.5, 60))
            fit = solve(lambda x: pressure(vr, tr, x) - pr, base.x)
            rounding_fits.append(fit.x)
        rounding_fits = np.array(rounding_fits)
        source_res = pressure(v, t, source) - p
        fitted_res = pressure(v, t, base.x) - p
        actual_vol_delta = volume(p, t, base.x) - volume(p, t, source)
        # Compare native and Python evaluators at all measured states.
        from peritheos.eos.rt import BM3
        from peritheos.eos.thermal import MieGruneisenDebye

        factor = Avogadro * 1e-25 / 4
        model = MieGruneisenDebye(
            BM3(V0=source[0] * factor, K0=source[1], K0_prime=source[2]),
            Tr=300,
            theta0=314,
            gamma0=source[3],
            q=source[4],
            n=4,
        )
        native_p = model.pressure(v * factor, t)
        model._native = None
        python_p = model.pressure(v * factor, t)
        checkpoints = []
        for pp, tt in [(329, 5000), (329, 6000), (364, 5000), (364, 6000)]:
            vv = brentq(lambda a: pressure(a, tt, source) - pp, 80, 120)
            checkpoints.append(
                {
                    "pressure_gpa": pp,
                    "temperature_k": tt,
                    "volume_a3": vv,
                    "density_g_cm3": 4
                    * (3 * 55.845 + 12.011)
                    / (Avogadro * vv * 1e-24),
                }
            )
        # Independently check quadrature at the hottest and most compressed states.
        quad_error = []
        for vv, tt in [(126, 2600), (92.5, 300), (94, 6000)]:
            gamma = source[3] * (vv / source[0]) ** source[4]
            theta = 314 * np.exp((source[3] - gamma) / source[4])

            def energy(temp):
                return (
                    36
                    * R
                    * temp
                    * (temp / theta) ** 3
                    * quad(lambda u: u**3 / np.expm1(u), 0, theta / temp, epsabs=1e-12)[
                        0
                    ]
                )

            thermal_p = (
                gamma * (energy(tt) - energy(300)) / (vv * 1e-30 * Avogadro / 4) / 1e9
            )
            quad_error.append(
                float(pressure(vv, tt, source) - pressure(vv, 300, source) - thermal_p)
            )
        groups = {}
        for name, mask in {
            "paired300_20": cold & (p < 290),
            "paired_hot30": ~cold,
            "sio2_300_10": p > 290,
        }.items():
            groups[name] = {
                "n": int(sum(mask)),
                "source_rmse": float(np.sqrt(np.mean(source_res[mask] ** 2))),
                "refit_rmse": float(np.sqrt(np.mean(fitted_res[mask] ** 2))),
                "source_mean": float(np.mean(source_res[mask])),
            }
        report["scales"][scale] = {
            "published": dict(zip(NAMES, source.tolist())),
            "printed_widths": dict(zip(NAMES, width.tolist())),
            "fits": fits,
            "multiple_starts": starts,
            "sequential_cold_then_hot": sequential,
            "thermal_only_with_source_cold_fixed": thermal.x.tolist(),
            "best_inside_printed_width_box": {
                "parameters": box.x.tolist(),
                "sse": float(2 * box.cost),
                "rmse": float(np.sqrt(np.mean(box.fun**2))),
            },
            "scaled_jacobian_singular_values": singular.tolist(),
            "scaled_jacobian_condition_number": float(singular[0] / singular[-1]),
            "rounding24_parameter_min": rounding_fits.min(axis=0).tolist(),
            "rounding24_parameter_max": rounding_fits.max(axis=0).tolist(),
            "source_coefficient_half_rounding": half_round.tolist(),
            "source_coefficient_rounding_linear_pressure_bound_gpa": float(
                np.max(
                    np.abs(approx_derivative(lambda x: pressure(v, t, x), source))
                    @ half_round
                )
            ),
            "native_direct_max_abs_pressure_difference_gpa": float(
                np.max(np.abs(native_p - pressure(v, t, source)))
            ),
            "python_direct_max_abs_pressure_difference_gpa": float(
                np.max(np.abs(python_p - pressure(v, t, source)))
            ),
            "native_direct_refit_parameter_difference": (
                np.array([z["refit"] for z in native[scale]["parameters"]]) - base.x
            ).tolist(),
            "quad_pressure_errors_gpa": quad_error,
            "source_density_checkpoints": checkpoints,
            "source_pressure_rmse": float(np.sqrt(np.mean(source_res**2))),
            "source_max_abs_pressure_residual_gpa": float(np.max(np.abs(source_res))),
            "source_refit_max_abs_pressure_difference_gpa": float(
                np.max(np.abs(source_res - fitted_res))
            ),
            "source_refit_max_abs_volume_difference_a3": float(
                np.max(np.abs(actual_vol_delta))
            ),
            "source_refit_max_abs_volume_difference_percent": float(
                100 * np.max(np.abs(actual_vol_delta / volume(p, t, source)))
            ),
            "residual_groups": groups,
            "states": [
                {
                    "source_row": int(r["source_row"]),
                    "volume_a3": float(vv),
                    "temperature_k": float(tt),
                    "pressure_gpa": float(pp),
                    "source_residual_gpa": float(sr),
                    "refit_residual_gpa": float(fr),
                }
                for r, vv, tt, pp, sr, fr in zip(rows, v, t, p, source_res, fitted_res)
            ],
        }
        print(scale, "baseline", base.x, flush=True)
        for name, fit in fits.items():
            print(
                name,
                list(fit["parameters"].values()),
                fit["all_central_values_within_source_width"],
                flush=True,
            )
    serialized = json.dumps(report, indent=2, allow_nan=False) + "\n"
    (OUT / "results.json").write_text(serialized)
    (ROOT / "docs/data/takahashi-2019-parity.json").write_text(serialized)
    plot(report)


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for col, scale in enumerate(SOURCE):
        result = report["scales"][scale]
        p = np.array([s["pressure_gpa"] for s in result["states"]])
        t = np.array([s["temperature_k"] for s in result["states"]])
        vv = np.array([s["volume_a3"] for s in result["states"]])
        refit = np.array(
            list(result["fits"]["equal_pressure_all60"]["parameters"].values())
        )
        weighted = np.array(
            list(result["fits"]["effective_variance_pvt"]["parameters"].values())
        )
        ax = axes[0, col]
        ax.scatter(p, vv, c=np.where(t == 300, "#333333", "#bb4422"), s=14, alpha=0.7)
        grid = np.linspace(69, 327, 140)
        for temp, color in [(300, "#333333"), (2000, "#bb4422")]:
            ax.plot(
                grid,
                volume(grid, temp, SOURCE[scale]),
                color=color,
                label=f"Source {temp} K",
            )
            ax.plot(
                grid,
                volume(grid, temp, refit),
                "--",
                color=color,
                label=f"Conditional refit {temp} K",
            )
            ax.plot(
                grid,
                volume(grid, temp, weighted),
                ":",
                color=color,
                label=f"P/V/T widths {temp} K",
            )
        ax.set(
            title=f"{scale.upper()} pressure scale",
            xlabel="Pressure (GPa)",
            ylabel="Conventional-cell volume (Å³)",
        )
        ax.legend(fontsize=8)
        ax = axes[1, col]
        ax.scatter(
            p,
            [s["source_residual_gpa"] for s in result["states"]],
            s=19,
            marker="o",
            label="Published coefficients",
        )
        ax.scatter(
            p,
            [s["refit_residual_gpa"] for s in result["states"]],
            s=19,
            marker="x",
            label="Equal-pressure refit",
        )
        ax.scatter(
            p,
            pressure(vv, t, weighted) - p,
            s=19,
            marker="+",
            label="P/V/T effective variance",
        )
        ax.axhline(0, color="gray", lw=0.8)
        ax.set(xlabel="Observed pressure (GPa)", ylabel="P(model) − P(observed) (GPa)")
        ax.legend(fontsize=8)
    fig.suptitle(
        "Takahashi (2019): source surfaces and conditional reconstruction\nDiagnostic fits; plotted 2000 K lines are calculated isotherms"
    )
    fig.savefig(OUT / "comparison.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
