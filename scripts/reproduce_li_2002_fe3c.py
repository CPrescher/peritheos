"""Repository-only Li (2002) Fe3C replay; no catalog or fitting-library edits.

Run: uv run python scripts/reproduce_li_2002_fe3c.py
All fits below are conditional diagnostic replays, not recovered author weights.
The implementation deliberately uses direct NumPy/SciPy, not Peritheos fitting.
"""

from __future__ import annotations

import csv
import hashlib
import json
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy import odr
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data/li-2002-fe3c"
OUTPUT = ROOT / "docs/data/li-2002-fe3c-reproduction.json"
PDF = ROOT / "output/li-2002-fe3c-audit/sources/li-2002.pdf"
SOURCE = np.array([174.0, 4.8])
WIDTHS = np.array([6.0, 0.8])
V0 = 155.28
AXIAL_SOURCE = {"a": 0.00184, "b": 0.00191, "c": 0.00200}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows():
    with (DATA / "table1-source-tokens.csv").open(newline="") as stream:
        raw = list(csv.DictReader(stream))
    normalized = []
    for row in raw:
        out = {k: row[k] for k in ("row_id", "kind", "footnote")}
        for field in ("pressure_gpa", "a_a", "b_a", "c_a", "volume_a3", "neon_a_a"):
            token_key = {
                "pressure_gpa": "pressure_error_token",
                "volume_a3": "volume_error_token",
                "neon_a_a": "neon_a_error_token",
            }.get(field, field[0] + "_error_token")
            val, err = row[field], row[token_key]
            out[field] = float(val) if val else None
            out[field + "_error_width"] = parse_width(val, err) if val else None
        out["calculated_abc_a3"] = out["a_a"] * out["b_a"] * out["c_a"]
        out["calculated_neon_cell_a3"] = (
            out["neon_a_a"] ** 3 if out["neon_a_a"] is not None else None
        )
        normalized.append(out)
    return raw, normalized


def parse_width(value, token):
    """Literal last-digit rule of Table 1 footnote b; no sigma assumption."""
    if not token:
        return None
    if token.startswith("("):
        place = Decimal(10) ** Decimal(value).as_tuple().exponent
        return float(Decimal(token[1:-1]) * place)
    return float(token)


def arrays():
    selected = [r for r in rows()[1] if r["kind"] == "measured"]
    return selected, [
        np.array([r[f] for r in selected])
        for f in (
            "pressure_gpa",
            "volume_a3",
            "pressure_gpa_error_width",
            "volume_a3_error_width",
        )
    ]


def pressure(volume, parameters=SOURCE, v0=V0):
    """Literal Eq. (2), GPa and any consistent volume units, at 300 K."""
    k, kp = parameters
    q = (v0 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k * (q**7 - q**5) * (1 + 0.75 * (kp - 4) * (q**2 - 1))


def dpdv(volume, parameters=SOURCE, v0=V0):
    k, kp = parameters
    volume = np.asarray(volume)
    q = (v0 / volume) ** (1 / 3)
    a = 0.75 * (kp - 4)
    return (
        -0.5
        * k
        * q
        / volume
        * ((7 * q**6 - 5 * q**4) * (1 + a * (q**2 - 1)) + (q**7 - q**5) * 2 * a * q)
    )


def volume(pressures, parameters=SOURCE, v0=V0):
    return np.array(
        [
            brentq(lambda v: pressure(v, parameters, v0) - p, 0.7 * v0, v0)
            for p in np.atleast_1d(pressures)
        ]
    )


def residuals(parameters, mode, p, v, sp, sv):
    if mode == "equal_pressure":
        return pressure(v, parameters) - p
    if mode == "pressure_width":
        return (pressure(v, parameters) - p) / sp
    if mode == "equal_volume":
        return volume(p, parameters) - v
    if mode == "volume_width":
        return (volume(p, parameters) - v) / sv
    if mode == "effective_variance":
        denom = np.sqrt(sp**2 + (dpdv(v, parameters) * sv) ** 2)
        return (pressure(v, parameters) - p) / denom
    raise ValueError(mode)


def curve_metrics(parameters, p, v):
    rp = pressure(v, parameters) - p
    rv = volume(p, parameters) - v
    grid = np.linspace(0, 30.5, 301)
    dv = volume(grid, parameters) - volume(grid)
    return {
        "rms_pressure_gpa": float(np.sqrt(np.mean(rp**2))),
        "max_abs_pressure_gpa": float(np.max(np.abs(rp))),
        "rms_volume_a3": float(np.sqrt(np.mean(rv**2))),
        "max_abs_volume_a3": float(np.max(np.abs(rv))),
        "max_abs_curve_volume_difference_from_source_a3": float(np.max(np.abs(dv))),
        "pressure_residuals_gpa": rp.tolist(),
        "volume_residuals_a3": rv.tolist(),
    }


def fit(mode):
    _, (p, v, sp, sv) = arrays()
    result = least_squares(
        residuals,
        SOURCE,
        args=(mode, p, v, sp, sv),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    cov = np.linalg.inv(result.jac.T @ result.jac)
    red = float(np.dot(result.fun, result.fun) / (len(p) - 2))
    widths = np.sqrt(np.diag(cov))
    return {
        "parameters": result.x.tolist(),
        "difference_from_source": (result.x - SOURCE).tolist(),
        "difference_in_printed_source_widths": ((result.x - SOURCE) / WIDTHS).tolist(),
        "both_central_values_inside_source_error_widths": bool(
            np.all(np.abs(result.x - SOURCE) <= WIDTHS)
        ),
        "conditional_absolute_covariance": cov.tolist(),
        "conditional_absolute_widths": widths.tolist(),
        "residual_scaled_widths": (widths * np.sqrt(red)).tolist(),
        "K0_K0_prime_correlation": float(cov[0, 1] / np.prod(widths)),
        "dimensionless_jacobian_condition_number": float(
            np.linalg.cond(result.jac * SOURCE)
        ),
        "sum_squared_objective_residuals": float(np.dot(result.fun, result.fun)),
        "objective_sum_squares_per_dof": red,
        "optimizer_success": bool(result.success),
        **curve_metrics(result.x, p, v),
    }


def orthogonal_fit():
    _, (p, v, sp, sv) = arrays()
    result = odr.ODR(
        odr.RealData(v, p, sx=sv, sy=sp),
        odr.Model(lambda b, x: pressure(x, b)),
        beta0=SOURCE,
        sstol=1e-12,
        partol=1e-12,
    ).run()
    return {
        "parameters": result.beta.tolist(),
        "conditional_absolute_covariance": result.cov_beta.tolist(),
        "conditional_absolute_widths": np.sqrt(np.diag(result.cov_beta)).tolist(),
        "residual_scaled_widths": result.sd_beta.tolist(),
        "objective_sum_squares_per_dof": float(result.res_var),
        "stop_reason": result.stopreason,
        **curve_metrics(result.beta, p, v),
    }


def linear_pressure_fit():
    """Independent algebra: solve K and K*(Kprime-4) by weighted linear LS."""
    _, (p, v, sp, _) = arrays()
    q = (V0 / v) ** (1 / 3)
    a = 1.5 * (q**7 - q**5)
    matrix = np.column_stack((a, a * 0.75 * (q**2 - 1)))
    k, b = np.linalg.lstsq(matrix / sp[:, None], p / sp, rcond=None)[0]
    return np.array([k, 4 + b / k])


def reference_pressure_sensitivity():
    _, (p, v, sp, _) = arrays()
    q = (V0 / v) ** (1 / 3)
    a = 1.5 * (q**7 - q**5)
    matrix = np.column_stack((a, a * 0.75 * (q * q - 1)))
    k, b = np.linalg.lstsq(matrix / sp[:, None], (p - 0.0001) / sp, rcond=None)[0]
    return {
        "P0_gpa": 0.0001,
        "shifted_pressure_width_parameters": [k, 4 + b / k],
        "change_from_literal_Eq2_fit": (
            np.array([k, 4 + b / k]) - linear_pressure_fit()
        ).tolist(),
    }


def axial_fits():
    _, norm = rows()
    selected, (p, _, sp, _) = arrays()
    anchor = next(r for r in norm if r["kind"] == "measured_anchor")
    output = {}
    for axis, beta in AXIAL_SOURCE.items():
        l0 = anchor[axis + "_a"]
        lengths = np.array([r[axis + "_a"] for r in selected])
        sl = np.array([r[axis + "_a_error_width"] for r in selected])
        fits = {}
        for mode in ("pressure_width", "effective_variance"):

            def objective(k):
                pp = pressure(lengths**3, [k[0], SOURCE[1]], l0**3)
                pl = dpdv(lengths**3, [k[0], SOURCE[1]], l0**3) * 3 * lengths**2
                denom = (
                    sp if mode == "pressure_width" else np.sqrt(sp**2 + (pl * sl) ** 2)
                )
                return (pp - p) / denom

            f = least_squares(
                objective, [1 / (3 * beta)], ftol=1e-12, xtol=1e-12, gtol=1e-12
            )
            k = f.x[0]
            width = np.sqrt(np.linalg.inv(f.jac.T @ f.jac)[0, 0]) / (3 * k**2)
            bfit = 1 / (3 * k)
            fits[mode] = {
                "beta_gpa_inverse": float(bfit),
                "conditional_absolute_width": float(width),
                "residual_scaled_width": float(
                    width * np.sqrt(np.dot(f.fun, f.fun) / 11)
                ),
                "difference_in_source_widths": float((bfit - beta) / 0.00002),
                "matches_source_central_value_at_printed_precision": bool(
                    abs(bfit - beta) <= 0.000005
                ),
            }
        output[axis] = {
            "source_beta": beta,
            "source_width": 0.00002,
            "fixed_l0_a": l0,
            "fixed_pseudo_volume_K0_prime": 4.8,
            "fits": fits,
        }
    return output


def graphical_checkpoints():
    raw = json.loads((DATA / "figure2-curve-readings.json").read_text())
    c = raw["calibration_pixels"]
    output = []
    for r in raw["readings"]:
        p = 30 * (r["x"] - c["pressure_0_x"]) / (c["pressure_30_x"] - c["pressure_0_x"])
        y = (r["y_min"] + r["y_max"]) / 2
        v = 135 + 20 * (c["volume_135_y"] - y) / (c["volume_135_y"] - c["volume_155_y"])
        calc = float(volume([p])[0])
        output.append(
            {
                "pressure_gpa": p,
                "read_volume_a3": v,
                "calculated_volume_a3": calc,
                "difference_a3": calc - v,
                "passes_reading_bound": abs(calc - v) <= raw["volume_reading_bound_a3"],
            }
        )
    return {"bound_a3": raw["volume_reading_bound_a3"], "checkpoints": output}


def reproduce():
    raw, norm = rows()
    selected, (p, v, sp, sv) = arrays()
    modes = (
        "equal_pressure",
        "pressure_width",
        "equal_volume",
        "volume_width",
        "effective_variance",
    )
    bm2 = least_squares(lambda k: (pressure(v, [k[0], 4]) - p) / sp, [180])
    return {
        "citation": "J. Li, H. K. Mao, Y. Fei, E. Gregoryanz, M. Eremets, C. S. Zha (2002), Phys Chem Minerals 29:166-169, doi:10.1007/s00269-001-0224-4",
        "scope": "Direct numerical source audit; the separate catalog mapping preserves original coefficients. This replay does not mutate the catalog or assess recommendations.",
        "verdict": {
            "volumetric_bm3_source_equation": "parity_against_five_independent_graphical_checkpoints",
            "volumetric_bm3_parameter_agreement": "conditional_parity_inside_printed_error_widths",
            "axial_bm3_parameter_agreement": "conditional_printed_precision_parity_for_all_three_axes",
            "uncertainty_reproduction": "conditional_absolute_widths_compatible_but_source_confidence_unknown",
            "author_workflow": "not_recovered",
            "calibration_replay": "unavailable_original_ruby_wavelengths_and_separate_pressure_arrays_missing",
            "catalog": "published_original_coefficients_integrated_as_cementite_li_2002_bm3_with_conditional_parameter_parity",
        },
        "source_parameters": {
            "V0_a3_fixed": V0,
            "V0_measured_error_width_a3": 0.11,
            "K0_gpa": 174,
            "K0_error_width_gpa": 6,
            "K0_prime": 4.8,
            "K0_prime_error_width": 0.8,
            "temperature_k": 300,
            "room_temperature_width_k": 1,
            "confidence_convention": None,
            "source_covariance": None,
        },
        "selection": {
            "all_table_rows": len(raw),
            "measured_rows": 13,
            "fit_row_ids": [r["row_id"] for r in selected],
            "reference_card_excluded": "T1-01 (PDF no. 35-0772)",
            "fixed_anchor": "T1-02; no fitted parameters obtain leverage from it",
            "all_12_nonambient_rows_used_in_every_volumetric_fit": True,
        },
        "normalized_rows": norm,
        "input_sha256": {
            str(x.relative_to(ROOT)): sha256(x)
            for x in sorted(DATA.glob("*"))
            if x.is_file()
        },
        "source_pdf_sha256_current": sha256(PDF) if PDF.exists() else None,
        "source_metrics": {
            **curve_metrics(SOURCE, p, v),
            "max_abs_pressure_residual_in_printed_widths": float(
                np.max(np.abs((pressure(v) - p) / sp))
            ),
            "effective_sum_squares": float(
                np.sum(residuals(SOURCE, "effective_variance", p, v, sp, sv) ** 2)
            ),
        },
        "fits": {m: fit(m) for m in modes},
        "orthogonal_distance_fit": orthogonal_fit(),
        "independent_linear_pressure_width_fit": linear_pressure_fit().tolist(),
        "one_bar_reference_sensitivity": reference_pressure_sensitivity(),
        "graphical_source_equation_check": graphical_checkpoints(),
        "axial_fits": axial_fits(),
        "bm2_diagnostic_not_a_printed_source_parameterization": {
            "K0_gpa": float(bm2.x[0]),
            "fixed_K0_prime": 4,
            **curve_metrics([bm2.x[0], 4], p, v),
        },
        "uncertainty_qualification": "Table widths are treated AS IF independent standard errors solely to calculate conditional covariance. Absolute inverse JtJ and residual-scaled widths are both retained; neither establishes source confidence, systematic/shared calibration covariance, or propagation of fixed V0 uncertainty. Equal-weight absolute widths have arbitrary unit-residual scale. Parameter-dependent effective variance uses no log-variance likelihood term. ODR adjusts both coordinates. No arbitrary row masks, scale recalibration, free V0, or replacement EOS were tried.",
        "software": {"numpy": np.__version__, "scipy": __import__("scipy").__version__},
    }


def main():
    result = reproduce()
    OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(OUTPUT.relative_to(ROOT))
    for name, f in result["fits"].items():
        print(name, f["parameters"], f["conditional_absolute_widths"])
    plot(result)


def plot(result):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, (p, v, sp, sv) = arrays()
    grid = np.linspace(0, 30.5, 250)
    fig, (ax, res) = plt.subplots(
        2, 1, figsize=(7.4, 7.4), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
    )
    ax.errorbar(
        p,
        v,
        xerr=sp,
        yerr=sv,
        fmt="o",
        ms=4,
        color="#334155",
        capsize=2,
        label="Table 1: measured (printed error widths)",
    )
    ax.plot(
        grid,
        volume(grid),
        color="#111827",
        lw=2,
        label="Li: 174 GPa, 4.8; fixed V0=155.28 A³",
    )
    styles = [
        ("pressure_width", "Pressure-width fit", "#0284c7", "--"),
        ("effective_variance", "Effective-variance fit", "#d97706", ":"),
        ("equal_pressure", "Equal-pressure fit", "#7c3aed", "-."),
    ]
    res.errorbar(p, pressure(v) - p, yerr=sp, fmt="o", ms=4, color="#334155", capsize=2)
    for key, label, color, style in styles:
        params = result["fits"][key]["parameters"]
        ax.plot(grid, volume(grid, params), color=color, ls=style, lw=1.5, label=label)
        res.plot(p, pressure(v, params) - p, color=color, ls=style, lw=1.4)
    ax.set_ylabel("Conventional-cell volume (Å³)")
    ax.set_title(
        "Li et al. (2002): Fe₃C at 300 K\nConditional refits of the original-scale printed observations",
        loc="left",
    )
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    res.axhline(0, color="0.5", lw=0.8)
    res.set_ylabel("Pmodel - Ptable (GPa)")
    res.set_xlabel("Pressure (GPa)")
    for panel in (ax, res):
        panel.grid(alpha=0.18)
        panel.spines[["top", "right"]].set_visible(False)
    fig.text(
        0.12,
        0.015,
        "Fixed V0 follows the source; weighting remains unconfirmed. No physical-validation claim.",
        fontsize=8,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    target = ROOT / "docs/literature-reproductions/figures/li-2002-fe3c-replay.png"
    target.parent.mkdir(exist_ok=True)
    fig.savefig(target, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
