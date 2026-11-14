"""Compare recovered pyrope observations with the unchanged Angel (2022) EOS.

Repository-only: ``python -m scripts.audit_pyrope_angel_2022_observations``.
This audits three recovered tables, not the original mixed-data refinement.
Only one instrumental volume scale is profiled for the Milani shape check;
no EOS coefficient, temperature, error, selection or relative weight is fitted.
The adiabatic conversion is audit-only, not a new public caloric API.
"""

from __future__ import annotations

import csv
import hashlib
import json

import numpy as np
from scipy.constants import Avogadro, R
from scipy.integrate import quad
from scipy.optimize import brentq

from scripts.audit_garnet_angel_2022 import FOLDER, independent_pressure, source

DATA = FOLDER / "pyrope"


def observations(filename):
    with (DATA / filename).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def properties(parameters, pressure, temperature, gas_constant=R):
    """Independent BM3 derivatives and EosFit q-compromise KS convention.

    EosFit's q-compromise conversion uses gamma(V)=gamma0*V/V0, despite
    constant theta. See the official console PROPS receipt and CrysFML
    Get_Grun_V. These are the source's operational approximations, not a
    claim of a consistent quasiharmonic thermodynamic potential.
    """
    p = parameters
    v = brentq(
        lambda v: independent_pressure(p, v, temperature, gas_constant) - pressure,
        0.65 * p[1],
        1.12 * p[1],
        xtol=1e-12,
    )
    eta = (p[1] / v) ** (1 / 3)
    a = 0.75 * (p[3] - 4)
    kt = (
        0.5
        * p[2]
        * eta
        * (
            (7 * eta**6 - 5 * eta**4) * (1 + a * (eta**2 - 1))
            + (eta**7 - eta**5) * 2 * a * eta
        )
    )
    y = p[11] / temperature
    cv = (
        9
        * p[13]
        * gas_constant
        / y**3
        * quad(lambda x: x**4 * np.exp(x) / np.expm1(x) ** 2 if x else 0, 0, y)[0]
    )
    slope = p[18] * cv / (1000 * p[1])
    alpha = slope / kt
    gamma = p[18] * v / p[1]
    return {
        "volume_cm3_mol": v,
        "kt_gpa": kt,
        "alpha_per_k": alpha,
        "gamma": gamma,
        "ks_gpa": kt * (1 + alpha * gamma * temperature),
    }


def summarize(residuals):
    r = np.asarray(residuals)
    return {
        "count": len(r),
        "mean": float(r.mean()),
        "rms": float(np.sqrt(np.mean(r**2))),
        "max_absolute": float(np.max(np.abs(r))),
    }


def hartwig_audit(parameters):
    """Audit literal measured-volume column; never repair source table cells."""
    rows = observations("hartwig-galkin-2021-table-3.csv")
    t = np.array([float(row["temperature_k"]) for row in rows])
    v = np.array([float(row["measured_molar_volume_cm3_mol"]) for row in rows])
    g = np.array(
        [
            properties(parameters, 0.0001, ti)["volume_cm3_mol"] / parameters[1]
            for ti in t
        ]
    )
    scale = float(np.dot(g, v) / np.dot(g, g))
    residual = v - scale * g
    checked = []
    for row, gi, ri in zip(rows, g, residual):
        a = float(row["lattice_a_angstrom"])
        # Printed last-place rounding only, not an invented measurement error.
        a_half_step = 0.00005
        vm_token = row["measured_molar_volume_cm3_mol_token"]
        v_half_step = 0.5 * 10 ** -len(vm_token.split(".")[1])
        converted = a**3 * Avogadro / 8e24
        conversion_rounding_bound = (
            (a + a_half_step) ** 3 - a**3
        ) * Avogadro / 8e24 + v_half_step
        calculated_column = float(row["source_calculated_molar_volume_cm3_mol"])
        delta = float(row["source_delta_volume_cm3_mol"])
        delta_token = row["source_delta_volume_cm3_mol_token"].replace("\u2009", "")
        delta_half_step = 0.5 * 10 ** -len(delta_token.split(".")[1])
        calculated_half_step = 0.5 * 10 ** -len(
            row["source_calculated_molar_volume_cm3_mol_token"].split(".")[1]
        )
        difference = float(row["measured_molar_volume_cm3_mol"]) - converted
        delta_discrepancy = (
            float(row["measured_molar_volume_cm3_mol"]) - calculated_column - delta
        )
        checked.append(
            {
                "source_row": int(row["source_row"]),
                "temperature_k": float(row["temperature_k"]),
                "measured_molar_volume_cm3_mol": float(vm_token),
                "calculated_on_profiled_scale_cm3_mol": scale * gi,
                "measured_minus_eos_cm3_mol": float(ri),
                "angel_initial_assigned_temperature_error_k": 0.01
                * float(row["temperature_k"]),
                "angel_initial_assigned_volume_error_cm3_mol": 0.03,
                "volume_from_printed_lattice_a_cm3_mol": converted,
                "measured_volume_minus_lattice_volume_cm3_mol": difference,
                "lattice_volume_discrepancy_exceeds_printed_rounding": abs(difference)
                > conversion_rounding_bound,
                "source_delta_discrepancy_cm3_mol": delta_discrepancy,
                "source_delta_discrepancy_exceeds_printed_rounding": abs(
                    delta_discrepancy
                )
                > v_half_step + calculated_half_step + delta_half_step,
            }
        )
    return {
        "protocol": "All 12 literal Vmeas rows at 1 bar, fixed source EOS, one constant molar-volume scale profiled by unweighted least squares (equivalent to the common initial 0.03 cm3/mol error). Neither T nor error multipliers fitted. Not the original 2022 normalization, effective-variance objective or increased weighting.",
        "profiled_molar_volume_at_298k_cm3_mol": scale,
        "residual_summary_cm3_mol": summarize(residual),
        "initial_assigned_volume_error_only_standardized_rms": float(
            np.sqrt(np.mean((residual / 0.03) ** 2))
        ),
        "temperature_range_k": [float(t.min()), float(t.max())],
        "source_checks": "Compare a^3 NA/Z using exact NA and Z=8; flag differences beyond printed-number rounding only. Compare Vmeas-Vcalc with the literal deltaV column. Flags establish internal source-table discrepancies under these conventions, not corrected measurements or a measurement-error/confidence estimate.",
        "qualification": "Publisher measured-volume column recovered; exact column/corrections used by Angel 2022 and increased Hartwig error/weight multiplier remain unavailable. Do not substitute lattice-derived volumes or calculated/interpolated columns, or infer the multiplier from apparent curve agreement.",
        "rows": checked,
    }


def audit():
    p, covariance = source("pyrope")
    milani = observations("milani-2015-table-a1.csv")
    t = np.array([float(row["temperature_k"]) for row in milani])
    v = np.array([float(row["unit_cell_volume_a3"]) for row in milani])
    sd = np.array([float(row["unit_cell_volume_sd_a3"]) for row in milani])
    g = np.array([properties(p, 0.0001, ti)["volume_cm3_mol"] / p[1] for ti in t])
    # Conditional diagonal volume-only objective; original normalization unknown.
    scale = float(np.sum(g * v / sd**2) / np.sum(g**2 / sd**2))
    vr = v - scale * g
    first, last = list(t).index(303), list(t).index(1073)
    observed_ratio = v[last] / v[first]
    ratio_error = observed_ratio * np.sqrt(
        (sd[last] / v[last]) ** 2 + (sd[first] / v[first]) ** 2
    )
    milani_rows = [
        {
            "source_row": int(row["source_row"]),
            "temperature_k": ti,
            "observed_cell_volume_a3": vi,
            "reported_volume_sd_a3": si,
            "calculated_volume_on_profiled_scale_a3": scale * gi,
            "observed_minus_calculated_a3": ri,
        }
        for row, ti, vi, si, gi, ri in zip(milani, t, v, sd, g, vr)
    ]
    zou_rows = []
    for row in observations("zou-2012b-table-i.csv"):
        ti, pressure = float(row["temperature_k"]), float(row["pressure_gpa"])
        ks = properties(p, pressure, ti)["ks_gpa"]
        dp, dt = 0.001, 0.01
        dkp = (
            properties(p, pressure + dp, ti)["ks_gpa"]
            - properties(p, pressure - dp, ti)["ks_gpa"]
        ) / (2 * dp)
        dkt = (
            properties(p, pressure, ti + dt)["ks_gpa"]
            - properties(p, pressure, ti - dt)["ks_gpa"]
        ) / (2 * dt)
        assigned_t_error = 0.02 * ti  # Angel (2022), supplementary S2 assignment.
        effective = float(
            np.sqrt(
                float(row["ks_gpa_error"]) ** 2
                + (dkp * float(row["pressure_gpa_error"])) ** 2
                + (dkt * assigned_t_error) ** 2
            )
        )
        zr = {
            "source_row": int(row["source_row"]),
            "temperature_k": ti,
            "pressure_gpa": pressure,
            "observed_ks_gpa": float(row["ks_gpa"]),
            "reported_ks_error_gpa": float(row["ks_gpa_error"]),
            "calculated_ks_gpa": ks,
            "observed_minus_calculated_ks_gpa": float(row["ks_gpa"]) - ks,
            "angel_assigned_temperature_error_k": assigned_t_error,
            "conditional_diagonal_effective_ks_error_gpa": effective,
            "volume_used_by_angel_2022": ti == 300,
        }
        if ti == 300:
            # S2 describes recovering V0 from V/(V/V0). Use the published ratio
            # directly, avoiding a falsely precise single reconstructed V0.
            normalized_v = p[1] * float(row["relative_volume"])
            zr["ratio_normalized_volume_cm3_mol"] = normalized_v
            zr["observed_minus_calculated_pressure_gpa"] = (
                pressure - independent_pressure(p, normalized_v, ti)
            )
        zou_rows.append(zr)
    receipt = json.loads((DATA / "eosfit-properties/checks.json").read_text())
    calculator_rows = []
    for row in receipt["states"]:
        calculated = properties(
            p, row["pressure_gpa"], row["temperature_k"], receipt["gas_constant"]
        )
        differences = {key: calculated[key] - row[key] for key in calculated}
        for key, difference in differences.items():
            if abs(difference) > receipt["absolute_tolerances"][key]:
                raise ValueError(
                    f"Official calculator discrepancy in {key}: {difference}"
                )
        calculator_rows.append({**row, "independent_minus_eosfit": differences})
    result = {
        "scope": "Partial observation audit; unchanged published EOS; no joint refit",
        "author_workflow_reproduced": False,
        "source_eos_sha256": hashlib.sha256(
            (FOLDER / "sources/pyrope-official.eos").read_bytes()
        ).hexdigest(),
        "data_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in (
                "milani-2015-table-a1.csv",
                "zou-2012b-table-i.csv",
                "hartwig-galkin-2021-table-3.csv",
            )
        },
        "sinogeikin_2000_scale": {
            "source_parameter_index_one_based": 53,
            "calculated_modulus_multiplier_to_observed_scale": p[53],
            "inverse_observed_multiplier_to_eos_scale": 1 / p[53],
            "covariance_marginal_width": float(np.sqrt(covariance[52, 52])),
            "article_rounded_value": 0.978,
            "supplement_S2_rounded_value": 0.977,
            "confidence": None,
        },
        "milani_expansion": {
            "protocol": "All 33 rows; 1 bar; source EOS coefficients fixed; one constant multiplicative volume calibration fitted by reported diagonal volume standard deviations only. No assigned T errors, anchor uncertainty or inter-row covariance. This profiles shape agreement, not the author's molar-volume normalization or joint pressure-residual objective.",
            "profiled_cell_volume_at_298k_on_observed_scale_a3": scale,
            "residual_summary_a3": summarize(vr),
            "volume_only_standardized_residual_rms": float(
                np.sqrt(np.mean((vr / sd) ** 2))
            ),
            "scale_independent_303_to_1073k_volume_ratio": {
                "observed": observed_ratio,
                "calculated": g[last] / g[first],
                "reported_volume_only_diagonal_error": ratio_error,
                "qualification": "A constant volume scale cancels. Error assumes independent reported volume standard deviations; T error, common calibration covariance and EOS parameter uncertainty are omitted. Not a confidence statement for the published joint fit.",
            },
            "rows": milani_rows,
        },
        "zou_moduli": {
            "protocol": "All 40 published KS evaluated at each reported P,T; no scale fitted. KS=KT*(1+alpha*gamma*T), gamma=gamma0*V/V0; fixed-theta Debye pressure. Volumes at high T are excluded from the 2022 fit and are never used to predict KS here.",
            "effective_errors": "Diagnostic propagation into KS: reported KS and P errors plus the 2 percent T errors assigned in Angel S2, treated as independent. No P-KS covariance recovered. Not the original EosFit pressure-residual objective, confidence statement, or fit statistic.",
            "residual_summary_gpa": summarize(
                [row["observed_minus_calculated_ks_gpa"] for row in zou_rows]
            ),
            "conditional_effective_error_standardized_rms": float(
                np.sqrt(
                    np.mean(
                        [
                            (
                                row["observed_minus_calculated_ks_gpa"]
                                / row["conditional_diagonal_effective_ks_error_gpa"]
                            )
                            ** 2
                            for row in zou_rows
                        ]
                    )
                )
            ),
            "room_temperature_ratio_normalized_pressure_residual_summary_gpa": summarize(
                [
                    row["observed_minus_calculated_pressure_gpa"]
                    for row in zou_rows
                    if row["volume_used_by_angel_2022"]
                ]
            ),
            "rows": zou_rows,
        },
        "hartwig_expansion": hartwig_audit(p),
        "official_calculator_checks": calculator_rows,
        "qualification": "Zou subset agrees at the scale of reported/assigned errors. Milani expansion shape has measurable residual tension on its reported volume precision alone. Hartwig's 12 measured-volume rows are recovered but its source table has internal discrepancies and the 2022 increased weight is unrecovered. These comparisons do not establish the complete joint-fit statistic or author-workflow parity. Original normalizations, source-used Hartwig values, weighting and other unrecovered inputs remain necessary for a source-protocol replay. No recommendation/default or catalog classification changed by this audit.",
    }
    return result


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from scripts.report_plot_themes import dark_svg

    with plt.rc_context({"svg.hashsalt": "pyrope-observation-audit", "font.size": 10}):
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
        rows = report["milani_expansion"]["rows"]
        axes[0].errorbar(
            [row["temperature_k"] for row in rows],
            [row["observed_minus_calculated_a3"] for row in rows],
            yerr=[row["reported_volume_sd_a3"] for row in rows],
            fmt="o",
            markersize=3,
            capsize=2,
            color="#2463a2",
        )
        axes[0].set_title("Milani (2015): 33 expansion volumes")
        axes[0].set_xlabel("Temperature (K)")
        axes[0].set_ylabel("Observed − calculated V (Å³/cell)")
        rows = report["zou_moduli"]["rows"]
        for row in rows:
            axes[1].plot(
                [row["temperature_k"]] * 2,
                [
                    row["observed_minus_calculated_ks_gpa"]
                    - row["reported_ks_error_gpa"],
                    row["observed_minus_calculated_ks_gpa"]
                    + row["reported_ks_error_gpa"],
                ],
                color="#aaaaaa",
                linewidth=0.7,
            )
        dots = axes[1].scatter(
            [row["temperature_k"] for row in rows],
            [row["observed_minus_calculated_ks_gpa"] for row in rows],
            c=[row["pressure_gpa"] for row in rows],
            cmap="viridis",
            s=22,
            zorder=3,
        )
        fig.colorbar(dots, ax=axes[1], label="Pressure (GPa)", pad=0.03)
        axes[1].set_title("Zou (2012b): 40 adiabatic moduli")
        axes[1].set_xlabel("Temperature (K)")
        axes[1].set_ylabel("Observed − calculated KS (GPa)")
        for ax in axes:
            ax.axhline(0, color="#777777", lw=0.7)
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(alpha=0.15)
        fig.text(
            0.5,
            0.015,
            "Published EOS fixed; one Milani volume scale profiled. Bars: reported V or KS errors only. No joint refit.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        path = DATA / "observation-residuals.svg"
        fig.savefig(path, metadata={"Date": None})
        (DATA / "observation-residuals-dark.svg").write_bytes(
            dark_svg(path.read_bytes())
        )
        plt.close(fig)


if __name__ == "__main__":
    report = audit()
    (DATA / "audit.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    plot(report)
