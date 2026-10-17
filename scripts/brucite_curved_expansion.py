"""Audit the ambient-data curved expansion fit with adopted Fukui moduli.

No optimization of the published moduli or temperature slopes is performed.
The four source-isotherm variants and the continuous published-slope variant
are deliberately assessed separately. --check verifies the saved report.
"""

from __future__ import annotations

import argparse
import copy
import json

import numpy as np
import scipy

from scripts.brucite_reproduce import (
    ROOT,
    bm3,
    comparison_volume,
    conditional_continuous_pressure,
    metrics,
    observations,
    source,
)
from scripts.brucite_reproduce import (
    reference_volume as published_reference_volume,
)
from scripts.brucite_v0_audit import ambient_linear_alpha_fits

REPORT = ROOT / "docs/data/brucite-curved-expansion.json"
PREFIX = "brucite_peritheos_fukui_2003_curved_expansion"
FIT_KEY = "unweighted_volume_free_300k_anchor"


def reference_volume(temperature, fit):
    dt = np.asarray(temperature) - 300
    return fit["V300_a3"] * np.exp(
        fit["alpha_300_per_k"] * dt + 0.5 * fit["dalpha_dT_per_k2"] * dt**2
    )


def covariance_at_temperature(temperature, fit):
    """Propagate the shared 3-parameter fit to a single slice's V0."""
    dt = temperature - 300
    volume = reference_volume(temperature, fit)
    derivative = np.array([volume / fit["V300_a3"], volume * dt, 0.5 * volume * dt**2])
    covariance = np.array(fit["parameter_covariance"]["matrix"])
    return float(derivative @ covariance @ derivative)


def reproduce():
    fit = ambient_linear_alpha_fits()[FIT_KEY]
    src = source()["sources"]["fukui_2003"]
    rows = observations("fukui_2003", "table4")
    curves = []
    continuous_residuals = []
    for parameters in src["parameters"]["isotherms"]:
        t = parameters["temperature_k"]
        selected = [row for row in rows if float(row["temperature_k"]) == t]
        v = np.array([float(row["cell_volume_a3"]) for row in selected])
        p = np.array([float(row["pressure_gpa"]) for row in selected])
        v0 = float(reference_volume(t, fit))
        residual = bm3(v, v0, parameters["K0_gpa"], parameters["K0_prime"]) - p
        dt = t - 300
        k = 41.8 + src["parameters"]["dK0_dT_gpa_per_k"]["value"] * dt
        kp = 6.6 + src["parameters"]["dK0_prime_dT_per_k"]["value"] * dt
        thermal_residual = bm3(v, v0, k, kp) - p
        continuous_residuals.extend(thermal_residual.tolist())
        curves.append(
            {
                "identifier": f"{PREFIX}_{t}k_bm3_refit",
                "temperature_k": t,
                "V0_a3": v0,
                "V0_conditional_variance_a6": covariance_at_temperature(t, fit),
                "source_rows": [int(row["source_row"]) for row in selected],
                "observations": len(selected),
                "published_moduli": parameters,
                "pressure_residuals_gpa": residual.tolist(),
                **metrics(residual),
                "continuous_variant": {
                    "K0_gpa": k,
                    "K0_prime": kp,
                    "pressure_residuals_gpa": thermal_residual.tolist(),
                    **metrics(thermal_residual),
                },
            }
        )
    return {
        "audit_date": "2026-10-16",
        "primary_source": {
            "doi": src["doi"],
            "access_url": src["url"],
            "pdf_sha256": src["pdf_sha256"],
        },
        "fit": fit,
        "fit_dataset": "brucite_fukui_2003_table1",
        "validation_dataset": "brucite_fukui_2003_table4",
        "fit_objective": "Unweighted sum of squared volume residuals [V300*exp(alpha300*(T-300)+0.5*dalpha_dT*(T-300)^2)-V_observed]^2, in A^6; all 11 printed ambient Table 1 rows",
        "fixed_temperature_k": 300,
        "covariance_convention": "SSE/(11-3)*inverse(J.T@J), transformed to physical parameters V300, alpha300, dalpha_dT; published moduli/slopes treated as fixed, their uncertainty/covariance and pressure-scale systematics excluded",
        "isotherms": curves,
        "continuous_variant": {
            "identifier": f"{PREFIX}_thermal_refit",
            "observations": 37,
            "K0_300_gpa": 41.8,
            "Kprime_300": 6.6,
            "dK_dT_gpa_per_k": src["parameters"]["dK0_dT_gpa_per_k"]["value"],
            "dKprime_dT_per_k": src["parameters"]["dK0_prime_dT_per_k"]["value"],
            **metrics(continuous_residuals),
            "qualification": "Printed slopes anchored to the printed 300 K moduli, explicitly adopted as this Peritheos model's convention. This continuous version does not pass exactly through the other three Table 6 moduli. Separately fitted original-author intercepts are not recovered; use the isotherm variants to retain each tabulated modulus exactly.",
        },
        "author_workflow_reproduced": False,
        "physical_validation": "Empirical volume representation only. High-temperature V0 is extrapolated from ambient observations ending at 592.9 K. Curvature is weakly constrained; no H-order, entropy, dehydration-equilibrium, or pressure-standard validation claimed.",
    }


def ledger_outcome(record):
    report = reproduce()
    fit = report["fit"]
    thermal = "thermal" in record
    if thermal:
        result = report["continuous_variant"]
        values = {
            "V0": fit["V300_a3"],
            "alpha0": fit["alpha_300_per_k"],
            "alpha1": fit["dalpha_dT_per_k2"],
        }
        order = ["V300_a3", "alpha_300_per_k", "dalpha_dT_per_k2"]
        errors = dict(
            zip(values, [fit["parameter_standard_errors"][key] for key in order])
        )
    else:
        result = next(
            r for r in report["isotherms"] if r["identifier"] == record["identifier"]
        )
        values = {"V0": result["V0_a3"]}
        errors = {"V0": np.sqrt(result["V0_conditional_variance_a6"])}
    stored = {
        **record["eos"]["parameters"],
        **record.get("thermal", {}).get("parameters", {}),
    }
    fixed = (
        {
            "K0": 41.8,
            "K0_prime": 6.6,
            "Tr": 300.0,
            "dK_dT": result["dK_dT_gpa_per_k"],
            "dKprime_dT": result["dKprime_dT_per_k"],
        }
        if thermal
        else {
            "K0": result["published_moduli"]["K0_gpa"],
            "K0_prime": result["published_moduli"]["K0_prime"],
        }
    )
    if any(stored[name] != value for name, value in fixed.items()):
        raise ValueError("Curved-expansion adopted published modulus/slope changed")
    if (
        thermal
        and record["thermal"].get("thermal_expansion_law")
        != "linear_reference_temperature"
    ):
        raise ValueError("Curved-expansion alpha reference convention changed")
    comparisons = [
        {
            "parameter": key,
            "published": stored[key],
            "published_error": float(errors[key]),
            "refit": value,
            "refit_error": float(errors[key]),
            "difference": value - stored[key],
            "relative_difference": abs((value - stored[key]) / stored[key]),
            "within_combined_2sigma": None,
            "similar": bool(np.isclose(value, stored[key], rtol=1e-7, atol=1e-15)),
        }
        for key, value in values.items()
    ]
    if not all(p["similar"] for p in comparisons):
        raise ValueError("Registered curved-expansion brucite refit is stale")
    return {
        "status": "parity",
        "parity_basis": "registered_peritheos_refit_reproduced",
        "reproduction_status": "ambient_expansion_refit_reproduced",
        "original_publication_reproduction_status": "not_reproduced",
        "original_author_workflow_parity": False,
        "phase_pure_fit": False,
        "dataset_identifiers": [report["fit_dataset"]],
        "validation_dataset_identifiers": [report["validation_dataset"]],
        "observations": 11,
        "validation_observations": result["observations"],
        "parameters": comparisons,
        "solver_success": True,
        "fixed_source_parameters": fixed,
        "objective": report["fit_objective"],
        "rmse_volume_a3": fit["volume_rmse_a3"],
        "rmse_gpa": result["pressure_rmse_gpa"],
        "accepted_pressure_gpa": record["validity"]["pressure_gpa"],
        "reproduction_report": "docs/data/brucite-curved-expansion.json",
        "reason": "Peritheos ambient expansion refit with published moduli adopted, not refitted. The continuous variant uses explicit 300 K anchors and the published temperature slopes; individual isotherms preserve Table 5/6 moduli exactly. Parity verifies this new fit, not original-author coefficients, covariance, or workflow.",
    }


def build_records(existing_records):
    """Add four fixed-source-modulus slices and one continuous slope variant."""
    report = reproduce()
    fit = report["fit"]
    ambient_covariance = copy.deepcopy(fit["parameter_covariance"])
    result = []
    for index in range(5):
        thermal = index == 4
        assessment = (
            report["continuous_variant"] if thermal else report["isotherms"][index]
        )
        record = copy.deepcopy(
            existing_records[4] if thermal else existing_records[index]
        )
        for key in ("catalog_access", "derived_from_record", "derivation"):
            record.pop(key, None)
        temperature = None if thermal else assessment["temperature_k"]
        record.update(
            identifier=assessment["identifier"],
            label="Peritheos curved expansion + Fukui published "
            + (
                "thermal slopes; refit" if thermal else f"{temperature} K moduli; refit"
            ),
            default=False,
            record_kind="refit",
            fit_datasets=[report["fit_dataset"]],
            fixed_parameters=["K0", "K0_prime"],
            notes="Peritheos unweighted fit to all 11 ambient observations; source moduli adopted unchanged. This is the purple expansion-law variant, not an author correction or original-fit reproduction. Shared expansion covariance correlates the four isotherm variants. "
            + report["physical_validation"],
        )
        record["eos"]["parameters"]["V0"] = (
            fit["V300_a3"] if thermal else assessment["V0_a3"]
        )
        variance = (
            fit["parameter_covariance"]["matrix"][0][0]
            if thermal
            else assessment["V0_conditional_variance_a6"]
        )
        record["parameter_errors"]["V0"] = float(np.sqrt(variance))
        record["parameter_covariance"] = {
            "parameter_order": ["rt_eos.V0", "alpha0", "alpha1"] if thermal else ["V0"],
            "matrix": ambient_covariance["matrix"] if thermal else [[variance]],
        }
        record["parameter_error_confidence"] = None
        record["parameter_provenance"].update(
            V0="Peritheos unweighted free-anchor alpha(T) fit to Fukui Table 1, all 11 rows; covariance propagated to this reference temperature. Not the published Table 2/3 volume law.",
            K0="Adopted unchanged from Fukui Table 5/6 at the listed isotherm; no modulus fitted by this record.",
            K0_prime="Adopted unchanged from Fukui Table 5/6 at the listed isotherm; source error convention retained.",
            determination_method="Ambient powder-XRD volume regression with adopted experimental compression moduli from the same publication; no shock, acoustic or theoretical parameters substituted.",
            parameter_errors="Expansion errors/covariance are residual-scaled local regression estimates. Source modulus widths are retained separately; prediction covariance conditions on those moduli/slopes and excludes their unknown covariance and calibration systematics.",
        )
        record["fit_provenance"] = {
            "software": {
                "name": "Peritheos scripts/brucite_v0_audit.py and brucite_curved_expansion.py; scipy.optimize.least_squares",
                "version": scipy.__version__,
            },
            "dataset": report["fit_dataset"],
            "selection": {
                "source_rows": list(range(1, 12)),
                "exclusions": [],
                "temperature_range_k": [299.3, 592.9],
            },
            "objective": report["fit_objective"],
            "refined_parameters": ["V0_300", "alpha_300", "dalpha_dT"],
            "fixed_parameters": ["Tr", "K0", "K0_prime"]
            + (["dK_dT", "dKprime_dT"] if thermal else []),
            "covariance_scaling": report["covariance_convention"],
            "shared_expansion_parameters": {
                key: fit[key]
                for key in ("V300_a3", "alpha_300_per_k", "dalpha_dT_per_k2")
            },
            "shared_expansion_parameter_covariance": ambient_covariance,
            "statistics": {
                "observations": 11,
                "refined_parameter_count": 3,
                "degrees_of_freedom": 8,
                "volume_rmse_a3": fit["volume_rmse_a3"],
            },
            "validation_dataset": report["validation_dataset"],
            "validation_source_rows": list(range(1, 38))
            if thermal
            else assessment["source_rows"],
            "validation_statistics": {
                key: assessment[key]
                for key in (
                    "observations",
                    "pressure_rmse_gpa",
                    "maximum_absolute_pressure_residual_gpa",
                )
            },
            "reproduction": {
                "report": "docs/data/brucite-curved-expansion.json",
                "author_workflow_reproduced": False,
            },
        }
        if thermal:
            record["thermal"]["thermal_expansion_law"] = "linear_reference_temperature"
            record["thermal"]["parameters"].update(
                alpha0=fit["alpha_300_per_k"], alpha1=fit["dalpha_dT_per_k2"]
            )
            record["thermal"]["parameter_errors"].update(
                alpha0=fit["parameter_standard_errors"]["alpha_300_per_k"],
                alpha1=fit["parameter_standard_errors"]["dalpha_dT_per_k2"],
            )
            record["thermal"]["fixed_parameters"] = ["Tr", "dK_dT", "dKprime_dT"]
            record["parameter_provenance"]["thermal"] = (
                "Peritheos alpha0 is alpha(300), alpha1=dalpha/dT; both fitted only to Table 1. Integrate alpha0+alpha1*(T-300) exponentially. dK/dT=-0.0319 GPa/K and dKprime/dT=0.00571 K^-1 retained from Fukui p. 513 with Table 5/6's 300 K anchors. Separately fitted original intercepts are not recovered."
            )
            record["parameter_provenance"]["temperature_ref"] = (
                "Explicit 300 K anchor for this new continuous Peritheos model; not an assertion of recovered author regression intercepts."
            )
            record["notes"] += " " + assessment["qualification"]
            record["validity"]["notes"].append(assessment["qualification"])
        else:
            record["validity"]["notes"].append(
                "Only this listed temperature is executable. V0 comes from the shared Peritheos curved expansion law; the source isotherm's moduli are unchanged."
            )
        validation = record["scientific_validation"]
        validation.update(
            status="primary_source_validated",
            note="New Peritheos ambient-expansion refit with source moduli adopted; primary Table 1 regression and independent Table 4 pressure residuals checked. Original published records/defaults preserved. No author-workflow reproduction or recovery of an unreported author expansion law claimed.",
            verified_fields=[
                "equation",
                "primary_observations",
                "units",
                "reference_state",
                "fixed_published_moduli",
                "independent_refit",
                "conditional_covariance",
                "validity",
            ],
            numerical_reproduction={
                "ambient_refit": fit,
                "pressure_validation": assessment,
                "author_workflow_reproduced": False,
            },
            usage_recommendation="Explicit selection of this nondefault refit. Use within its represented pressure/temperature domain for volume calculations. Curvature is weakly constrained; high-temperature V0 is extrapolated. Covariance excludes source modulus/slopes and pressure-calibration uncertainty.",
        )
        validation["primary_source_check"]["finding"] = (
            "Fukui Table 1 supplies all 11 fitted ambient volumes; Tables 5/6 and p. 513 supply adopted moduli/slopes. Table 4 supplies independent pressure benchmarks. The curved expansion coefficients are Peritheos results, not published author coefficients."
        )
        result.append(record)
    return result


def write_comparison_plot(report):
    """Compare the continuous model with both original-source conventions.

    Residual statistics include every printed Table 4 row, including points
    outside the continuous catalog record's common pressure interval.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    rows = observations("fukui_2003", "table4")
    card = json.loads((ROOT / "peritheos/data/materials/brucite.eosmat").read_text())
    record = next(
        r
        for r in card["eos_records"]
        if r["identifier"] == report["continuous_variant"]["identifier"]
    )
    pressure_min, pressure_max = record["validity"]["pressure_gpa"]
    colors = {"new": "#7639a5", "published": "#c64c40", "old": "#77838b"}
    fig = plt.figure(figsize=(13.2, 11.5), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.075, right=0.98, bottom=0.17, top=0.79, hspace=0.31, wspace=0.18
    )
    residual_axes = []
    largest_residual = 0.0
    comparison = []
    for index, assessment in enumerate(report["isotherms"]):
        temperature = assessment["temperature_k"]
        selected = [r for r in rows if float(r["temperature_k"]) == temperature]
        volume = np.array([float(r["cell_volume_a3"]) for r in selected])
        pressure = np.array([float(r["pressure_gpa"]) for r in selected])
        volume_error = np.array([float(r["cell_volume_error_a3"]) for r in selected])
        pressure_error = np.array([float(r["pressure_error_gpa"]) for r in selected])
        source_moduli = assessment["published_moduli"]
        continuous = assessment["continuous_variant"]

        def curves(v):
            return {
                "new": bm3(
                    v,
                    reference_volume(temperature, report["fit"]),
                    continuous["K0_gpa"],
                    continuous["K0_prime"],
                ),
                "published": bm3(
                    v,
                    published_reference_volume(temperature),
                    source_moduli["K0_gpa"],
                    source_moduli["K0_prime"],
                ),
                "old": conditional_continuous_pressure(v, temperature),
            }

        predictions = curves(volume)
        residuals = {key: value - pressure for key, value in predictions.items()}
        parameters = {
            "new": (
                reference_volume(temperature, report["fit"]),
                continuous["K0_gpa"],
                continuous["K0_prime"],
            ),
            "published": (
                published_reference_volume(temperature),
                source_moduli["K0_gpa"],
                source_moduli["K0_prime"],
            ),
            "old": (
                published_reference_volume(temperature),
                continuous["K0_gpa"],
                continuous["K0_prime"],
            ),
        }
        volume_residuals = {}
        for key, values in parameters.items():
            predicted_volume = comparison_volume(pressure, *values)
            np.testing.assert_allclose(
                bm3(predicted_volume, *values), pressure, atol=1e-10
            )
            volume_residuals[key] = predicted_volume - volume
        np.testing.assert_allclose(
            residuals["new"], continuous["pressure_residuals_gpa"], atol=1e-12
        )
        comparison.append(
            {
                "temperature_k": temperature,
                "observations": len(selected),
                "new_continuous": metrics(residuals["new"]),
                "published_isotherm": metrics(residuals["published"]),
                "original_continuous_reconstruction": metrics(residuals["old"]),
                "volume_rmse_a3": {
                    key: float(np.sqrt(np.mean(value**2)))
                    for key, value in volume_residuals.items()
                },
            }
        )
        inner = outer[index // 2, index % 2].subgridspec(
            2, 1, height_ratios=[2.4, 1], hspace=0.06
        )
        ax = fig.add_subplot(inner[0])
        residual_ax = fig.add_subplot(inner[1], sharex=ax)
        residual_axes.append(residual_ax)
        pressure_padding = float(np.ptp(pressure)) * 0.035
        grid = np.linspace(
            pressure.min() - pressure_padding, pressure.max() + pressure_padding, 500
        )
        lines = {
            key: comparison_volume(grid, *values) for key, values in parameters.items()
        }
        for key, linestyle in (("old", "--"), ("published", "-"), ("new", "-")):
            values = lines[key]
            if key == "new":
                inside = (grid >= pressure_min) & (grid <= pressure_max)
                ax.plot(grid, values, color=colors[key], lw=2, linestyle=":")
                ax.plot(grid, np.where(inside, values, np.nan), color=colors[key], lw=2)
            else:
                ax.plot(grid, values, color=colors[key], lw=1.7, linestyle=linestyle)
            residual_ax.errorbar(
                pressure,
                volume_residuals[key],
                xerr=pressure_error,
                yerr=volume_error,
                color=colors[key],
                fmt={"old": "^", "published": "s", "new": "o"}[key],
                ms=4,
                markerfacecolor="white" if key != "new" else colors[key],
                elinewidth=0.75,
                capsize=2,
                alpha=0.9,
            )
            largest_residual = max(
                largest_residual,
                float(np.max(abs(volume_residuals[key]) + volume_error)),
            )
        ax.errorbar(
            pressure,
            volume,
            xerr=pressure_error,
            yerr=volume_error,
            fmt="o",
            color="#242424",
            ms=4.2,
            capsize=2,
            elinewidth=0.9,
            zorder=6,
        )
        ax.set_title(
            f"{temperature} K  ·  {len(selected)} observations", loc="left", fontsize=12
        )
        ax.set_ylabel("Unit-cell volume (Å³)")
        ax.tick_params(labelbottom=False)
        residual_ax.axhline(0, color="#555555", lw=0.8, zorder=0)
        residual_ax.set_xlabel("Pressure (GPa)")
        residual_ax.set_ylabel("ΔV (Å³)")
        residual_ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        pressure_padding = float(np.ptp(pressure)) * 0.035
        ax.set_xlim(
            pressure.min() - pressure_padding, pressure.max() + pressure_padding
        )
        for axis in (ax, residual_ax):
            axis.grid(alpha=0.18)
            axis.spines[["top", "right"]].set_visible(False)
    for axis in residual_axes:
        axis.set_ylim(-largest_residual * 1.07, largest_residual * 1.07)
    fig.suptitle(
        "Brucite: continuous expansion refit versus Fukui et al. (2003)",
        y=0.97,
        fontsize=17,
    )
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="o",
                color="#242424",
                linestyle="none",
                label="Published Table 4 data",
            ),
            Line2D(
                [],
                [],
                color=colors["new"],
                lw=2,
                label="New continuous refit · curved expansion + published slopes",
            ),
            Line2D(
                [],
                [],
                color=colors["published"],
                lw=1.7,
                label="Published isotherms · original expansion + tabulated moduli",
            ),
            Line2D(
                [],
                [],
                color=colors["old"],
                lw=1.7,
                linestyle="--",
                label="Original continuous reconstruction · original expansion + published slopes",
            ),
        ],
        loc="upper left",
        bbox_to_anchor=(0.075, 0.937),
        frameon=False,
        fontsize=10,
    )
    fig.text(
        0.075,
        0.063,
        "ΔV = calculated − observed volume at observed pressure; bars show reported P and V errors, excluding model uncertainty.\n"
        f"Purple dotted portions lie outside the continuous catalog interval ({pressure_min:g}–{pressure_max:g} GPa).\n"
        "Both continuous curves use published slopes anchored at the 300 K moduli; original fitted intercepts were not reported.\n"
        "Curved expansion was fitted only to ambient data (299.3–592.9 K); higher-temperature reference volumes are extrapolated.",
        fontsize=9,
        va="center",
        linespacing=1.45,
    )
    stem = ROOT / "docs/data/brucite-continuous-comparison"
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    return comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.check:
        if REPORT.read_text() != payload:
            raise SystemExit("Brucite curved-expansion report is stale")
    else:
        REPORT.write_text(payload)
    print("Brucite curved-expansion report is current.")
    if args.plot:
        print(json.dumps(write_comparison_plot(report), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
