"""Audit the original Hu (2019) fit and retain its full experimental shock path.

Repository maintenance: ``uv run python -m scripts.hu_2019_fe3c`` refreshes
the numerical report; ``--integrate`` appends the record and complete dataset
once. Published coefficients, all six shots and the porous precursor are kept.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from jsonschema import Draft202012Validator

from peritheos import eosmat_schema, validate_eosmat_document
from peritheos._eosmat_format import format_eosmat
from scripts.standardize_reconstruction_metadata import assessment_for_record

ROOT = Path(__file__).resolve().parents[1]
IDENTIFIER = "cementite_nonmagnetic_hu_2019_linear_hugoniot"
DATASET = "cementite_hu_2019_table1_hugoniot"
REPORT = "docs/literature-reproductions/hu-2019-fe3c.md"
NUMERICAL = "docs/data/hu-2019-fe3c.json"
TOKENS = ROOT / "peritheos/data/datasets/cementite-hu-2019-table1-tokens.csv"
CARD = ROOT / "peritheos/data/materials/cementite_nonmagnetic.eosmat"
RHO0 = 7.525
MOLAR_MASS = 179.546
V0 = 4 * MOLAR_MASS * 1e24 / (RHO0 * 6.02214076e23)
PUBLISHED = np.array([4.627, 1.614])
WIDTHS = np.array([0.073, 0.028])
P_LIMITS = [85.6, 245.0]
CAVEAT = (
    "Possible melting at higher pressures; density data do not resolve the phase state."
)
REFERENCE = {
    "authors": [
        "Hu, X.",
        "Fei, Y.",
        "Yang, J.",
        "Cai, Y.",
        "Ye, S.",
        "Qi, M.",
        "Liu, F.",
        "Zhang, M.",
    ],
    "year": 2019,
    "title": "Phase Stability and Thermal Equation of State of Iron Carbide Fe3C to 245 GPa",
    "source": "Geophysical Research Letters",
    "volume": "46",
    "issue": "20",
    "locator": "11018-11024",
    "doi": "10.1029/2019GL084545",
}
QUALIFICATION = (
    "Published six-shot linear Us-up relation retained, with conditional coefficient "
    "agreement within both reported error widths under an equal-Us-weight diagnostic. "
    "Source confidence, covariance and exact fitting weights are unspecified; this "
    "is not recovered author-workflow or uncertainty parity. The published empirical "
    "principal Hugoniot is retained over its full experimental range, 85.6-245 GPa. "
    "Possible melting at higher pressures; density data do not resolve the phase state. "
    "All six shots remain data and source-fit inputs. Mean porous "
    "precursor density 7.525(0.010) g/cm3 is explicitly adopted as the common path "
    "normalization; individual measured initial densities remain unchanged. "
    "V0 is a mass-equivalent bulk volume including approximately 2% porosity, not "
    "a diffraction lattice volume or the Sata nonmagnetic V0. The 300 K reference "
    "is conventional precursor normalization, not measured shock temperature. "
    "The separate thermal reduction and the 7.982 g/cm3 Gruneisen reference "
    "are not installed as this Hugoniot."
)
RECOMMENDATION = (
    "Retain the published experimental Hugoniot over 85.6-245 GPa for the porous precursor. "
    + CAVEAT
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_token(token):
    """Table 1 parentheses contain absolute errors, not last-digit integers."""
    value, sep, width = token.partition("(")
    return float(value), float(width.rstrip(")")) if sep else None


def observations():
    with TOKENS.open(newline="", encoding="utf-8") as f:
        raw = list(csv.DictReader(f))
    rows = []
    for r in raw:
        d = {"shot": int(r["shot"]), "impactor_code": 0 if r["impactor"] == "Cu" else 1}
        for key in (
            "initial_density_g_cm3",
            "impact_velocity_km_s",
            "shock_velocity_km_s",
            "particle_velocity_km_s",
            "density_g_cm3",
            "pressure_gpa",
            "calculated_hugoniot_temperature_k",
        ):
            d[key], d[key + "_error_width"] = parse_token(r[key])
        d["initial_density_g_cm3_error_width"] = 0.010
        d["solid_subset"] = int(d["pressure_gpa"] <= 171.6)
        rows.append(d)
    return raw, rows


def shock_velocity(up, parameters=PUBLISHED):
    return parameters[0] + parameters[1] * np.asarray(up)


def up_from_pressure(p):
    c0, slope = PUBLISHED
    return float(2 * p / (RHO0 * (np.sqrt(c0**2 + 4 * slope * p / RHO0) + c0)))


def up_domain():
    # Outward rounding avoids excluding the inclusive pressure endpoints after
    # independent inversion, without implying extra experimental precision.
    lower, upper = [up_from_pressure(p) for p in P_LIMITS]
    observed_up = [r["particle_velocity_km_s"] for r in observations()[1]]
    lower, upper = min(lower, min(observed_up)), max(upper, max(observed_up))
    return [float(np.floor(lower * 1e10) / 1e10), float(np.ceil(upper * 1e10) / 1e10)]


def reproduce():
    _, rows = observations()
    up = np.array([r["particle_velocity_km_s"] for r in rows])
    us = np.array([r["shock_velocity_km_s"] for r in rows])
    widths = np.array([r["shock_velocity_km_s_error_width"] for r in rows])
    a = np.column_stack([np.ones(len(up)), up])
    coefficients = np.linalg.lstsq(a, us, rcond=None)[0]
    residual = a @ coefficients - us
    covariance = np.linalg.inv(a.T @ a) * np.sum(residual**2) / (len(up) - 2)
    source_residual = shock_velocity(up) - us
    states = []
    for r in rows:
        rho, s, u = (
            r["initial_density_g_cm3"],
            r["shock_velocity_km_s"],
            r["particle_velocity_km_s"],
        )
        states.append(
            {
                "shot": r["shot"],
                "replayed_pressure_gpa": rho * s * u,
                "pressure_difference_gpa": rho * s * u - r["pressure_gpa"],
                "replayed_density_g_cm3": rho * s / (s - u),
                "density_difference_g_cm3": rho * s / (s - u) - r["density_g_cm3"],
            }
        )
    return {
        "record_identifier": IDENTIFIER,
        "published": {
            "c0_km_s": 4.627,
            "c0_error_width_km_s": 0.073,
            "s": 1.614,
            "s_error_width": 0.028,
            "confidence": None,
            "covariance": None,
        },
        "source_equation_check": {
            "source": "All six Table 1 Us-up states, independently deposited in PANGAEA",
            "residual_sign": "Us(model)-Us(source), km/s",
            "residuals_km_s": source_residual.tolist(),
            "rms_km_s": float(np.sqrt(np.mean(source_residual**2))),
            "max_abs_residual_over_reported_Us_width": float(
                max(abs(source_residual) / widths)
            ),
            "all_within_reported_Us_width": bool(
                np.all(abs(source_residual) <= widths)
            ),
            "qualification": "Descriptive width comparison at printed up; source confidence and Us-up covariance are unknown. Not a formal sigma or chi-square test.",
        },
        "independent_diagnostic_refit": {
            "source_shots": [r["shot"] for r in rows],
            "objective": "Equal vertical Us residuals on all six shots; source objective/weights unknown",
            "parameters_c0_s": coefficients.tolist(),
            "residual_scaled_standard_errors": np.sqrt(np.diag(covariance)).tolist(),
            "covariance": covariance.tolist(),
            "covariance_convention": "RSS/(6-2) times inverse normal matrix; conditional on printed up treated as exact. No source confidence inference or initial-state uncertainty propagation.",
            "residuals_km_s": residual.tolist(),
            "rms_km_s": float(np.sqrt(np.mean(residual**2))),
            "within_published_coefficient_widths": (
                abs(coefficients - PUBLISHED) <= WIDTHS
            ).tolist(),
            "author_workflow_reproduced": False,
        },
        "rankine_hugoniot_replay": {
            "rows": states,
            "max_abs_pressure_difference_gpa": float(
                max(abs(s["pressure_difference_gpa"]) for s in states)
            ),
            "max_abs_density_difference_g_cm3": float(
                max(abs(s["density_difference_g_cm3"]) for s in states)
            ),
            "notes": "Each shot uses its own printed initial density; calculated source P/rho are not independent measurements. Printed precision limits exact replay.",
        },
        "operational_path": {
            "precursor_density_g_cm3": RHO0,
            "density_error_width_g_cm3": 0.010,
            "bulk_mass_equivalent_V0_a3": V0,
            "formula_units": 4,
            "formula_molar_mass_g_mol": MOLAR_MASS,
            "pressure_gpa": P_LIMITS,
            "particle_velocity_km_s": up_domain(),
            "historical_lower_pressure_group_shots": [
                r["shot"] for r in rows if r["solid_subset"]
            ],
            "all_source_fit_shots": [r["shot"] for r in rows],
            "interpretation": "Full empirical experimental shock path; " + CAVEAT,
            "initial_temperature_k": 300,
            "initial_temperature_basis": "Conventional ambient normalization; not a reported numeric precursor measurement or a constant shock temperature",
        },
        "qualification": QUALIFICATION,
    }


def ledger_outcome(record):
    result = reproduce()
    refit = result["independent_diagnostic_refit"]
    if not result["source_equation_check"]["all_within_reported_Us_width"] or not all(
        refit["within_published_coefficient_widths"]
    ):
        raise AssertionError(
            "Hu source checks changed; review rather than reclassify silently"
        )
    parameters = []
    for i, name in enumerate(("c0", "s")):
        value = refit["parameters_c0_s"][i]
        parameters.append(
            {
                "parameter": name,
                "published": float(PUBLISHED[i]),
                "published_error": float(WIDTHS[i]),
                "refit": value,
                "refit_error": refit["residual_scaled_standard_errors"][i],
                "difference": value - float(PUBLISHED[i]),
                "relative_difference": abs(value / float(PUBLISHED[i]) - 1),
                "within_published_error_width": True,
                "similar": True,
                "within_combined_2sigma": None,
            }
        )
    return {
        "status": "parity",
        "parity_basis": "within_reported_parameter_errors",
        "comparison_scope": "parameter_agreement",
        "conditional": True,
        "author_workflow_reproduced": False,
        "original_author_workflow_parity": False,
        "dataset_identifiers": [DATASET],
        "fit_kind": "linear_us_up_hugoniot",
        "objective": "equal_shock_velocity_residuals",
        "weighting": "Equal Us weights; printed Up treated as exact only for this diagnostic. Both quoted velocity widths are retained in the plot, not used as statistical sigmas.",
        "selection": "All six original Table 1 shots over the full experimental Hugoniot, 85.6-245 GPa",
        "columns": {
            "particle_velocity": "particle_velocity_km_s",
            "shock_velocity": "shock_velocity_km_s",
            "particle_velocity_sigma": "particle_velocity_km_s_error_width",
            "shock_velocity_sigma": "shock_velocity_km_s_error_width",
        },
        "absolute_sigma": False,
        "free_parameters": ["c0", "s"],
        "fixed_parameters": ["V0", "rho0", "P0"],
        "parameters": parameters,
        "observations": 6,
        "degrees_of_freedom": 4,
        "rmse_shock_velocity_km_s": refit["rms_km_s"],
        "published_rmse_shock_velocity_km_s": result["source_equation_check"][
            "rms_km_s"
        ],
        "parameter_covariance": refit["covariance"],
        "covariance_convention": refit["covariance_convention"],
        "reduced_chi_square": None,
        "solver_success": True,
        "qualification": QUALIFICATION,
        "reproduction_report": REPORT,
    }


def catalog_dataset():
    raw, rows = observations()
    columns = []
    for name in rows[0]:
        if name in {"shot", "impactor_code", "solid_subset"}:
            columns.append(
                {"name": name, "quantity": "index", "unit": "1", "role": "flag"}
            )
            continue
        base = name.removesuffix("_error_width")
        quantity, unit = (
            ("density", "g/cm^3")
            if "density" in name
            else ("pressure", "GPa")
            if base == "pressure_gpa"
            else ("temperature", "K")
            if "temperature" in name
            else ("velocity", "km/s")
        )
        col = {
            "name": name,
            "quantity": quantity,
            "unit": unit,
            "role": "uncertainty" if name.endswith("_error_width") else "value",
        }
        if name.endswith("_error_width"):
            col["of"] = base
        columns.append(col)
    return {
        "identifier": DATASET,
        "kind": "shock_hugoniot",
        "description": "Complete six-shot Table 1, with per-shot porous initial densities, Cu/Ta impactors (codes 0/1), impact velocity, Us/up, derived P/rho and model-calculated temperatures. All error tokens retained. The contributor solid_subset flag preserves the former lower-pressure grouping for provenance only; it does not limit the full experimental path or identify phase states.",
        "reference": copy.deepcopy(REFERENCE),
        "source_location": "Final Table 1 and footnote; official PANGAEA Sheet1!A3:I8",
        "row_count": 6,
        "columns": columns,
        "rows": [[r[c["name"]] for c in columns] for r in rows],
        "source_tokens": raw,
        "used_by_eos_records": [IDENTIFIER],
        "uncertainty_convention": "Parentheses give explicit absolute uncertainty widths, not last-digit counts. Confidence and coordinate covariance unspecified. Initial density width 0.010 g/cm3 from Table 1 footnote; impact velocity accuracy 0.2% is not labelled a standard deviation. Missing impact-velocity error remains null.",
        "quantity_origins": {
            "initial_density_g_cm3": "Archimedes measurement, separately for each shot",
            "impact_velocity_km_s": "Measured magnetic flyer velocity; quoted accuracy 0.2%",
            "shock_velocity_km_s": "Electric-pin shock measurement",
            "particle_velocity_km_s": "Source-reduced shock coordinate; full impactor/reduction settings not recovered",
            "density_g_cm3": "Calculated Rankine-Hugoniot mass conservation",
            "pressure_gpa": "Calculated Rankine-Hugoniot momentum conservation",
            "calculated_hugoniot_temperature_k": "Source thermal-model output from Eq. 3; NOT measured temperature. Error considers fitted C0/S uncertainties.",
        },
        "provenance": {
            "source_manifest": "datasets/cementite-hu-2019-source.json",
            "source_tokens_resource": "datasets/" + TOKENS.name,
            "source_tokens_sha256": digest(TOKENS),
            "source_dataset_doi": "10.1594/PANGAEA.905474",
        },
        "license": "CC-BY-4.0",
        "source_attribution": "Fei, Yingwei (2019): Shockwave hugoniot data for Fe3C [dataset]. PANGAEA, doi:10.1594/PANGAEA.905474; supplement to Hu et al. (2019), doi:10.1029/2019GL084545.",
    }


def catalog_record():
    record = {
        "identifier": IDENTIFIER,
        "label": "Hu et al. (2019), porous-precursor experimental Hugoniot, 85.6-245 GPa",
        "model": "LinearUsUpHugoniot",
        "record_kind": "published",
        "determination_method": "experimental",
        "equation_kind": "hugoniot",
        "default_for": "hugoniot",
        "temperature_ref": 300,
        "parameters": {
            "V0": {"value": V0, "error": None, "fixed": True},
            "rho0": {"value": RHO0, "error": 0.010, "fixed": True},
            "c0": {"value": 4.627, "error": 0.073},
            "s": {"value": 1.614, "error": 0.028},
            "P0": {"value": 0, "error": None, "fixed": True},
        },
        "parameter_error_confidence": None,
        "experimental_pressure_range_gpa": P_LIMITS,
        "experimental_temperature_range_k": [300, 300],
        "pressure_range_status": "reported_exactly",
        "validity": {
            "pressure_gpa": P_LIMITS,
            "temperature_k": [300, 300],
            "notes": [
                CAVEAT,
                "Temperature bounds denote conventional precursor normalization only. Actual shock temperatures vary and are source-calculated, not measured.",
                "Approximately 2% porous precursor. Bulk mass-equivalent volumes must not be used as dense-crystal diffraction lattice volumes.",
            ],
        },
        "loading_path": "principal",
        "branch_kind": "transformed",
        "initial_state": {
            "phase": "ferromagnetic Pnma cementite, approximately 2% porous",
            "material_identifier": "cementite",
            "temperature_k": 300,
            "temperature_basis": "Conventional normalization of ambient precursor; exact numeric temperature not reported",
            "pressure_gpa": 0,
            "density_g_cm3": RHO0,
            "density_error_width_g_cm3": 0.010,
            "density_basis": "Explicitly adopted source-average Archimedes bulk density; per-shot densities retained, no porosity correction",
        },
        "volume_basis": {
            "kind": "formula_units",
            "formula_units": 4,
            "molar_mass_g_mol": MOLAR_MASS,
            "notes": "Operational four-formula mass basis including pore volume, not the Pnma crystal lattice. Standard Fe=55.845 and C=12.011 g/mol; exact Avogadro constant.",
        },
        "branch_domain": {
            "particle_velocity_km_s": up_domain(),
            "kind": "experimental_coverage",
            "boundary_status": "inferred",
            "notes": [
                "Domain encloses all six measured particle velocities and the source pressure extrema inverted with the adopted mean bulk density. Bounds are rounded outwards to 10 decimals for inclusive inversion. This is experimental path coverage, not a phase-stability boundary."
            ],
        },
        "source_fit_domain": {
            "source_shots": [1, 2, 3, 4, 5, 6],
            "pressure_gpa": [85.6, 245],
            "particle_velocity_km_s": [1.585, 3.286],
            "calculated_temperature_k": [1190, 5530],
            "installed_shots": [1, 2, 3, 4, 5, 6],
        },
        "path_interpretation": "Empirical experimental shock path; " + CAVEAT,
        "scope_history": [
            {
                "date": "2026-11-02",
                "superseded_pressure_restriction_gpa": [85.6, 171.6],
                "current_pressure_gpa": P_LIMITS,
                "reason": "The former solid-only restriction was a catalog interpretation, not a source-fit limit. Full published experimental coverage retained with the possible-melting caveat; source values and diagnostic refit unchanged.",
            }
        ],
        "fit_datasets": [DATASET],
        "reference": copy.deepcopy(REFERENCE),
        "notes": QUALIFICATION,
        "pressure_calibration": {
            "status": "not_applicable",
            "audit_date": "2026-11-01",
            "methods": [
                {
                    "kind": "shock_wave",
                    "reference": copy.deepcopy(REFERENCE),
                    "source_location": "Methods and Results; Table 1",
                    "notes": "P=rho0*Us*up, using individual shot densities for source-row replay. Impact velocity and Us measured; Cu/Ta impactors retained. Full source reduction of up from the impact configuration is not recovered. No static MgO/NaCl pressure scale is substituted.",
                }
            ],
            "recalculation": {
                "status": "not_applicable",
                "notes": "No static pressure recalibration. Algebraic per-shot Rankine-Hugoniot P/rho replay is retained separately; full author raw timing/impactor reduction unavailable.",
            },
        },
        "range_provenance": {
            "pressure": "All six Table 1 shots span 85.6-245 GPa and enter the published fit. Section 3.1.2 discusses possible melting at higher pressures; it does not impose a 171.6 GPa fit limit.",
            "temperature": "Ambient precursor normalized to 300 K; no constant-T shock path or measured 1190-5530 K PVT dataset implied.",
        },
        "evidence": {
            "parameter_provenance": {
                "determination_method": "Six experimental shock Us-up pairs and source sample bulk density; no theoretical thermal or transition-energy coefficient enters the installed path.",
                "c0": "Abstract/Results and Figure 1: 4.627 +/-0.073 km/s",
                "s": "Same original linear relation: 1.614 +/-0.028",
                "rho0": "Methods: mean sample Archimedes density 7.525(0.010) g/cm3, explicitly adopted here as common path normalization. Source Table 1 shot densities vary 7.478-7.546 and are preserved separately.",
                "V0": "Derived 4*M/(rho0*NA)*1e24 A3; four Fe3C formula units, M=179.546 g/mol. Includes pore volume; not Sata V0=149.456 A3 or measured diffraction V0.",
                "P0": "Ambient pressure approximated by zero in source shock conservation equations",
                "structure": "Source precursor XRD: Pnma Z=4; shock response is source-assigned nonmagnetic, with possible melting at higher pressures. Existing independent Prescher structure is representative solid evidence, not in-situ Hu shock diffraction or proof that the full path remains solid.",
            },
            "scientific_validation": {
                "status": "primary_source_validated",
                "audit_date": "2026-11-01",
                "report_path": REPORT,
                "note": QUALIFICATION,
                "usage_recommendation": "qualified_use",
                "recommendation_reason": RECOMMENDATION,
                "recommendation_date": "2026-11-02",
                "recommendation_audit": REPORT
                + "#full-experimental-range-and-possible-melting",
                "warning": CAVEAT,
                "primary_source_check": {
                    "access_url": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL084545",
                    "doi": REFERENCE["doi"],
                    "locations": [
                        "Methods, sample density/porosity and starting-phase XRD",
                        "Results, Us-up and conservation equations",
                        "Table 1 and footnote",
                        "Section 3.1.2, possible melting",
                        "Official PANGAEA.905474 workbook, Sheet1!A1:I8",
                    ],
                    "finding": QUALIFICATION,
                },
                "verified_fields": [
                    "original_linear_shock_coefficients_and_errors",
                    "all_six_source_shots_and_reported_errors",
                    "source_equation_against_velocity_observations",
                    "conditional_parameter_agreement",
                    "per_shot_Rankine_Hugoniot_replay",
                    "porous_initial_state_and_operational_mass_basis",
                    "full_experimental_shock_path_and_phase_caveat",
                ],
                "primary_data_check": {
                    "status": "bundled",
                    "dataset_identifiers": [DATASET],
                    "finding": "All six deposited Table 1 shots match the publisher table; all source widths, impactors, per-shot initial densities and calculated temperatures retained. All six belong to the installed experimental Hugoniot and its source-fit comparison.",
                },
                "unresolved": [
                    "Source coefficient/error confidence, covariance and exact fitting weights unspecified.",
                    "No in-situ diffraction or direct temperature measurement during shock; possible high-pressure melting is not resolved.",
                    "Exact porous-to-dense correction and full impedance/timing reduction are not supplied; none is invented.",
                    "Separate thermal model and its 7.982 g/cm3 reference, gamma law and transformation energy are not represented by this linear path.",
                ],
            },
            "reproduction": {
                "author_workflow_reproduced": False,
                "source_equation_check": reproduce()["source_equation_check"],
                "rankine_hugoniot_replay": reproduce()["rankine_hugoniot_replay"],
            },
        },
    }
    outcome = {
        "record_identifier": IDENTIFIER,
        "material_identifier": "cementite_nonmagnetic",
        "usage_recommendation": "qualified_use",
        "recommendation_reason": RECOMMENDATION,
        **ledger_outcome(record),
    }
    outcome["qualification"] += " " + RECOMMENDATION
    record["evidence"]["reproduction"]["assessment"] = assessment_for_record(
        record, outcome
    )
    return record


def integrate():
    original = json.loads(CARD.read_text())
    document = copy.deepcopy(original)
    if any(r["identifier"] == IDENTIFIER for r in document["eos_records"]):
        raise ValueError("Hu record already installed; no replacement authorized")
    document["eos_records"].append(catalog_record())
    document["datasets"].append(catalog_dataset())
    validate_eosmat_document(document)
    Draft202012Validator(eosmat_schema(4)).validate(document)
    rendered = format_eosmat(document)
    restored = json.loads(rendered)
    assert restored["eos_records"][:-1] == original["eos_records"]
    assert restored["datasets"][:-1] == original["datasets"]
    assert all(
        restored[k] == v
        for k, v in original.items()
        if k not in {"eos_records", "datasets"}
    )
    CARD.write_text(rendered)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrate", action="store_true")
    args = parser.parse_args()
    (ROOT / NUMERICAL).write_text(
        json.dumps(reproduce(), indent=2, allow_nan=False) + "\n"
    )
    if args.integrate:
        integrate()
