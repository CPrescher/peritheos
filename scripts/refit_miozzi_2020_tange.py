"""Refit Miozzi Fe observations on the Tange 2009 MgO scale.

The default regenerates evidence and an attributed derived observation table.
--register additionally writes the explicitly independent, nondefault refit into
iron.eosmat. Original source tables and published coefficients are never edited.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
from functools import lru_cache

import numpy as np
import scipy
from scipy.optimize import least_squares

from scripts.reconstruct_miozzi_2020_iron import (
    NA,
    ROOT,
    bm3,
    conditional_covariance,
    debye_energy,
    derivative,
    fe_pressure,
    metrics,
    read_data,
)

RECORD_ID = "iron_miozzi_2020_tange_2009_bm3_mgd_refit"
DATASET_ID = "iron_miozzi_2020_tange_2009_vinet_pvt"
CALIBRATION_ID = "mgo_b1_tange_2009_vinet"
REPORT = ROOT / "docs/data/miozzi-2020-tange-refit.json"
CSV_PATH = ROOT / "peritheos/data/datasets/iron-miozzi-2020-tange-2009-vinet-pvt.csv"
IRON_PATH = ROOT / "peritheos/data/materials/iron.eosmat"
NAMES = ["V0", "K0", "K0_prime", "gamma0", "q"]
SCALES = np.array([23, 100, 5, 1, 1])
LOWER = np.array([20, 20, 1, 0.1, -5])
UPPER = np.array([25, 400, 10, 5, 10])
TANGE = {
    "vinet": {
        "V0": 74.698,
        "K0": 160.63,
        "K0_prime": 4.367,
        "theta0": 761,
        "gamma0": 1.442,
        "a": 0.138,
        "b": 5.4,
    },
    "bm3": {
        "V0": 74.698,
        "K0": 160.64,
        "K0_prime": 4.221,
        "theta0": 761,
        "gamma0": 1.431,
        "a": 0.29,
        "b": 3.5,
    },
}
QUALIFICATION = (
    "Independent Peritheos equal-pressure-weight refit of all 131 Miozzi (2020) "
    "Fe observations; 116 pressures recalculated from paired MgO V,T using "
    "Tange (2009) Fit3-Vinet and 15 He pressures unchanged. This is not "
    "Miozzi's published EOS or a reproduction/correction of the original fit. "
    "V0,K0,K0_prime,gamma0,q fitted jointly; theta0=420 K, Tr=300 K, n=1 "
    "and two Fe atoms/cell fixed. Errors/covariance are conditional local "
    "standard errors using RSS/(131-5), assuming independent equal-variance "
    "pressure residuals. Calibrant systematics, shared/run correlations, "
    "predictor errors, fixed-theta uncertainty and model discrepancy are not "
    "included. Measured coverage is a marginal envelope, not a rectangular "
    "stability domain or validation for core-condition extrapolation."
)


def tange_pressure(volume, temperature, model="vinet"):
    """Independent Table 4 / Eqs. 2-5 calculation, MgO n=2 and Z=4."""
    m = TANGE[model]
    v, t = np.broadcast_arrays(volume, temperature)
    ratio = v / m["V0"]
    if model == "bm3":
        cold = bm3(v, m["V0"], m["K0"], m["K0_prime"])
    else:
        eta = ratio ** (1 / 3)
        cold = (
            3
            * m["K0"]
            * (1 - eta)
            / eta**2
            * np.exp(1.5 * (m["K0_prime"] - 1) * (1 - eta))
        )
    gamma = m["gamma0"] * (1 + m["a"] * (ratio ** m["b"] - 1))
    theta = m["theta0"] * np.exp(
        -m["gamma0"]
        * ((1 - m["a"]) * np.log(ratio) + m["a"] / m["b"] * (ratio ** m["b"] - 1))
    )
    return (
        cold
        + gamma
        * (debye_energy(theta, t, 2) - debye_energy(theta, 300, 2))
        / (v * NA / 4e24)
        / 1000
    )


def target_pressures(data, model):
    p = data["pressure_gpa"].copy()
    mask = data["mgo"]
    p[mask] = tange_pressure(
        data["mgo_volume_a3"][mask], data["temperature_k"][mask], model
    )
    return p


def target_sigma(data, parameters, model):
    v, t = data["volume_a3"], data["temperature_k"]
    fv = derivative(lambda x: fe_pressure(x, t, parameters), v, 1e-4)
    ft = derivative(lambda x: fe_pressure(v, x, parameters), t, 1e-2)
    mask = data["mgo"]
    mv, mt = data["mgo_volume_a3"][mask], t[mask]
    pv = derivative(lambda x: tange_pressure(x, mt, model), mv, 1e-4)
    pt = derivative(lambda x: tange_pressure(mv, x, model), mt, 1e-2)
    variance = (fv * data["volume_error_a3"]) ** 2
    variance[mask] += (pv * data["mgo_volume_error_a3"][mask]) ** 2
    variance[mask] += ((ft[mask] - pt) * data["temperature_error_k"][mask]) ** 2
    variance[~mask] += (
        data["pressure_error_gpa"][~mask] ** 2
        + (ft[~mask] * data["temperature_error_k"][~mask]) ** 2
    )
    return np.sqrt(variance)


def fit(
    data,
    model="vinet",
    weighting="equal",
    initial=None,
    fixed_v0=None,
    selection=None,
    fixed_k0=None,
):
    p = target_pressures(data, model)
    v, t = data["volume_a3"], data["temperature_k"]
    free = np.array([True] * 5)
    if fixed_v0 is not None:
        free[0] = False
    if fixed_k0 is not None:
        free[1] = False
    start = np.array(
        [22.81, 129, 6.24, 1.11, 0.3] if initial is None else initial, dtype=float
    )
    if fixed_v0 is not None:
        start[0] = fixed_v0
    if fixed_k0 is not None:
        start[1] = fixed_k0
    scales = SCALES[free]
    z = start[free] / scales

    def full(values):
        result = start.copy()
        result[free] = values * scales
        return result

    select = np.ones(len(p), dtype=bool) if selection is None else selection
    errors = (
        target_sigma(data, full(z), model)
        if weighting == "effective_errors"
        else np.ones_like(p)
    )
    usable = select & np.isfinite(errors) & (errors > 0)
    converged = False
    for cycle in range(100):
        if weighting == "effective_errors":
            errors = target_sigma(data, full(z), model)
        result = least_squares(
            lambda values: (
                (fe_pressure(v, t, full(values))[usable] - p[usable]) / errors[usable]
            ),
            z,
            bounds=(LOWER[free] / scales, UPPER[free] / scales),
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
            max_nfev=3000,
        )
        change = float(np.max(np.abs(result.x - z)))
        z = result.x
        if weighting == "equal" or change < 1e-7:
            converged = True
            break
    if not converged:
        raise ValueError("Effective-weight iterations did not converge")
    covariance, dof, reduced, factor, condition = conditional_covariance(
        result,
        scales,
        int(np.sum(usable)),
        "effective_errors" if weighting == "effective_errors" else "unweighted",
    )
    fitted = full(z)
    residual = fe_pressure(v, t, fitted) - p
    return {
        "parameters": dict(zip(NAMES, fitted.tolist())),
        "standard_errors": dict(
            zip(np.array(NAMES)[free], np.sqrt(np.diag(covariance)).tolist())
        ),
        "free_parameter_names": np.array(NAMES)[free].tolist(),
        "covariance": covariance.tolist(),
        "degrees_of_freedom": dof,
        "reduced_objective": reduced,
        "covariance_scale": factor,
        "jacobian_condition_number": condition,
        "solver_success": bool(result.success),
        "active_bounds": result.active_mask.tolist(),
        "weight_cycles": cycle + 1,
        "scaled_parameter_change": change,
        "fit_row_metrics": metrics(residual[usable]),
        "selection_metrics": metrics(residual[select]),
        "excluded_combined_row_numbers": (
            np.flatnonzero(select & ~usable) + 1
        ).tolist(),
    }


def derived_rows(original, data):
    p = target_pressures(data, "vinet")
    rows = []
    for i, row in enumerate(original):
        output = {
            "source_medium": row["medium"],
            "source_row": row["source_row"],
            "source_page": row["source_page"],
            "pressure_gpa": f"{p[i]:.12g}" if data["mgo"][i] else row["pressure_gpa"],
            "printed_pressure_gpa": row["pressure_gpa"],
            "printed_pressure_error_gpa": row["pressure_error_gpa"],
            **{
                name: row.get(name, "")
                for name in (
                    "temperature_k",
                    "temperature_error_k",
                    "volume_a3",
                    "volume_error_a3",
                    "mgo_volume_a3",
                    "mgo_volume_error_a3",
                )
            },
        }
        if data["mgo"][i]:
            mv, t = data["mgo_volume_a3"][i], data["temperature_k"][i]
            pv = derivative(lambda x: tange_pressure(x, t), mv, 1e-4)
            pt = derivative(lambda x: tange_pressure(mv, x), t, 1e-2)
            sigma = np.sqrt(
                (pv * data["mgo_volume_error_a3"][i]) ** 2
                + (pt * data["temperature_error_k"][i]) ** 2
            )
            output["conditional_pressure_error_gpa"] = (
                f"{sigma:.12g}" if np.isfinite(sigma) else ""
            )
            output["pressure_temperature_covariance_gpa_k"] = (
                f"{pt * data['temperature_error_k'][i] ** 2:.12g}"
                if np.isfinite(data["temperature_error_k"][i])
                else ""
            )
        else:
            output["conditional_pressure_error_gpa"] = row["pressure_error_gpa"]
            output["pressure_temperature_covariance_gpa_k"] = (
                ""  # Unavailable, not known zero.
            )
        rows.append(output)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return rows, stream.getvalue()


@lru_cache(maxsize=1)
def reproduce():
    from peritheos import get_eos_record, get_material_document

    original, data, hashes = read_data()
    calibration = next(
        r
        for r in get_material_document("mgo")["eos_records"]
        if r["identifier"] == CALIBRATION_ID
    )
    for name, value in TANGE["vinet"].items():
        stored = calibration["eos"]["parameters"].get(
            name, calibration["thermal"]["parameters"].get(name)
        )
        if value != stored:
            raise ValueError(f"Independent calibration parameter differs: {name}")
    mask = data["mgo"]
    native = get_eos_record(CALIBRATION_ID).pressure(
        data["mgo_volume_a3"][mask], data["temperature_k"][mask]
    )
    delta = native - target_pressures(data, "vinet")[mask]
    if np.max(np.abs(delta)) > 1e-9:
        raise ValueError("Independent/native Tange calibration disagreement")
    primary = fit(data)
    multistart = [
        fit(data, initial=start)
        for start in (
            [22.5, 160, 5.5, 2, 1],
            [23, 100, 7, 1.2, -0.5],
            [22.8, 130, 6.2, 3, 3],
        )
    ]
    for check in multistart:
        if not np.allclose(
            list(check["parameters"].values()),
            list(primary["parameters"].values()),
            rtol=2e-6,
            atol=2e-6,
        ):
            raise ValueError("Starting-point dependence in preferred fit")
    _, table = derived_rows(original, data)
    p = target_pressures(data, "vinet")
    return {
        "format": "peritheos.miozzi-2020-tange-refit",
        "format_version": 1,
        "record_identifier": RECORD_ID,
        "dataset_identifier": DATASET_ID,
        "calibration_record_identifier": CALIBRATION_ID,
        "calibration_parameters": TANGE,
        "calibration_record_sha256": hashlib.sha256(
            json.dumps(calibration, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "original_csv_sha256": hashes,
        "derived_csv_sha256": hashlib.sha256(table.encode()).hexdigest(),
        "fixed_fe_parameters": {
            "theta0": 420,
            "Tr": 300,
            "n": 1,
            "formula_units_per_cell": 2,
        },
        "fit_bounds": {
            name: [float(lo), float(hi)] for name, lo, hi in zip(NAMES, LOWER, UPPER)
        },
        "primary": primary,
        "multistart": multistart,
        "sensitivity": {
            "tange_bm3_calibration": fit(data, "bm3"),
            "shared_temperature_effective_errors": fit(
                data, weighting="effective_errors"
            ),
            "published_V0_fixed": fit(data, fixed_v0=22.81),
            "K0_160_gpa_fixed": fit(data, fixed_k0=160),
            "heated_only_fixed_V0": fit(
                data, fixed_v0=22.81, selection=data["temperature_k"] > 300
            ),
        },
        "observed_ranges": {
            "pressure_gpa": [float(p.min()), float(p.max())],
            "temperature_k": [
                float(data["temperature_k"].min()),
                float(data["temperature_k"].max()),
            ],
            "volume_a3": [
                float(data["volume_a3"].min()),
                float(data["volume_a3"].max()),
            ],
        },
        "pressure_recalculation_metrics": metrics((p - data["pressure_gpa"])[mask]),
        "independent_native_calibration_max_difference_gpa": float(
            np.max(np.abs(delta))
        ),
        "qualification": QUALIFICATION,
        "sensitivity_qualification": "Sensitivity fits are alternatives, not an uncertainty distribution or separately registered EOS. Effective errors include shared MgO/Fe T covariance, require complete coordinate errors, and assume printed widths are standard deviations; other coordinate/inter-row covariance and calibration uncertainty are unavailable.",
    }


def catalog_record(report):
    from peritheos import get_material_document

    original = next(
        r
        for r in get_material_document("iron")["eos_records"]
        if r["identifier"] == "iron_miozzi_2020_bm3_mgd"
    )
    calibration = next(
        r
        for r in get_material_document("mgo")["eos_records"]
        if r["identifier"] == CALIBRATION_ID
    )
    f = report["primary"]
    parameters, errors = f["parameters"], f["standard_errors"]
    reference = copy.deepcopy(original["reference"])
    document = "docs/literature-reproductions/miozzi-2020-iron.md"
    record = {
        "identifier": RECORD_ID,
        "label": "Peritheos refit of Miozzi (2020) hcp-Fe — Tange (2009) MgO pressures, equal weights",
        "record_kind": "refit",
        "determination_method": "experimental",
        "default": False,
        "reference": reference,
        "supporting_references": [copy.deepcopy(calibration["reference"])],
        "eos": {
            "type": "BM3",
            "model": "birch_murnaghan_3",
            "parameters": {name: parameters[name] for name in NAMES[:3]},
        },
        "parameter_errors": {name: errors[name] for name in NAMES[:3]},
        "parameter_error_confidence": None,
        "parameter_covariance": {
            "parameter_order": [
                "rt_eos.V0",
                "rt_eos.K0",
                "rt_eos.K0_prime",
                "gamma0",
                "q",
            ],
            "matrix": f["covariance"],
        },
        "fixed_parameters": [],
        "temperature_ref": 300,
        "reference_pressure_gpa": 0,
        "thermal": {
            "type": "MieGruneisenDebye",
            "model": "mie_gruneisen_debye",
            "debye_temperature_law": "integrated_gruneisen",
            "parameters": {
                "Tr": 300,
                "theta0": 420,
                "gamma0": parameters["gamma0"],
                "q": parameters["q"],
                "n": 1,
            },
            "parameter_errors": {
                "Tr": None,
                "theta0": None,
                "gamma0": errors["gamma0"],
                "q": errors["q"],
                "n": None,
            },
            "fixed_parameters": ["Tr", "theta0", "n"],
        },
        "experimental_pressure_range_gpa": report["observed_ranges"]["pressure_gpa"],
        "experimental_temperature_range_k": report["observed_ranges"]["temperature_k"],
        "pressure_range_status": "reference_parameterization",
        "range_provenance": {
            "notes": "Actual extrema after explicit MgO pressure recalculation; not the supplied pressure bounds or a rectangular stability domain."
        },
        "fit_datasets": [DATASET_ID],
        "notes": QUALIFICATION,
        "validity": {
            "pressure_gpa": report["observed_ranges"]["pressure_gpa"],
            "temperature_k": report["observed_ranges"]["temperature_k"],
            "notes": [
                "Marginal observation ranges only; zero-pressure coefficients extrapolated from compressed hcp Fe.",
                "No independent experimental validation or validated core extrapolation.",
            ],
        },
        "parameter_provenance": {
            "determination_method": "Independent Peritheos regression of Miozzi experimental observations with Tange recalibration; reference DOI attributes data, not authorship of these coefficients.",
            **{
                name: "Jointly fitted to all 131 recalibrated/source Fe observations; conditional equal-pressure-weight standard error and full covariance retained."
                for name in NAMES
            },
            "theta0": "420 K fixed, adopted from Miozzi Section 3.3; its uncertainty is not propagated.",
            "volume_basis": "Conventional hcp cell with two Fe atoms; n=1 per Fe formula, molar volume=V_cell*N_A/(2e24).",
            "uncertainties": QUALIFICATION,
        },
        "fit_provenance": {
            "software": {
                "name": "Peritheos independent SciPy least_squares audit",
                "version": scipy.__version__,
                "version_date": "2026-10-02",
            },
            "dataset": DATASET_ID,
            "selection": {
                "predicate": "All 15 He and 116 MgO rows; no exclusions or replacements of original coordinates",
                "included_rows": 131,
                "excluded_rows": 0,
            },
            "objective": "sum_i[(P_Fe(V_Fe_i,T_i)-P_target_i)^2] with equal pressure weights. P_target is Tange Fit3-Vinet MgO pressure for MgO rows and unchanged source pressure for He rows.",
            "refined_parameters": [
                "rt_eos.V0",
                "rt_eos.K0",
                "rt_eos.K0_prime",
                "gamma0",
                "q",
            ],
            "fixed_parameters": ["theta0", "Tr", "n"],
            "covariance_scaling": "SVD inverse of J.T J multiplied by RSS/(131-5); conditional local standard errors, not a complete uncertainty budget.",
            "statistics": {
                "observations": 131,
                "degrees_of_freedom": 126,
                "pressure_rmse_gpa": f["fit_row_metrics"]["rmse_gpa"],
                "solver_success": True,
                "jacobian_condition_number": f["jacobian_condition_number"],
            },
            "reproduction": {
                "script": "scripts/refit_miozzi_2020_tange.py",
                "report": REPORT.relative_to(ROOT).as_posix(),
                "original_csv_sha256": report["original_csv_sha256"],
                "calibration_record": CALIBRATION_ID,
                "calibration_record_sha256": report["calibration_record_sha256"],
                "derived_csv_sha256": report["derived_csv_sha256"],
            },
        },
        "pressure_calibration": {
            "status": "partially_resolved",
            "audit_date": "2026-10-02",
            "methods": [
                {
                    "kind": "equation_of_state",
                    "reference_eos_record": CALIBRATION_ID,
                    "reference": copy.deepcopy(calibration["reference"]),
                    "source_location": "Tange Table 4 Fit3-Vinet; paired MgO V,T in Miozzi supplementary table",
                    "scope": "All 116 MgO-series pressures explicitly recalculated by Peritheos; not the calibration selected in Miozzi's paper.",
                },
                *copy.deepcopy(original["pressure_calibration"]["methods"][1:]),
            ],
            "recalculation": {
                "status": "ready",
                "notes": "MgO pressures explicitly recalculated using the exact bundled Tange model; derived rows and code bundled. He source pressures retained; their raw optical readings and combination remain unavailable.",
            },
            "notes": "Analyst-selected Tange calibration, not inferred authors' settings. He provenance remains partial. No source pressure column overwritten.",
        },
        "scientific_validation": {
            "status": "primary_source_validated",
            "audit_date": "2026-10-02",
            "refit_status": "conditional",
            "reproduction_status": "independent_refit_reproduced",
            "original_publication_reproduction_status": "not_reproduced",
            "note": "Validation concerns source-observation traceability, exact Tange calculation and numerical reproduction of this separate Peritheos fit. It does not reproduce Miozzi's published fit or establish independent predictive validity. "
            + QUALIFICATION,
            "verified_fields": [
                "phase",
                "source_observations",
                "equation",
                "units",
                "reference_state",
                "fixed_parameters",
                "pressure_recalculation",
                "refit_parameters",
                "refit_covariance",
                "numerical_reproduction",
            ],
            "primary_source_check": {
                "doi": reference["doi"],
                "access_url": "https://doi.org/" + reference["doi"],
                "locations": [
                    "Miozzi official supplementary Fe-He and Fe-MgO observations",
                    "Miozzi Section 3.3 adopted theta0",
                    "Tange Table 4 Fit3-Vinet and Eqs. 2-5",
                ],
                "finding": "Data and fixed assumptions are primary-source traced; fitted coefficients and pressure calibration choice belong to the independent Peritheos refit.",
            },
            "primary_data_check": {
                "status": "bundled",
                "dataset_identifiers": [DATASET_ID],
                "source_locations": [
                    "All 131 Miozzi supplementary rows, retained in derived Tange input table"
                ],
                "finding": "All source rows preserved with per-medium/source-row identity, original P/errors, paired calibrant volumes, conditional recalculated P/errors and available P-T covariance. Equal weights use every central coordinate.",
            },
            "independent_numerical_check": {
                "script": "scripts/refit_miozzi_2020_tange.py",
                "report": REPORT.relative_to(ROOT).as_posix(),
                "native_calibration_max_difference_gpa": report[
                    "independent_native_calibration_max_difference_gpa"
                ],
            },
        },
        "reproduction": {
            "documentation": document,
            "summary_status": "independent_refit_reproduced",
            "original_publication_reproduction_status": "not_reproduced",
            "display_label": "Independent Peritheos refit — Tange MgO scale",
            "reason": QUALIFICATION,
        },
    }
    original_rows, data, _ = read_data()
    rows, table = derived_rows(original_rows, data)
    columns = []
    quantities = {
        "pressure_gpa": ("pressure", "GPa"),
        "printed_pressure_gpa": ("source_pressure", "GPa"),
        "printed_pressure_error_gpa": ("source_pressure", "GPa"),
        "conditional_pressure_error_gpa": ("pressure", "GPa"),
        "temperature_k": ("temperature", "K"),
        "temperature_error_k": ("temperature", "K"),
        "volume_a3": ("unit_cell_volume", "angstrom^3/conventional_unit_cell"),
        "volume_error_a3": ("unit_cell_volume", "angstrom^3/conventional_unit_cell"),
        "mgo_volume_a3": (
            "calibrant_unit_cell_volume",
            "angstrom^3/conventional_unit_cell",
        ),
        "mgo_volume_error_a3": (
            "calibrant_unit_cell_volume",
            "angstrom^3/conventional_unit_cell",
        ),
        "pressure_temperature_covariance_gpa_k": (
            "pressure_temperature_covariance",
            "GPa*K",
        ),
    }
    errors_of = {
        "printed_pressure_error_gpa": "printed_pressure_gpa",
        "conditional_pressure_error_gpa": "pressure_gpa",
        "temperature_error_k": "temperature_k",
        "volume_error_a3": "volume_a3",
        "mgo_volume_error_a3": "mgo_volume_a3",
    }
    for name in rows[0]:
        quantity, unit = quantities.get(name, (name, "dimensionless"))
        column = {
            "name": name,
            "quantity": quantity,
            "unit": unit,
            "role": "uncertainty"
            if name in errors_of
            else "value"
            if name in quantities
            else "flag",
        }
        if name in errors_of:
            column["of"] = errors_of[name]
        columns.append(column)
    dataset = {
        "identifier": DATASET_ID,
        "kind": "pressure_volume_temperature",
        "description": "Independent Peritheos Tange-calibrated inputs from all 131 Miozzi Fe observations; originals retained separately.",
        "reference": reference,
        "source_location": "All rows in the official Miozzi supplementary Fe-He and Fe-MgO tables; derived Tange pressures by scripts/refit_miozzi_2020_tange.py",
        "source_url": "https://www.mdpi.com/2075-163X/10/2/100/s1",
        "license": "Derived from Miozzi (2020) supplementary observations; attribution retained, article CC BY 4.0. Source PDFs are not redistributed.",
        "columns": columns,
        "resource": {
            "path": CSV_PATH.relative_to(ROOT / "peritheos/data").as_posix(),
            "sha256": report["derived_csv_sha256"],
            "media_type": "text/csv",
        },
        "used_by_eos_records": [RECORD_ID],
        "source_lineage": [
            {
                "dataset": "iron_miozzi_2020_he_pvt"
                if "-he-" in path
                else "iron_miozzi_2020_mgo_pvt",
                "resource_path": path,
                "sha256": digest,
            }
            for path, digest in report["original_csv_sha256"].items()
        ],
        "pressure_reconstruction": {
            "scope": "source_medium == mgo",
            "calibration_record": CALIBRATION_ID,
            "equation": "P_target=P_Tange(V_MgO,T); same measured T as Fe",
            "source_pressure_column": "printed_pressure_gpa",
            "target_pressure_column": "pressure_gpa",
            "unchanged_scope": "source_medium == he",
            "calibration_record_sha256": report["calibration_record_sha256"],
        },
        "uncertainty": {
            "type": "Original coordinate errors plus conditional MgO pressure propagation; source confidence unspecified",
            "notes": "Equal-weight refit does not use coordinate errors. Derived MgO pressure error omits calibration uncertainty; P-T covariance retained. He P-T covariance unknown and blank. Missing errors stay blank, original zeros remain zero.",
        },
        "notes": "The canonical pressure_gpa is the explicit fit coordinate: recalculated for MgO and unchanged for He. printed_pressure_gpa always retains the source value. This is a derived table, not an author-supplied Tange-pressure column. The typed pressure_reductions contract currently supports optical-scale transformations only; EOS recalculation provenance is retained in pressure_reconstruction and the record metadata.",
    }
    return record, dataset, table


def ledger_outcome(record):
    report = reproduce()
    f = report["primary"]
    stored = {**record["eos"]["parameters"], **record["thermal"]["parameters"]}
    errors = {**record["parameter_errors"], **record["thermal"]["parameter_errors"]}
    comparisons = []
    for name in NAMES:
        p, value = stored[name], f["parameters"][name]
        close = bool(np.isclose(p, value, rtol=2e-6, atol=2e-6))
        if not close:
            raise ValueError(f"Registered Tange refit is stale: {name}")
        comparisons.append(
            {
                "parameter": name,
                "published": p,
                "refit": value,
                "difference": value - p,
                "relative_difference": abs(value - p) / abs(p),
                "published_error": errors[name],
                "refit_error": f["standard_errors"][name],
                "within_combined_2sigma": close,
                "within_reported_error": close,
                "similar": close,
            }
        )
    return {
        "status": "parity",
        "parity_basis": "registered_peritheos_refit_reproduced",
        "fit_kind": "stored_refit_reproduction",
        "dataset_identifiers": [DATASET_ID],
        "observations": 131,
        "degrees_of_freedom": 126,
        "parameters": comparisons,
        "free_parameters": NAMES,
        "solver_success": True,
        "selection": "All 131 source observations; 116 Tange Fit3-Vinet MgO pressures plus 15 unchanged He pressures",
        "objective": "equal_weight_pressure_residuals",
        "absolute_sigma": False,
        "rmse_gpa": f["fit_row_metrics"]["rmse_gpa"],
        "qualification": QUALIFICATION,
        "original_publication_reproduction_status": "not_reproduced",
        "reproduction_status": "independent_refit_reproduced",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--register", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.register and args.check:
        parser.error("Choose --register or --check")
    report = reproduce()
    record, dataset, table = catalog_record(report)
    if args.check:
        saved = json.loads(REPORT.read_text(encoding="utf-8"))
        from scripts.check_numerical_archive import check_csv

        check_csv(
            CSV_PATH.read_text(encoding="utf-8"),
            table,
            {
                "pressure_gpa",
                "conditional_pressure_error_gpa",
                "pressure_temperature_covariance_gpa_k",
            },
        )
        if (
            saved["derived_csv_sha256"]
            != hashlib.sha256(CSV_PATH.read_bytes()).hexdigest()
        ):
            raise SystemExit("Derived Tange table is stale")
        if (
            saved["calibration_record_sha256"] != report["calibration_record_sha256"]
            or saved["original_csv_sha256"] != report["original_csv_sha256"]
        ):
            raise SystemExit("Source/calibration fingerprint is stale")
        # JSON serialization sorts report keys; compare by parameter names.
        if any(
            not np.isclose(
                saved["primary"]["parameters"][name],
                report["primary"]["parameters"][name],
                rtol=2e-6,
                atol=2e-6,
            )
            for name in NAMES
        ):
            raise SystemExit("Stored Tange parameters are stale")
        if not np.allclose(
            saved["primary"]["covariance"],
            report["primary"]["covariance"],
            rtol=1e-4,
            atol=1e-6,
        ):
            raise SystemExit("Stored Tange covariance is stale")
        iron = json.loads(IRON_PATH.read_text(encoding="utf-8"))
        registered = next(
            r for r in iron["eos_records"] if r["identifier"] == RECORD_ID
        )
        ledger_outcome(registered)
        if not np.allclose(
            registered["parameter_covariance"]["matrix"],
            report["primary"]["covariance"],
            rtol=1e-4,
            atol=1e-6,
        ):
            raise SystemExit("Registered Tange covariance is stale")
        registered_dataset = next(
            d for d in iron["datasets"] if d["identifier"] == DATASET_ID
        )
        if registered_dataset["resource"]["sha256"] != report["derived_csv_sha256"]:
            raise SystemExit("Registered Tange dataset fingerprint is stale")
        return
    REPORT.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    CSV_PATH.write_text(table, encoding="utf-8")
    if args.register:
        iron = json.loads(IRON_PATH.read_text(encoding="utf-8"))
        iron["eos_records"] = [
            r for r in iron["eos_records"] if r["identifier"] != RECORD_ID
        ] + [record]
        iron["datasets"] = [
            d for d in iron["datasets"] if d["identifier"] != DATASET_ID
        ] + [dataset]
        IRON_PATH.write_text(
            json.dumps(iron, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
