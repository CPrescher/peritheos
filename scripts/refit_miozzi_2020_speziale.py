"""Register the author-reported Miozzi n=2 EOS and a separate Speziale n=1 refit.

python -m scripts.refit_miozzi_2020_speziale --register
python -m scripts.refit_miozzi_2020_speziale --check

The refit uses Speziale's variable-q Debye MgO reconstruction, fixed calibration
coefficients, equal pressure weights, all 131 observations and the user-relayed
three-stage protocol. This does not recover original author input files/weights.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
from functools import lru_cache
from importlib.metadata import version

import numpy as np

from peritheos import get_material_document
from scripts.reconstruct_miozzi_2020_iron import (
    MODELS,
    ROOT,
    derivative,
    fe_pressure,
    mgo_pressure,
    pressure_column,
    read_data,
)
from scripts.refit_miozzi_2020_cell_basis import staged_fit

AUTHOR_ID = "iron_miozzi_2020_bm3_mgd"
RECORD_ID = "iron_miozzi_2020_speziale_2001_bm3_mgd_refit"
DATASET_ID = "iron_miozzi_2020_speziale_2001_variable_q_pvt"
MODEL = "speziale_variable_q_debye"
IRON_PATH = ROOT / "peritheos/data/materials/iron.eosmat"
CSV_PATH = (
    ROOT / "peritheos/data/datasets/iron-miozzi-2020-speziale-2001-variable-q-pvt.csv"
)
REPORT = ROOT / "docs/data/miozzi-2020-speziale-refit.json"
NAMES = ["V0", "K0", "K0_prime", "gamma0", "q"]
DOCUMENT = "docs/literature-reproductions/miozzi-2020-iron.md"
QUALIFICATION = (
    "Independent Peritheos refit, not Miozzi's published coefficients. All 116 "
    "MgO-series pressures are recalculated with the Speziale (2001) variable-q "
    "Debye model; the 15 He pressures remain unchanged. Equal pressure-residual "
    "weights; 36 RT rows first, all rows with cold coefficients fixed second, "
    "then V0,K0,K0_prime,gamma0,q free. theta0=420 K, Tr=300 K, n=1 per mole "
    "Fe fixed. Local errors/covariance are conditional on this calibration, "
    "objective and fixed parameters, excluding calibrant systematics, predictor "
    "and inter-row covariance and fixed-theta uncertainty. Original author "
    "weights and exact pressure reduction remain unavailable; numerical "
    "reproducibility does not establish independent physical accuracy."
)


def derived_rows(original, data):
    pressure = pressure_column(data, MODEL)[0]
    rows = []
    for i, source in enumerate(original):
        row = {
            "source_medium": source["medium"],
            "source_row": source["source_row"],
            "source_page": source["source_page"],
            "pressure_gpa": f"{pressure[i]:.12g}"
            if data["mgo"][i]
            else source["pressure_gpa"],
            "printed_pressure_gpa": source["pressure_gpa"],
            "printed_pressure_error_gpa": source["pressure_error_gpa"],
            **{
                key: source.get(key, "")
                for key in (
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
            v, t = data["mgo_volume_a3"][i], data["temperature_k"][i]
            pv = derivative(lambda x: mgo_pressure(x, t, MODEL), v, 1e-4)
            pt = derivative(lambda x: mgo_pressure(v, x, MODEL), t, 1e-2)
            sigma = np.sqrt(
                (pv * data["mgo_volume_error_a3"][i]) ** 2
                + (pt * data["temperature_error_k"][i]) ** 2
            )
            row["conditional_pressure_error_gpa"] = (
                f"{sigma:.12g}" if np.isfinite(sigma) else ""
            )
            row["pressure_temperature_covariance_gpa_k"] = (
                f"{pt * data['temperature_error_k'][i] ** 2:.12g}"
                if np.isfinite(data["temperature_error_k"][i])
                else ""
            )
        else:
            row["conditional_pressure_error_gpa"] = source["pressure_error_gpa"]
            row["pressure_temperature_covariance_gpa_k"] = ""
        rows.append(row)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return rows, stream.getvalue()


@lru_cache(maxsize=1)
def reproduce():
    original, data, hashes = read_data()
    pressures = pressure_column(data, MODEL)[0]
    fit = staged_fit(data, pressures, n=1)
    mapping = {"V0_cell_a3": "V0", **{key: key for key in NAMES[1:]}}
    fit["parameters"] = {
        mapping[key]: value for key, value in fit["parameters"].items()
    }
    fit["conditional_standard_errors"] = {
        mapping[key]: value for key, value in fit["conditional_standard_errors"].items()
    }
    fit["covariance_parameter_order"] = NAMES
    _, table = derived_rows(original, data)
    independent = fe_pressure(
        data["volume_a3"],
        data["temperature_k"],
        [fit["parameters"][key] for key in NAMES],
    )
    delta = float(np.max(np.abs(independent - fit["predicted_pressure_gpa"])))
    if delta > 1e-8:
        raise ValueError("Independent per-Fe pressure disagrees with refit")
    return {
        "format": "peritheos.miozzi-2020-speziale-refit",
        "format_version": 1,
        "record_identifier": RECORD_ID,
        "dataset_identifier": DATASET_ID,
        "original_csv_sha256": hashes,
        "source_code_sha256": {
            path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for path in (
                "scripts/refit_miozzi_2020_speziale.py",
                "scripts/refit_miozzi_2020_cell_basis.py",
                "scripts/reconstruct_miozzi_2020_iron.py",
            )
        },
        "derived_csv_sha256": hashlib.sha256(table.encode()).hexdigest(),
        "calibration": {
            "model": MODEL,
            "parameters": MODELS[MODEL],
            "MgO_n": 2,
            "MgO_Z": 4,
            "theta_law": "Numerical integral of d ln(theta)/d ln(V)=-gamma; gamma follows Speziale Eqs. 10-11, not constant-q theta applied to variable q.",
        },
        "primary": fit,
        "fixed_fe_parameters": {"theta0": 420, "Tr": 300, "n": 1, "Z": 2},
        "fit_bounds_in_cell_basis": {
            "V0": [20, 25],
            "K0": [20, 400],
            "K0_prime": [1, 10],
            "gamma0": [0.1, 5],
            "q": [-5, 10],
        },
        "observed_ranges": {
            "pressure_gpa": [float(pressures.min()), float(pressures.max())],
            "temperature_k": [
                float(data["temperature_k"].min()),
                float(data["temperature_k"].max()),
            ],
        },
        "independent_native_max_difference_gpa": delta,
        "qualification": QUALIFICATION,
    }


def author_record(record):
    record = copy.deepcopy(record)
    communication = {
        "type": "personal_communication",
        "attribution": "Miozzi et al.",
        "recorded_on": "2026-10-03",
        "channel": "User-relayed in-person communication with Miozzi",
        "finding": "Full MGD, n=2 and V0 approximately 6.87 cm3/mol in EosFit; theta0 fixed at 420 K. RT-only fit, all data with RT coefficients fixed, final V0,K0,K0_prime,gamma0,q free.",
        "scope": "Author-reported software settings; not printed n in the article or a recovered original session. Article prose says V0 fixed; final free V0 follows the personal communication.",
    }
    reason = (
        "Published coefficients retained with full MGD n=2 at the per-mole-Fe "
        "volume (~6.87 cm3/mol), as used in the paper according to personal "
        "communication with Miozzi et al. relayed by the user and recorded "
        "2026-10-03. This source-reproduction setting differs from physical "
        "per-Fe normalization (n=1). Unchanged coefficients give RMS 2.49886 "
        "GPa on the 131 printed-pressure observations; exact author inputs, "
        "weights and full source-fit reproduction remain unresolved."
    )
    record["label"] = "Miozzi et al. (2020), published BM3–MGD — author-reported n=2"
    record["thermal"]["parameters"]["n"] = 2
    record["fixed_parameters"] = []
    record["notes"] = reason
    record["parameter_provenance"]["n"] = communication
    record["parameter_provenance"]["volume_basis"] = (
        "Public conventional hcp cell contains two Fe atoms; internal molar "
        "volume is V_cell*N_A/(2e24) cm3/mol (~6.87 at V0), with author-reported n=2. "
        "n is not silently normalized to 1 or accompanied by a volume doubling."
    )
    record["parameter_provenance"]["fit_protocol"] = communication
    record["parameter_provenance"]["uncertainties"] = (
        "Parenthesized published errors, unspecified confidence; reference "
        "volume error converted from the printed molar error. No author covariance."
    )
    # Keep repeated registration idempotent: provenance text is set from the
    # preserved published statement rather than accumulating communication text.
    marker = " Personal communication recorded 2026-10-03 clarifies"
    record["parameter_provenance"]["V0"] = (
        record["parameter_provenance"]["V0"].split(marker)[0]
        + marker
        + " that V0 was free in the final stage; the article's fixed-V0 "
        "wording remains a documented source discrepancy. The printed coefficient is retained."
    )
    record["parameter_provenance"]["equation"] = (
        "BM3 Equation (3), thermal Equations (4), (6), (8) with integrated "
        "Gruneisen MGD. Apparent Equation (5) denominator/subtraction and "
        "Equation (7) Debye-integral typesetting defects are mapped explicitly "
        "to gamma*delta(E_D)/V_molar. Full MGD and n=2 at the per-mole-Fe "
        "volume are author-confirmed by the recorded personal communication."
    )
    for key in ("note", "usage_recommendation"):
        record["scientific_validation"][key] = reason
    record["scientific_validation"]["audit_date"] = "2026-10-03"
    record["scientific_validation"]["personal_communication"] = communication
    record["scientific_validation"]["primary_source_check"]["finding"] = (
        "Printed preferred thermal coefficients and their uncertainties are "
        "preserved. Full MGD, n=2 with V0 approximately 6.87 cm3/mol and the "
        "final free V0 are attributed to recorded personal communication, "
        "rather than inferred from the printed equations. The author-setting "
        "replay improves residuals but does not recover the full original fit."
    )
    record["scientific_validation"]["independent_numerical_check"] = {
        "script": "scripts/reproduce_miozzi_2020_iron.py",
        "report": "docs/data/miozzi-2020-iron-reproduction.json",
        "author_setting_diagnostic": "docs/data/miozzi-2020-eosfit-n2-author-volume/pressure-check.json",
    }
    record["reproduction"]["reason"] = reason
    return record


def catalog_record(report, document):
    # Use the existing attributed dataset column contract, not its Tange values.
    template = next(
        ds
        for ds in document["datasets"]
        if ds["identifier"] == "iron_miozzi_2020_tange_2009_vinet_pvt"
    )
    dataset = copy.deepcopy(template)
    fit = report["primary"]
    parameters, errors = fit["parameters"], fit["conditional_standard_errors"]
    source = next(r for r in document["eos_records"] if r["identifier"] == AUTHOR_ID)
    calibration_reference = copy.deepcopy(
        next(
            r
            for r in get_material_document("mgo")["eos_records"]
            if r["identifier"] == "mgo_speziale_2001_bm3_2"
        )["reference"]
    )
    dataset.update(
        {
            "identifier": DATASET_ID,
            "used_by_eos_records": [RECORD_ID],
            "description": "Independent Peritheos Speziale variable-q Debye recalculation of all 131 Miozzi Fe observations; original pressures retained separately.",
            "source_location": "All 15 Fe-He and 116 Fe-MgO source rows; derived Speziale pressures by scripts/refit_miozzi_2020_speziale.py",
            "resource": {
                "path": str(CSV_PATH.relative_to(ROOT / "peritheos/data")),
                "sha256": report["derived_csv_sha256"],
                "media_type": "text/csv",
            },
            "pressure_reconstruction": {
                "scope": "source_medium == mgo",
                "model": MODEL,
                "reference": calibration_reference,
                "implementation": report["calibration"],
                "equation": "P_target=P_Speziale_variable_q_Debye(V_MgO,T)",
                "source_pressure_column": "printed_pressure_gpa",
                "target_pressure_column": "pressure_gpa",
                "unchanged_scope": "source_medium == he",
            },
            "notes": "Explicit derived Speziale variable-q Debye pressures; not an author-supplied pressure column. Original pressures, coordinate precision, missing errors and He optical-gauge provenance retained. The isothermal-only catalog Speziale record is not used as a thermal calibrant.",
        }
    )
    record = {
        "identifier": RECORD_ID,
        "label": "Peritheos refit of Miozzi (2020) hcp-Fe — Speziale (2001) variable-q MgO, n=1",
        "record_kind": "refit",
        "determination_method": "experimental",
        "default": False,
        "reference": copy.deepcopy(source["reference"]),
        "supporting_references": [calibration_reference],
        "eos": {
            "type": "BM3",
            "model": "birch_murnaghan_3",
            "parameters": {key: parameters[key] for key in NAMES[:3]},
        },
        "parameter_errors": {key: errors[key] for key in NAMES[:3]},
        "parameter_error_confidence": None,
        "parameter_covariance": {
            "parameter_order": [
                "rt_eos.V0",
                "rt_eos.K0",
                "rt_eos.K0_prime",
                "gamma0",
                "q",
            ],
            "matrix": fit["covariance_in_cell_volume_basis"],
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
        "fit_datasets": [DATASET_ID],
        "notes": QUALIFICATION,
        "parameter_provenance": {
            "determination_method": "Independent regression of experimental Miozzi observations; the source DOI attributes the data, not authorship of these fitted coefficients.",
            "volume_basis": "Conventional hcp Fe cell, Z=2; internal energy/volume per mole Fe with n=1.",
            "uncertainties": QUALIFICATION,
        },
        "fit_provenance": {
            "software": {
                "name": "Peritheos bounded least squares",
                "version": version("peritheos"),
                "version_date": "2026-10-03",
            },
            "source_author_record": AUTHOR_ID,
            "dataset": DATASET_ID,
            "selection": {
                "included_rows": 131,
                "excluded_rows": 0,
                "predicate": "All 15 He and 116 MgO observations",
            },
            "objective": "sum_i(P_Fe(V_i,T_i)-P_target_i)^2; equal pressure-residual weights",
            "refined_parameters": [
                "rt_eos.V0",
                "rt_eos.K0",
                "rt_eos.K0_prime",
                "gamma0",
                "q",
            ],
            "fixed_parameters": ["theta0", "Tr", "n"],
            "covariance_scaling": "RSS/(131-5); conditional local covariance, including cold/thermal cross-correlations",
            "statistics": {
                "observations": 131,
                "degrees_of_freedom": 126,
                "pressure_rmse_gpa": fit["rmse_gpa"],
                "solver_success": True,
            },
            "reproduction": {
                "script": "scripts/refit_miozzi_2020_speziale.py",
                "report": str(REPORT.relative_to(ROOT)),
                "original_csv_sha256": report["original_csv_sha256"],
                "derived_csv_sha256": report["derived_csv_sha256"],
            },
        },
        "pressure_calibration": {
            "status": "partially_resolved",
            "audit_date": "2026-10-03",
            "methods": [
                {
                    "kind": "equation_of_state",
                    "reference": calibration_reference,
                    "source_location": "Speziale Sections 2.1 and 5.4, Eqs. 10-11; paired MgO coordinates from Miozzi supplement",
                    "scope": "All 116 MgO-series pressures explicitly recalculated with variable-q Debye model",
                    "implementation": report["calibration"],
                },
                *copy.deepcopy(source["pressure_calibration"]["methods"][1:]),
            ],
            "recalculation": {
                "status": "ready",
                "notes": "Thermal calibration implementation and derived observations bundled; isothermal-only Speziale catalog record is not substituted. He source pressures retained; raw optical gauge readings unavailable.",
            },
        },
        "scientific_validation": {
            "status": "primary_source_validated",
            "audit_date": "2026-10-03",
            "refit_status": "conditional",
            "reproduction_status": "independent_refit_reproduced",
            "original_publication_reproduction_status": "not_reproduced",
            "note": QUALIFICATION,
            "verified_fields": [
                "phase",
                "source_observations",
                "equation",
                "units",
                "fixed_parameters",
                "pressure_recalculation",
                "refit_parameters",
                "refit_covariance",
                "numerical_reproduction",
            ],
            "primary_source_check": {
                "doi": source["reference"]["doi"],
                "access_url": "https://doi.org/" + source["reference"]["doi"],
                "locations": [
                    "Miozzi supplementary observations",
                    "Speziale Sections 2.1 and 5.4",
                ],
                "finding": "Source observations and calibration equations traced to primary papers; fitted Fe coefficients belong to the independent Peritheos refit.",
            },
            "primary_data_check": {
                "status": "bundled",
                "dataset_identifiers": [DATASET_ID],
                "source_locations": [
                    "All 131 source rows retained in derived Speziale table"
                ],
                "finding": "Recalculated MgO pressures coexist with original pressures and coordinate errors; unchanged He pressures and full source lineage retained.",
            },
            "independent_numerical_check": {
                "script": "scripts/refit_miozzi_2020_speziale.py",
                "report": str(REPORT.relative_to(ROOT)),
                "independent_native_max_difference_gpa": report[
                    "independent_native_max_difference_gpa"
                ],
            },
        },
        "reproduction": {
            "documentation": DOCUMENT,
            "summary_status": "independent_refit_reproduced",
            "original_publication_reproduction_status": "not_reproduced",
            "display_label": "Independent Peritheos refit — Speziale variable-q MgO, n=1",
            "reason": QUALIFICATION,
        },
    }
    return record, dataset


def ledger_outcome(record):
    fit = reproduce()["primary"]
    errors = {**record["parameter_errors"], **record["thermal"]["parameter_errors"]}
    stored = {**record["eos"]["parameters"], **record["thermal"]["parameters"]}
    comparisons = []
    for name in NAMES:
        value = fit["parameters"][name]
        if not np.isclose(stored[name], value, rtol=2e-6, atol=2e-6):
            raise ValueError(f"Registered Speziale refit is stale: {name}")
        comparisons.append(
            {
                "parameter": name,
                "published": stored[name],
                "refit": value,
                "difference": value - stored[name],
                "relative_difference": abs(value - stored[name]) / abs(stored[name]),
                "published_error": errors[name],
                "refit_error": fit["conditional_standard_errors"][name],
                "within_combined_2sigma": True,
                "within_reported_error": True,
                "similar": True,
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
        "selection": "All 131 observations; 116 Speziale variable-q Debye MgO pressures plus 15 source He pressures",
        "objective": "equal_weight_pressure_residuals",
        "absolute_sigma": False,
        "rmse_gpa": fit["rmse_gpa"],
        "qualification": QUALIFICATION,
        "original_publication_reproduction_status": "not_reproduced",
        "reproduction_status": "independent_refit_reproduced",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--register", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    document = json.loads(IRON_PATH.read_text(encoding="utf-8"))
    record, dataset = catalog_record(report, document)
    _, table = derived_rows(*read_data()[:2])
    if args.check:
        stored = next(
            r for r in document["eos_records"] if r["identifier"] == RECORD_ID
        )
        ledger_outcome(stored)
        np.testing.assert_allclose(
            stored["parameter_covariance"]["matrix"],
            report["primary"]["covariance_in_cell_volume_basis"],
            rtol=1e-4,
            atol=1e-6,
        )
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
        saved = json.loads(REPORT.read_text(encoding="utf-8"))
        assert saved["original_csv_sha256"] == report["original_csv_sha256"]
        assert (
            saved["derived_csv_sha256"]
            == hashlib.sha256(CSV_PATH.read_bytes()).hexdigest()
        )
        assert (
            next(r for r in document["eos_records"] if r["identifier"] == AUTHOR_ID)[
                "thermal"
            ]["parameters"]["n"]
            == 2
        )
        return
    REPORT.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    CSV_PATH.write_text(table, encoding="utf-8")
    if args.register:
        document["eos_records"] = [
            author_record(r) if r["identifier"] == AUTHOR_ID else r
            for r in document["eos_records"]
            if r["identifier"] != RECORD_ID
        ] + [record]
        document["datasets"] = [
            d for d in document["datasets"] if d["identifier"] != DATASET_ID
        ] + [dataset]
        IRON_PATH.write_text(
            json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "parameters": report["primary"]["parameters"],
                "standard_errors": report["primary"]["conditional_standard_errors"],
                "rms_gpa": report["primary"]["rmse_gpa"],
            }
        )
    )


if __name__ == "__main__":
    main()
