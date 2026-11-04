"""Staged Hirose/Fei Au refit on the published, author-confirmed Speziale EOS.

python -m scripts.refit_hirose_fei_gold_speziale --register
python -m scripts.refit_hirose_fei_gold_speziale --check

Adds only its own optional gold record and derived dataset. Published records,
original tables, defaults and shared generated ledgers/manifests are preserved.
Fei's own EOS work remains paused; his 2004 measurements are source inputs here.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
import re
from importlib.metadata import version
from pathlib import Path

import numpy as np
import scipy
from jsonschema import Draft202012Validator
from scipy.optimize import least_squares

from peritheos import eosmat_schema, get_eos_record, get_material_document
from scripts.check_numerical_archive import check_csv
from scripts.reproduce_fei_2007_gold import check_saved
from scripts.reproduce_hirose_2008_gold import source_pressure
from scripts.reproduce_speziale_2001_mgo import pressure as mgo_pressure

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
RECORD_ID = "gold_hirose_2008_fei_2004_speziale_2001_bm3_refit"
DATASET_ID = "gold_hirose_2008_fei_2004_speziale_2001_recalculated_pvt"
MGO_ID = "mgo_speziale_2001_bm3_2"
CSV = DATA / "gold-hirose-2008-fei-2004-speziale-2001-recalculated-pvt.csv"
REPORT = ROOT / "docs/data/hirose-2008-fei-2004-gold-speziale-2001-refit.json"
MATERIAL = ROOT / "peritheos/data/materials/gold.eosmat"
NAMES = ["beta1", "beta2", "beta3", "kprime_log_coefficient"]
INITIAL = np.array([1.03e-6, 3.95e-10, 1.61e-13, 3.61e-4])
SCALES = np.array([1e-6, 1e-10, 1e-13, 1e-4])
NOTE = (
    "Independent Peritheos staged refit of Hirose (2008) and Fei (2004) Au-MgO "
    "measurements. All 47 pressure targets are recalculated with the published, "
    "author-confirmed Speziale (2001) variable-q Debye EOS. Nine 300 K Hirose "
    "rows determine K0-prime with V0=67.85 A^3 and K0=167 GPa fixed; K0-prime "
    "is then fixed while 26 Fei and 12 hot Hirose rows determine beta1,beta2,"
    "beta3 and the logarithmic Kprime coefficient. Ambient expansivity remains "
    "source-fixed. Equal pressure-residual weights are an analyst choice. "
    "This is not a Hirose/Fei author fit or a Fei EOS refit. Stage-conditional "
    "curvature widths/covariance are archived only; catalog parameter errors "
    "remain unset because predictor, calibrant, fixed-anchor and inter-stage "
    "uncertainty/covariance are not propagated. Numerical reproducibility does "
    "not establish absolute pressure accuracy or independent predictive validity."
)


def gold_pressure(volume, temperature, parameters, kprime):
    """Independent Hirose fit #2 expression with the newly determined cold K'."""
    v, t = np.broadcast_arrays(volume, temperature)
    b1, b2, b3, a = parameters
    v0 = 67.85 * np.exp(3.824e-5 * (t - 300) + 1.499e-8 * (t**2 - 300**2) / 2)
    k0 = 1 / (1 / 167 + b1 * (t - 300) + b2 * (t**2 - 300**2) + b3 * (t**3 - 300**3))
    kp = kprime + a * (t - 300) * np.log(t / 300)
    x = (v0 / v) ** (1 / 3)
    return 1.5 * k0 * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def derived_rows():
    rows, hashes = [], {}
    for source, name in (
        ("Hirose2008", "gold-hirose-2008-table1.csv"),
        ("Fei2004", "gold-fei-2004-table1.csv"),
    ):
        path = DATA / name
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open() as stream:
            original = list(csv.DictReader(stream))
        for number, row in enumerate(original, 1):
            # Preserve source numeric strings and missing widths; calculated
            # pressures are in their own explicitly derived coordinate.
            rows.append(
                {
                    "source": source,
                    "source_row": number,
                    "run": row["run"],
                    "fit_stage": "cold"
                    if float(row["temperature_k"]) == 300
                    else "thermal",
                    "volume_a3": row["volume_a3"],
                    "volume_uncertainty_a3": row["volume_uncertainty_a3"],
                    "temperature_k": row["temperature_k"],
                    "mgo_volume_a3": row["mgo_volume_a3"],
                    "reported_pressure_gpa": row["pressure_gpa"],
                    "reported_pressure_width_gpa": row["pressure_gpa_uncertainty"],
                    "pressure_gpa": float(
                        mgo_pressure(
                            float(row["mgo_volume_a3"]), float(row["temperature_k"])
                        )
                    ),
                }
            )
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return rows, hashes, stream.getvalue()


def stage_report(fit, scales):
    covariance = (
        np.linalg.inv(fit.jac.T @ fit.jac)
        * (fit.fun @ fit.fun)
        / (len(fit.fun) - len(fit.x))
    )
    covariance *= scales[:, None] * scales[None, :]
    return {
        "observations": len(fit.fun),
        "degrees_of_freedom": len(fit.fun) - len(fit.x),
        "parameters": (fit.x * scales).tolist(),
        "conditional_standard_errors": np.sqrt(np.diag(covariance)).tolist(),
        "conditional_covariance": covariance.tolist(),
        "scaled_jacobian_condition_number": float(np.linalg.cond(fit.jac)),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        "pressure_mae_gpa": float(np.mean(np.abs(fit.fun))),
        "solver_success": bool(fit.success),
    }


def reproduce():
    rows, hashes, table = derived_rows()
    v, vm, t, p = (
        np.array([float(r[k]) for r in rows])
        for k in ("volume_a3", "mgo_volume_a3", "temperature_k", "pressure_gpa")
    )
    cold, hot = t == 300, t > 300
    cold_fit = least_squares(
        lambda z: source_pressure(v[cold], t[cold], fit=0, coefficients=z) - p[cold],
        [5.58],
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    kp = float(cold_fit.x[0])
    hot_fit = least_squares(
        lambda z: gold_pressure(v[hot], t[hot], z * SCALES, kp) - p[hot],
        INITIAL / SCALES,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=4000,
    )
    assert cold_fit.success and hot_fit.success
    assert (
        np.max(np.abs(gold_pressure(v, t, INITIAL, 5.58) - source_pressure(v, t)))
        < 1e-9
    )
    native = get_eos_record(MGO_ID).pressure(vm, t)
    assert np.max(np.abs(native - p)) < 1e-9
    predicted = gold_pressure(v, t, hot_fit.x * SCALES, kp)
    delta = predicted - p
    report = {
        "record_identifier": RECORD_ID,
        "dataset_identifier": DATASET_ID,
        "source_hashes": hashes,
        "derived_csv_sha256": hashlib.sha256(table.encode()).hexdigest(),
        "implementation_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "software": {
            "name": "Peritheos / scipy.optimize.least_squares",
            "version": version("peritheos"),
            "scipy": scipy.__version__,
            "numpy": np.__version__,
            "version_date": "2026-10-12",
        },
        "calibration": {
            "reference_eos_record": MGO_ID,
            "gamma_law": "Speziale Equation 11",
            "theta_law": "integral of d ln(theta) = -gamma(V) d ln(V)",
            "reference_temperature_k": 300,
            "native_equation_max_difference_gpa": float(np.max(np.abs(native - p))),
        },
        "cold_stage": stage_report(cold_fit, np.ones(1)),
        "thermal_stage": stage_report(hot_fit, SCALES),
        "thermal_parameter_names": NAMES,
        "fixed_anchors": {
            "V0": 67.85,
            "K0": 167,
            "Tr": 300,
            "alpha0": 3.824e-5,
            "alpha1": 1.499e-8,
        },
        "observed_ranges": {
            "pressure_gpa": [float(p.min()), float(p.max())],
            "temperature_k": [float(t.min()), float(t.max())],
            "hot_pressure_gpa": [float(p[hot].min()), float(p[hot].max())],
            "hot_temperature_k": [float(t[hot].min()), float(t[hot].max())],
        },
        "residuals_gpa": delta.tolist(),
        "selection": "All nine cold Hirose rows first, then all 26 Fei2004 plus 12 hot Hirose2008 rows; 8.64 GPa Fei source row retained; repeated Hirose2006 rows not added",
        "objective": "sum of squared pressure residuals with equal weights at each stage",
        "uncertainty_scope": NOTE,
        "original_author_fit_reproduced": False,
        "fei_eos_work_paused": True,
    }
    return report, table


def catalog_objects(report):
    document = get_material_document("gold")
    source = next(
        r
        for r in document["eos_records"]
        if r["identifier"] == "gold_hirose_2008_bm3_fit2"
    )
    fei = next(
        d for d in document["datasets"] if d["identifier"] == "gold_fei_2004_table1"
    )
    # The raw reference mapping also preserves its source locator and DOI.
    mgo = next(
        r
        for r in get_material_document("mgo")["eos_records"]
        if r["identifier"] == MGO_ID
    )["reference"]
    p = dict(zip(NAMES, report["thermal_stage"]["parameters"]))
    record = {
        "identifier": RECORD_ID,
        "label": "Peritheos refit: Hirose (2008) + Fei (2004), Au on author-confirmed Speziale (2001) MgO",
        "record_kind": "refit",
        "determination_method": "experimental",
        "equation_kind": "thermal",
        "default": False,
        "reference": copy.deepcopy(source["reference"]),
        "supporting_references": [copy.deepcopy(fei["reference"]), copy.deepcopy(mgo)],
        "volume_basis": copy.deepcopy(source["volume_basis"]),
        "eos": {
            "type": "BM3",
            "model": "birch_murnaghan_3",
            "parameters": {
                "V0": 67.85,
                "K0": 167,
                "K0_prime": report["cold_stage"]["parameters"][0],
            },
        },
        "parameter_errors": {"V0": None, "K0": None, "K0_prime": None},
        "parameter_error_confidence": None,
        "parameter_covariance": None,
        "fixed_parameters": ["V0", "K0"],
        "temperature_ref": 300,
        "thermal": {
            "type": "AlphaKT",
            "model": "thermal_reference_state",
            "thermal_expansion_law": "linear_temperature",
            "bulk_modulus_law": "reciprocal_cubic",
            "parameters": {
                "Tr": 300,
                "alpha0": 3.824e-5,
                "alpha1": 1.499e-8,
                "dK_dT": 0,
                **p,
            },
            "parameter_errors": {
                k: None for k in ("Tr", "alpha0", "alpha1", "dK_dT", *NAMES)
            },
            "fixed_parameters": ["Tr", "alpha0", "alpha1", "dK_dT"],
        },
        "experimental_pressure_range_gpa": report["observed_ranges"]["pressure_gpa"],
        "experimental_temperature_range_k": report["observed_ranges"]["temperature_k"],
        "pressure_range_status": "reference_parameterization",
        "validity": {
            "pressure_gpa": report["observed_ranges"]["pressure_gpa"],
            "temperature_k": report["observed_ranges"]["temperature_k"],
            "notes": [
                "Marginal observed coverage only, not a validated P-T rectangle. Hot coverage: "
                + str(report["observed_ranges"]["hot_pressure_gpa"])
                + " GPa and "
                + str(report["observed_ranges"]["hot_temperature_k"])
                + " K. No 3000 K extrapolation recommendation inherited."
            ],
        },
        "fit_datasets": [DATASET_ID],
        "notes": NOTE,
        "parameter_provenance": {
            "determination_method": "Independent Peritheos regression; source DOI attributes data and functional form, not authorship of these coefficients.",
            "cold_stage": "Kprime from nine recalculated 300 K observations, then fixed in thermal regression",
            "thermal_stage": "Four thermal coefficients from 38 recalculated hot observations; ambient expansivity source-fixed",
            "uncertainties": NOTE,
        },
        "fit_provenance": {
            "software": report["software"],
            "dataset": DATASET_ID,
            "selection": {
                "included_rows": 47,
                "cold_rows": 9,
                "thermal_rows": 38,
                "predicate": report["selection"],
            },
            "objective": report["objective"],
            "refined_parameters": ["rt_eos.K0_prime", *NAMES],
            "fixed_parameters": [
                "rt_eos.V0",
                "rt_eos.K0",
                "Tr",
                "alpha0",
                "alpha1",
                "dK_dT",
            ],
            "statistics": {
                "observations": 47,
                "thermal_pressure_rmse_gpa": report["thermal_stage"][
                    "pressure_rmse_gpa"
                ],
                "thermal_pressure_mae_gpa": report["thermal_stage"]["pressure_mae_gpa"],
                "solver_success": True,
            },
            "reproduction": {
                "script": "scripts/refit_hirose_fei_gold_speziale.py",
                "report": REPORT.relative_to(ROOT).as_posix(),
                "source_hashes": report["source_hashes"],
                "derived_csv_sha256": report["derived_csv_sha256"],
            },
        },
        "pressure_calibration": {
            "status": "resolved",
            "audit_date": "2026-10-12",
            "methods": [
                {
                    "kind": "equation_of_state",
                    "reference": copy.deepcopy(mgo),
                    "reference_eos_record": MGO_ID,
                    "scope": "All 47 pressure targets explicitly recalculated from paired MgO volumes and source T",
                    "source_location": "Speziale Eq. 11 plus author-confirmed theta integral; derived CSV and refit report",
                }
            ],
            "recalculation": {
                "status": "ready",
                "notes": "Published Speziale model is the explicitly chosen new calibration. Historical Fei/Hirose reducer remains unresolved and is not claimed reproduced.",
            },
        },
        "scientific_validation": {
            "status": "primary_source_validated",
            "audit_date": "2026-10-12",
            "refit_status": "conditional",
            "reproduction_status": "independent_refit_reproduced",
            "original_publication_reproduction_status": "not_reproduced",
            "note": NOTE,
            "verified_fields": [
                "source_observations",
                "equation",
                "units",
                "pressure_recalculation",
                "fixed_parameters",
                "refit_parameters",
            ],
            "primary_source_check": {
                "access_url": "https://doi.org/10.1016/j.pepi.2008.03.002",
                "doi": "10.1016/j.pepi.2008.03.002",
                "locations": [
                    "Hirose (2008), Table 1 and Section 3.1 / Table 2 fit #2",
                    "Fei (2004), Table 1",
                    "Speziale (2001), Equation (11) and recorded author confirmation",
                    "Peritheos staged regression and derived-data report",
                ],
                "finding": NOTE,
            },
            "primary_data_check": {
                "status": "bundled",
                "dataset_identifiers": [DATASET_ID],
                "finding": "All 47 source observations are retained with paired MgO volumes; the pressure targets are explicitly recalculated, not historical author pressures.",
            },
            "independent_numerical_check": {
                "script": "scripts/refit_hirose_fei_gold_speziale.py",
                "report": REPORT.relative_to(ROOT).as_posix(),
            },
        },
    }
    rows = list(csv.DictReader(io.StringIO(derived_rows()[2])))
    quantities = {
        "volume_a3": ("volume", "angstrom^3", "value"),
        "volume_uncertainty_a3": ("volume", "angstrom^3", "uncertainty"),
        "mgo_volume_a3": ("calibrant_volume", "angstrom^3", "value"),
        "temperature_k": ("temperature", "K", "value"),
        "reported_pressure_gpa": ("comparison_pressure", "GPa", "value"),
        "reported_pressure_width_gpa": ("comparison_pressure", "GPa", "uncertainty"),
        "pressure_gpa": ("pressure", "GPa", "value"),
    }
    columns = []
    for key in rows[0]:
        quantity, unit, role = quantities.get(key, (key, "1", "flag"))
        col = {"name": key, "quantity": quantity, "unit": unit, "role": role}
        if key == "volume_uncertainty_a3":
            col["of"] = "volume_a3"
        if key == "reported_pressure_width_gpa":
            col["of"] = "reported_pressure_gpa"
        columns.append(col)
    dataset = {
        "identifier": DATASET_ID,
        "kind": "recalculated_pressure_volume_temperature",
        "reference": copy.deepcopy(source["reference"]),
        "supporting_references": [copy.deepcopy(fei["reference"]), copy.deepcopy(mgo)],
        "description": "Nine cold and twelve hot Hirose observations plus 26 hot Fei observations; pressure recalculated with published Speziale MgO",
        "source_location": "Hirose2008 Table1 and Fei2004 Table1; source identifiers and row numbers preserved",
        "columns": columns,
        "used_by_eos_records": [RECORD_ID],
        "resource": {
            "path": CSV.relative_to(ROOT / "peritheos/data").as_posix(),
            "sha256": report["derived_csv_sha256"],
            "media_type": "text/csv",
        },
        "pressure_reconstruction": {
            "reference_eos_record": MGO_ID,
            "source_pressure_column": "reported_pressure_gpa",
            "target_pressure_column": "pressure_gpa",
        },
        "uncertainty": {
            "notes": "Original pressure widths refer to reported pressures only; no uncertainty on recalculated pressure or missing Au-volume width is invented. Source confidence conventions and covariance are unknown."
        },
        "provenance": {
            "source_datasets": ["gold_hirose_2008_table1", "gold_fei_2004_table1"],
            "source_hashes": report["source_hashes"],
        },
        "notes": "Third-party source rights remain with their owners; no source-rights transfer or source CC0 claim. Fei EOS work remains paused; source measurements are used for the separate combined refit.",
    }
    return record, dataset


def append_object(text, key, obj):
    """Add/update this refit's member without reserializing other members."""
    match = re.search(r'\n "' + re.escape(key) + r'": \[', text)
    assert match, key
    entries, end = json.JSONDecoder().raw_decode(text, match.end() - 1)
    if any(e["identifier"] == obj["identifier"] for e in entries):
        index = match.end()
        for _ in entries:
            while text[index].isspace() or text[index] == ",":
                index += 1
            entry, stop = json.JSONDecoder().raw_decode(text, index)
            if entry["identifier"] == obj["identifier"]:
                encoded = json.dumps(obj, indent=1).replace("\n", "\n  ")
                return text[:index] + encoded + text[stop:]
            index = stop
    encoded = "\n".join("  " + line for line in json.dumps(obj, indent=1).splitlines())
    return text[: end - 1].rstrip() + ",\n" + encoded + "\n " + text[end - 1 :]


def ledger_outcome(record):
    """Reproduce this staged refit, without treating it as source-author parity."""
    report, table = reproduce()
    archived = json.loads(REPORT.read_text(encoding="utf-8"))
    check_archives(archived, report, table)
    values = {
        "rt_eos.K0_prime": report["cold_stage"]["parameters"][0],
        **dict(zip(NAMES, report["thermal_stage"]["parameters"])),
    }
    stored = {
        "rt_eos.K0_prime": record["eos"]["parameters"]["K0_prime"],
        **{name: record["thermal"]["parameters"][name] for name in NAMES},
    }
    comparisons = [
        {
            "parameter": name,
            "published": stored[name],
            "published_error": None,
            "refit": value,
            "refit_error": None,
            "relative_difference": abs(value - stored[name]) / abs(stored[name]),
            "within_combined_2sigma": None,
            "similar": bool(np.isclose(stored[name], value, rtol=1e-4, atol=1e-16)),
        }
        for name, value in values.items()
    ]
    assert all(item["similar"] for item in comparisons)
    return {
        "status": "parity",
        "parity_basis": "registered_peritheos_refit_reproduced",
        "reproduction_status": "independent_refit_reproduced",
        "original_publication_reproduction_status": "not_reproduced",
        "dataset_identifiers": [DATASET_ID],
        "observations": 47,
        "fit_kind": "staged_cold_then_thermal_bm3_reference_state",
        "parameters": comparisons,
        "objective": report["objective"],
        "rmse_gpa": report["thermal_stage"]["pressure_rmse_gpa"],
        "solver_success": True,
        "qualification": NOTE,
        "reason": NOTE,
        # Keep the archived run's versions and hashes with its diagnostics.
        # The independently replayed coefficients are compared above.
        "reproduction": archived,
    }


def check_archives(saved, current, table):
    """Verify archived bytes exactly and recalculated pressure columns numerically."""
    assert hashlib.sha256(CSV.read_bytes()).hexdigest() == saved["derived_csv_sha256"]
    assert hashlib.sha256(table.encode()).hexdigest() == current["derived_csv_sha256"]
    # Versions and implementation hashes describe the archived run. A replay
    # may use another supported environment or a revised verification script.
    provenance = {"derived_csv_sha256", "implementation_sha256", "software"}
    assert (
        hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        == current["implementation_sha256"]
    )
    assert saved["software"].keys() == current["software"].keys()
    assert saved["software"]["name"] == current["software"]["name"]
    check_saved(
        {k: v for k, v in saved.items() if k not in provenance},
        {k: v for k, v in current.items() if k not in provenance},
    )
    check_csv(CSV.read_text(encoding="utf-8"), table, {"pressure_gpa"}, atol=1e-8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--register", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report, table = reproduce()
    if args.check:
        check_archives(json.loads(REPORT.read_text(encoding="utf-8")), report, table)
        record = get_eos_record(RECORD_ID)
        rows = list(csv.DictReader(io.StringIO(table)))
        predicted = np.array(
            [
                record.pressure(float(r["volume_a3"]), float(r["temperature_k"]))
                for r in rows
            ]
        )
        p = np.array([float(r["pressure_gpa"]) for r in rows])
        np.testing.assert_allclose(predicted - p, report["residuals_gpa"], atol=1e-8)
    else:
        REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        CSV.write_text(table)
        if args.register:
            before = MATERIAL.read_text(encoding="utf-8")
            record, dataset = catalog_objects(report)
            after = append_object(
                append_object(before, "eos_records", record), "datasets", dataset
            )
            Draft202012Validator(eosmat_schema()).validate(json.loads(after))
            old, new = json.loads(before), json.loads(after)
            for key, own_id in (("eos_records", RECORD_ID), ("datasets", DATASET_ID)):
                assert [r for r in old[key] if r["identifier"] != own_id] == [
                    r for r in new[key] if r["identifier"] != own_id
                ], "Registration changed an unrelated record/dataset"
            assert MATERIAL.read_text(encoding="utf-8") == before, (
                "Concurrent gold edit; retry"
            )
            MATERIAL.write_text(after)
    print(
        f"{RECORD_ID}: cold Kprime={report['cold_stage']['parameters'][0]:.8f}; thermal MAE={report['thermal_stage']['pressure_mae_gpa']:.6f} GPa"
    )


if __name__ == "__main__":
    main()
