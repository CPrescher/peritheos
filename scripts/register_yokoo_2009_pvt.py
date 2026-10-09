"""Register bounded Au/Pt Yokoo Table-output pressure reconstructions."""

from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np

from peritheos import get_eos_record
from scripts.fit_yokoo_2009_gold_thermal import reconstruct as gold_reconstruction
from scripts.fit_yokoo_2009_platinum_thermal import (
    reconstruct as platinum_reconstruction,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
RECORDS = {
    "gold": "gold_yokoo_2009_pvt_reconstruction",
    "platinum": "platinum_yokoo_2009_pvt_reconstruction",
}
OUTPUT = ROOT / "docs/data/yokoo-2009-pvt-library-validation.json"


def reports():
    return {"gold": gold_reconstruction(), "platinum": platinum_reconstruction()}


def record_definition(metal, document, report):
    original = next(
        r
        for r in document["eos_records"]
        if r["identifier"] == f"{metal}_yokoo_2009_vinet_300k"
    )
    fit = report["fits"]["primary"]
    p = fit["parameters"]
    with (
        ROOT
        / "peritheos/data/datasets/tsuchiya-kawamura-2002-table1-electronic-pressure.csv"
    ).open() as stream:
        electronic = list(csv.DictReader(stream))
    column = (
        "gold_corrected_electronic_pressure_gpa"
        if metal == "gold"
        else "platinum_electronic_pressure_gpa"
    )
    table = "III" if metal == "gold" else "V"
    grid = f"{metal}_yokoo_2009_table{3 if metal == 'gold' else 5}_isochores"
    correction_t = fit.get("residual_temperature_k", [0, 3000])
    correction_p = fit.get("residual_pressure_gpa", [0, 0])
    notes = (
        f"Diagnostic reconstruction of Yokoo (2009) Table {table} calculated PVT outputs; published analytical PVT remains unreproduced. Not a validated replacement for the source pressure standard. "
        "Cold BM3 V0 is Vc at 0 K; ambient normalization is Vc/cold_volume_ratio. Phonon and electronic pressures are absolute relative to 0 K. "
        "Electronic pressure is assumed volume independent and interpolated linearly; no extrapolation. Source first-liquid annotations are preserved and excluded from the primary fit. "
        "Numerical bounds do not guarantee solid Au/Pt throughout the rectangular domain; use the source phase markers. "
        + (
            "Pt has an additional derived residual-pressure correction, linear between source temperatures, retained separately from the original electronic table. Its physical origin is unresolved. "
            if metal == "platinum"
            else ""
        )
        + "No parameter errors or covariance are inferred from model-output residuals. This pressure model supplies no caloric potential."
    )
    return {
        "identifier": RECORDS[metal],
        "label": f"Yokoo et al. (2009), {metal.title()} PVT — diagnostic Table {table} reconstruction",
        "default": False,
        "catalog_access": "explicit_selection",
        "reference": copy.deepcopy(original["reference"]),
        "record_kind": "derived",
        "equation_kind": "thermal",
        "determination_method": "hybrid",
        "eos": {
            "type": "BM3",
            "model": "birch_murnaghan_3",
            "parameters": {
                "V0": report["cold_volume_cell_a3"],
                "K0": p["k0_gpa"],
                "K0_prime": p["k0_prime"],
            },
        },
        "parameter_errors": {"V0": None, "K0": None, "K0_prime": None},
        "parameter_error_confidence": None,
        "parameter_covariance": None,
        "fixed_parameters": ["K0"] if metal == "gold" else ["K0", "K0_prime"],
        "thermal": {
            "fixed_parameters": ["Tr", "theta0", "n"]
            if metal == "gold"
            else ["Tr", "theta0", "gamma0", "n"],
            "type": "AsymptoticDebyeTabulatedPressure",
            "model": "asymptotic_debye_tabulated_pressure",
            "parameters": {
                "Tr": 300.0,
                "theta0": p["theta0_k"],
                "gamma0": p["gamma0"],
                "a": p["a"],
                "b": p["b"],
                "n": 1.0,
                "cold_volume_ratio": p["vc_over_v0"],
            },
            "parameter_errors": {
                k: None
                for k in ["Tr", "theta0", "gamma0", "a", "b", "n", "cold_volume_ratio"]
            },
            "configuration": {
                "electronic_temperature_k": [
                    float(r["temperature_k"]) for r in electronic
                ],
                "electronic_pressure_gpa": [float(r[column]) for r in electronic],
                "interpolation": "linear",
                "volume_ratio_range": [0.6, 1.0],
                "temperature_range_k": [0, 3000],
                "residual_temperature_k": correction_t,
                "residual_pressure_gpa": correction_p,
            },
        },
        "validity": {
            "temperature_k": [0, 3000],
            "volume_ratio": [0.6, 1.0],
            "notes": [
                "Bounds of pressure reconstruction, not independently measured or guaranteed solid-phase validity. First-liquid markers excluded from fitting; missing source cells are not fabricated observations."
            ],
        },
        "derivation": {
            "diagnostic_only": True,
            "published_analytical_pvt_reproduced": False,
            "pressure_convention_audit": "docs/data/yokoo-2009-pressure-conventions.json",
            "source_kind": "published_table",
            "source_identifier": grid,
            "source_version": "Yokoo et al. PRB 80, 104114 (2009)",
            "method": report["conventions"]["objective"],
            "software": {
                "script": f"scripts/fit_yokoo_2009_{metal}_thermal.py",
                "report": f"docs/data/yokoo-2009-{metal}-thermal-reconstruction.json",
            },
            "input_sha256": report["input_sha256"],
            "parameter_comparison": report["parameter_comparison"],
            "source_phase_markers": [
                r
                for r in report["states"]
                if r["source_phase_annotation"] != "unmarked"
            ],
            "experimental_observations": 0,
            "author_fit_reproduced": False,
        },
        "parameter_provenance": {
            "V0": "Inferred cold Vc from model outputs, distinct from ambient volume normalization",
            "thermal": f"Derived Table {table} pressure reconstruction; published and reconstructed coefficients separately archived",
            "electronic_pressure": "Tsuchiya-Kawamura (2002) Table I, Au-only 2003 erratum applied; explicit linear interpolation",
            "residual_pressure": "Empirical model-output residual correction; not identified as electronic pressure"
            if metal == "platinum"
            else "Zero; no additional residual correction",
        },
        "pressure_calibration": {
            "status": "not_applicable",
            "methods": copy.deepcopy(original["pressure_calibration"]["methods"]),
            "recalculation": {
                "status": "not_applicable",
                "notes": "Reconstruction of source-calculated pressures, not an experimental recalibration",
            },
            "audit_date": "2026-10-09",
        },
        "scientific_validation": {
            "status": "not_reproduced",
            "audit_date": "2026-10-09",
            "note": "Diagnostic table-output reconstruction. Its numerical implementation, inversions and canonical interchange are checked, but the published analytical Yokoo PVT parameterization is not reproduced. Au uses altered phonon coefficients; Pt additionally uses an empirical residual-pressure correction. These are not validated replacements for the published pressure standard. See the pressure-convention audit; source caloric data and original optimization are not required to evaluate PVT.",
            "verified_fields": [
                "equation",
                "parameters",
                "units",
                "reference_state",
                "published_uncertainties",
                "validity",
                "source_phase_annotations",
            ],
            "primary_source_check": {
                **copy.deepcopy(
                    original["scientific_validation"]["primary_source_check"]
                ),
                "finding": f"Equations 1-12 and Table {table} supply the pressure reconstruction inputs; numerical output agreement is validated separately from original author-fit reproduction.",
            },
            "primary_data_check": {
                "status": "parameterization_only",
                "audit_date": "2026-10-09",
                "dataset_identifiers": [
                    grid,
                    "tsuchiya_kawamura_2002_table1_electronic_pressure",
                ],
                "finding": "Derived source outputs, not primary experimental observations. Exact author optimization is not reproduced.",
                "reproduction_resource": f"docs/data/yokoo-2009-{metal}-thermal-reconstruction.json",
            },
        },
        "notes": notes,
    }


def register():
    for metal, report in reports().items():
        path = MATERIALS / f"{metal}.eosmat"
        document = json.loads(path.read_text())
        record = record_definition(metal, document, report)
        existing = next(
            (
                i
                for i, r in enumerate(document["eos_records"])
                if r["identifier"] == record["identifier"]
            ),
            None,
        )
        if existing is None:
            document["eos_records"].append(record)
        else:
            document["eos_records"][existing] = record
        for dataset in document["datasets"]:
            if (
                dataset["identifier"]
                in record["scientific_validation"]["primary_data_check"][
                    "dataset_identifiers"
                ]
            ):
                if record["identifier"] not in dataset["used_by_eos_records"]:
                    dataset["used_by_eos_records"].append(record["identifier"])
        path.write_text(
            json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n"
        )
    audit_path = ROOT / "peritheos/data/primary-source-audit.json"
    audit = json.loads(audit_path.read_text())
    for metal in RECORDS:
        document = json.loads((MATERIALS / f"{metal}.eosmat").read_text())
        record = next(
            r for r in document["eos_records"] if r["identifier"] == RECORDS[metal]
        )
        check = record["scientific_validation"]
        entry = {
            "material": document["identifier"],
            "file": f"{metal}.eosmat",
            "record": record["identifier"],
            "label": record["label"],
            "doi": record["reference"]["doi"].lower(),
            "status": check["status"],
            "note": check["note"],
            "primary_source_check": check["primary_source_check"],
        }
        audit["records"] = [
            r for r in audit["records"] if r["record"] != record["identifier"]
        ] + [entry]
    audit["records"].sort(key=lambda r: (r["material"], r["record"]))
    audit["summary"] = {
        "records": len(audit["records"]),
        **dict(sorted(Counter(r["status"] for r in audit["records"]).items())),
    }
    audit_path.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    documents = [json.loads(path.read_text()) for path in MATERIALS.glob("*.eosmat")]
    records = [r for d in documents for r in d["eos_records"]]
    path = MATERIALS / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest["eos_records"] = len(records)
    manifest["scientific_validation"]["counts"] = dict(
        sorted(Counter(r["scientific_validation"]["status"] for r in records).items())
    )
    manifest["pressure_calibration"]["status_counts"] = dict(
        sorted(Counter(r["pressure_calibration"]["status"] for r in records).items())
    )
    manifest["pressure_calibration"]["recalculation_counts"] = dict(
        sorted(
            Counter(
                r["pressure_calibration"]["recalculation"]["status"] for r in records
            ).items()
        )
    )
    path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )


def validate_library():
    results = {}
    for metal, report in reports().items():
        record = get_eos_record(RECORDS[metal])
        rows = report["states"]
        volume = np.array(
            [
                float(r["volume_ratio"]) * report["conventions"]["ambient_cell_a3"]
                for r in rows
            ]
        )
        temperature = np.array([float(r["temperature_k"]) for r in rows])
        oracle = np.array([r["model_pressure_gpa"] for r in rows])
        pressure = record.pressure(volume, temperature)
        from peritheos import Material, get_material_document

        material = Material.from_eosmat(
            get_material_document(metal),
            record_identifiers=[record.identifier],
            require_primary_validation=False,
        )
        restored = Material.from_eosmat(
            material.to_eosmat(), require_primary_validation=False
        ).eos_records[0]
        from peritheos.errors import MaterialError

        try:
            Material.from_dict(material.to_snapshot_dict())
        except MaterialError:
            snapshot_rejected = True
        else:
            raise AssertionError(
                "An unvalidated composition must not import as a validated snapshot"
            )
        results[metal] = {
            "record_identifier": record.identifier,
            "kind": "derived_table_output_pressure_reconstruction",
            "states": len(rows),
            "source_fit": report["fits"]["primary"]["fit_states"],
            "volume_holdout": report["volume_interpolation_holdout"]["held_out_states"],
            "library_vs_independent_pressure_max_gpa": float(
                np.max(abs(pressure - oracle))
            ),
            "volume_inversion_max_cell_a3": float(
                np.max(abs(record.volume(pressure, temperature) - volume))
            ),
            "temperature_inversion_max_k": float(
                np.max(
                    abs(
                        record.eos.temperature(pressure, volume * record.volume_scale)
                        - temperature
                    )
                )
            ),
            "eosmat_roundtrip_max_gpa": float(
                np.max(abs(restored.pressure(volume, temperature) - pressure))
            ),
            "snapshot_import_rejected_for_unvalidated_composition": snapshot_rejected,
            "scientific_validation_status": record.scientific_validation_status,
            "published_analytical_pvt_reproduced": False,
            "ambient_volume_cell_a3": record.reference_volume,
            "cold_volume_cell_a3": report["cold_volume_cell_a3"],
            "original_author_fit_reproduced": False,
            "independent_experimental_validation": False,
        }
    return {
        "scope": "Numerical verification of diagnostic table-output reconstructions; published analytical PVT reproduction remains unresolved",
        "metals": results,
    }


def ledger_outcome(record):
    metal = next(m for m, k in RECORDS.items() if k == record["identifier"])
    report = reports()[metal]
    return {
        "status": "source_reconstruction",
        "reconstruction_kind": "derived_table_output_pressure_reconstruction",
        "observations": 0,
        "derived_source_states": len(report["states"]),
        "dataset_identifiers": record["scientific_validation"]["primary_data_check"][
            "dataset_identifiers"
        ],
        "rmse_gpa": report["fits"]["primary"]["fit_states"]["rmse_gpa"],
        "author_fit_reproduced": False,
        "published_analytical_pvt_reproduced": False,
        "diagnostic_only": True,
        "report": f"docs/data/yokoo-2009-{metal}-thermal-reconstruction.json",
        "reason": "Numerically checked diagnostic reconstruction of published table outputs. The published analytical PVT parameterization remains unreproduced: Au uses altered phonon coefficients and Pt adds an empirical pressure correction. This source-output status does not validate the published Yokoo pressure standard. The separate pressure-convention audit rules out cold-reference offsets and tests last-digit rounding; no caloric or original-fit inputs are required for PVT evaluation.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not args.check:
        register()
    report = validate_library()
    if args.check:
        from scripts.fit_yokoo_2009_gold_thermal import check_reconstruction

        check_reconstruction(json.loads(OUTPUT.read_text()), report)
    else:
        OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
