"""Integrate source-specific Dewaele 2008 records without changing old numbers.

Run after recover_dewaele_2008_metals.py; repository maintenance only.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json

from peritheos._eosmat_format import format_eosmat
from peritheos.eosmat import validate_eosmat_document
from scripts.reproduce_dewaele_2008_metals import (
    DATA,
    NAMES,
    REPORT,
    ROOT,
    Z,
    reproduce,
    source,
    table4_specs,
    table6_specs,
)
from scripts.standardize_reconstruction_metadata import assessment_for_record

DATE = "2026-10-23"
REF = {
    "authors": ["Dewaele", "Torrent", "Loubeyre", "Mezouar"],
    "year": 2008,
    "title": "Compression curves of transition metals in the Mbar range: Experiments and projector augmented-wave calculations",
    "source": "Physical Review B",
    "volume": "78",
    "locator": "104102",
    "doi": "10.1103/PhysRevB.78.104102",
}


def experimental_observation_view(material, records):
    """Reviewed Table IV coordinates; retain source rows outside the fit mask."""
    iron = material == "fe"
    view = {
        "identifier": f"{material}_dewaele_2008_table4_source_view",
        "eos_records": records,
        "evidence_kind": "experimental",
        "notes": (
            "Published revised-pressure column on the 2004 input scale; "
            "parenthesized helium-EOS pressures remain as printed. Atomic "
            "volumes are converted to the declared conventional cell. No "
            "rowwise errors or missing calibration readings are synthesized."
        ),
        "coordinates": {
            "pressure": {
                "kind": "column",
                "column": "pressure_revised_gpa"
                if iron
                else "ruby_pressure_dewaele_gpa",
            },
            "volume": {
                "kind": "column",
                "column": "atomic_volume_a3",
                "unit": "angstrom^3/conventional_unit_cell",
                "factor": Z[material],
                "reference": "Publication crystal phase and conventional-cell atom count.",
            },
        },
        "phase_reviewed": True,
        "fit_selection": {
            "notes": "Unweighted pressure-residual Table IV replay; V0 fixed, with source-specific fixed ultrasonic K0 alternatives."
        },
    }
    if iron:
        view["row_selection"] = [{"column": "phase", "values": ["bcc"]}]
    else:
        view["fit_selection"]["flag_column"] = "fit_included"
        view["fit_selection"]["notes"] += (
            " The flag excludes only run-2 Ag/Mo rows above 100 GPa, as stated in Section III."
        )
    return view


def dataset(identifier, filename, records, description, location, note):
    path = DATA / filename
    with path.open() as stream:
        names = next(csv.reader(stream))
    columns = []
    for name in names:
        if "pressure" in name and name.endswith("gpa"):
            quantity, unit, role = "pressure", "GPa", "value"
        elif name == "atomic_volume_a3":
            quantity, unit, role = "atomic_volume", "angstrom^3/atom", "value"
        elif name.endswith("angstrom"):
            quantity, unit, role = "lattice_parameter", "angstrom", "value"
        else:
            quantity, unit, role = name, "1", "flag"
        columns.append({"name": name, "quantity": quantity, "unit": unit, "role": role})
    metadata = {
        "identifier": identifier,
        "kind": "literature_reconstruction_diagnostic"
        if "figure4" in identifier
        else "pressure_volume_lattice_parameters_with_calibrant",
        "description": description,
        "source_location": location,
        "source_url": "https://doi.org/10.1103/PhysRevB.78.104102",
        "reference": REF,
        "columns": columns,
        "resource": {
            "path": f"datasets/{filename}",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "media_type": "text/csv",
        },
        "used_by_eos_records": records,
        "notes": note,
        "license": "Factual numerical transcription and vector-coordinate recovery with attribution. No explicit source data licence was found; source article and figures retain third-party rights. No licence grant over the publisher PDF is asserted.",
    }
    if identifier.endswith("_dewaele_2008_source_table"):
        material = identifier.removesuffix("_dewaele_2008_source_table")
        metadata["observation_views"] = [
            experimental_observation_view(material, records)
        ]
    return metadata


def base_record(material, rid, parameters, location, temperature):
    return {
        "identifier": rid,
        "label": "Dewaele et al. (2008), " + location,
        "model": "Vinet",
        "temperature_ref": temperature,
        "parameters": {
            k: {"value": float(p) * (Z[material] if k == "V0" else 1), "error": None}
            for k, p in zip(NAMES, parameters)
        },
        "parameter_error_confidence": None,
        "parameter_covariance": None,
        "record_kind": "published",
        "equation_kind": "isothermal",
        "reference": REF,
        "evidence": {
            "parameter_provenance": {
                k: f"{location}; atomic V0 multiplied by Z={Z[material]}"
                if k == "V0"
                else location
                for k in NAMES
            },
            "source_lineage": [
                {
                    "role": "equation and published coefficients",
                    "doi": REF["doi"],
                    "citation": f"Equation (1); {location}",
                }
            ],
        },
    }


def audit(record, outcome, location, status):
    finding = outcome.get("qualification", outcome.get("reason"))
    check = {
        "status": status,
        "audit_date": DATE,
        "note": finding,
        "verified_fields": [
            "equation",
            "parameters",
            "units",
            "reference_state",
            "phase",
            "published_uncertainties",
            "fixed_parameters",
        ],
        "primary_source_check": {
            "doi": REF["doi"],
            "access_url": "https://doi.org/" + REF["doi"],
            "locations": ["Equation (1)", location],
            "source_pdf_sha256": source()["source_pdf_sha256"],
            "finding": finding,
        },
        "primary_data_check": {
            "status": "bundled"
            if "table4" in record["identifier"]
            else "parameterization_only",
            "dataset_identifiers": outcome["dataset_identifiers"],
            "finding": finding,
        },
    }
    if status == "not_reproduced":
        record["catalog_access"] = "explicit_selection"
        check.update(
            usage_recommendation="not_recommended",
            reproduction_status="not_reproduced",
            recommendation_reason=finding,
            recommendation_date=DATE,
            recommendation_audit="docs/literature-reproductions/dewaele-2008-metals.md",
        )
    if status == "deferred":
        check["unresolved"] = finding
    record["evidence"]["scientific_validation"] = check
    assessment_row = {
        **outcome,
        **{
            key: check[key]
            for key in (
                "usage_recommendation",
                "recommendation_reason",
                "recommendation_date",
                "recommendation_audit",
            )
            if key in check
        },
    }
    record["evidence"]["reproduction"] = {
        "assessment": assessment_for_record(record, assessment_row)
    }
    if "curve_validation" in outcome:
        curve = outcome["curve_validation"]
        checks = curve["checks"]
        record["evidence"]["reproduction"]["source_equation_check"] = {
            "passed": curve["passed"],
            "source_curve_indices": sorted({c["source_curve_index"] for c in checks}),
            "checkpoint_count": len(checks),
            "max_absolute_volume_difference_a3_per_atom": max(
                (abs(c["volume_difference_a3"]) for c in checks), default=None
            ),
            "source_thermal_pressure_shift_gpa": curve[
                "source_thermal_pressure_shift_gpa"
            ],
            "method": curve["method"],
            "independent_paw_refit_available": False,
            "report": "docs/data/dewaele-2008-metals-reproduction.json",
        }


def experiment_record(material, rid, parameters, errors, fixed, outcome):
    constrained = "K0" in fixed
    location = f"Table IV, {'fourth' if constrained else 'third'} column, p. 104102-5"
    r = base_record(material, rid, parameters, location, 300.0)
    for name, error in zip(NAMES, errors):
        r["parameters"][name]["error"] = error
        if name in fixed:
            r["parameters"][name]["fixed"] = True
    r["parameter_error_confidence"] = 0.95
    r["determination_method"] = "experimental"
    r["evidence"]["parameter_provenance"]["determination_method"] = (
        "DAC observations on the 2004 revised ruby scale; V0 fixed from separate 0-5 GPa reduction"
        + ("; K0 fixed to adopted ultrasonic value." if constrained else ".")
    )
    r["fit_datasets"] = outcome["dataset_identifiers"]
    r["experimental_pressure_range_gpa"] = outcome["observed_pressure_range_gpa"]
    r["pressure_range_status"] = "reported_qualitatively"
    r["label"] = (
        f"Dewaele et al. (2008), Table IV Vinet, {'US K0=' + str(parameters[1]) + ' GPa fixed' if constrained else 'K0 fitted'}"
    )
    r["notes"] = (
        "Published Table IV coefficients, distinct from the later 2019 all-free-parameter recalibration. V0 is fixed after fitting low-pressure (0-5 GPa) data, not refitted in the full compression regression. "
        + (
            "K0 is adopted from ultrasound and fixed; its adopted uncertainty is not a fitted parameter error. "
            if constrained
            else "K0 and K0-prime are fitted. "
        )
        + "All printed rows are retained; run-2 Ag/Mo rows with revised pressure >100 GPa are excluded as stated in Section III. Parenthetical errors in Table IV are 95% widths; no covariance is given. The new 2008 ruby upgrade is an output of the later cross-calibration discussion, not the input pressure scale of these fits."
    )
    r["pressure_calibration"] = {
        "status": "partially_resolved"
        if material in {"nickel", "zinc_hcp"}
        else "resolved",
        "methods": [
            {
                "kind": "ruby_fluorescence",
                "reference_calibration_record": "ruby_dewaele_2004",
                "source_location": "Section IV, Tables II-IV: P_R-prime from reference 3",
            }
        ],
        "recalculation": {
            "status": "missing_calibrant_observations"
            if material in {"nickel", "zinc_hcp"}
            else "ready",
            "notes": "Revised-scale pressures are tabulated. Ruby-only pressures can be inverted to the corrected wavelength ratio. For parenthesized run-6 rows the pressure came from helium: raw He volumes are absent, so an exact mixed-gauge recalibration is unavailable. No wavelength or helium volume is synthesized.",
        },
        "audit_date": DATE,
    }
    if material in {"nickel", "zinc_hcp"}:
        r["pressure_calibration"]["methods"].append(
            {
                "kind": "equation_of_state",
                "reference": {
                    "authors": [
                        "Loubeyre",
                        "LeToullec",
                        "Hanfland",
                        "Ulivi",
                        "Datchi",
                        "Hausermann",
                    ],
                    "year": 1998,
                    "source": "Physical Review B",
                    "doi": "10.1103/PhysRevB.57.10403",
                },
                "source_location": "Table II caption, parenthesized run-6 rows; He Vinet V0=24.804 A3/atom, K0=0.165 GPa, K0-prime=7.483",
                "notes": "Exact source coefficients retained here as provenance; helium standard not registered as a material EOS. No raw paired He volumes are printed.",
            }
        )
    audit(r, outcome, location, "primary_source_validated")
    return r


def theory_record(material, rid, parameters, functional, state, outcome):
    location = f"Table VI, {material} {functional}, {'parenthesized static-lattice' if state == 'static' else 'unparenthesized 298 K'} values, p. 104102-10"
    r = base_record(
        material, rid, parameters, location, 0.0 if state == "static" else 298.0
    )
    r["label"] = (
        f"Dewaele et al. (2008), PAW {functional.upper()}, {'static lattice' if state == 'static' else '298 K'}, with semicore"
    )
    r["determination_method"] = "theoretical" if state == "static" else "hybrid"
    r["pressure_range_status"] = "theoretical"
    r["notes"] = (
        "Published Vinet fit to PAW calculations with semicore states (Tables V-VI). The parenthesized Table VI numbers are static-lattice coefficients, not errors. The unparenthesized values are independently tabulated fits after addition of the estimated constant 298 K thermal pressure, including zero-point motion; this record is one isotherm, not a continuous thermal EOS. Missing parameter errors and covariance remain unknown. Raw DFT fitting points are not printed; linked Figure 4 vertices are graphical checks, not fitting observations. The paper's experimental comparison rows do not define additional theoretical inputs or new experiments."
    )
    r["evidence"]["parameter_provenance"]["determination_method"] = (
        "PAW static-lattice calculation"
        + (
            "; source 298 K pressure correction uses adopted literature Debye estimates, hence hybrid."
            if state == "298k"
            else "."
        )
    )
    r["pressure_calibration"] = {
        "status": "not_applicable",
        "methods": [
            {
                "kind": "ab_initio",
                "source_location": "Sections V-VI, Tables V-VI",
                "scope": "PAW with semicore states; " + functional,
            }
        ],
        "recalculation": {
            "status": "not_applicable",
            "notes": "No experimental pressure gauge defines the PAW pressure. The 298 K correction is separately stated in Section V.",
        },
        "audit_date": DATE,
    }
    r["source_thermal_pressure_298k_gpa"] = source()["thermal_pressure_298k_gpa"][
        material
    ]
    passed = outcome["curve_validation"]["passed"]
    status = (
        "primary_source_validated"
        if passed
        else "not_reproduced"
        if outcome["observations"]
        else "deferred"
    )
    audit(r, outcome, location, status)
    return r


def correct_legacy(document, material):
    if material not in {"silver", "nickel"}:
        return
    records = {r["identifier"]: r for r in document["eos_records"]}
    old = records[f"{material}_dewaele_2008_vinet_1"]
    mao = records[f"{material}_dewaele_2019_mao_vinet"]
    dor = records[f"{material}_dewaele_2019_dor_vinet"]
    original_parameters = copy.deepcopy(old["parameters"])
    migration = copy.deepcopy(
        old["evidence"]["scientific_validation"].get("migration_source")
    )
    old.update(
        reference=copy.deepcopy(mao["reference"]),
        pressure_calibration=copy.deepcopy(mao["pressure_calibration"]),
        pressure_scale_variant_of=mao["identifier"],
        fit_datasets=copy.deepcopy(mao["fit_datasets"]),
        evidence=copy.deepcopy(mao["evidence"]),
    )
    old["parameters"] = original_parameters
    old["label"] = "Dewaele (2019), Mao-scale Vinet; historical 2008 identifier"
    old["notes"] = (
        "Compatibility lookup: this identifier historically attributed later 2019 Mao-scale coefficients to the 2008 paper and compared them with revised-scale observations. Its numerical coefficients and uncertainties are preserved exactly. Table IV of the 2008 paper is represented separately by the *_dewaele_2008_table4_* records. The source attribution, calibration and audit are corrected to the matching 2019 Table 1 Mao row; this is not an additional experiment."
    )
    old["evidence"]["scientific_validation"]["audit_date"] = DATE
    if migration:
        old["evidence"]["scientific_validation"]["migration_source"] = migration
    mao.pop("pressure_scale_variant_of", None)
    dor["pressure_scale_variant_of"] = mao["identifier"]
    for ds in document.get("datasets", []):
        if (
            ds["identifier"] == f"{material}_dewaele_2008_table2_compression"
            and material == "silver"
        ):
            ds["notes"] = (
                "Both printed pressure-scale columns are preserved. Table IV (2008) uses the revised 2004 scale, fixed V0 and an explicit run-2 exclusion; Table 1 (2019) refits the data with V0 free on the named Mao or Dorogokupets scale. Do not infer one universal fit coordinate from the dataset identifier. Common atomic-volume uncertainty 0.01 A3/atom; conditional pressure precision 0.05 GPa at 1 GPa to 2 GPa at 150 GPa, without rowwise errors."
            )


def main():
    # Initial integration predates the user's decision to separate PAW research
    # records. Do not let this historical command restore them to the database.
    from scripts.move_dewaele_2008_paw_to_research import COLLECTION

    if COLLECTION.exists():
        raise RuntimeError(
            "PAW records now live in the research collection. "
            "Use python -m scripts.update_dewaele_2008_assessments to refresh "
            "the 13 catalog and 54 research records separately."
        )
    report = reproduce()
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    documents = {
        m: json.loads((ROOT / f"peritheos/data/materials/{m}.eosmat").read_text())
        for m in Z
    }
    additions = {m: [] for m in Z}
    for m, rid, p, e, fixed in table4_specs():
        additions[m].append(
            experiment_record(m, rid, p, e, fixed, report["records"][rid])
        )
    for m, rid, p, f, state in table6_specs():
        additions[m].append(theory_record(m, rid, p, f, state, report["records"][rid]))
    with (DATA / "dewaele-2008-table2-complete.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    for m, document in documents.items():
        correct_legacy(document, m)
        owned_ids = {r["identifier"] for r in additions[m]}
        document["eos_records"] = [
            r for r in document["eos_records"] if r["identifier"] not in owned_ids
        ] + additions[m]
        if m in {"silver", "nickel"}:
            rid = f"{m}_dewaele_2008_table4_vinet"
            document["eos_records"].sort(key=lambda r: r["identifier"] != rid)
        datasets = []
        if m in source()["tables"]["IV"]:
            filename = (
                "dewaele-2008-table3-fe.csv"
                if m == "fe"
                else f"{m}-dewaele-2008-source-table2.csv"
            )
            if m != "fe":
                with (DATA / filename).open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                    writer.writeheader()
                    writer.writerows(r for r in rows if r["material"] == m)
            datasets.append(
                dataset(
                    f"{m}_dewaele_2008_source_table",
                    filename,
                    [
                        r["identifier"]
                        for r in additions[m]
                        if "table4" in r["identifier"]
                    ],
                    "Complete source table rows, including excluded observations and gauge flags; selection belongs to each fit.",
                    "Table III" if m == "fe" else "Table II",
                    "Source tokens, run identities, reported lattice data and mixed-gauge markings retained. Common Table II volume uncertainty 0.01 A3/atom is reported without an explicit confidence convention; pressure precision endpoints do not supply rowwise errors. No rowwise uncertainty array is synthesized.",
                )
            )
        datasets.append(
            dataset(
                f"{m}_dewaele_2008_figure4",
                f"{m}-dewaele-2008-figure4-vector.csv",
                [r["identifier"] for r in additions[m] if "paw" in r["identifier"]],
                "Figure 4 with-semicore static-lattice vector vertices; material and functional columns select the applicable source curve.",
                "Figure 4, p. 104102-9; Section V thermal-pressure estimates",
                "Graphical model-curve checkpoints, not raw PAW observations and not fit inputs. Volume and pressure reading halfwidths are diagnostic bounds, not standard deviations. The source's explicit constant PTH is applied only when checking the distinct 298 K coefficient set.",
            )
        )
        dsids = {d["identifier"] for d in datasets}
        document["datasets"] = [
            d for d in document.get("datasets", []) if d["identifier"] not in dsids
        ] + datasets
        validate_eosmat_document(document)
        (ROOT / f"peritheos/data/materials/{m}.eosmat").write_text(
            format_eosmat(document)
        )
    print(
        f"Integrated {sum(map(len, additions.values()))} source records; retained historical Ag/Ni coefficients."
    )


if __name__ == "__main__":
    main()
