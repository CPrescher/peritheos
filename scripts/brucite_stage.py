"""Build published-range brucite records and preserve separate audit snapshots."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json

import numpy as np

from scripts.brucite_reproduce import (
    COLUMNS,
    DATA,
    REPORT,
    ROOT,
    STAGED,
    TABLES,
    parse_number,
    reference_volume,
    source,
)
from scripts.eosmat_io import render_material


def staged_document():
    original = source()
    fukui = original["sources"]["fukui_2003"]
    structure = original["sources"]["structure"]
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    reference = {
        "authors": ["Fukui", "Ohtaka", "Suzuki", "Funakoshi"],
        "year": 2003,
        "title": fukui["title"],
        "source": "Physics and Chemistry of Minerals",
        "volume": "30",
        "locator": "511-516",
        "doi": fukui["doi"],
    }
    ids = [
        f"brucite_fukui_2003_{s['temperature_k']}k_bm3"
        for s in fukui["parameters"]["isotherms"]
    ]
    datasets = []
    for owner, table, filename in TABLES:
        ref = reference
        if owner == "xia_1998":
            ref = {
                "authors": ["Xia", "Weidner", "Zhao"],
                "year": 1998,
                "title": "Equation of state of brucite: Single-crystal Brillouin spectroscopy study and polycrystalline pressure-volume-temperature measurement",
                "source": "American Mineralogist",
                "volume": "83",
                "locator": "68-74",
                "doi": original["sources"][owner]["doi"],
            }
        columns = []
        for name in COLUMNS:
            unit = "1"
            if "temperature" in name:
                unit = "K"
            elif "pressure" in name:
                unit = "GPa"
            elif "lattice" in name:
                unit = "angstrom"
            elif "volume" in name:
                unit = "angstrom^3/conventional_unit_cell"
            quantity = name
            for prefix in ("temperature", "pressure", "cell_volume"):
                if name.startswith(prefix):
                    quantity = "volume" if prefix == "cell_volume" else prefix
            if name.startswith("lattice_a"):
                quantity = "lattice_a"
            if name.startswith("lattice_c"):
                quantity = "lattice_c"
            column = {
                "name": name,
                "quantity": quantity,
                "unit": unit,
                "role": "value",
            }
            if "error" in name:
                column["role"] = (
                    "standard_deviation" if owner == "fukui_2003" else "uncertainty"
                )
                column["of"] = name.replace("_error", "")
            columns.append(column)
        datasets.append(
            {
                "identifier": f"brucite_{owner}_{table}",
                "kind": "pressure_volume_temperature_lattice_parameters",
                "description": f"Every printed row of {owner} {table}; original errors retained.",
                "reference": ref,
                "source_location": original["sources"][owner][table]["location"],
                "source_url": original["sources"][owner]["url"],
                "columns": columns,
                "resource": {
                    "path": f"datasets/{filename}",
                    "media_type": "text/csv",
                    "sha256": report["datasets"][filename]["sha256"],
                },
                "used_by_eos_records": ids if owner == "fukui_2003" else [],
                "row_count": report["datasets"][filename]["rows"],
                "uncertainty": {
                    "notes": original["sources"][owner][table]["uncertainties"]
                },
                "license": "Source-data licence unspecified; contributor transcription/normalization/metadata rights only: CC0-1.0. See datasets/brucite.LICENSE.md.",
                "rights": original["rights"],
                "notes": "Blank means not reported or not applicable; it is never an invented zero. Table 1 has no numerical pressure column. Xia rows are recovered evidence for a held alternative, not inputs to the Fukui fits. Raw printed tokens are in datasets/brucite-source.json.",
            }
        )
    records = []
    for index, (params, assessment) in enumerate(
        zip(fukui["parameters"]["isotherms"], report["isotherms"])
    ):
        temperature = params["temperature_k"]
        records.append(
            {
                "identifier": ids[index],
                "label": f"Fukui et al. (2003), {temperature} K lattice-average BM3; staged, phase scope unresolved",
                "default": index == 0,
                "record_kind": "published",
                "equation_kind": "isothermal",
                "determination_method": "experimental",
                "reference": reference,
                "eos": {
                    "type": "BM3",
                    "model": "birch_murnaghan_3",
                    "parameters": {
                        "V0": float(reference_volume(temperature)),
                        "K0": params["K0_gpa"],
                        "K0_prime": params["K0_prime"],
                    },
                },
                "temperature_ref": temperature,
                "parameter_errors": {
                    "V0": 0.012 if temperature == 300 else None,
                    "K0": params["K0_error_gpa"],
                    "K0_prime": params["K0_prime_error"],
                },
                "parameter_error_confidence": None,
                "parameter_covariance": None,
                "fixed_parameters": ["V0"] if temperature != 300 else [],
                "parameter_provenance": {
                    "determination_method": "Powder X-ray isothermal compression; V0(T) from independent ambient expansion measurements. No acoustic or shock coefficients substituted.",
                    "equation": "Equation (2), p. 513: ordinary third-order Birch-Murnaghan.",
                    "V0": "Table 2: V0(300)=40.917(12) A^3. For high-T slices, p. 513 states V0 was fixed from ambient expansion using Eq. (1) and Table 3 constant alpha=6.69(18)e-5 K^-1. V0(T)=40.917 exp[6.69e-5(T-300)]. Static compression V0 treatment is not separately stated; Table 2 anchoring is explicit here and subject to integration review. Derived high-T V0 errors are missing because covariance is unreported, not zero.",
                    "K0": params["location"],
                    "K0_prime": params["location"],
                    "temperature_ref": "Exact listed isotherm; no interpolation or temperature extrapolation between records.",
                },
                "experimental_pressure_range_gpa": assessment["pressure_range_gpa"],
                "experimental_temperature_range_k": [temperature, temperature],
                "pressure_range_status": "reported_exactly",
                "validity": {
                    "pressure_gpa": [0, 22.56]
                    if temperature == 300
                    else assessment["pressure_range_gpa"],
                    "temperature_k": [temperature, temperature],
                    "notes": [
                        "Review-only numerical domain, not accepted thermodynamic phase stability. Records remain pending and require explicit opt-in to execute.",
                        "Historical staged source-equation snapshot. Current catalog accepts the layered lattice-average EOS across the source range, with H-order and source-normalization qualifications recorded separately.",
                        "No state at another temperature is validated by an isothermal record. No continuous EOS crosses the H transition or dehydration boundary.",
                    ],
                },
                "fit_datasets": [
                    "brucite_fukui_2003_table1",
                    "brucite_fukui_2003_table4",
                ],
                "fit_provenance": {
                    "software": {
                        "name": "not reported for EOS fitting",
                        "version": "not reported",
                    },
                    "dataset": "brucite_fukui_2003_table4",
                    "selection": {
                        "temperature_k": temperature,
                        "source_rows": assessment["source_rows"],
                        "source_exclusions": "not reported; every printed same-T row retained in diagnostics",
                    },
                    "objective": "not reported; diagnostics compare three explicit objectives",
                    "source_weights": "not reported",
                    "source_covariance": "not reported",
                    "source_software": "not reported for EOS fitting",
                    "refined_parameters": ["K0", "K0_prime"],
                    "fixed_parameters": ["V0"] if temperature != 300 else [],
                    "statistics": {"observations": assessment["observations"]},
                    "diagnostic_refits": assessment["diagnostic_refits"],
                    "diagnostic_refits_are_records": False,
                },
                "pressure_calibration": {
                    "status": "resolved",
                    "audit_date": "2026-10-16",
                    "methods": [
                        {
                            "kind": "equation_of_state",
                            "material": "Au",
                            "reference_eos_record": "gold_anderson_1989_bm3_1",
                            "reference": {
                                "authors": ["Anderson", "Isaak", "Yamamoto"],
                                "year": 1989,
                                "title": "Anharmonicity and the equation of state for gold",
                                "source": "Journal of Applied Physics",
                                "volume": "65",
                                "locator": "1534-1543",
                                "doi": "10.1063/1.342969",
                            },
                            "source_location": "High-pressure experiments, p. 512; reference p. 516",
                            "scope": "Every Table 4 pressure; Au mixed with brucite 1:10 by weight",
                        }
                    ],
                    "recalculation": {
                        "status": "missing_calibrant_observations",
                        "notes": "No paired Au lattice parameters or volumes tabulated; pressure errors alone cannot reconstruct them. No inverse-EOS synthetic marker observations are created.",
                    },
                },
                "experimental_configuration": {
                    "sample": "Synthetic commercial Mg(OH)2 reagent, 99.9% purity, Iwatani Chemical Ltd.; Au mixture under pressure.",
                    "apparatus": "SPEED-1500 multi-anvil; BL04B1 SPring-8 energy-dispersive diffraction; MgO medium, graphite sample capsule.",
                    "path": "Compressed, annealed at 900 K for over 10 minutes, collected on cooling.",
                    "temperature": fukui["pressure_calibration"]["temperature"],
                },
                "scientific_validation": {
                    "status": "pending_primary_source_check",
                    "audit_date": "2026-10-16",
                    "note": "Source equation and printed coefficients audited; retained as staged research candidate. Direct primary structural verification is resolved using He (2023). Broad phase-resolved acceptance and exact original fit workflow remain unresolved."
                    + (
                        " At 873 K the printed curve has an unexplained 0.486 GPa RMS offset and a maximum residual of 6.71 combined observational standard errors; conditional coefficient agreement does not resolve this."
                        if temperature == 873
                        else ""
                    ),
                    "primary_source_check": {
                        "doi": reference["doi"],
                        "access_url": fukui["url"],
                        "source_equation_validated": True,
                        "locations": ["Eqs. 1-2 and Tables 1-6, pp. 512-515"],
                    },
                    "primary_data_check": {
                        "status": "bundled",
                        "dataset_identifiers": [
                            "brucite_fukui_2003_table1",
                            "brucite_fukui_2003_table4",
                        ],
                        "finding": "All 48 sample observations, including negative reported pressure, repeated temperatures and errors preserved; original Au readings unavailable.",
                    },
                    "numerical_reproduction": {
                        "source_equation_validated": True,
                        "published_curve": assessment["published_curve"],
                        "conditional_coefficient_agreement": "Pressure-sigma diagnostic fits agree with printed K0 and Kprime within each reported error; objective itself is not source-confirmed.",
                        "author_workflow_reproduced": False,
                    },
                    "audit_report": "docs/literature-reproductions/brucite.md",
                },
                "notes": "Preserves source K0 and Kprime; high-T V0 derived transparently from the same source. No accepted independent refit. Table 5 prints 300 K Kprime error 0.3, Table 6 prints 0.4; Table 5 convention retained. Unreported coefficient covariance prevents treating propagated parameter errors as an independent validation statistic. Not a continuous thermal record and not an H-order/dehydration phase model.",
            }
        )
    a = float(parse_number(structure["cell"]["a_angstrom"])[0])
    c = float(parse_number(structure["cell"]["c_angstrom"])[0])
    return {
        "format": "peritheos.material",
        "format_version": 3,
        "identifier": "brucite",
        "name": "Brucite",
        "formula": "Mg(OH)2",
        "phase": "Layered brucite; ambient ordered-H P -3 m 1 structural reference",
        "aliases": ["Mg(OH)2", "magnesium hydroxide brucite"],
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        },
        "cell_contents": "One Mg(OH)2 formula unit, five atoms, in the primitive hexagonal setting of trigonal P -3 m 1. Conventional and primitive cells coincide; all EOS volumes use Z=1.",
        "symmetry": "TRIGONAL",
        "space_group": "P -3 m 1",
        "space_group_number": 164,
        "formula_units_per_cell": 1,
        "lattice": {"a": a, "b": a, "c": c, "alpha": 90, "beta": 90, "gamma": 120},
        "atom_sites": [
            {
                "label": "Mg1",
                "element": "Mg",
                "wyckoff": "1a",
                "x": 0,
                "y": 0,
                "z": 0,
                "occupancy": 1,
            },
            {
                "label": "O1",
                "element": "O",
                "wyckoff": "2d",
                "x": 1 / 3,
                "y": 2 / 3,
                "z": 0.2189,
                "occupancy": 1,
            },
            {
                "label": "H1",
                "element": "H",
                "wyckoff": "2d",
                "x": 1 / 3,
                "y": 2 / 3,
                "z": 0.4108,
                "occupancy": 1,
            },
        ],
        "peaks": [],
        "source": {
            "structure_reference": structure["reference"],
            "structure_source_kind": structure["source_kind"],
            "structure_access_url": structure["access_url"],
            "structure_dataset": "datasets/brucite-he-2023-table43.cif",
            "structure_dataset_sha256": hashlib.sha256(
                (DATA / "brucite-he-2023-table43.cif").read_bytes()
            ).hexdigest(),
            "structure_metric_volume_a3": float(np.sqrt(3) / 2 * a**2 * c),
            "structure_printed_volume_a3": 40.9437,
            "structure_errors": {
                "a_angstrom": 0.00003,
                "c_angstrom": 0.00005,
                "O_z": 0.0002,
                "H_z": 0.0004,
            },
            "structure_refinement": structure["refinement"],
            "structure_atomic_displacements_a2": structure["atomic_displacements_a2"],
            "structure_normalization": structure["normalization"],
            "structure_cell_vs_eos": "Ambient structural cell has a different sample/reference from Fukui V0=40.917(12) A^3; no rescaling applied.",
            "structure_acceptance_limit": structure["acceptance_limit"],
            "source_manifest": "datasets/brucite-source.json",
        },
        "integration_status": "Historical staged source-equation snapshot. The separate current catalog card includes published-range lattice-average curves and independently attributed table-data refits.",
        "notes": "Four independently listed isothermal source curves; no thermal model inferred by freezing Kprime. 102 factual observations are recovered. All source-parameter records remain pending. No first-order transition, high-pressure H-order phase, or dehydration product is combined into an accepted EOS. Scalar volume scaling is not an anisotropic high-pressure structural refinement.",
        "datasets": datasets,
        "eos_records": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    for path, document in (
        (STAGED, staged_document()),
        (ROOT / "peritheos/data/materials/brucite.eosmat", catalog_document()),
        (STAGED.with_name("brucite-diagnostic-archive.eosmat"), audit_document()),
        (
            STAGED.with_name("brucite-conditional-thermal.eosmat"),
            conditional_thermal_document(),
        ),
    ):
        payload = (
            render_material(document)
            if path == ROOT / "peritheos/data/materials/brucite.eosmat"
            else json.dumps(document, indent=2, allow_nan=False) + "\n"
        )
        if args.check:
            if path.read_text(encoding="utf-8") != payload:
                raise SystemExit(f"Staged brucite card is outdated: {path}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload, encoding="utf-8")
        print(path.relative_to(ROOT))


def conditional_thermal_document():
    """Add an opt-in audit surface without reclassifying it as a source fit."""
    document = staged_document()
    record = copy.deepcopy(document["eos_records"][0])
    record.update(
        identifier="brucite_fukui_2003_conditional_thermal",
        label="Fukui (2003) slopes with assumed 300 K intercepts; conditional audit surface",
        default=False,
        record_kind="diagnostic",
        equation_kind="thermal",
        derived_from_record="brucite_fukui_2003_300k_bm3",
        fit_datasets=[],
        experimental_pressure_range_gpa=[-0.06, 22.56],
        experimental_temperature_range_k=[300, 873],
        thermal={
            "type": "AlphaKT",
            "model": "thermal_reference_state",
            "thermal_expansion_law": "constant",
            "reference_volume_law": "integrated_expansivity",
            "parameters": {
                "Tr": 300.0,
                "alpha0": 6.69e-5,
                "dK_dT": -0.0319,
                "dKprime_dT": 0.00571,
            },
            "parameter_errors": {
                "Tr": None,
                "alpha0": 0.18e-5,
                "dK_dT": 0.0029,
                "dKprime_dT": 0.00045,
            },
            "fixed_parameters": ["Tr"],
        },
        validity={
            "pressure_gpa": [5.95, 15.48],
            "temperature_k": [300, 873],
            "notes": [
                "Conservative research evaluation domain only; no validated phase-stability rectangle is established.",
                "Source temperature-slope intercepts are unreported. This diagnostic assumes the separately printed 300 K slice supplies them; it is not an exact published surface or independent refit.",
            ],
        },
        notes="Conditional audit surface enabled by dKprime_dT support. Uses printed source temperature slopes and assumed 300 K intercepts. No regression performed; no new accepted phase-spanning EOS. Broad phase scope and source observation holds remain; the ambient primary structure check is resolved.",
    )
    record.pop("fit_provenance")
    record["parameter_provenance"]["thermal"] = (
        "Fukui p. 513 reports dK/dT=-0.0319(29) GPa/K and dKprime/dT=0.00571(45) K^-1. "
        "The 300 K values K0=41.8 GPa and Kprime=6.6 are adopted as assumed intercepts; the author regression's separately fitted intercepts are unreported."
    )
    record["scientific_validation"].update(
        note="Diagnostic conditional surface, not a validated published thermal parameterization. Runtime mapping and the ambient primary structure are verified; source line intercepts, broad phase scope and isotherm residual questions remain unresolved.",
        numerical_reproduction={
            "independent_equation_checkpoint": {
                "volume_a3": 34,
                "temperature_k": 873,
                "pressure_gpa": 14.051454515344124,
            },
            "source_intercepts_reproduced": False,
            "author_workflow_reproduced": False,
        },
    )
    document["eos_records"].append(record)
    for dataset in document["datasets"]:
        if "fukui" in dataset["identifier"]:
            dataset["used_by_eos_records"].append(record["identifier"])
    document["notes"] += (
        " This variant also includes the explicitly diagnostic conditional thermal surface; the original four-isotherm staged card remains separate."
    )
    return document


def audit_document():
    """Preserve published coefficients and historical refits as audit evidence."""
    document = staged_document()
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    phase_note = (
        "Layered brucite lattice-average empirical volume EOS, including H rearrangement; "
        "not a hydrogen-order-resolved thermodynamic fit or a phase-boundary model. "
        "Fukui describes a higher-order H change; Sun (2025) XRD supports the layered "
        "average lattice beyond this range. Ambient atomic sites are a reference only."
    )
    records = document["eos_records"]
    for record, assessment in zip(records, report["isotherms"]):
        temperature = assessment["temperature_k"]
        record["label"] = f"Fukui et al. (2003), {temperature} K lattice-average BM3"
        record["validity"] = {
            "pressure_gpa": [0, 22.56]
            if temperature == 300
            else assessment["pressure_range_gpa"],
            "temperature_k": [temperature, temperature],
            "notes": [
                phase_note,
                "Measured pressure span at this temperature; no temperature interpolation by an isothermal record.",
            ],
        }
        validation = record["scientific_validation"]
        validation.update(
            status="primary_source_validated"
            if temperature != 873
            else "not_reproduced",
            note="Primary equation, coefficients, units, uncertainties and source PVT observations directly checked. Published parameters preserved; exact author regression is not reproduced. "
            + phase_note,
            verified_fields=[
                "equation",
                "parameters",
                "units",
                "reference_state",
                "published_uncertainties",
                "validity",
                "source_phase_annotations",
            ],
        )
        validation["primary_source_check"]["finding"] = (
            "BM3 and published moduli transcribed from Tables 5-6; V0 from Eq. (1) and Tables 2-3. "
            "Independent source observations checked; author weights and covariance unreported."
        )
        record["notes"] = (
            phase_note
            + " Published coefficients unchanged. Original workflow and calibration re-reduction remain unavailable."
        )
        if temperature == 873:
            record["catalog_access"] = "explicit_selection"
            validation["note"] += (
                " Printed Table 4/6 and Fig. 6(c) normalization are inconsistent: the plotted markers "
                "imply V0 about 42.70 A^3 rather than the printed expansion's 42.516 A^3. "
                "Original source curve has 0.486 GPa RMS residual; select the independent table refit for data-based use."
            )
            validation["numerical_reproduction"]["figure6c_audit"] = report[
                "figure6c_audit"
            ]
    thermal = copy.deepcopy(conditional_thermal_document()["eos_records"][-1])
    thermal.update(
        identifier="brucite_fukui_2003_300k_anchored_thermal",
        label="Fukui (2003), published slopes with 300 K anchors; qualified reconstruction",
        catalog_access="explicit_selection",
        validity={
            "pressure_gpa": [5.95, 15.48],
            "temperature_k": [300, 873],
            "notes": [
                phase_note,
                "Common pressure interval of the four measured isotherms; individual source isotherms retain their larger measured ranges. No ambient high-temperature or dehydration extrapolation accepted.",
                "Literal source-defined 300 K anchors; independently regressed intercepts and exact original surface unrecovered.",
            ],
        },
        notes="Published slopes and 300 K anchors preserved. Original independently regressed intercepts are unrecovered; not a refit. "
        + phase_note,
    )
    thermal["scientific_validation"].update(
        status="not_reproduced",
        note="Qualified literal reconstruction, not recovery of the author's fitted continuous surface. Source normalization inconsistency remains recorded; separately attributed Peritheos fits are available.",
        usage_recommendation="Explicit source-comparison selection only within the common measured pressure interval. Prefer the independent thermal refit for reproducible table-data interpolation; neither model is a phase-equilibrium calculation or pressure standard.",
    )
    records.append(thermal)
    for fit_key, identifier, temperature in (
        ("isotherm_873k", "brucite_peritheos_fukui_2003_873k_bm3_refit", 873),
        ("continuous_thermal", "brucite_peritheos_fukui_2003_thermal_refit", None),
    ):
        fit = report["peritheos_refits"][fit_key]
        record = copy.deepcopy(records[3] if temperature else thermal)
        record.pop("catalog_access", None)
        # These regressions depend on source observations and the explicit
        # reference-volume law, not on execution of either qualified curve.
        record.pop("derived_from_record", None)
        record.update(
            identifier=identifier,
            label="Peritheos unweighted fit to Fukui (2003) Table 4, "
            + ("873 K BM3" if temperature else "continuous thermal BM3"),
            default=False,
            record_kind="refit",
            fixed_parameters=["V0"],
            derivation={
                "source_kind": "published_table",
                "source_identifier": "brucite_fukui_2003_table4",
                "method": report["peritheos_refits"]["objective"],
            },
            fit_datasets=["brucite_fukui_2003_table1", "brucite_fukui_2003_table4"],
            parameter_covariance=fit["parameter_covariance"],
            notes="New Peritheos regression to printed observations, not an author fit or correction of published coefficients. No figure curve or digitized marker is fitted. "
            + phase_note,
        )
        record["eos"]["parameters"].update(
            {k: fit["parameters"][k] for k in ("K0", "K0_prime")}
        )
        record["parameter_errors"].update(
            {k: fit["parameter_errors"][k] for k in ("K0", "K0_prime")}
        )
        record["parameter_provenance"]["determination_method"] = (
            "Peritheos independent regression to primary powder-XRD PVT observations; "
            "published expansion law supplies fixed reference volumes."
        )
        record["parameter_provenance"]["parameter_errors"] = (
            "One conditional OLS standard error from residual-scaled covariance. "
            "No exact finite-sample confidence interval is asserted; fixed-input "
            "uncertainty and calibrant systematics are excluded."
        )
        for key in ("K0", "K0_prime"):
            record["parameter_provenance"][key] = (
                "Peritheos independent unweighted pressure-residual regression to Fukui Table 4, with printed source V0 and alpha fixed; not a published coefficient."
            )
        if not temperature:
            record["thermal"]["fixed_parameters"] = ["Tr", "alpha0"]
            record["parameter_covariance"] = copy.deepcopy(fit["parameter_covariance"])
            record["parameter_covariance"]["parameter_order"] = [
                "rt_eos.K0",
                "rt_eos.K0_prime",
                "dK_dT",
                "dKprime_dT",
            ]
            record["parameter_provenance"]["temperature_ref"] = (
                "300 K anchors for the jointly fitted linear K0(T)/Kprime(T) laws; "
                "interpolation over 300-873 K."
            )
            for key in ("dK_dT", "dKprime_dT"):
                record["thermal"]["parameters"][key] = fit["parameters"][key]
                record["thermal"]["parameter_errors"][key] = fit["parameter_errors"][
                    key
                ]
            record["parameter_provenance"]["thermal"] = (
                "Peritheos jointly regressed K0(300), Kprime(300) and both linear slopes on all 37 Table 4 observations. Source V0(300) and alpha fixed. No source slope or intercept reproduction claimed."
            )
        record["fit_provenance"] = {
            "software": {
                "name": "Peritheos scripts/brucite_reproduce.py; scipy.optimize.least_squares",
                "version": report["software"]["scipy"],
            },
            "dataset": "brucite_fukui_2003_table4",
            "selection": {
                "source_rows": report["isotherms"][3]["source_rows"]
                if temperature
                else list(range(1, 38)),
                "exclusions": [],
            },
            "objective": report["peritheos_refits"]["objective"],
            "refined_parameters": list(fit["parameters"]),
            "fixed_parameters": ["V0", "alpha0", "Tr"],
            "covariance_scaling": fit["covariance_scaling"],
            "statistics": {
                "observations": 8 if temperature else 37,
                **{
                    k: fit[k]
                    for k in (
                        "degrees_of_freedom",
                        "pressure_rmse_gpa",
                        "maximum_absolute_pressure_residual_gpa",
                    )
                },
            },
            "reproduction": {
                "report": "docs/data/brucite-reproduction.json",
                "author_workflow_reproduced": False,
            },
        }
        record["scientific_validation"].update(
            status="primary_source_validated",
            note="Independent Peritheos refit: primary observations and fixed reference-volume law verified, optimizer and covariance reproducible. Coefficients belong to Peritheos, not Fukui. "
            + phase_note,
            numerical_reproduction=fit,
            usage_recommendation="Use within the represented measured/interpolation domain. Conditional model-fit covariance excludes source V0/alpha uncertainty and calibrant systematics; residual scatter is not pressure-standard accuracy.",
        )
        record["scientific_validation"]["primary_source_check"]["finding"] = (
            "All selected printed Table 4 observations retained with original errors; source Eq. (1)/Tables 2-3 supply fixed V0(T). Optimized coefficients are explicitly attributed to Peritheos."
        )
        if not temperature:
            record["notes"] += (
                " Volume interpolation only: thermal derivatives are not physically validated. "
                "The fit predicts weak negative expansivity near the upper P/T corner "
                "(about -4.50e-6 K^-1 at 15.48 GPa/873 K); this is a model consequence, "
                "not an observed anomaly or evidence of a phase transition."
            )
            record["fit_provenance"]["physical_validation"] = {
                "thermal_derivatives_validated": False,
                "corner_diagnostic": {
                    "pressure_gpa": 15.48,
                    "temperature_k": 873,
                    "volumetric_expansivity_per_k_approx": -4.50e-6,
                },
            }
        records.append(record)
    from scripts.brucite_curved_expansion import build_records

    records.extend(build_records(records))
    installed = [r["identifier"] for r in records]
    for dataset in document["datasets"]:
        dataset["used_by_eos_records"] = (
            installed if "fukui" in dataset["identifier"] else []
        )
    document["phase"] = (
        "Layered brucite; lattice-average volume EOS including H rearrangement"
    )
    document["integration_status"] = (
        "Published-range lattice-average isotherms and independent table-data refits installed; original 873 K curve and source-slope reconstruction remain explicitly qualified."
    )
    document["notes"] = (
        phase_note
        + " Direct ambient average F-free hydrogenated structure from He (2023), Table 4.3. Source structure and EOS samples differ; no rescaling. All 102 Fukui/Xia observations preserved."
    )
    # The source manifest retains the historical investigation state. This
    # current curation decision supersedes its earlier integration proposal.
    document["phase_scope"] = copy.deepcopy(source()["catalog_disposition"])
    document["phase_scope"].update(
        held_isotherms_k=[873],
        thermal_status="Not recommended: published thermal-EOS parity not achieved. All refits are historical diagnostics only.",
        source_873k_status="Not reproduced; unchanged coefficients archived outside executable selection.",
        thermal_interpolation_pressure_gpa=None,
        qualified_thermal_temperature_k=None,
    )
    document["phase_references"] = [
        source()["sources"][key]
        for key in ("duffy_1995_phase", "wu_2026_phase", "sun_2025_phase")
    ]
    document["curation_disposition"] = {
        "date": "2026-10-17",
        "published_thermal_eos_parity": "not_achieved",
        "thermal_eos_recommendation": "not_recommended",
        "reason": "The printed expansion law, 873 K data/figure normalization and continuous source reconstruction are not mutually reproduced. No source correction or validated replacement has been established.",
        "research_scope": "Further model development is outside the Peritheos curation scope. Existing refits are historical diagnostics only.",
        "retained_catalog_isotherms_k": [300, 473, 673],
    }
    for record in records:
        validation = record["scientific_validation"]
        record.pop("catalog_access", None)
        if (
            validation["status"] == "primary_source_validated"
            and record.get("record_kind") != "refit"
        ):
            validation["usage_recommendation"] = (
                "Source-checked isotherm at its stated temperature only. This does not validate the full published thermal EOS or an interpolation between temperatures."
            )
            continue
        is_refit = record.get("record_kind") == "refit"
        validation.update(
            status="deferred" if is_refit else "not_reproduced",
            audit_date="2026-10-17",
            usage_recommendation="Not recommended. Historical audit evidence only; excluded from the executable catalog.",
            note="Diagnostic Peritheos refit only; numerical replay does not establish published-EOS parity or a validated replacement."
            if is_refit
            else "Published-EOS parity not achieved. The 873 K normalization and continuous source reconstruction remain unresolved; published coefficients are preserved without correction.",
        )
        record["notes"] += (
            " Current disposition: not recommended; historical audit evidence only."
        )
    document["integration_status"] = (
        "Audit archive: full thermal EOS not recommended; unresolved source curves and all Peritheos refits withheld from the executable catalog."
    )
    return document


def catalog_document():
    """Bundle only the three individually source-checked published isotherms."""
    document = audit_document()
    document["eos_records"] = [
        r
        for r in document["eos_records"]
        if r["scientific_validation"]["status"] == "primary_source_validated"
    ]
    installed = [r["identifier"] for r in document["eos_records"]]
    for dataset in document["datasets"]:
        dataset["used_by_eos_records"] = (
            installed if "fukui" in dataset["identifier"] else []
        )
    document["integration_status"] = (
        "Only published 300/473/673 K isotherms installed. Full thermal EOS not recommended: published-data parity not achieved. Original 873 K/continuous coefficients and all refits retained in the diagnostic archive."
    )
    return document


if __name__ == "__main__":
    main()
