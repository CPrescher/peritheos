#!/usr/bin/env python3
"""Build 200 metal EOS records dominated by experimental compression fits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any

from build_deltacodesdft_metal_eos import EXISTING_PHASES, NEW_PHASES

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"
DATASET = (
    ROOT / "peritheos" / "data" / "datasets" / "experimental-metal-eos-batch-200.csv"
)
DORFMAN_REFIT = ROOT / "docs" / "data" / "dorfman-2012-cocompression-refit.json"
DORFMAN_DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "dorfman-2012-tables-s1-s6-cocompression.csv"
)
DORFMAN_DATASET_ID = "dorfman_2012_tables_s1_s6_cocompression"
AUDIT_DATE = "2026-09-06"
DORFMAN_AUDIT_DATE = "2026-09-07"
CM3_MOL_TO_A3_FORMULA = 1.6605390671738466

DELTA_EXPERIMENT = """element,V0,K0,K0p
Li,20.4119306779,13.0656516899,3.510
Na,37.1751953558,7.9126957362,4.130
Mg,22.6471798130,38.6053690056,4.800
Al,16.2686134104,77.1357712089,4.450
K,72.1993569140,3.7509702193,4.089
Ca,42.9453933306,15.8522348531,3.100
Sc,24.7695611521,44.5429241793,2.800
Ti,17.4972491019,107.3532076344,3.400
V,13.8102371627,165.8391206631,4.135
Cr,11.8218181990,204.5647000914,6.895
Mn,11.9680120678,174.7008209027,6.600
Fe,11.6393714028,175.1086987963,4.600
Co,10.9563191311,198.3952571603,4.260
Ni,10.8077248858,192.4580551235,4.000
Cu,11.6468691265,144.2785999228,4.880
Zn,14.8589080560,64.6681064605,4.400
Rb,89.1556853499,3.5526213033,3.885
Sr,55.6036124476,11.9842858226,2.485
Y,32.9504356012,37.3135136902,2.200
Zr,23.1775566622,84.2834017272,2.575
Nb,17.9728074489,173.1796015761,4.015
Mo,15.5081291945,276.2110604310,3.980
Ru,13.4494319310,335.5206725035,6.610
Rh,13.5713064752,277.1283701334,4.500
Pd,14.5649181140,187.1936487566,5.000
Ag,16.8502671224,105.7066536338,4.725
Cd,21.0072915257,50.6575327406,4.900
In,25.7032021592,44.7057199504,5.350
Sn,33.9673514422,42.8175727666,4.000
Cs,110.3240829143,2.3449684759,3.790
Ba,62.2874792710,10.5948206033,2.430
Hf,22.2975348853,110.6674189102,3.950
Ta,17.9333133566,202.6879051772,3.750
W,15.7959730837,327.4742143494,4.320
Re,14.6193984266,380.7610091505,5.410
Os,13.8455627662,424.6083092316,4.500
Ir,14.0640862983,362.2372836709,4.830
Pt,15.0221694263,285.5122223387,5.180
Au,16.8216099724,182.0053319987,6.400
Tl,27.9937800663,37.4420913632,3.000
Pb,29.8632199952,46.3406063802,5.335
Bi,35.1255842981,32.0779415682,2.400
"""

SUN_LOW = """solid,pmax,VN_K0,VN_K0p,BN_K0,BN_K0p,MRS3_K0,MRS3_K0p,SMS3_K0,SMS3_K0p,SMS4_K0,SMS4_K0p
Cu,34,140.95,4.798,141.46,4.572,141.47,4.691,141.47,4.711,141.38,4.722
Mo,60,269.51,3.857,270.55,3.742,267.24,4.077,267.82,4.039,267.32,4.027
W,70,314.99,3.802,313.00,3.833,312.82,3.954,312.40,3.993,312.73,3.949
Zn,16,60.204,5.680,60.615,5.475,60.237,5.715,60.252,5.712,60.298,5.679
Ag,28,105.50,5.564,105.65,5.537,105.32,5.742,105.61,5.695,105.39,5.692
Pt,70,281.42,4.879,282.09,4.954,281.27,5.115,281.18,5.132,281.36,5.100
Ti,22,97.788,3.842,98.005,3.582,97.766,3.680,97.868,3.716,97.806,3.705
Ta,48,197.60,4.012,197.75,3.706,197.68,3.810,197.64,3.809,197.20,3.842
Au,46,184.07,5.029,184.64,4.983,184.74,5.118,184.30,5.165,184.49,5.125
Pd,70,193.37,5.565,196.03,5.016,195.18,5.227,194.94,5.250,195.25,5.209
Zr,20,95.913,2.339,95.614,2.527,95.254,2.691,95.285,2.719,95.173,2.740
Cr,46,190.63,4.586,190.49,4.851,191.54,4.858,190.49,4.979,191.22,4.887
Co,46,196.40,4.145,196.75,4.175,195.97,4.359,196.27,4.359,195.86,4.352
Ni,46,188.32,4.614,188.94,4.610,188.59,4.752,188.35,4.803,188.72,4.730
Nb,38,167.28,3.447,169.99,3.704,168.00,4.093,168.31,4.142,167.90,4.104
Cd,14,50.959,4.704,50.798,5.423,50.767,5.554,50.704,5.592,50.726,5.553
Al,18,74.762,5.036,75.802,4.766,75.416,4.945,75.406,4.956,75.516,4.902
Th,12,52.062,3.920,52.181,4.067,52.112,4.171,52.054,4.215,52.124,4.164
V,36,159.80,3.449,159.99,3.510,159.66,3.626,159.50,3.670,159.82,3.614
In,1,40.064,5.017,40.228,4.913,39.727,5.347,39.946,5.280,39.830,5.289
Be,26,119.97,3.434,120.42,3.429,120.25,3.530,120.16,3.555,120.24,3.520
Pb,12,44.026,4.996,44.223,5.066,44.005,5.277,44.066,5.256,44.117,5.217
Sn,12,43.636,5.338,43.957,5.114,43.841,5.283,43.826,5.296,43.847,5.266
Mg,8,34.527,4.049,34.669,3.834,34.614,3.9626,34.586,3.998,34.612,3.957
Ca,4,19.480,2.532,19.397,2.478,19.336,2.579,19.299,2.632,19.336,2.585
Tl,18,35.268,5.551,35.570,5.404,35.362,5.663,35.336,5.697,35.434,5.613
"""

SUN_V0 = {
    row["solid"]: float(row["v0_cm3_mol"])
    for row in csv.DictReader(
        io.StringIO(
            (DATASET.parent / "sun-2010-table2-morse-eos-parameters.csv").read_text()
        )
    )
}
SUN_MODELS = {
    "VN": ("Vinet", "vinet"),
    "MRS3": ("Morse3", "morse_3"),
    "SMS3": ("SunMorse3", "sun_morse_3"),
    "SMS4": ("SunMorse4", "sun_morse_4"),
    "BN": ("Baonza", "baonza"),
}
SUN_BAONZA_PRIORITY = {"Cu", "Mo", "W", "Ag", "Pt", "Ta", "Au"}

DEWAELE_2019 = """solid,material,V0_mao,K0_mao,K0p_mao,eV0_mao,eK0_mao,eK0p_mao,V0_dor,K0_dor,K0p_dor,pmin,pmax
Au,gold,16.983,166.4,5.47,0.020,2.0,0.06,16.986,163.4,6.04,0,131
Pt,platinum,15.099,273.4,4.83,0.025,2.5,0.08,15.098,270.8,5.50,0,95
Cu,copper,11.810,135.3,4.91,0.015,1.5,0.06,11.81,133.1,5.38,0,155
Ta,tantalum,18.020,197.9,3.17,0.018,3.7,0.10,18.019,196.1,3.64,0,90
Al,aluminum,16.573,76.32,4.16,0.019,1.5,0.06,16.584,74.2,4.52,0,155
W,tungsten,15.862,298.3,3.82,0.016,4.1,0.11,15.858,298.6,4.37,0,155
Co,cobalt_hcp,11.077,197.0,3.85,0.012,3.2,0.20,11.077,194.85,4.36,0,66
Ag,silver,17.070,100.2,5.70,0.016,1.6,0.09,17.088,96.6,6.22,0,124
Mo,molybdenum,15.569,270.3,3.34,0.021,3.9,0.12,15.567,269.3,3.87,0,124
Ni,nickel,10.954,177.5,4.83,0.018,2.4,0.09,10.952,176.0,5.322,0,157
Zn,zinc_hcp,15.147,64.3,5.30,0.019,1.2,0.10,15.155,62.2,5.705,0,157
Be,beryllium_hcp,8.133,115.2,2.94,0.005,1.1,0.05,8.134,113.4,3.29,0,95
Pb,lead_hcp,28.063,71.8,4.40,0.055,2.1,0.08,28.058,70.0,4.77,13,131
Re,rhenium,14.737,350.5,3.98,0.020,8.0,0.17,14.734,350.5,4.62,0,144
Fe,iron,11.209,164.5,4.96,0.050,7.9,0.16,11.177,168.4,5.33,17,204
"""

DEWAELE_2004 = """solid,material,scale,V0,K0,K0p,eK0,eK0p,fixed
Au,gold,mao_ruby,16.962,172.5,5.40,1.4,0.08,
Pt,platinum,mao_ruby,15.095,275.3,4.78,2.0,0.08,
Ta,tantalum,mao_ruby,18.035,198.2,3.07,3.1,0.14,
W,tungsten,mao_ruby,15.862,298.3,3.81,3.6,0.10,
Cu,copper,mao_ruby,11.810,135.1,4.91,1.1,0.05,
Al,aluminum,mao_ruby,16.573,76.3,4.16,1.1,0.05,
Pt,platinum,revised_ruby,15.095,277,5.08,,0.02,K0
Ta,tantalum,revised_ruby,18.035,194,3.52,,0.03,K0
"""

DORFMAN_2012 = """solid,material,V0,K0_fixed,K0p_fixed,eK0p_fixed,K0_free,eK0_free,K0p_free,eK0p_free,pmax
Au,gold,67.85,167,5.88,0.02,167,3,5.84,0.02,259
Mo,molybdenum,31.17,261,4.19,0.02,271,3,3.89,0.09,207
Pt,platinum,60.38,277,5.43,0.02,280,2,5.29,0.07,226
"""

THERMAL_2025 = """fit,type,model,V0_cm3_mol,K0,K0p,eV0,eK0,eK0p,theta0,etheta0,gamma0,egamma0,q,eq,rmse
1,BM3,birch_murnaghan_3,6.756,174.7,4.790,0.010,1.7,0.014,1209,73,2.86,0.10,0.84,0.05,4.50
2,BM3,birch_murnaghan_3,6.753,175.1,4.787,,0.5,0.009,1205,73,2.86,0.10,0.84,0.04,4.50
5,Vinet,vinet,6.753,151.6,5.845,,0.5,0.010,960,69,3.51,0.13,1.28,0.05,5.08
"""

ELEMENT_TO_MATERIAL = {
    **EXISTING_PHASES,
    **{element: values[0] for element, values in NEW_PHASES.items()},
    "Mn": "manganese_alpha",
}


def rows(text: str) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(text.strip())))


def reference(
    authors: list[str], year: int, title: str, source: str, doi: str
) -> dict[str, Any]:
    return {
        "authors": authors,
        "year": year,
        "title": title,
        "source": source,
        "doi": doi,
    }


def dorfman_dataset(
    material: str, dor_ref: dict[str, Any], dor_refit: dict[str, Any]
) -> dict[str, Any]:
    record_ids = [
        f"{material}_dorfman_2012_tange_mgo_k0_{mode}_vinet"
        for mode in ("fixed", "free")
    ]
    columns = [
        ("observation_index", "observation_index", "dimensionless", "flag"),
        ("source_table", "source_table", "dimensionless", "flag"),
        ("source_pdf_page", "source_page", "dimensionless", "flag"),
        ("run", "experiment_identifier", "dimensionless", "flag"),
        ("row_in_run", "measurement_sequence", "dimensionless", "flag"),
        ("material", "material_identifier", "dimensionless", "flag"),
        ("volume_source_token", "source_token", "dimensionless", "flag"),
        (
            "volume_a3_conventional_cell",
            "volume",
            "angstrom^3/conventional_unit_cell",
            "value",
        ),
        (
            "volume_standard_error_a3",
            "volume",
            "angstrom^3/conventional_unit_cell",
            "standard_error",
        ),
        (
            "reported_pressure_source_token",
            "source_token",
            "dimensionless",
            "flag",
        ),
        ("reported_pressure_gpa", "pressure", "GPa", "value"),
        ("single_peak_mgo", "single_peak_indicator", "dimensionless", "flag"),
    ]
    column_metadata = [
        {"name": name, "quantity": quantity, "unit": unit, "role": role}
        for name, quantity, unit, role in columns
    ]
    column_metadata[8]["of"] = "volume_a3_conventional_cell"
    return {
        "identifier": DORFMAN_DATASET_ID,
        "kind": "simultaneous_unit_cell_volumes",
        "description": (
            "All 368 non-missing conventional-cell volume measurements from 165 "
            "aligned observations in auxiliary Tables S1-S6."
        ),
        "reference": dor_ref,
        "source_location": "Auxiliary Tables S1-S6",
        "source_url": dor_refit["source"]["supporting_url"],
        "license": "CC0-1.0",
        "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
        "license_scope": (
            "Peritheos-created factual CSV transcription, normalization, column naming, "
            "and arrangement, solely to the extent contributors hold rights; excludes "
            "the article, publisher PDF, and all third-party rights."
        ),
        "columns": column_metadata,
        "resource": {
            "path": "datasets/dorfman-2012-tables-s1-s6-cocompression.csv",
            "sha256": hashlib.sha256(DORFMAN_DATASET.read_bytes()).hexdigest(),
            "media_type": "text/csv",
        },
        "provenance": {
            "type": "lossless_text_table_transcription",
            "official_supporting_pdf_sha256": dor_refit["source"]["supporting_sha256"],
            "extractor": "Poppler pdftotext -tsv, pages 1-5",
        },
        "used_by_eos_records": record_ids,
        "notes": (
            "Exact printed volume and pressure tokens are retained beside normalized "
            "numeric values. Reported pressures are source-derived EOS outputs used only "
            "for transcription checks, not fit observations. The three starred MgO "
            "single-peak measurements are retained. Run AN012 is listed in article Table "
            "1 but absent from Tables S1-S6. The CC0 dedication is scoped by the adjacent "
            "LICENSE file and does not relicense the source article or PDF."
        ),
    }


def validated(
    doi: str, locations: list[str], finding: str, access_url: str
) -> dict[str, Any]:
    return {
        "status": "primary_source_validated",
        "note": "Equation, printed coefficients, units, fixed parameters, and reported fitting envelope checked directly against the primary source.",
        "audit_date": AUDIT_DATE,
        "verified_fields": [
            "equation",
            "parameters",
            "units",
            "reference_state",
            "phase",
            "published_uncertainties",
            "validity",
        ],
        "primary_source_check": {
            "access_url": access_url,
            "doi": doi,
            "locations": locations,
            "finding": finding,
        },
        "primary_data_check": {
            "status": "parameterization_only",
            "audit_date": AUDIT_DATE,
            "source_locations": locations,
            "finding": "Published coefficients are executable; this import does not claim a new refit of row-level observations.",
        },
    }


def calibration(status: str, description: str) -> dict[str, Any]:
    if "ruby" in description:
        kind = "ruby_fluorescence"
    elif "mgo_reference" in description:
        kind = "equation_of_state"
    elif description == "compiled_experimental_reference":
        kind = "ambient_pressure"
    else:
        kind = "other"
    method: dict[str, Any] = {"kind": kind, "source_location": description}
    if kind == "ruby_fluorescence":
        method["reference_calibration_record"] = (
            "ruby_dorogokupets_oganov_2007"
            if description.startswith("dor")
            else "ruby_dewaele_2004"
            if description.startswith("revised")
            else "ruby_mao_1986"
        )
    if kind == "equation_of_state":
        method["reference"] = "Tange et al. (2009) MgO pressure scale"
    return {
        "status": status,
        "methods": [method],
        "recalculation": {
            "status": "missing_calibrant_observations"
            if status == "resolved"
            else "not_possible",
            "notes": "A coefficient-level import cannot reconstruct unreported row-level calibration observations.",
        },
        "audit_date": AUDIT_DATE,
    }


def base_record(
    identifier: str,
    label: str,
    ref: dict[str, Any],
    eos_type: str,
    eos_model: str,
    v0: float,
    k0: float,
    k0p: float,
    pmin: float,
    pmax: float,
    notes: str,
    doi: str,
    locations: list[str],
    finding: str,
    access_url: str,
    pressure_calibration: dict[str, Any],
    errors: dict[str, float | None] | None = None,
    fixed: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "identifier": identifier,
        "label": label,
        "record_kind": "published",
        "equation_kind": "isothermal",
        "reference": ref,
        "eos": {
            "type": eos_type,
            "model": eos_model,
            "parameters": {"V0": v0, "K0": k0, "K0_prime": k0p},
        },
        "parameter_errors": errors or {"V0": None, "K0": None, "K0_prime": None},
        "parameter_error_confidence": None,
        "parameter_covariance": None,
        "fixed_parameters": fixed or [],
        "temperature_ref": 300.0,
        "experimental_pressure_range_gpa": [pmin, pmax],
        "pressure_range_status": "reported_exactly",
        "validity": {
            "pressure_gpa": [pmin, pmax],
            "temperature_k": [300.0, 300.0],
            "notes": [
                "Published fit/data envelope, not a phase-stability or extrapolation guarantee."
            ],
        },
        "notes": notes,
        "pressure_calibration": pressure_calibration,
        "parameter_provenance": {
            "V0": locations[-1],
            "K0": locations[-1],
            "K0_prime": locations[-1],
        },
        "scientific_validation": validated(doi, locations, finding, access_url),
        "source_lineage": [
            {
                "role": "scientific equation and coefficients",
                "citation": f"{ref['title']}, {', '.join(locations)}",
                "doi": doi,
            }
        ],
    }


def load_documents() -> dict[str, dict[str, Any]]:
    return {
        path.stem: json.loads(path.read_text()) for path in MATERIALS.glob("*.eosmat")
    }


def ensure_lead_hcp(documents: dict[str, dict[str, Any]]) -> None:
    if "lead_hcp" in documents:
        return
    documents["lead_hcp"] = {
        "format": "peritheos.material",
        "format_version": 3,
        "identifier": "lead_hcp",
        "name": "Lead (hcp)",
        "formula": "Pb",
        "phase": "hexagonal close-packed lead, high-pressure phase",
        "cell_contents": "Two Pb atoms in the conventional hcp cell; lattice-shape parameters are not inferred from the EOS table.",
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        },
        "formula_units_per_cell": 2,
        "atom_sites": [],
        "peaks": [],
        "notes": "Phase-specific card created for the hcp-Pb compression fits of Dewaele (2019).",
        "eos_records": [],
    }


def add(
    records_by_material: dict[str, list[dict[str, Any]]],
    material: str,
    record: dict[str, Any],
) -> None:
    records_by_material.setdefault(material, []).append(record)


def build_records(
    documents: dict[str, dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    delta_ref = reference(
        ["Lejaeghere", "Van Speybroeck", "Van Oost", "Cottenier"],
        2014,
        "Error estimates for solid-state density-functional theory predictions: an overview by means of the ground-state elemental crystals",
        "Critical Reviews in Solid State and Materials Sciences",
        "10.1080/10408436.2013.772503",
    )
    for row in rows(DELTA_EXPERIMENT):
        element = row["element"]
        material = ELEMENT_TO_MATERIAL[element]
        z = float(documents[material]["formula_units_per_cell"])
        rec = base_record(
            f"{material}_lejaeghere_2014_experimental_corrected_bm3",
            f"Lejaeghere et al. (2014), {element} experimental-corrected BM3",
            delta_ref,
            "BM3",
            "birch_murnaghan_3",
            float(row["V0"]) * z,
            float(row["K0"]),
            float(row["K0p"]),
            0.0,
            0.0,
            "Experimental reference parameters corrected for zero-point motion and room-temperature effects. Only the reference-state coefficients are published; the zero-width pressure envelope prevents implying an unreported high-pressure validation range.",
            delta_ref["doi"],
            [
                "official Delta-project archive exp.txt",
                "reference-data methodology in Sections III-IV",
            ],
            "The official archive labels exp.txt as experimental data with zero-point and thermal corrections and prints per-atom V0, B0, and B0-prime.",
            "https://archive.materialscloud.org/record/2023.133",
            calibration("unresolved", "compiled_experimental_reference"),
        )
        rec["temperature_ref"] = 0.0
        rec["validity"]["temperature_k"] = [0.0, 0.0]
        add(result, material, rec)

    sun_ref = reference(
        ["Sun", "Wu", "Guo", "Cai"],
        2010,
        "Two Universal Equations of State for Solids",
        "Zeitschrift fuer Naturforschung A",
        "10.1515/zna-2010-1-202",
    )
    for row in rows(SUN_LOW):
        solid = row["solid"]
        material = f"{solid.lower()}_sun_2010_legacy"
        v0 = SUN_V0[solid] * CM3_MOL_TO_A3_FORMULA
        prefixes = ["VN", "MRS3", "SMS3", "SMS4"] + (
            ["BN"] if solid in SUN_BAONZA_PRIORITY else []
        )
        for prefix in prefixes:
            eos_type, eos_model = SUN_MODELS[prefix]
            record = base_record(
                f"{solid.lower()}_sun_2010_low_{prefix.lower()}",
                f"Sun et al. (2010), {solid} low-pressure {prefix}",
                sun_ref,
                eos_type,
                eos_model,
                v0,
                float(row[f"{prefix}_K0"]),
                float(row[f"{prefix}_K0p"]),
                0.0,
                float(row["pmax"]),
                "Low-pressure refit of the source's legacy experimental compression compilation. The phase and individual pressure calibrations are unresolved, so this record remains on the existing phase-unresolved legacy card.",
                sun_ref["doi"],
                ["equations (3)-(15)", "Table 3"],
                "Table 3 prints the low-pressure range and Vinet, Baonza, MRS3, SMS3, and SMS4 coefficients; V0 is the experimental molar volume from Table 2.",
                "https://doi.org/10.1515/zna-2010-1-202",
                calibration("unresolved", "legacy_compiled_compression"),
            )
            add(result, material, record)

    d19_ref = reference(
        ["Dewaele"],
        2019,
        "Equations of State of Simple Solids (Including Pb, NaCl and LiF) Compressed in Helium or Neon in the Mbar Range",
        "Minerals",
        "10.3390/min9110684",
    )
    for row in rows(DEWAELE_2019):
        material = row["material"]
        z = float(documents[material]["formula_units_per_cell"])
        for scale in ("mao", "dor"):
            # Table 1 states that the unprinted uncertainties for the
            # Dorogokupets reduction are the same as those in the Mao column.
            errors = {
                "V0": float(row["eV0_mao"]) * z,
                "K0": float(row["eK0_mao"]),
                "K0_prime": float(row["eK0p_mao"]),
            }
            rec = base_record(
                f"{material}_dewaele_2019_{scale}_vinet",
                f"Dewaele (2019), {row['solid']} {scale.upper()}-scale Vinet",
                d19_ref,
                "Vinet",
                "vinet",
                float(row[f"V0_{scale}"]) * z,
                float(row[f"K0_{scale}"]),
                float(row[f"K0p_{scale}"]),
                float(row["pmin"]),
                float(row["pmax"]),
                "Vinet fit to static diamond-anvil-cell X-ray compression data in helium or neon. The two records preserve the published Mao-1986 and Dorogokupets-2007 ruby-scale sensitivity; they are alternative reductions of the same data, not independent experiments.",
                d19_ref["doi"],
                ["Equation (1)", "Table 1"],
                "Table 1 prints both pressure-scale reductions, 95% fit errors for the Mao reduction, pressure medium, gauge, and exact data domain.",
                "https://www.mdpi.com/2075-163X/9/11/684",
                calibration("partially_resolved", "linked_tungsten_scale")
                if row["solid"] in {"Re", "Fe"}
                else calibration("resolved", f"{scale}_ruby_scale"),
                errors,
            )
            rec["parameter_error_confidence"] = 0.95
            add(result, material, rec)

    d04_ref = reference(
        ["Dewaele", "Loubeyre", "Mezouar"],
        2004,
        "Equations of state of six metals above 94 GPa",
        "Physical Review B",
        "10.1103/PhysRevB.70.094112",
    )
    for row in rows(DEWAELE_2004):
        material = row["material"]
        z = float(documents[material]["formula_units_per_cell"])
        errors = {
            "V0": None,
            "K0": float(row["eK0"]) if row["eK0"] else None,
            "K0_prime": float(row["eK0p"]) if row["eK0p"] else None,
        }
        fixed = ["V0"] + ([row["fixed"]] if row["fixed"] else [])
        rec = base_record(
            f"{material}_dewaele_2004_{row['scale']}_vinet",
            f"Dewaele et al. (2004), {row['solid']} {row['scale'].replace('_', ' ')} Vinet",
            d04_ref,
            "Vinet",
            "vinet",
            float(row["V0"]) * z,
            float(row["K0"]),
            float(row["K0p"]),
            0.0,
            144.0,
            "Static synchrotron X-ray compression fit. This record fills a pressure-scale variant absent from the existing catalog; the shared 0-144 GPa envelope is conservative relative to the six source data series.",
            d04_ref["doi"],
            ["Equation (2)", "Table II"],
            "Table II prints the Vinet coefficients for both the classical and revised ruby calibrations, with 95% fitting errors and fixed coefficients in bold.",
            "https://doi.org/10.1103/PhysRevB.70.094112",
            calibration("resolved", row["scale"]),
            errors,
            fixed,
        )
        rec["parameter_error_confidence"] = 0.95
        add(result, material, rec)

    dor_ref = reference(
        ["Dorfman", "Prakapenka", "Meng", "Duffy"],
        2012,
        "Intercomparison of pressure standards (Au, Pt, Mo, MgO, NaCl and Ne) to 2.5 Mbar",
        "Journal of Geophysical Research",
        "10.1029/2012JB009292",
    )
    dor_refit = json.loads(DORFMAN_REFIT.read_text(encoding="utf-8"))
    tange_ref = reference(
        ["Tange", "Nishihara", "Tsuchiya"],
        2009,
        "Unified analyses for P-V-T equation of state of MgO: A solution for pressure-scale problems in high P-T experiments",
        "Journal of Geophysical Research",
        "10.1029/2008JB005813",
    )
    correction_ref = reference(
        ["Dorfman", "Prakapenka", "Meng", "Duffy"],
        2012,
        "Correction to Intercomparison of pressure standards (Au, Pt, Mo, MgO, NaCl and Ne) to 2.5 Mbar",
        "Journal of Geophysical Research",
        "10.1029/2012JB009800",
    )
    for row in rows(DORFMAN_2012):
        material = row["material"]
        for mode in ("fixed", "free"):
            fixed = ["V0"] + (["K0"] if mode == "fixed" else [])
            errors = {
                "V0": None,
                "K0": None if mode == "fixed" else float(row["eK0_free"]),
                "K0_prime": float(row[f"eK0p_{mode}"]),
            }
            rec = base_record(
                f"{material}_dorfman_2012_tange_mgo_k0_{mode}_vinet",
                f"Dorfman et al. (2012), {row['solid']} K0-{mode} Vinet",
                dor_ref,
                "Vinet",
                "vinet",
                float(row["V0"]),
                float(row[f"K0_{mode}"]),
                float(row[f"K0p_{mode}"]),
                0.0,
                float(row["pmax"]),
                "Simultaneous co-compression pressure-standard fit using the Tange et al. (2009) MgO EOS as reference. Fixed- and free-K0 solutions are correlated alternatives, not independent experiments.",
                dor_ref["doi"],
                ["Equation (2)", "Table 2"],
                "Table 2 prints conventional-cell V0 and fixed/free-K0 Vinet solutions for Au, Mo, and Pt.",
                "https://doi.org/10.1029/2012JB009292",
                calibration("resolved", "tange_2009_mgo_reference"),
                errors,
                fixed,
            )
            refit = dor_refit["fits"][mode]
            comparison = refit["published_parameter_comparison"][row["solid"]]
            source_range = dor_refit["selection"][
                "printed_pressure_range_gpa_by_material"
            ][row["solid"]]
            rec["pressure_range_status"] = "reference_parameterization"
            rec["validity"]["notes"] = [
                "The 0-Pmax range follows the published Table 2 pressure-standard parameterization; it is not the observation envelope.",
                f"The available auxiliary rows containing {row['solid']} span {source_range[0]}-{source_range[1]} GPa in the source's printed derived pressures.",
            ]
            rec["notes"] = (
                "Published simultaneous co-compression pressure-standard fit using "
                "the Tange et al. (2009) MgO Fit 3 Vinet EOS as the fixed 300 K "
                "anchor. Fixed- and free-K0 solutions are correlated alternatives, "
                "not independent experiments. The coefficients remain the printed "
                "Table 2 values: the independent source-row refit is a diagnostic "
                "and does not overwrite them."
            )
            rec["fit_datasets"] = [DORFMAN_DATASET_ID]
            rec["pressure_calibration"] = {
                "status": "resolved",
                "methods": [
                    {
                        "kind": "equation_of_state",
                        "reference": tange_ref,
                        "reference_eos_record": "mgo_b1_tange_2009_vinet",
                        "source_location": (
                            "Dorfman et al. Section 3.2 and Table 2; Tange et al. "
                            "(2009) Table 4, Fit 3 Vinet"
                        ),
                        "scope": "Fixed 300 K MgO anchor for the coupled fit",
                        "notes": "V0=74.698 A^3, K0=160.6 GPa, and K0-prime=4.37.",
                    }
                ],
                "recalculation": {
                    "status": "ready",
                    "notes": (
                        "A lossless, checksummed CSV transcription of the paired "
                        "conventional-cell volumes is bundled. The original publisher "
                        "PDF is not redistributed. The reproducer can verify and "
                        "regenerate the CSV from that exact source file."
                    ),
                },
                "audit_date": DORFMAN_AUDIT_DATE,
            }
            rec["scientific_validation"] = {
                "status": "primary_source_validated",
                "note": (
                    "Equation, coefficients, simultaneous-fit objective, fixed MgO "
                    "anchor, auxiliary paired-volume rows, correction, and fit "
                    "reproducibility were checked directly. The available source "
                    "rows do not reproduce the published coefficients."
                ),
                "audit_date": DORFMAN_AUDIT_DATE,
                "verified_fields": [
                    "equation",
                    "parameters",
                    "units",
                    "reference_state",
                    "phase",
                    "published_uncertainties",
                    "validity",
                    "pressure_calibration",
                    "fit_scope",
                    "primary_data",
                    "fit_reproducibility",
                ],
                "primary_source_check": {
                    "access_url": "https://doi.org/10.1029/2012JB009292",
                    "doi": dor_ref["doi"],
                    "locations": [
                        "Table 1",
                        "Equations (2)-(3)",
                        "Section 3.2",
                        "Table 2",
                        "auxiliary Tables S1-S6",
                    ],
                    "finding": (
                        "The article defines a coupled pressure-difference objective "
                        "and Table 2 prints the two correlated Vinet solutions. The "
                        "official auxiliary PDF provides 165 aligned volume rows; "
                        "AN012 is listed in Table 1 but absent from that PDF."
                    ),
                },
                "primary_data_check": {
                    "status": "bundled",
                    "audit_date": DORFMAN_AUDIT_DATE,
                    "source_locations": ["auxiliary Tables S1-S6"],
                    "finding": (
                        "The checksummed official PDF yields 165 rows, 368 volume "
                        "values, and 241 simultaneous material pairs. A lossless CSV "
                        "transcription is bundled under a scoped CC0 dedication. A "
                        "literal Equation (3) refit converges but does not recover the "
                        "printed Table 2 coefficients."
                    ),
                    "dataset_identifiers": [DORFMAN_DATASET_ID],
                    "resource": (
                        "peritheos/data/datasets/"
                        "dorfman-2012-tables-s1-s6-cocompression.csv"
                    ),
                    "reproduction_resource": (
                        "docs/data/dorfman-2012-cocompression-refit.json"
                    ),
                    "source_sha256": dor_refit["source"]["supporting_sha256"],
                    "observation_rows": dor_refit["selection"]["observation_rows"],
                    "pair_count": dor_refit["selection"]["pair_count"],
                    "available_observation_pressure_range_gpa": source_range,
                    "missing_run": dor_refit["selection"]["missing_run"],
                    "fit_mode": mode,
                    "objective": refit["objective"],
                    "published_objective": refit["published_objective"],
                    "record_parameter_comparison": comparison,
                },
            }
            rec["source_lineage"] = [
                {
                    "role": "simultaneous-fit equation and published coefficients",
                    "citation": (
                        f"{dor_ref['title']}, Equations (2)-(3), Section 3.2, and Table 2"
                    ),
                    "doi": dor_ref["doi"],
                },
                {
                    "role": "paired conventional-cell volume observations",
                    "citation": "Dorfman et al. (2012), auxiliary Tables S1-S6",
                    "doi": dor_ref["doi"],
                    "sha256": dor_refit["source"]["supporting_sha256"],
                    "access_url": dor_refit["source"]["supporting_url"],
                },
                {
                    "role": "fixed MgO pressure anchor",
                    "citation": f"{tange_ref['title']}, Table 4 Fit 3 Vinet",
                    "doi": tange_ref["doi"],
                },
                {
                    "role": "published correction (Figure S3 caption only)",
                    "citation": correction_ref["title"],
                    "doi": correction_ref["doi"],
                },
            ]
            add(result, material, rec)

    thermal_ref = reference(
        ["Zhang", "Zhang", "Kuang", "Xiong"],
        2025,
        "Equation of State Parameters of hcp-Fe Up to Super-Earth Interior Conditions",
        "Crystals",
        "10.3390/cryst15030221",
    )
    z = float(documents["iron"]["formula_units_per_cell"])
    for row in rows(THERMAL_2025):
        fit = row["fit"]
        v0 = float(row["V0_cm3_mol"]) * CM3_MOL_TO_A3_FORMULA * z
        errors = {
            "V0": float(row["eV0"]) * CM3_MOL_TO_A3_FORMULA * z if row["eV0"] else None,
            "K0": float(row["eK0"]),
            "K0_prime": float(row["eK0p"]),
        }
        rec = base_record(
            f"iron_zhang_2025_fit{fit}_{row['model']}_mgd",
            f"Zhang et al. (2025), hcp-Fe Fit #{fit} {row['type']}+MGD",
            thermal_ref,
            row["type"],
            row["model"],
            v0,
            float(row["K0"]),
            float(row["K0p"]),
            0.0,
            1374.0,
            f"Published mixed-data thermal EOS Fit #{fit}; RMSE={row['rmse']} GPa. The fit combines static and dynamic experiments with theoretical data, so it is not labeled purely experimental. Fit #1 is the authors' preferred relaxed-V0 solution.",
            thermal_ref["doi"],
            ["Equations (1)-(6)", "Table 2"],
            "Table 2 prints the complete BM3/Vinet plus Mie-Gruneisen-Debye coefficients, uncertainties, RMSE, and fit strategy.",
            "https://www.mdpi.com/2073-4352/15/3/221",
            calibration("unresolved", "mixed_static_dynamic_and_theoretical_dataset"),
            errors,
            ["V0"] if not row["eV0"] else [],
        )
        rec["equation_kind"] = "thermal"
        rec["validity"] = {
            "pressure_gpa": [0.0, 1374.0],
            "temperature_k": [300.0, 12000.0],
            "notes": [
                "Envelope of the compiled mixed P-V-T dataset; not a phase-stability boundary."
            ],
        }
        rec["thermal"] = {
            "type": "MieGruneisenDebye",
            "model": "mie_gruneisen_debye",
            "parameters": {
                "Tr": 300.0,
                "theta0": float(row["theta0"]),
                "gamma0": float(row["gamma0"]),
                "q": float(row["q"]),
                "n": 1.0,
            },
            "parameter_errors": {
                "theta0": float(row["etheta0"]),
                "gamma0": float(row["egamma0"]),
                "q": float(row["eq"]),
                "Tr": None,
                "n": None,
            },
            "debye_temperature_law": "integrated_gruneisen",
        }
        rec["parameter_provenance"].update(
            {
                "theta0": "Table 2",
                "gamma0": "Table 2",
                "q": "Table 2",
                "n": "one atom per Fe formula unit",
            }
        )
        add(result, "iron", rec)

    total = sum(map(len, result.values()))
    if total != 200:
        raise AssertionError(
            f"curated batch must contain exactly 200 records, found {total}"
        )
    return result


def write_dataset(records_by_material: dict[str, list[dict[str, Any]]]) -> None:
    fields = [
        "family",
        "material_identifier",
        "record_identifier",
        "formula",
        "model",
        "V0",
        "K0",
        "K0_prime",
        "pressure_min_gpa",
        "pressure_max_gpa",
        "temperature_min_k",
        "temperature_max_k",
        "reference_doi",
    ]
    DATASET.parent.mkdir(parents=True, exist_ok=True)
    with DATASET.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for material, records in sorted(records_by_material.items()):
            for record in records:
                identifier = record["identifier"]
                family = (
                    "thermal_mixed"
                    if "zhang_2025" in identifier
                    else "experimental_corrected"
                    if "lejaeghere_2014" in identifier
                    else "experimental_fit"
                )
                validity = record["validity"]
                params = record["eos"]["parameters"]
                writer.writerow(
                    {
                        "family": family,
                        "material_identifier": material,
                        "record_identifier": identifier,
                        "formula": record["label"].split(",")[1].strip().split()[0],
                        "model": record["eos"]["model"],
                        "V0": params["V0"],
                        "K0": params["K0"],
                        "K0_prime": params["K0_prime"],
                        "pressure_min_gpa": validity["pressure_gpa"][0],
                        "pressure_max_gpa": validity["pressure_gpa"][1],
                        "temperature_min_k": validity["temperature_k"][0],
                        "temperature_max_k": validity["temperature_k"][1],
                        "reference_doi": record["reference"]["doi"],
                    }
                )


def apply_records(
    documents: dict[str, dict[str, Any]],
    records_by_material: dict[str, list[dict[str, Any]]],
) -> None:
    identifiers = {
        record["identifier"]
        for records in records_by_material.values()
        for record in records
    }
    for material, additions in records_by_material.items():
        document = documents[material]
        document["eos_records"] = [
            record
            for record in document["eos_records"]
            if record["identifier"] not in identifiers
        ] + additions
        if material in {"gold", "molybdenum", "platinum"}:
            dor_ref = next(
                record["reference"]
                for record in additions
                if "_dorfman_2012_tange_mgo_k0_" in record["identifier"]
            )
            dor_refit = json.loads(DORFMAN_REFIT.read_text(encoding="utf-8"))
            datasets = [
                dataset
                for dataset in document.get("datasets", [])
                if dataset["identifier"] != DORFMAN_DATASET_ID
            ]
            document["datasets"] = datasets + [
                dorfman_dataset(material, dor_ref, dor_refit)
            ]
        (MATERIALS / f"{material}.eosmat").write_text(
            json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    print(
        f"wrote {len(identifiers)} records across {len(records_by_material)} material documents"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apply", action="store_true", help="update bundled material documents"
    )
    args = parser.parse_args()
    documents = load_documents()
    ensure_lead_hcp(documents)
    records = build_records(documents)
    write_dataset(records)
    if args.apply:
        apply_records(documents, records)
    else:
        print(
            f"validated {sum(map(len, records.values()))} records; pass --apply to write material documents"
        )


if __name__ == "__main__":
    main()
