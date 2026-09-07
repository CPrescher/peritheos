#!/usr/bin/env python3
"""Build the Zhang et al. (2025) hcp-Fe fit dataset and EOS records.

The source workbook is an official CC-BY-4.0 supplementary file.  It is not
vendored: pass its path with ``--source-workbook``.  The SHA-256 check below
prevents silently ingesting a different workbook revision.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
MATERIAL = ROOT / "peritheos" / "data" / "materials" / "iron.eosmat"
DATASET = (
    ROOT / "peritheos" / "data" / "datasets" / "iron-zhang-2025-tables-s1-s3-s4-pvt.csv"
)
DATASET_ID = "iron_zhang_2025_tables_s1_s3_s4_pvt"
ARTICLE_DOI = "10.3390/cryst15030221"
ARTICLE_URL = (
    "https://mdpi-res.com/d_attachment/crystals/crystals-15-00221/"
    "article_deploy/crystals-15-00221.pdf"
)
SUPPLEMENT_URL = (
    "https://mdpi-res.com/d_attachment/crystals/crystals-15-00221/"
    "article_deploy/crystals-15-00221-s001.zip"
)
SOURCE_WORKBOOK = "crystals-3484051-supplementary.xlsx"
SOURCE_WORKBOOK_SHA256 = (
    "a3654047a5cabf03633c9ce280a264fa2c7bfae1c86d4e4c04028e715a0cd4a2"
)
SOURCE_ZIP_SHA256 = "518ed6f4a17e5dee22b69bfe395fe160e81bf03ea3e9090a8573e251a77ab5a8"
SOURCE_PDF_SHA256 = "f7da2f0446bd76d35b25f363e7f4ee4266db130636b7ef907ddd459af9cc1fe8"
AUDIT_DATE = "2026-09-07"
CM3_MOL_TO_A3_FORMULA = 1.6605390671738466
FORMULA_UNITS_PER_HCP_CELL = 2
VOLUME_SCALE = CM3_MOL_TO_A3_FORMULA * FORMULA_UNITS_PER_HCP_CELL

FIELDS = [
    "source_table",
    "source_row",
    "source_reference_as_published",
    "data_class",
    "dynamic_path",
    "pressure_gpa_raw",
    "pressure_uncertainty_gpa_raw",
    "pressure_calibrant",
    "pressure_gpa_fit",
    "pressure_uncertainty_gpa_fit",
    "pressure_scale",
    "pressure_correction_applied",
    "molar_volume_cm3_mol",
    "molar_volume_uncertainty_cm3_mol",
    "volume_a3_conventional_cell",
    "volume_uncertainty_a3_conventional_cell",
    "temperature_k",
    "temperature_uncertainty_k",
    "temperature_provenance",
    "fit1_observed_minus_model_gpa",
    "fit2_observed_minus_model_gpa",
    "fit5_observed_minus_model_gpa",
    "source_cell_note",
    "used_in_fit",
]

FIT_SPECS = {
    "1": {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "V0": 6.756,
        "eV0": 0.010,
        "K0": 174.7,
        "eK0": 1.7,
        "K0p": 4.790,
        "eK0p": 0.014,
        "theta0": 1209.0,
        "etheta0": 73.0,
        "gamma0": 2.86,
        "egamma0": 0.10,
        "q": 0.84,
        "eq": 0.05,
        "rmse": 4.504302,
    },
    "2": {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "V0": 6.753,
        "eV0": None,
        "K0": 175.1,
        "eK0": 0.5,
        "K0p": 4.787,
        "eK0p": 0.009,
        "theta0": 1205.0,
        "etheta0": 73.0,
        "gamma0": 2.86,
        "egamma0": 0.10,
        "q": 0.84,
        "eq": 0.04,
        "rmse": 4.504408,
    },
    "5": {
        "type": "Vinet",
        "model": "vinet",
        "V0": 6.753,
        "eV0": None,
        "K0": 151.6,
        "eK0": 0.5,
        "K0p": 5.845,
        "eK0p": 0.010,
        "theta0": 960.0,
        "etheta0": 69.0,
        "gamma0": 3.51,
        "egamma0": 0.13,
        "q": 1.28,
        "eq": 0.05,
        "rmse": 5.084659,
    },
}

SOURCE_COVARIANCE = {
    "1": [
        [0.00010006, -0.016735, 0.00011166, 0.10853, -9.8975e-05, -0.00012682],
        [-0.016735, 3.0208, -0.022456, -13.222, 0.0031549, 0.014549],
        [0.00011166, -0.022456, 0.0001994, 0.21017, 0.00034559, 4.7259e-05],
        [0.10853, -13.222, 0.21017, 5395.9, 4.759, 0.52581],
        [-9.8975e-05, 0.0031549, 0.00034559, 4.759, 0.010316, 0.0036313],
        [-0.00012682, 0.014549, 4.7259e-05, 0.52581, 0.0036313, 0.0020943],
    ],
    "2": [
        [0, 0, 0, 0, 0, 0],
        [0, 0.22244, -0.0037781, 4.9088, -0.013392, -0.0066644],
        [0, -0.0037781, 7.4552e-05, 0.090054, 0.00045511, 0.0001888],
        [0, 4.9088, 0.090054, 5293.4, 4.8789, 0.67663],
        [0, -0.013392, 0.00045511, 4.8789, 0.010191, 0.0035072],
        [0, -0.0066644, 0.0001888, 0.67663, 0.0035072, 0.0019348],
    ],
    "5": [
        [0, 0, 0, 0, 0, 0],
        [0, 0.22988, -0.004717, 6.2192, -0.013427, -0.0061313],
        [0, -0.004717, 0.00011013, 0.048069, 0.0006019, 0.00022626],
        [0, 6.2192, 0.048069, 4807.1, 5.991, 0.69938],
        [0, -0.013427, 0.0006019, 5.991, 0.016374, 0.0048411],
        [0, -0.0061313, 0.00022626, 0.69938, 0.0048411, 0.0024021],
    ],
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _number(value: Any, *, cell: str) -> tuple[float, str]:
    note = ""
    if isinstance(value, str):
        note = f"{cell} is numeric text in the source workbook; parsed as a number"
    try:
        return float(value), note
    except (TypeError, ValueError) as error:
        raise ValueError(f"expected a number at {cell}, found {value!r}") from error


def _column_index(coordinate: str) -> int:
    letters = re.match(r"[A-Z]+", coordinate)
    if letters is None:
        raise ValueError(f"invalid XLSX coordinate {coordinate!r}")
    result = 0
    for character in letters.group():
        result = result * 26 + ord(character) - ord("A") + 1
    return result


def _xlsx_sheets(
    workbook_path: Path, widths: dict[str, int]
) -> dict[str, list[tuple[int, tuple[Any, ...]]]]:
    """Read selected worksheets with only the standard library.

    The source workbook contains formatting through Excel's maximum column, so
    this reader deliberately ignores cells beyond each table's data width.
    """
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    package_rel_ns = "http://schemas.openxmlformats.org/package/2006/relationships"
    with zipfile.ZipFile(workbook_path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{main_ns}}}si"):
                shared_strings.append(
                    "".join(node.text or "" for node in item.iter(f"{{{main_ns}}}t"))
                )

        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            relationship.attrib["Id"]: relationship.attrib["Target"]
            for relationship in relationships.findall(
                f"{{{package_rel_ns}}}Relationship"
            )
        }
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet_paths = {}
        for sheet in workbook.findall(f".//{{{main_ns}}}sheet"):
            name = sheet.attrib["name"]
            if name not in widths:
                continue
            target = targets[sheet.attrib[f"{{{rel_ns}}}id"]]
            sheet_paths[name] = "xl/" + target.lstrip("/").removeprefix("xl/")

        result: dict[str, list[tuple[int, tuple[Any, ...]]]] = {}
        for name, width in widths.items():
            root = ET.fromstring(archive.read(sheet_paths[name]))
            rows = []
            for row in root.findall(f".//{{{main_ns}}}row"):
                row_number = int(row.attrib["r"])
                if row_number < 3:
                    continue
                values: list[Any] = [None] * width
                for cell in row.findall(f"{{{main_ns}}}c"):
                    column = _column_index(cell.attrib["r"])
                    if column > width:
                        continue
                    value_node = cell.find(f"{{{main_ns}}}v")
                    cell_type = cell.attrib.get("t")
                    if cell_type == "inlineStr":
                        value = "".join(
                            node.text or "" for node in cell.iter(f"{{{main_ns}}}t")
                        )
                    elif value_node is None or value_node.text is None:
                        value = None
                    elif cell_type == "s":
                        value = shared_strings[int(value_node.text)]
                    elif cell_type in {"str", "e"}:
                        value = value_node.text
                    elif cell_type == "b":
                        value = value_node.text == "1"
                    else:
                        value = float(value_node.text)
                    values[column - 1] = value
                if values[0] is not None:
                    rows.append((row_number, tuple(values)))
            result[name] = rows
    return result


def extract_rows(workbook_path: Path) -> list[dict[str, Any]]:
    if _sha256(workbook_path) != SOURCE_WORKBOOK_SHA256:
        raise ValueError(f"unexpected SHA-256 for {workbook_path}")
    result: list[dict[str, Any]] = []
    table_specs = {
        "Table S1": (16, "static_experiment"),
        "Table S3": (13, "dynamic_experiment"),
        "Table S4": (12, "ab_initio"),
    }
    sheets = _xlsx_sheets(
        workbook_path, {table: width for table, (width, _) in table_specs.items()}
    )
    for table, (width, data_class) in table_specs.items():
        for source_row, values in sheets[table]:
            reference = str(values[0])
            pressure, pressure_note = _number(values[1], cell=f"{table}!B{source_row}")
            volume, volume_note = _number(values[3], cell=f"{table}!D{source_row}")
            temperature, temperature_note = _number(
                values[5], cell=f"{table}!F{source_row}"
            )
            pressure_error = float(values[2] or 0.0)
            volume_error = float(values[4] or 0.0)
            temperature_error = float(values[6] or 0.0)
            dynamic_path = ""
            calibrant = ""
            pressure_scale = ""
            corrected = 0
            fit_pressure = pressure
            fit_pressure_error = pressure_error
            if table == "Table S1":
                calibrant = str(values[7] or "")
                if values[8] is not None:
                    fit_pressure = float(values[8])
                    fit_pressure_error = float(values[9])
                    pressure_scale = str(values[10])
                    corrected = 1
                residuals = (values[11], values[12], values[15])
                temperature_provenance = "experimental_reported"
            elif table == "Table S3":
                dynamic_path = str(values[7])
                residuals = (values[8], values[9], values[12])
                if reference.startswith("Huang et al."):
                    temperature_provenance = "model_calculated_huang_2022"
                elif dynamic_path == "Ramp":
                    temperature_provenance = "theory_assigned_zhuang_2021_isentrope"
                else:
                    temperature_provenance = "theory_assigned_zhuang_2021_hugoniot"
            else:
                residuals = (values[7], values[8], values[11])
                temperature_provenance = "ab_initio_calculated"
            result.append(
                {
                    "source_table": table,
                    "source_row": source_row,
                    "source_reference_as_published": reference,
                    "data_class": data_class,
                    "dynamic_path": dynamic_path,
                    "pressure_gpa_raw": pressure,
                    "pressure_uncertainty_gpa_raw": pressure_error,
                    "pressure_calibrant": calibrant,
                    "pressure_gpa_fit": fit_pressure,
                    "pressure_uncertainty_gpa_fit": fit_pressure_error,
                    "pressure_scale": pressure_scale,
                    "pressure_correction_applied": corrected,
                    "molar_volume_cm3_mol": volume,
                    "molar_volume_uncertainty_cm3_mol": volume_error,
                    "volume_a3_conventional_cell": volume * VOLUME_SCALE,
                    "volume_uncertainty_a3_conventional_cell": volume_error
                    * VOLUME_SCALE,
                    "temperature_k": temperature,
                    "temperature_uncertainty_k": temperature_error,
                    "temperature_provenance": temperature_provenance,
                    "fit1_observed_minus_model_gpa": float(residuals[0]),
                    "fit2_observed_minus_model_gpa": float(residuals[1]),
                    "fit5_observed_minus_model_gpa": float(residuals[2]),
                    "source_cell_note": "; ".join(
                        note
                        for note in (pressure_note, volume_note, temperature_note)
                        if note
                    ),
                    "used_in_fit": 1,
                }
            )
    counts = {
        table: sum(row["source_table"] == table for row in result)
        for table in table_specs
    }
    if counts != {"Table S1": 1078, "Table S3": 132, "Table S4": 103}:
        raise AssertionError(f"unexpected selected-row counts: {counts}")
    return result


def write_rows(rows: list[dict[str, Any]]) -> str:
    DATASET.parent.mkdir(parents=True, exist_ok=True)
    with DATASET.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return _sha256(DATASET)


def _reference() -> dict[str, Any]:
    return {
        "authors": ["Zhang", "Zhang", "Kuang", "Xiong"],
        "year": 2025,
        "title": "Equation of State Parameters of hcp-Fe Up to Super-Earth Interior Conditions",
        "source": "Crystals",
        "volume": "15",
        "locator": "221",
        "doi": ARTICLE_DOI,
    }


def _covariance(fit: str) -> dict[str, Any]:
    source = SOURCE_COVARIANCE[fit]
    scales = [VOLUME_SCALE, 1.0, 1.0, 1.0, 1.0, 1.0]
    matrix = [
        [value * scales[row] * scales[column] for column, value in enumerate(values)]
        for row, values in enumerate(source)
    ]
    return {
        "parameter_order": [
            "rt_eos.V0",
            "rt_eos.K0",
            "rt_eos.K0_prime",
            "theta0",
            "gamma0",
            "q",
        ],
        "matrix": matrix,
    }


def build_record(fit: str) -> dict[str, Any]:
    spec = FIT_SPECS[fit]
    v0 = spec["V0"] * VOLUME_SCALE
    fixed = ["V0"] if spec["eV0"] is None else []
    identifier = f"iron_zhang_2025_fit{fit}_{spec['model']}_mgd"
    preferred = fit == "1"
    return {
        "identifier": identifier,
        "label": f"Zhang et al. (2025), hcp-Fe Fit #{fit} {spec['type']}+MGD",
        "record_kind": "published",
        "equation_kind": "thermal",
        "reference": _reference(),
        "volume_basis": {
            "kind": "formula_units",
            "formula_units": FORMULA_UNITS_PER_HCP_CELL,
            "molar_mass_g_mol": 55.845,
        },
        "eos": {
            "type": spec["type"],
            "model": spec["model"],
            "parameters": {"V0": v0, "K0": spec["K0"], "K0_prime": spec["K0p"]},
        },
        "parameter_errors": {
            "V0": None if spec["eV0"] is None else spec["eV0"] * VOLUME_SCALE,
            "K0": spec["eK0"],
            "K0_prime": spec["eK0p"],
        },
        "parameter_error_confidence": None,
        "parameter_covariance": _covariance(fit),
        "fixed_parameters": fixed,
        "temperature_ref": 300.0,
        "experimental_pressure_range_gpa": [10.7, 1374.0],
        "pressure_range_status": "reported_exactly",
        "fit_datasets": [DATASET_ID],
        "validity": {
            "pressure_gpa": [10.7, 1374.0],
            "temperature_k": [298.0, 12000.0],
            "notes": [
                "Envelope of the 1,313-row mixed fit dataset, not a phase-stability boundary.",
                "Static Table S1 reaches 4695 K although Methods section 2.1 says 3400 K; the deposited rows and residuals demonstrate that the higher-temperature rows were fitted.",
            ],
        },
        "notes": (
            f"Published mixed-input thermal EOS Fit #{fit}; Supplementary Table S6 "
            f"reports RMSE={spec['rmse']} GPa. It fits 1,078 static experimental, "
            "132 dynamic P-V, and 103 ab-initio rows. Dynamic temperatures are "
            "calculated/model-assigned rather than measured. "
            + (
                "This is the authors' preferred relaxed-V0 solution."
                if preferred
                else "Fit #1 is the authors' preferred relaxed-V0 solution."
            )
        ),
        "pressure_calibration": {
            "status": "partially_resolved",
            "methods": [
                {
                    "kind": "equation_of_state",
                    "reference": {
                        "authors": [
                            "Fei",
                            "Ricolleau",
                            "Frank",
                            "Mibe",
                            "Shen",
                            "Prakapenka",
                        ],
                        "year": 2007,
                        "source": "PNAS",
                        "doi": "10.1073/pnas.0609013104",
                    },
                    "source_location": "Supplementary Table S1 pressure-scale column",
                    "scope": "499 static rows recalculated to Fei et al. (2007)",
                },
                {
                    "kind": "other",
                    "reference": {
                        "authors": ["Ye", "Prakapenka", "Meng", "Shim"],
                        "year": 2018,
                        "source": "High Pressure Research",
                        "doi": "10.1080/08957959.2018.1493477",
                    },
                    "source_location": "Supplementary Table S1 pressure-scale column",
                    "scope": "142 static rows recalculated to Ye et al. (2018)",
                    "notes": "The workbook names Ye et al. (2018), whereas Methods section 2.1 cites Fei et al. (2007) and Dorfman et al. (2012).",
                },
                {
                    "kind": "other",
                    "source_location": "Supplementary Table S1 blank pressure-scale cells",
                    "scope": "437 static rows retain the compiled source pressure",
                    "notes": "The sheet says these were already internally consistent or lacked a compatible internally consistent calibrant scale.",
                },
                {
                    "kind": "shock_wave",
                    "source_location": "Supplementary Table S3",
                    "scope": "72 Hugoniot and 60 ramp-compression P-V rows",
                },
                {
                    "kind": "ab_initio",
                    "source_location": "Supplementary Table S4",
                    "scope": "103 calculated P-V-T rows; no experimental gauge",
                },
            ],
            "recalculation": {
                "status": "missing_calibrant_observations",
                "notes": "Table S1 deposits the final corrected pressures and errors but not the row-wise calibrant volumes/temperatures needed to reproduce the 641 recalculations upstream.",
            },
            "audit_date": AUDIT_DATE,
        },
        "parameter_provenance": {
            "V0": "Table 2/S6 molar V0 converted exactly to a two-atom hcp conventional-cell volume",
            "K0": "Table 2 and Supplementary Table S6",
            "K0_prime": "Table 2 and Supplementary Table S6",
            "theta0": "Table 2 and Supplementary Table S6",
            "gamma0": "Table 2 and Supplementary Table S6",
            "q": "Table 2 and Supplementary Table S6",
            "n": "one atom per Fe formula unit",
            "covariance": "Supplementary Table S5; V0 row and column converted to the stored cell-volume basis",
        },
        "fit_provenance": {
            "software": {"name": "EosFit-GUI", "version": "not reported"},
            "dataset": DATASET_ID,
            "selection": {
                "predicate": "all numeric rows in Supplementary Tables S1, S3, and S4",
                "included_rows": 1313,
                "included_static_experimental_rows": 1078,
                "included_dynamic_pv_rows": 132,
                "included_ab_initio_rows": 103,
                "excluded_rows": 325,
                "exclusion": "all Table S2 self-calibrated hcp-Fe rows",
            },
            "objective": "Unweighted least squares of pressure residuals; supplement defines residual as observed/calculated pressure minus fitted pressure.",
            "weights": [],
            "refined_parameters": [
                name
                for name in ("V0", "K0", "K0_prime", "theta0", "gamma0", "q")
                if name not in fixed
            ],
            "fixed_parameters": fixed + ["Tr=300 K", "n=1", "R=8.314 J mol^-1 K^-1"],
            "statistics": {
                "observations": 1313,
                "published_rmse_gpa": spec["rmse"],
            },
            "reproduction": {
                "script": "scripts/reproduce_zhang_2025_hcp_iron.py",
                "solver": "scipy.optimize.least_squares",
                "scope_status": "final_input_parity_upstream_reduction_partial",
                "reproduced_scope": [
                    "all 1,313 deposited final fit inputs",
                    "published residual columns",
                    "Fits 1, 2, and 5 coefficients",
                    "Supplementary Table S5 covariance",
                ],
                "unresolved_upstream_steps": [
                    "641 static-pressure recalculations lack row-wise calibrant observations",
                    "the generation/interpolation procedure for 122 Zhuang-derived dynamic temperatures is not documented",
                ],
                "objective_inference": "An independently computed unweighted Jacobian covariance matches Table S5 within 0.6%; no row weights are stated by the source.",
            },
        },
        "scientific_validation": {
            "status": "primary_source_validated",
            "note": "Equations, selected rows, coefficients, covariance, units, and residual convention checked against the article and official supplement; source types are retained per row.",
            "audit_date": AUDIT_DATE,
            "verified_fields": [
                "equation",
                "parameters",
                "units",
                "reference_state",
                "phase",
                "published_uncertainties",
                "covariance",
                "fit_dataset",
                "validity",
            ],
            "primary_source_check": {
                "access_url": ARTICLE_URL,
                "doi": ARTICLE_DOI,
                "locations": [
                    "Methods sections 2.1-2.2",
                    "Equations (1)-(6)",
                    "Table 2",
                    "Supplementary Tables S1-S6",
                ],
                "finding": "The official supplement supplies every selected fit row, all five residual columns, covariance matrices, exact RMSE values, and the excluded self-calibrated Table S2 rows.",
            },
            "primary_data_check": {
                "status": "bundled",
                "audit_date": AUDIT_DATE,
                "dataset_identifiers": [DATASET_ID],
                "source_locations": ["Supplementary Tables S1, S3, and S4"],
                "finding": "All 1,313 published fit inputs are bundled with experimental, dynamic, ab-initio, calibration, and temperature-provenance labels; an independent six/five-parameter refit recovers Fits 1, 2, and 5.",
            },
        },
        "source_lineage": [
            {
                "role": "scientific equations, selection rules, and coefficients",
                "citation": "Zhang et al. (2025), Methods 2.1-2.2, Equations (1)-(6), Table 2",
                "doi": ARTICLE_DOI,
            },
            {
                "role": "selected P-V-T inputs, residuals, covariance, and fit comparison",
                "citation": "Zhang et al. (2025), Supplementary Tables S1-S6",
                "doi": ARTICLE_DOI,
            },
            {
                "role": "model basis for assigned Hugoniot and ramp temperatures",
                "citation": "Zhuang et al. (2021), hcp-Fe Hugoniot and isentrope calculations",
                "doi": "10.1103/PhysRevB.103.144102",
            },
            {
                "role": "calculated temperatures for ten laser-shock rows",
                "citation": "Huang et al. (2022), laser-driven shock compression of iron",
                "doi": "10.1038/s41467-022-28255-2",
            },
        ],
        "thermal": {
            "type": "MieGruneisenDebye",
            "model": "mie_gruneisen_debye",
            "parameters": {
                "Tr": 300.0,
                "theta0": spec["theta0"],
                "gamma0": spec["gamma0"],
                "q": spec["q"],
                "n": 1.0,
            },
            "parameter_errors": {
                "theta0": spec["etheta0"],
                "gamma0": spec["egamma0"],
                "q": spec["eq"],
                "Tr": None,
                "n": None,
            },
            "fixed_parameters": ["Tr", "n"],
            "debye_temperature_law": "integrated_gruneisen",
            "gas_constant_j_mol_k": 8.314,
        },
    }


def build_dataset(resource_sha256: str) -> dict[str, Any]:
    value_columns = {
        "pressure_gpa_raw": ("pressure", "GPa", "value"),
        "pressure_uncertainty_gpa_raw": ("pressure", "GPa", "uncertainty"),
        "pressure_gpa_fit": ("pressure", "GPa", "value"),
        "pressure_uncertainty_gpa_fit": ("pressure", "GPa", "uncertainty"),
        "molar_volume_cm3_mol": ("molar_volume", "cm^3/mol", "value"),
        "molar_volume_uncertainty_cm3_mol": ("molar_volume", "cm^3/mol", "uncertainty"),
        "volume_a3_conventional_cell": (
            "volume",
            "angstrom^3/conventional_unit_cell",
            "value",
        ),
        "volume_uncertainty_a3_conventional_cell": (
            "volume",
            "angstrom^3/conventional_unit_cell",
            "uncertainty",
        ),
        "temperature_k": ("temperature", "K", "value"),
        "temperature_uncertainty_k": ("temperature", "K", "uncertainty"),
        "fit1_observed_minus_model_gpa": ("pressure_residual", "GPa", "value"),
        "fit2_observed_minus_model_gpa": ("pressure_residual", "GPa", "value"),
        "fit5_observed_minus_model_gpa": ("pressure_residual", "GPa", "value"),
    }
    columns = []
    for name in FIELDS:
        quantity, unit, role = value_columns.get(name, (name, "dimensionless", "flag"))
        column: dict[str, Any] = {
            "name": name,
            "quantity": quantity,
            "unit": unit,
            "role": role,
        }
        if name.endswith("uncertainty_gpa_raw"):
            column["of"] = "pressure_gpa_raw"
        elif name.endswith("uncertainty_gpa_fit"):
            column["of"] = "pressure_gpa_fit"
        elif name == "molar_volume_uncertainty_cm3_mol":
            column["of"] = "molar_volume_cm3_mol"
        elif name == "volume_uncertainty_a3_conventional_cell":
            column["of"] = "volume_a3_conventional_cell"
        elif name == "temperature_uncertainty_k":
            column["of"] = "temperature_k"
        columns.append(column)
    return {
        "identifier": DATASET_ID,
        "kind": "pressure_volume_temperature",
        "description": "The exact 1,313-row mixed-input regression dataset for Zhang et al. Fits 1, 2, and 5: 1,078 static experimental rows, 132 dynamic experimental P-V rows with calculated/model-assigned temperatures, and 103 ab-initio rows.",
        "reference": _reference(),
        "source_location": "Supplementary Tables S1, S3, and S4",
        "source_url": SUPPLEMENT_URL,
        "license": "CC BY 4.0",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "columns": columns,
        "resource": {
            "path": f"datasets/{DATASET.name}",
            "sha256": resource_sha256,
            "media_type": "text/csv",
        },
        "used_by_eos_records": [
            f"iron_zhang_2025_fit{fit}_{FIT_SPECS[fit]['model']}_mgd"
            for fit in ("1", "2", "5")
        ],
        "notes": "Table S2's 325 hcp-Fe-self-calibrated static rows were excluded by the authors and are not bundled as fit inputs. Source spelling, including 'Zhaung' in 38 Table S4 rows, is preserved. Table S1 row 461 stores its volume as numeric text and is explicitly coerced. Source residuals use observed/calculated minus model pressure.",
        "provenance": {
            "type": "transcribed_from_official_supplement",
            "method": "Selected every numeric row from S1, S3, and S4; used corrected S1 pressure/error where present, otherwise the compiled pressure/error; converted molar volume to a two-atom hcp cell using the exact Avogadro constant.",
            "source_article_pdf_sha256": SOURCE_PDF_SHA256,
            "source_supplement_zip_sha256": SOURCE_ZIP_SHA256,
            "source_workbook_name": SOURCE_WORKBOOK,
            "source_workbook_sha256": SOURCE_WORKBOOK_SHA256,
            "quality_control": {
                "rows": 1313,
                "table_s1_rows": 1078,
                "table_s3_rows": 132,
                "table_s4_rows": 103,
                "table_s2_excluded_rows": 325,
                "s1_corrected_pressure_rows": 641,
                "s1_retained_pressure_rows": 437,
                "inferred_dynamic_temperature_relations": {
                    "status": "exact_relationship_in_deposited_values_not_published_methodology",
                    "zhuang_hugoniot_rows": 62,
                    "zhuang_hugoniot_temperature_k": "0.01491*P_gpa^2 + 21.77176*P_gpa - 395.42764",
                    "zhuang_hugoniot_max_abs_reconstruction_k": 2.3e-12,
                    "zhuang_ramp_rows": 60,
                    "zhuang_ramp_temperature_k": "4.22e-7*P_gpa^3 - 0.00155*P_gpa^2 + 2.58948*P_gpa + 920.817",
                    "zhuang_ramp_max_abs_reconstruction_k": 4.4e-12,
                },
            },
            "transcribed_on": AUDIT_DATE,
        },
    }


def apply(resource_sha256: str) -> None:
    document = json.loads(MATERIAL.read_text(encoding="utf-8"))
    identifiers = {
        f"iron_zhang_2025_fit{fit}_{FIT_SPECS[fit]['model']}_mgd" for fit in FIT_SPECS
    }
    document["eos_records"] = [
        record
        for record in document["eos_records"]
        if record["identifier"] not in identifiers
    ] + [build_record(fit) for fit in ("1", "2", "5")]
    document["datasets"] = [
        dataset
        for dataset in document.get("datasets", [])
        if dataset["identifier"] != DATASET_ID
    ] + [build_dataset(resource_sha256)]
    MATERIAL.write_text(
        json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-workbook", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    rows = extract_rows(args.source_workbook)
    resource_sha256 = write_rows(rows)
    if args.apply:
        apply(resource_sha256)
    print(f"wrote {len(rows)} rows to {DATASET} ({resource_sha256})")


if __name__ == "__main__":
    main()
