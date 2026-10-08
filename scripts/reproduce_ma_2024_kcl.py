#!/usr/bin/env python3
"""Audit the final Ma, Sumita and Murakami (2024) KCl data and joint fit.

The optional source extraction reads cached workbook values with the standard
library, verifies the author's file checksum, and preserves all S1-S4 cells.
No modeled pressure column or grid is used as a fit observation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
from scipy.constants import Avogadro, Boltzmann, R, hbar
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
RECORD_ID = "kcl_b2_ma_2024_bm3_mgd"
SOURCE_SHA256 = "e4484174250276f092a8abd0e3b7a4f76d9a7a2cb853b76cf5ca7d54da057b84"
SOURCE_URL = (
    "https://data.mendeley.com/public-files/datasets/7svmv9hvft/files/"
    "51aa0d16-dad5-49b2-a4c5-be70a64488e9/file_downloaded"
)
SOURCE = DATA / "ma_2024_sources/data-tables.xlsx"
SUPPLEMENT = SOURCE.parent / "supporting-information.docx"
SUPPLEMENT_SHA256 = "427318e7e0a536d7f7635da07ec305fa2bb809832f457cfeefcafc68918a6c01"
PUBLISHED = np.array([32.48, 21.33, 4.836])
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
TABLES = {
    "acoustic": (
        "kcl-ma-2024-table-s1-acoustic.csv",
        [
            "source_row",
            "run",
            "diamond_raman_shift_cm1",
            "diamond_pressure_gpa",
            "diamond_pressure_error_gpa",
            "molar_volume_cm3_mol",
            "molar_volume_error_cm3_mol",
            "vp_km_s",
            "vp_error_km_s",
            "vs_km_s",
            "vs_error_km_s",
            "kt_gpa",
            "kt_error_gpa",
            "theta_d_k",
            "theta_d_error_k",
            "model_pressure_gpa",
            "model_pressure_error_gpa",
        ],
    ),
    "walker": (
        "kcl-ma-2024-table-s1-walker-recalibrated.csv",
        [
            "source_row",
            "run",
            "molar_volume_cm3_mol",
            "source_temperature_k",
            "walker_pressure_gpa",
            "walker_pressure_error_gpa",
            "matsui_300k_pressure_gpa",
            "matsui_pressure_error_gpa",
            "model_pressure_gpa",
            "model_pressure_error_gpa",
        ],
    ),
    "cold_grid": (
        "kcl-ma-2024-table-s3-model-grid.csv",
        [
            "source_cell",
            "molar_volume_cm3_mol",
            "model_pressure_gpa",
            "model_pressure_error_gpa",
        ],
    ),
    "thermal_grid": (
        "kcl-ma-2024-table-s4-model-grid.csv",
        ["source_cell", "temperature_k", "model_thermal_pressure_gpa"],
    ),
}


def workbook_cells(path: Path) -> dict[str, dict[str, str | float]]:
    """Return exact cached source values, retaining table layout coordinates."""
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("This is not the final 7svmv9hvft.1 author workbook")
    with ZipFile(path) as archive:
        strings = [
            "".join(element.itertext())
            for element in ET.fromstring(archive.read("xl/sharedStrings.xml"))
        ]
        sheets = {}
        for index in range(1, 5):
            cells = {}
            for cell in ET.fromstring(
                archive.read(f"xl/worksheets/sheet{index}.xml")
            ).findall(".//s:c", NS):
                value = cell.find("s:v", NS)
                if value is not None and value.text is not None:
                    cells[cell.attrib["r"]] = (
                        strings[int(value.text)]
                        if cell.get("t") == "s"
                        else float(value.text)
                    )
            sheets[f"Table S{index}"] = cells
    return sheets


def extract(path: Path = SOURCE) -> dict[str, list[list]]:
    """Extract observations, reductions, and modeled grids into distinct files."""
    sheets = workbook_cells(path)
    s = sheets["Table S1"]
    rows = {
        "acoustic": [
            [i, *[s[f"{c}{i}"] for c in "ABCDEFGHIJKLMNOP"]] for i in range(3, 14)
        ],
        "walker": [[i, *[s[f"{c}{i}"] for c in "ABCEGHKNP"]] for i in range(15, 23)],
        "cold_grid": [],
        "thermal_grid": [],
    }
    for key, sheet, columns in [
        ("cold_grid", "Table S3", ["ABC", "DEF", "GHI", "JKL"]),
        ("thermal_grid", "Table S4", ["AB", "CD", "EF", "GH", "IJ", "KL", "MN", "OP"]),
    ]:
        for group in columns:
            for i in range(3, 62):
                if isinstance(sheets[sheet].get(f"{group[0]}{i}"), float):
                    rows[key].append(
                        [f"{group[0]}{i}", *[sheets[sheet][f"{c}{i}"] for c in group]]
                    )
    for key, table in rows.items():
        with (DATA / TABLES[key][0]).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(TABLES[key][1])
            writer.writerows(table)
    # Keep all text, comparison parameters, blank-layout information and notes
    # recoverable from the original workbook; no subset overwrites the source.
    (SOURCE.parent / "cached-cells.json").write_text(
        json.dumps(sheets, indent=2) + "\n"
    )
    return rows


def load(key: str) -> list[dict]:
    with (DATA / TABLES[key][0]).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    strings = {"source_cell", "run"}
    return [
        {k: v if k in strings else float(v) for k, v in row.items()} for row in rows
    ]


def publisher_supplement_check(path=SUPPLEMENT):
    """Compare unchanged publisher tables with cached final-deposit cells.

    Half a printed last digit is allowed for numerical rounding. Differences
    remain provenance findings; neither original is overwritten or repaired.
    """
    if hashlib.sha256(path.read_bytes()).hexdigest() != SUPPLEMENT_SHA256:
        raise ValueError("This is not the recovered publisher supplement")
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
        tables = root.find("w:body", ns).findall("w:tbl", ns)
        embedded = [
            name for name in archive.namelist() if name.startswith("word/embeddings/")
        ]
    rows = [
        [
            [
                "".join(n.text or "" for n in c.iter() if n.tag.endswith("}t")).strip()
                for c in row.findall("w:tc", ns)
            ]
            for row in table.findall("w:tr", ns)
        ]
        for table in tables
    ]
    if [len(table) for table in rows] != [22, 13, 48, 50]:
        raise ValueError("Unexpected publisher S1-S4 table layout")
    sheets = workbook_cells(SOURCE)
    mismatches, checked = [], 0
    for index, table in enumerate(rows):
        sheet = f"Table S{index + 1}"
        for i, row in enumerate(table):
            if i == 0:
                continue
            # S1 has a second header and footnote; S2 comparisons have prose.
            if index == 0 and i in {12, 21}:
                continue
            if index == 1 and i > 2:
                continue  # Only the two Ma parameterizations; other studies reorder.
            columns = (
                range(3, 12) if index == 1 else range(1 if index == 0 else 0, len(row))
            )
            workbook_columns = (
                "ABCDEFGHIJKLMNOP" if index != 0 or i < 12 else "ABCEGHKNP"
            )
            for j in columns:
                if j >= len(row) or not row[j]:
                    continue
                cell = f"{workbook_columns[j]}{i + 2}"
                source = sheets[sheet].get(cell)
                if source is None:
                    continue
                printed = row[j]
                try:
                    value = Decimal(printed)
                except InvalidOperation:
                    if index != 1 or not isinstance(source, str):
                        continue
                    differs = printed.replace(" ", "") != source.replace(" ", "")
                else:
                    if not isinstance(source, float):
                        continue
                    half_digit = Decimal(5).scaleb(value.as_tuple().exponent - 1)
                    differs = abs(Decimal(str(source)) - value) > half_digit + Decimal(
                        "1e-12"
                    )
                checked += 1
                if differs:
                    mismatches.append(
                        {
                            "table": sheet,
                            "workbook_cell": cell,
                            "publisher_text": printed,
                            "deposited_value": source,
                        }
                    )
    return {
        "sha256": SUPPLEMENT_SHA256,
        "table_rows_including_headers": [len(t) for t in rows],
        "checked_cells": checked,
        "differences_beyond_printed_rounding": mismatches,
        "embedded_files": embedded,
        "qualification": "Publisher S1/S3/S4 numerical cells and the two Ma S2 parameter rows compared with unchanged final workbook. Other-study S2 rows reorder and are not compared by position. Text S1 and Figures S1-S3 were visually inspected; no recalibration algorithm, regression objective, weights or numerical covariance are supplied.",
    }


def matsui_nacl_pressure(cell_volume_a3, temperature=300.0):
    """Matsui 2012 Table 2 BM4+MGD, conventional NaCl cell Z=4."""
    v0, k0, kp, kpp = 179.425, 23.7, 5.14, -0.392
    f = 0.5 * ((v0 / cell_volume_a3) ** (2.0 / 3.0) - 1.0)
    aa = 1.5 * (kp - 4.0)
    bb = (9.0 * k0 * kpp + 9.0 * kp**2 - 63.0 * kp + 143.0) / 6.0
    cold = 3.0 * f * k0 * (1.0 + 2.0 * f) ** 2.5 * (1.0 + aa * f + bb * f**2)
    molar = cell_volume_a3 * Avogadro * 1e-24 / 4.0
    return cold + thermal_pressure(
        molar, temperature, v0 * Avogadro * 1e-24 / 4.0, 1.56, 279.0, 0.96
    )


def pressure(volume, parameters=PUBLISHED):
    """Independent molar-volume BM3, Equation 16, in GPa."""
    v0, k0, kp = parameters
    ratio = np.asarray(volume) / v0
    f = ratio ** (-2.0 / 3.0) - 1.0
    return 1.5 * k0 * f * (1.0 + 0.75 * (kp - 4.0) * f) * ratio ** (-5.0 / 3.0)


def bulk_modulus(volume, parameters=PUBLISHED):
    """Independent KT(V), Equation 15, in GPa."""
    v0, k0, kp = parameters
    ratio = np.asarray(volume) / v0
    f = 1.0 - ratio ** (-2.0 / 3.0)
    return (
        k0
        * (1.0 + 0.5 * (5.0 - 3.0 * kp) * f + 27.0 / 8.0 * (kp - 4.0) * f**2)
        * ratio ** (-5.0 / 3.0)
    )


def debye_function(x):
    return (
        3.0
        * quad(lambda t: t**3 / np.expm1(t) if t else 0.0, 0.0, x, epsabs=1e-11)[0]
        / x**3
    )


def thermal_pressure(volume, temperature, v0=32.48, gamma0=1.92, theta0=251.0, q=1.0):
    """Equations 13-14 and 18-20, molar cm3 volume, two atoms, 300 K baseline."""
    gamma = gamma0 * (volume / v0) ** q
    theta = theta0 * np.exp(gamma0 / q * (1.0 - (volume / v0) ** q))
    energy_difference = (
        6.0
        * R
        * (
            temperature * debye_function(theta / temperature)
            - 300.0 * debye_function(theta / 300.0)
        )
    )
    # J/cm3 = MPa. The zero-point term cancels at fixed volume.
    return gamma / volume * energy_difference * 1e-3


def acoustic_reduction(row):
    """Check the printed Eq 9 and the atom-count convention used by source S1."""
    volume = row["molar_volume_cm3_mol"]
    vp, vs = row["vp_km_s"] * 1000.0, row["vs_km_s"] * 1000.0
    rho = 0.0745513 / (volume * 1e-6)  # assumed KCl kg/mol; source mass not deposited
    ks = rho * (vp**2 - 4.0 / 3.0 * vs**2) * 1e-9
    vd = (3.0 / (2.0 * vs**-3 + vp**-3)) ** (1.0 / 3.0)
    printed_theta = (
        hbar
        / (2.0 * Boltzmann)
        * (6.0 * np.pi**2 * Avogadro / (volume * 1e-6)) ** (1.0 / 3.0)
        * vd
    )
    theta = printed_theta * 2.0 * 2.0 ** (1.0 / 3.0)
    x = theta / 300.0
    cv = 6.0 * R * (4.0 * debye_function(x) - 3.0 * x / np.expm1(x))
    gamma = 1.92 * volume / 32.48
    kt = ks - gamma**2 * cv * 300.0 / volume * 1e-3
    return {
        "ks_gpa": float(ks),
        "printed_eq9_theta_d_k": float(printed_theta),
        "two_atom_theta_d_k": float(theta),
        "kt_gpa": float(kt),
        "source_theta_difference_k": float(theta - row["theta_d_k"]),
        "source_kt_difference_gpa": float(kt - row["kt_gpa"]),
    }


def walker_input_check():
    """Trace Ma's eight reductions to original Walker calibrant observations.

    Independently apply Matsui Table 2 BM4+MGD at the measured temperature,
    then use Walker's alphaKT to bring KCl pressure to 300 K. This explicit
    diagnostic does not silently replace Ma's deposited reduced input.
    """
    with (DATA / "kcl-walker-2002-table2-pvt.csv").open(newline="") as stream:
        raw = [
            r for r in csv.DictReader(stream) if float(r["temperature_celsius"]) <= 24.0
        ]
    results = []
    for row in load("walker"):
        candidates = [
            r
            for r in raw
            if abs(float(r["pressure_kbar"]) / 10.0 - row["walker_pressure_gpa"]) < 1e-8
            and float(r["temperature_celsius"]) + 273.15 == row["source_temperature_k"]
            and round(float(r["b2_kcl_cell_volume_a3"]) * Avogadro * 1e-24, 4)
            == row["molar_volume_cm3_mol"]
        ]
        if len(candidates) != 1:
            raise ValueError(
                "Walker pressure/temperature/volume row mapping is ambiguous"
            )
        original = candidates[0]
        cell = float(original["nacl_lattice_a_angstrom"]) ** 3
        lattice = float(original["nacl_lattice_a_angstrom"])
        p300 = matsui_nacl_pressure(cell)
        temperature = float(original["temperature_celsius"]) + 273.15
        pt = matsui_nacl_pressure(cell, temperature) - p300
        kcl_increment = 0.00275 * (300.0 - temperature)
        corrected = p300 + pt + kcl_increment
        # This inversion is a conditional diagnostic, not recovered source data.
        implied_lattice = brentq(
            lambda a: matsui_nacl_pressure(a**3) - row["matsui_300k_pressure_gpa"],
            5.0,
            5.7,
        )
        rounding_bound = max(
            abs(matsui_nacl_pressure((lattice + delta) ** 3) - p300)
            for delta in [-0.00005, 0.00005]
        )
        results.append(
            {
                "source_row": row["source_row"],
                "nacl_file": original["nacl_file"],
                "kcl_file": original["kcl_file"],
                "nacl_lattice_a_angstrom": float(original["nacl_lattice_a_angstrom"]),
                "nacl_lattice_a_esd_angstrom": float(
                    original["nacl_lattice_a_esd_angstrom"]
                ),
                "b2_kcl_cell_volume_a3": float(original["b2_kcl_cell_volume_a3"]),
                "b2_kcl_cell_volume_esd_a3": float(
                    original["b2_kcl_cell_volume_esd_a3"]
                ),
                "kcl_molar_volume_from_original_cell_cm3_mol": float(
                    original["b2_kcl_cell_volume_a3"]
                )
                * Avogadro
                * 1e-24,
                "matsui_300k_at_original_nacl_volume_gpa": float(p300),
                "source_temperature_k": temperature,
                "nacl_temperature_increment_gpa": float(pt),
                "kcl_temperature_increment_to_300k_gpa": kcl_increment,
                "printed_lattice_half_digit_pressure_bound_gpa": float(rounding_bound),
                "implied_nacl_lattice_at_300k_angstrom": implied_lattice,
                "implied_lattice_shift_at_300k_angstrom": implied_lattice - lattice,
                "matsui_at_measured_temperature_gpa": float(p300 + pt),
                "kcl_300k_pressure_diagnostic_gpa": float(corrected),
                "ma_deposited_recalibrated_pressure_gpa": row[
                    "matsui_300k_pressure_gpa"
                ],
                "difference_gpa": float(corrected - row["matsui_300k_pressure_gpa"]),
            }
        )
    return {
        "rows": results,
        "matsui_source_doi": "10.2138/am.2012.4136",
        "matsui_source_url": "https://rruff.info/doclib/am/vol97/AM97_1670.pdf",
        "method": "Matsui (2012) Equations 2 and 4-11, Table 2; conventional NaCl cell Z=4; Walker alphaKT=0.00275 GPa/K correction to 300 K.",
        "status": "deposited_inputs_preserved_upstream_recalculation_not_reproduced",
        "qualification": "Original Walker lattice observations do not reproduce Ma's deposited recalibrated pressures. Half a printed lattice digit cannot explain the difference; the implied lattice is conditional on direct Matsui evaluation at 300 K, not an inferred measurement. The actual correction prescription and upstream inputs remain unavailable. Preserve the deposited reductions in the source joint fit.",
    }


def joint_fit(mode="reported_errors"):
    """Fit all 11 KT-V and eight recalibrated P-V rows; never fit model grids.

    Source weighting/confidence/covariance are not deposited. Three transparent
    sensitivity objectives are offered, not an assertion of the authors' code.
    """
    a, b = load("acoustic"), load("walker")
    va, ka, ea, ev = (
        np.array([r[c] for r in a])
        for c in [
            "molar_volume_cm3_mol",
            "kt_gpa",
            "kt_error_gpa",
            "molar_volume_error_cm3_mol",
        ]
    )
    vb, pb, eb = (
        np.array([r[c] for r in b])
        for c in [
            "molar_volume_cm3_mol",
            "matsui_300k_pressure_gpa",
            "matsui_pressure_error_gpa",
        ]
    )
    if mode not in {"reported_errors", "unweighted", "effective_errors"}:
        raise ValueError("unknown joint-fit mode")

    def residual(p):
        ra = bulk_modulus(va, p) - ka
        rb = pressure(vb, p) - pb
        if mode != "unweighted":
            scale = ea
            if mode == "effective_errors":
                h = 1e-4
                derivative = (bulk_modulus(va + h, p) - bulk_modulus(va - h, p)) / (
                    2.0 * h
                )
                scale = np.sqrt(ea**2 + (derivative * ev) ** 2)
            ra, rb = ra / scale, rb / eb
        return np.r_[ra, rb]

    fitted = least_squares(
        residual,
        PUBLISHED,
        bounds=([30, 1, 1], [36, 60, 10]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    covariance = (
        np.linalg.inv(fitted.jac.T @ fitted.jac) * np.sum(fitted.fun**2) / (19 - 3)
    )
    return {
        "mode": mode,
        "success": bool(fitted.success),
        "observations": 19,
        "acoustic_rows": 11,
        "walker_rows": 8,
        "parameters_molar": dict(
            zip(["V0_cm3_mol", "K0_gpa", "K0_prime"], fitted.x.tolist())
        ),
        "difference_from_published": (fitted.x - PUBLISHED).tolist(),
        "residual_sum_squares": float(np.sum(fitted.fun**2)),
        "kt_rmse_gpa": float(np.sqrt(np.mean((bulk_modulus(va, fitted.x) - ka) ** 2))),
        "pv_rmse_gpa": float(np.sqrt(np.mean((pressure(vb, fitted.x) - pb) ** 2))),
        "diagnostic_covariance": covariance.tolist(),
        "diagnostic_errors": np.sqrt(np.diag(covariance)).tolist(),
        "covariance_provenance": "Residual-scaled inverse JtJ of this diagnostic objective, not the source covariance or independent physical errors.",
    }


def debye_fit(weighted=True):
    """Diagnostic source Eq 13 fit with q=1, using only acoustic Debye inputs."""
    rows = load("acoustic")
    volumes = np.array([r["molar_volume_cm3_mol"] for r in rows])
    theta = np.array([r["theta_d_k"] for r in rows])
    scale = np.array([r["theta_d_error_k"] for r in rows]) if weighted else 1.0
    fitted = least_squares(
        lambda p: (p[0] * np.exp(-p[1] * volumes) - theta) / scale, [1716, 0.0592]
    )
    aa, bb = fitted.x
    return {
        "weighted_by_reported_errors": weighted,
        "q_fixed": 1.0,
        "A_k": float(aa),
        "B_mol_cm3": float(bb),
        "gamma0": float(bb * 32.48),
        "theta0_k": float(aa * np.exp(-bb * 32.48)),
        "qualification": "Diagnostic Eq 13 fit to all 11 reduced acoustic Debye values; source weights and covariance unknown. Uses the published V0 to map A,B to gamma0,theta0.",
    }


def ledger_outcome(record):
    """Return the dedicated acoustic/Walker outcome for the catalog ledger."""
    audit = reproduce()
    fitted = audit["fits"]["reported_errors"]
    values = fitted["parameters_molar"]
    converted = [
        values["V0_cm3_mol"] / (Avogadro * 1e-24),
        values["K0_gpa"],
        values["K0_prime"],
    ]
    comparison = []
    for i, name in enumerate(["V0", "K0", "K0_prime"]):
        source = record["eos"]["parameters"][name]
        comparison.append(
            {
                "parameter": name,
                "published": source,
                "published_error": record["parameter_errors"][name],
                "refit": converted[i],
                "refit_error": None,
                "difference": converted[i] - source,
                "relative_difference": abs(converted[i] - source) / abs(source),
                "within_combined_2sigma": None,
                "similar": abs(converted[i] - source)
                <= record["parameter_errors"][name],
            }
        )
    return {
        "status": "similar"
        if all(p["similar"] for p in comparison)
        else "parity_not_achieved",
        "dataset_identifiers": record["fit_datasets"],
        "observations": 19,
        "selection": audit["selection"],
        "fit_kind": "joint_acoustic_kt_volume_and_recalibrated_walker_pv",
        "objective": "KT and pressure residuals divided by deposited row error widths; all 19 rows, V0/K0/K0_prime free",
        "absolute_sigma": False,
        "parameters": comparison,
        "free_parameters": ["V0", "K0", "K0_prime"],
        "rmse_gpa": fitted["pv_rmse_gpa"],
        "kt_rmse_gpa": fitted["kt_rmse_gpa"],
        "solver_success": fitted["success"],
        "solver_message": "dedicated Ma acoustic/Walker diagnostic completed",
        "qualification": "All three BM3 coefficients lie within their published error widths, but unknown source weights/error confidence/covariance prevent strict statistical parity. Gamma0/theta0 receive a separate Eq 13 acoustic diagnostic; high-temperature use is modeled. Direct upstream Walker recalibration remains unresolved. See the [Ma audit](literature-reproductions/ma-2024-kcl.md).",
        "source_grid_validation": audit["source_grid_validation"],
        "debye_temperature_fits": audit["debye_temperature_fits"],
        "upstream_recalibration_status": audit["walker_input_check"]["status"],
    }


def reproduce():
    cold, thermal = load("cold_grid"), load("thermal_grid")
    reductions = [{"run": r["run"], **acoustic_reduction(r)} for r in load("acoustic")]
    v = np.array([r["molar_volume_cm3_mol"] for r in cold])
    p = np.array([r["model_pressure_gpa"] for r in cold])
    precision = least_squares(
        lambda coefficients: pressure(v, coefficients) - p,
        PUBLISHED,
        xtol=1e-12,
        gtol=1e-12,
        ftol=1e-12,
    )
    # First-order P(V) variance has only five independent volume powers for
    # three BM3 coefficients; six covariance entries cannot be recovered.
    jac = np.column_stack(
        [
            np.imag(pressure(v, precision.x.astype(complex) + 1e-20j * direction))
            / 1e-20
            for direction in np.eye(3)
        ]
    )
    variance_design = np.column_stack(
        [
            jac[:, 0] ** 2,
            jac[:, 1] ** 2,
            jac[:, 2] ** 2,
            2 * jac[:, 0] * jac[:, 1],
            2 * jac[:, 0] * jac[:, 2],
            2 * jac[:, 1] * jac[:, 2],
        ]
    )
    temperatures = np.array([r["temperature_k"] for r in thermal])
    thermal_precision = least_squares(
        lambda coefficients: (
            np.array(
                [
                    thermal_pressure(
                        32.48, t, gamma0=coefficients[0], theta0=coefficients[1]
                    )
                    for t in temperatures
                ]
            )
            - np.array([r["model_thermal_pressure_gpa"] for r in thermal])
        ),
        [1.92, 251.0],
        xtol=1e-12,
        gtol=1e-12,
        ftol=1e-12,
    )
    return {
        "record_identifier": RECORD_ID,
        "source": {
            "article_doi": "10.1029/2024JB028819",
            "dataset_doi": "10.17632/7svmv9hvft.1",
            "url": SOURCE_URL,
            "sha256": SOURCE_SHA256,
            "license": "CC-BY-4.0",
            "excluded_older_deposit": "10.17632/mws6hnp49j.1",
        },
        "selection": "11 measured acoustic-volume rows reduced to KT plus eight Matsui-2012-recalibrated Walker-2002 300 K P-V rows; S1 modeled pressures and S3-S4 grids excluded from all fits.",
        "fits": {
            mode: joint_fit(mode)
            for mode in ["reported_errors", "effective_errors", "unweighted"]
        },
        "debye_temperature_fits": {
            "reported_errors": debye_fit(),
            "unweighted": debye_fit(False),
        },
        "source_grid_validation": {
            "cold_rows": len(cold),
            "thermal_rows": len(thermal),
            "cold_max_abs_error_gpa": max(
                abs(
                    float(pressure(r["molar_volume_cm3_mol"])) - r["model_pressure_gpa"]
                )
                for r in cold
            ),
            "thermal_max_abs_error_gpa": max(
                abs(
                    thermal_pressure(32.48, r["temperature_k"])
                    - r["model_thermal_pressure_gpa"]
                )
                for r in thermal
            ),
            "cold_max_error_over_reported_grid_error": max(
                abs(
                    float(pressure(r["molar_volume_cm3_mol"])) - r["model_pressure_gpa"]
                )
                / r["model_pressure_error_gpa"]
                for r in cold
            ),
            "grid_precision_diagnostic": {
                "cold_parameters_molar": precision.x.tolist(),
                "cold_max_abs_residual_gpa": float(np.max(np.abs(precision.fun))),
                "thermal_gamma0_theta0": thermal_precision.x.tolist(),
                "thermal_max_abs_residual_gpa": float(
                    np.max(np.abs(thermal_precision.fun))
                ),
                "meaning": "Reverse reconstruction of source grid precision only, not a refit of observations, source covariance, or selectable EOS. Published rounded coefficients remain authoritative.",
            },
            "qualification": "Grids are independently deposited benchmarks, not fit observations. All inferred grid coefficients round to the published Table 1 values; their precision explains the small grid discrepancy.",
        },
        "acoustic_reductions": reductions,
        "publisher_supplement_check": publisher_supplement_check(),
        "cold_grid_covariance_identifiability": {
            "symmetric_covariance_entries": 6,
            "first_order_variance_design_rank": int(
                np.linalg.matrix_rank(variance_design)
            ),
            "qualification": "S3 pressure-error envelopes alone cannot uniquely recover the six BM3 covariance entries under first-order propagation; no covariance is inferred or attached to the source EOS.",
        },
        "walker_input_check": walker_input_check(),
        "limitations": [
            "Source regression objective, weights, confidence of coefficient errors and full covariance are unavailable; diagnostic objectives do not recreate source covariance.",
            "Only room-temperature acoustic measurements to 85 GPa constrain the source; high-temperature use is a quasi-harmonic model extension, and higher pressures extrapolate the acoustic coverage.",
            "Equation 9 prints hbar/(2 kB) and Avogadro number without the two-atom factor. Source S1 values instead agree with hbar/kB and the two-atom number density; printed Equation 9 is not used to replace published theta0.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract-source", type=Path)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "docs/data/ma-2024-kcl-reproduction.json"
    )
    args = parser.parse_args()
    if args.extract_source:
        extract(args.extract_source)
    report = reproduce()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "grids": report["source_grid_validation"],
                "joint_fits": {
                    k: v["parameters_molar"] for k, v in report["fits"].items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
