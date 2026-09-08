#!/usr/bin/env python3
"""Audit the CaSiO3 literature-data refit reported by Fu et al. (2023).

The source tables are intentionally not bundled.  Supply a ``pdftotext -layout``
rendering of Sun et al. (2016) and the official Gréaux et al. (2019) source-data
workbook.  This script validates their checksums, reconstructs the row selection
visible in Fu Figure S3, and explores explicit weighting choices for Fu equations
17, 18, and 23-28.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

SUN_PDF_SHA256 = "10d8f6389b59fd7dcb94e349cbac0390e04a2ec55a815b8d238fd97797fccb2b"
GRE_AUX_SHA256 = "cde0964709863c3b487a054597fd3dd549738aaa9dc3f0c637a996932944479a"
MOLAR_MASS_G_MOL = 116.162
AVOGADRO = 6.02214076e23
R = 8.31446261815324
REFERENCE_TEMPERATURE_K = 300.0
DEBYE_TEMPERATURE_K = 1000.0
ATOM_COUNT = 5.0
PUBLISHED = np.array([45.4, 248.0, 126.0, 1.6, 1.42, 2.65, 1.54])
PARAMETERS = ("V0", "K0", "mu0", "mu0_prime", "gamma0", "q", "eta_s0")
REGISTERED_REFIT_RECORD = (
    "ca_perovskite_fu_2023_candidate_data_unweighted_bm3_mgd_refit"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parenthetical_sigma(value: str, digits: str) -> float:
    decimals = len(value.partition(".")[2])
    return float(digits) * 10.0 ** (-decimals)


def parse_sun_table(text_path: Path) -> list[dict[str, float]]:
    """Parse the three-across Table 1 layout from ``pdftotext -layout`` output."""
    text = text_path.read_text(encoding="utf-8")
    start = text.index("Table 1. Experimental P-V-T Data of CaSiO3 Perovskite")
    end = text.index("\n\n\n\n                  and", start)
    table = text[start:end]
    pattern = re.compile(
        r"(?P<p>\d+(?:\.\d+)?)\s+\((?P<ps>\d+)\)\s+"
        r"(?P<t>\d{4})\s+(?P<v>\d+(?:\.\d+)?)\s+\((?P<vs>\d+)\)"
    )
    rows = []
    for match in pattern.finditer(table):
        p_text, v_text = match["p"], match["v"]
        rows.append(
            {
                "pressure_gpa": float(p_text),
                "pressure_sigma_gpa": _parenthetical_sigma(p_text, match["ps"]),
                "temperature_k": float(match["t"]),
                "volume_a3": float(v_text),
                "volume_sigma_a3": _parenthetical_sigma(v_text, match["vs"]),
            }
        )
    if len(rows) != 144:
        raise ValueError(f"expected 144 Sun Table 1 rows, found {len(rows)}")
    return rows


def _xlsx_rows(path: Path, sheet_name: str) -> list[dict[str, object]]:
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rel_namespace = {
        "r": "http://schemas.openxmlformats.org/package/2006/relationships"
    }
    with zipfile.ZipFile(path) as archive:
        strings_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = [
            "".join(node.text or "" for node in item.findall(".//x:t", namespace))
            for item in strings_root.findall("x:si", namespace)
        ]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            relation.attrib["Id"]: relation.attrib["Target"]
            for relation in relationships.findall("r:Relationship", rel_namespace)
        }
        sheet = next(
            item
            for item in workbook.findall("x:sheets/x:sheet", namespace)
            if item.attrib["name"] == sheet_name
        )
        relation_id = sheet.attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        ]
        target = targets[relation_id].lstrip("/")
        if not target.startswith("xl/"):
            target = f"xl/{target}"
        root = ET.fromstring(archive.read(target))
        rows = []
        for row in root.findall("x:sheetData/x:row", namespace):
            values: dict[str, object] = {}
            for cell in row.findall("x:c", namespace):
                column = re.match(r"[A-Z]+", cell.attrib["r"])[0]
                value = cell.find("x:v", namespace)
                if value is None:
                    continue
                raw: object = value.text
                if cell.attrib.get("t") == "s":
                    raw = strings[int(str(raw))]
                else:
                    raw = float(str(raw))
                values[column] = raw
            rows.append(values)
    return rows


def parse_greaux_table(path: Path) -> list[dict[str, float]]:
    """Read the 34 cubic high-temperature rows in official sheet ``Figure 3b``."""
    rows = []
    for row in _xlsx_rows(path, "Figure 3b"):
        if not all(column in row for column in "CDEFGHIJKLM"):
            continue
        if not isinstance(row["G"], float):
            continue
        rows.append(
            {
                "pressure_nacl_gpa": float(row["C"]),
                "pressure_nacl_sigma_gpa": float(row["D"]),
                "pressure_finite_strain_gpa": float(row["E"]),
                "pressure_finite_strain_sigma_gpa": float(row["F"]),
                "temperature_k": float(row["G"]),
                "vp_km_s": float(row["H"]),
                "vp_sigma_km_s": float(row["I"]),
                "vs_km_s": float(row["J"]),
                "vs_sigma_km_s": float(row["K"]),
                "density_g_cm3": float(row["L"]),
                "density_sigma_g_cm3": float(row["M"]),
            }
        )
    if len(rows) != 34:
        raise ValueError(f"expected 34 Gréaux Figure 3b rows, found {len(rows)}")
    return rows


def _debye_energy(
    volume: float, temperature: float, v0: float, gamma0: float, q: float
) -> float:
    gamma = gamma0 * (volume / v0) ** q
    theta = DEBYE_TEMPERATURE_K * math.exp((gamma0 - gamma) / q)
    upper = theta / temperature
    integral = quad(lambda x: x**3 / math.expm1(x), 0.0, upper)[0]
    return 9.0 * ATOM_COUNT * R * temperature * (temperature / theta) ** 3 * integral


def _debye_cv(
    volume: float, temperature: float, v0: float, gamma0: float, q: float
) -> float:
    gamma = gamma0 * (volume / v0) ** q
    theta = DEBYE_TEMPERATURE_K * math.exp((gamma0 - gamma) / q)
    upper = theta / temperature
    integral = quad(lambda x: x**4 * math.exp(x) / math.expm1(x) ** 2, 0.0, upper)[0]
    return 9.0 * ATOM_COUNT * R * (temperature / theta) ** 3 * integral


def _predictions(
    parameters: np.ndarray,
    volume: np.ndarray,
    temperature: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    v0, k0, mu0, mu0_prime, gamma0, q, eta_s0 = parameters
    f = 0.5 * ((v0 / volume) ** (2.0 / 3.0) - 1.0)
    scale = (1.0 + 2.0 * f) ** 2.5
    pressure_300 = 3.0 * k0 * f * (1.0 + 2.0 * f) ** 2.5
    bulk_300 = scale * (k0 + (12.0 * k0 - 5.0 * k0) * f)
    shear_300 = scale * (
        mu0
        + (3.0 * k0 * mu0_prime - 5.0 * mu0) * f
        + (6.0 * k0 * mu0_prime - 24.0 * k0 - 14.0 * mu0 + 18.0 * k0) * f**2
    )
    energy = np.array(
        [_debye_energy(v, t, v0, gamma0, q) for v, t in zip(volume, temperature)]
    )
    energy_ref = np.array(
        [_debye_energy(v, REFERENCE_TEMPERATURE_K, v0, gamma0, q) for v in volume]
    )
    cv_t = np.array(
        [_debye_cv(v, t, v0, gamma0, q) * t for v, t in zip(volume, temperature)]
    )
    cv_t_ref = np.array(
        [
            _debye_cv(v, REFERENCE_TEMPERATURE_K, v0, gamma0, q)
            * REFERENCE_TEMPERATURE_K
            for v in volume
        ]
    )
    gamma = gamma0 * (volume / v0) ** q
    energy_density = (energy - energy_ref) / (volume * AVOGADRO * 1.0e-21)
    cv_t_density = (cv_t - cv_t_ref) / (volume * AVOGADRO * 1.0e-21)
    pressure = pressure_300 + gamma * energy_density
    bulk = (
        bulk_300 + (gamma + 1.0 - q) * gamma * energy_density - gamma**2 * cv_t_density
    )
    shear = shear_300 - eta_s0 * energy_density
    return pressure, bulk, shear


def _arrays(
    sun_rows: list[dict[str, float]], greaux_rows: list[dict[str, float]]
) -> dict[str, np.ndarray]:
    sun = [row for row in sun_rows if row["temperature_k"] <= 2200.0]
    g = greaux_rows
    rho = np.array([row["density_g_cm3"] for row in g])
    rho_sigma = np.array([row["density_sigma_g_cm3"] for row in g])
    vp = np.array([row["vp_km_s"] for row in g])
    vp_sigma = np.array([row["vp_sigma_km_s"] for row in g])
    vs = np.array([row["vs_km_s"] for row in g])
    vs_sigma = np.array([row["vs_sigma_km_s"] for row in g])
    bulk = rho * (vp**2 - 4.0 * vs**2 / 3.0)
    shear = rho * vs**2
    bulk_sigma = np.sqrt(
        ((vp**2 - 4.0 * vs**2 / 3.0) * rho_sigma) ** 2
        + (2.0 * rho * vp * vp_sigma) ** 2
        + (8.0 * rho * vs * vs_sigma / 3.0) ** 2
    )
    shear_sigma = np.sqrt((vs**2 * rho_sigma) ** 2 + (2.0 * rho * vs * vs_sigma) ** 2)
    return {
        "sun_v": np.array([row["volume_a3"] for row in sun]),
        "sun_t": np.array([row["temperature_k"] for row in sun]),
        "sun_p": np.array([row["pressure_gpa"] for row in sun]),
        "sun_ps": np.array([row["pressure_sigma_gpa"] for row in sun]),
        "sun_vs": np.array([row["volume_sigma_a3"] for row in sun]),
        "g_v": MOLAR_MASS_G_MOL / (0.602214076 * rho),
        "g_t": np.array([row["temperature_k"] for row in g]),
        "g_p": np.array([row["pressure_nacl_gpa"] for row in g]),
        "g_ps": np.array([row["pressure_nacl_sigma_gpa"] for row in g]),
        "g_vs": MOLAR_MASS_G_MOL / (0.602214076 * rho**2) * rho_sigma,
        "g_k": bulk,
        "g_ks": bulk_sigma,
        "g_mu": shear,
        "g_mus": shear_sigma,
    }


def fit_variant(arrays: dict[str, np.ndarray], weighting: str) -> dict[str, object]:
    def residual(parameters: np.ndarray) -> np.ndarray:
        sun_p, _, _ = _predictions(parameters, arrays["sun_v"], arrays["sun_t"])
        gp, gk, gmu = _predictions(parameters, arrays["g_v"], arrays["g_t"])
        groups = [
            sun_p - arrays["sun_p"],
            gp - arrays["g_p"],
            gk - arrays["g_k"],
            gmu - arrays["g_mu"],
        ]
        if weighting == "reported_sigma":
            sigmas = [arrays["sun_ps"], arrays["g_ps"], arrays["g_ks"], arrays["g_mus"]]
            groups = [values / sigma for values, sigma in zip(groups, sigmas)]
        elif weighting == "propagated_pv_sigma_no_temperature":
            step = 1.0e-4
            sun_plus = _predictions(
                parameters, arrays["sun_v"] + step, arrays["sun_t"]
            )[0]
            sun_minus = _predictions(
                parameters, arrays["sun_v"] - step, arrays["sun_t"]
            )[0]
            g_plus = _predictions(parameters, arrays["g_v"] + step, arrays["g_t"])[0]
            g_minus = _predictions(parameters, arrays["g_v"] - step, arrays["g_t"])[0]
            sun_effective_sigma = np.hypot(
                arrays["sun_ps"],
                (sun_plus - sun_minus) * arrays["sun_vs"] / (2.0 * step),
            )
            g_effective_sigma = np.hypot(
                arrays["g_ps"],
                (g_plus - g_minus) * arrays["g_vs"] / (2.0 * step),
            )
            groups = [
                groups[0] / sun_effective_sigma,
                groups[1] / g_effective_sigma,
                groups[2] / arrays["g_ks"],
                groups[3] / arrays["g_mus"],
            ]
        elif weighting == "equal_observable_groups":
            groups = [values / math.sqrt(values.size) for values in groups]
        elif weighting != "unweighted_absolute_gpa":
            raise ValueError(weighting)
        return np.concatenate(groups)

    fit = least_squares(
        residual,
        PUBLISHED,
        bounds=(
            [40.0, 100.0, 50.0, 0.0, 0.1, -5.0, 0.0],
            [50.0, 400.0, 250.0, 5.0, 4.0, 10.0, 5.0],
        ),
        x_scale="jac",
        max_nfev=5000,
    )
    return {
        "objective": weighting,
        "parameters": {name: float(value) for name, value in zip(PARAMETERS, fit.x)},
        "residual_sum_squares": float(np.sum(fit.fun**2)),
        "success": bool(fit.success),
        "message": fit.message,
    }


def audit(
    sun_text: Path, greaux_xlsx: Path, sun_pdf: Path | None = None
) -> dict[str, object]:
    if _sha256(greaux_xlsx) != GRE_AUX_SHA256:
        raise ValueError(
            "Gréaux workbook checksum does not match the audited official file"
        )
    if sun_pdf is not None and _sha256(sun_pdf) != SUN_PDF_SHA256:
        raise ValueError("Sun PDF checksum does not match the audited author copy")
    sun_rows = parse_sun_table(sun_text)
    greaux_rows = parse_greaux_table(greaux_xlsx)
    arrays = _arrays(sun_rows, greaux_rows)
    variants = [
        fit_variant(arrays, weighting)
        for weighting in (
            "unweighted_absolute_gpa",
            "reported_sigma",
            "propagated_pv_sigma_no_temperature",
            "equal_observable_groups",
        )
    ]
    return {
        "format": "peritheos.fu-2023-casio3-refit-audit",
        "format_version": 1,
        "audit_date": "2026-09-08",
        "sources": {
            "fu_2023_official_deposit": {
                "url": "http://www.minsocam.org/MSA/AmMin/TOC/2023/Apr2023_data/AM-23-48435.zip",
                "sha256": "12121bb769862e781e8c232c85f2f3caae83d1a6d37f6e18d76ba6fde61b9486",
                "contents": ["8435_supp1.pdf", "8435fu.cif"],
                "row_table_present": False,
            },
            "sun_2016_table1": {
                "doi": "10.1002/2016JB013062",
                "url": "https://www.jsg.utexas.edu/lin/files/SunLowerMantleEoSJGR2016.pdf",
                "pdf_sha256": SUN_PDF_SHA256,
                "role": "candidate fit observations",
                "redistributed": False,
                "reason": "no reusable table-data license identified",
            },
            "greaux_2019_source_data": {
                "doi": "10.1038/s41586-018-0816-5",
                "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-018-0816-5/MediaObjects/41586_2018_816_MOESM1_ESM.xlsx",
                "xlsx_sha256": GRE_AUX_SHA256,
                "role": "candidate fit observations",
                "redistributed": False,
                "reason": "Springer Nature copyright; no reusable data license identified",
            },
            "kawai_tsuchiya_2015": {
                "doi": "10.1002/2015GL063446",
                "role": "analytical comparison curve; not counted as fit observations",
            },
            "li_2006": {
                "doi": "10.1016/j.pepi.2005.12.006",
                "role": "analytical comparison curve; not counted as fit observations",
            },
            "thomson_2019": {
                "doi": "10.1038/s41586-019-1483-x",
                "role": "explicitly excluded by Fu Figure S3/Text S3",
            },
        },
        "source_rows": {
            "sun_2016_table1_total": len(sun_rows),
            "sun_2016_figure_s3_temperature_subset": int(arrays["sun_p"].size),
            "sun_2016_excluded_above_2200_k": len(sun_rows) - int(arrays["sun_p"].size),
            "greaux_2019_figure3b_cubic": len(greaux_rows),
            "greaux_2019_figure3a_tetragonal_excluded": 13,
        },
        "published_parameters": dict(zip(PARAMETERS, map(float, PUBLISHED))),
        "fit_observation_count": int(arrays["sun_p"].size + arrays["g_p"].size),
        "fit_output_count": int(arrays["sun_p"].size + 3 * arrays["g_p"].size),
        "observed_ranges": {
            "pressure_gpa": [
                float(min(arrays["sun_p"].min(), arrays["g_p"].min())),
                float(max(arrays["sun_p"].max(), arrays["g_p"].max())),
            ],
            "temperature_k": [
                float(min(arrays["sun_t"].min(), arrays["g_t"].min())),
                float(max(arrays["sun_t"].max(), arrays["g_t"].max())),
            ],
            "volume_a3_per_formula_unit": [
                float(min(arrays["sun_v"].min(), arrays["g_v"].min())),
                float(max(arrays["sun_v"].max(), arrays["g_v"].max())),
            ],
        },
        "temperature_row_counts": {
            "sun_2016_all": {
                str(int(temperature)): sum(
                    row["temperature_k"] == temperature for row in sun_rows
                )
                for temperature in sorted({row["temperature_k"] for row in sun_rows})
            },
            "greaux_2019_cubic": {
                str(int(temperature)): sum(
                    row["temperature_k"] == temperature for row in greaux_rows
                )
                for temperature in sorted({row["temperature_k"] for row in greaux_rows})
            },
        },
        "weighting_disclosure": "not published by Fu et al.",
        "registered_refit": {
            "record_identifier": REGISTERED_REFIT_RECORD,
            "objective": "unweighted_absolute_gpa",
            "status": "opt-in Peritheos refit; not a reproduction of Fu's undisclosed regression protocol",
        },
        "model_assumptions": {
            "reference_isotherm": "Fu equation 17 BM3 with K0_prime fixed at 4",
            "shear_reference": "Fu equation 18",
            "thermal_terms": "Fu equations 23-28 with the dimensionally standard Debye heat capacity",
            "greaux_pressure": "PNaCl, the measured pressure column; PFS is a derived finite-strain pressure diagnostic",
            "sun_selection": "all Table 1 rows at 1200-2200 K, matching the Figure S3 legend",
        },
        "sensitivity_fits": variants,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sun-text", type=Path, required=True)
    parser.add_argument("--sun-pdf", type=Path)
    parser.add_argument("--greaux-xlsx", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.sun_text, args.greaux_xlsx, args.sun_pdf)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
