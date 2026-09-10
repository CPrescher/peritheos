#!/usr/bin/env python3
"""Independent Fei (2016) equation reconstruction and partial-data diagnostics."""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
DOI = "10.1002/2016GL069456"
NA = 6.02214076e23
MASS = 55.845
R = 8.31446261815324
RHO0 = 8.2695
V0 = 2 * MASS / (NA * 1e-24 * RHO0)
RECORDS = (
    "iron_fei_2016_bm3",
    "iron_fei_2016_free_density_bm3",
    "iron_fei_2016_bm3_debye_quadratic",
)
STEMS = {
    "s1": "grl54635-sup-0002-2016gl069456-ts01.xlsx",
    "s2": "grl54635-sup-0003-2016gl069456-ts02.xlsx",
}
COLUMNS = {
    "s1": [
        "sample_name",
        "pressure_gpa",
        "a_angstrom",
        "a_error_angstrom",
        "c_angstrom",
        "c_error_angstrom",
        "volume_a3_conventional_cell",
        "density_g_cm3",
        "marker_a_angstrom",
        "marker_a_error_angstrom",
        "pressure_standard",
    ],
    "s2": [
        "sample_name",
        "temperature_k",
        "temperature_error_k",
        "pressure_gpa",
        "a_angstrom",
        "a_error_angstrom",
        "c_angstrom",
        "c_error_angstrom",
        "volume_a3_conventional_cell",
        "density_g_cm3",
        "thermal_pressure_gpa",
    ],
}


def csv_path(table):
    return DATA / f"iron-fei-2016-table-{table}.csv"


def normalized_csv(path, table):
    """Read the official sheet's cached numbers verbatim, including blank errors.

    XLSX XML extraction avoids recalculation or editing of the user's workbook.
    Numeric rows are identified by the published sample-name prefix.
    """
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as archive:
        strings = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = ["".join(s.itertext()) for s in strings]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["source_row", *COLUMNS[table]])
    for row in sheet.findall(".//x:sheetData/x:row", ns):
        values = [""] * 11
        for cell in row.findall("x:c", ns):
            column = "".join(c for c in cell.attrib["r"] if c.isalpha())
            value = cell.findtext("x:v", "", ns)
            if cell.attrib.get("t") == "s":
                value = strings[int(value)]
            if len(column) == 1 and "A" <= column <= "K":
                values[ord(column) - ord("A")] = value
        if values[0].startswith("c1"):
            writer.writerow([row.attrib["r"], *values])
    return output.getvalue()


def rows(table):
    with csv_path(table).open() as stream:
        return list(csv.DictReader(stream))


def bm3(volume, v0=V0, k0=172.7, kp=4.79):
    eta = (v0 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def thermal_pressure(volume, temperature, beta0=0.07, k=1.34):
    """Equation (2) evaluated independently on the conventional two-Fe cell."""
    ratio = volume / V0
    gamma = 1.74 * ratio**0.78
    theta = 422 * np.exp(1.74 / 0.78 * (1 - ratio**0.78))

    def energy(t):
        y = theta / t
        integral = quad(lambda x: x**3 / np.expm1(x), 0, y, epsabs=1e-11, epsrel=1e-11)[
            0
        ]
        return 9 * R * t * integral / y**3

    molar_volume_si = volume * NA * 1e-30 / 2
    vibrational = gamma * (energy(temperature) - energy(300)) / molar_volume_si / 1e9
    density_si = RHO0 * 1000 / ratio
    electronic = 2 * beta0 * ratio**k * density_si * (temperature**2 - 300**2) / 2e9
    return vibrational + electronic


def reproduce():
    static, thermal = rows("s1"), rows("s2")
    v = np.array([float(r["volume_a3_conventional_cell"]) for r in static])
    p = np.array([float(r["pressure_gpa"]) for r in static])
    fits = {}
    for free_density in (False, True):
        start = [V0, 172.7, 4.79] if free_density else [172.7, 4.79]

        def residual(x):
            return bm3(v, *x) - p if free_density else bm3(v, V0, *x) - p

        fit = least_squares(residual, start, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        params = dict(
            zip(
                ["V0", "K0", "K0_prime"] if free_density else ["K0", "K0_prime"],
                fit.x.tolist(),
            )
        )
        source = (
            (2 * MASS / (NA * 1e-24 * 8.3602), 191.44, 4.52)
            if free_density
            else (V0, 172.7, 4.79)
        )
        fits["free_density" if free_density else "fixed_density"] = {
            "observations": len(v),
            "parameters": params,
            "published_rmse_gpa": float(np.sqrt(np.mean((bm3(v, *source) - p) ** 2))),
            "refit_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        }
    tv = np.array([float(r["volume_a3_conventional_cell"]) for r in thermal])
    tt = np.array([float(r["temperature_k"]) for r in thermal])
    tp = np.array([float(r["thermal_pressure_gpa"]) for r in thermal])
    predicted = np.array([thermal_pressure(vv, t) for vv, t in zip(tv, tt)])
    full = np.array([float(r["pressure_gpa"]) for r in thermal])
    # Only static thermal rows are available. Fit A with m fixed; fitting the
    # source's beta0,k pair would falsely imply recovery of its shock constraints.
    unit = np.array(
        [
            thermal_pressure(vv, t, 1.0) - thermal_pressure(vv, t, 0.0)
            for vv, t in zip(tv, tt)
        ]
    )
    debye = predicted - 0.07 * unit
    beta_fit = float(unit @ (tp - debye) / (unit @ unit))
    return {
        "doi": DOI,
        "scope": "Partial supplementary-data diagnostics, not the published combined fits. Prior-study 300 K rows and shock P-V-T reductions/weights are not supplied by these supplements. No diagnostic coefficients replace published records.",
        "reference_volume_a3": V0,
        "static": fits,
        "thermal": {
            "observations": len(tv),
            "published_pth_rmse_gpa": float(np.sqrt(np.mean((predicted - tp) ** 2))),
            "published_total_pressure_rmse_gpa": float(
                np.sqrt(np.mean((bm3(tv) + predicted - full) ** 2))
            ),
            "s2_pth_vs_p_minus_bm3_max_abs_gpa": float(
                np.max(np.abs(tp - (full - bm3(tv))))
            ),
            "partial_static_only_beta0_k_fixed": beta_fit,
            "partial_refit_pth_rmse_gpa": float(
                np.sqrt(np.mean((debye + beta_fit * unit - tp) ** 2))
            ),
            "first_row_predicted_pth_gpa": float(predicted[0]),
            "last_row_predicted_pth_gpa": float(predicted[-1]),
        },
    }


def ledger_outcome(record):
    result = reproduce()
    return {
        "status": "not_refittable",
        "dataset_identifiers": record["fit_datasets"],
        "reason": result["scope"]
        + " Independent equation and partial-data checks: [Fei 2016 audit](literature-reproductions/fei-2016-iron.md).",
        "partial_validation": result["thermal"]
        if "thermal" in record
        else result["static"][
            "free_density"
            if "free_density" in record["identifier"]
            else "fixed_density"
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.source_directory:
        for table, filename in STEMS.items():
            content = normalized_csv(args.source_directory / filename, table)
            if args.write:
                csv_path(table).write_text(content)
            else:
                assert content == csv_path(table).read_text(), filename
    result = reproduce()
    if args.write:
        path = ROOT / "docs/data/fei-2016-iron-reproduction.json"
        path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
