"""Build phase-unresolved material cards for Sun et al. (2010).

This is a deterministic catalog-generation helper.  The selection order follows
the source's own model ranking: all SMS4 curves, then SMS3, then MRS3.  Brass is
never emitted because the source does not define its composition.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from importlib import resources
from pathlib import Path

AVOGADRO_CM3_TO_A3 = 1.6605390671738466
DOI = "10.1515/zna-2010-1-202"
PDF_SHA256 = "57b6c4ff8b59dcc3734aa9507de0f62610839363a8e5ae790ed8f71e6fa66db4"
TABLE_SHA256 = "e71613962b7df52559b8633546a7513b902ce1b5d8d4a7666e0476464123ca68"
MOLAR_MASS = {
    "H2": 2.01588,
    "Cu": 63.546,
    "Mo": 95.95,
    "W": 183.84,
    "Zn": 65.38,
    "Ag": 107.8682,
    "Pt": 195.084,
    "Ti": 47.867,
    "Ta": 180.94788,
    "Au": 196.96657,
    "Pd": 106.42,
    "Zr": 91.224,
    "Cr": 51.9961,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Al2O3": 101.960077,
    "Nb": 92.90637,
    "Cd": 112.414,
    "Al": 26.9815385,
    "Th": 232.0377,
    "V": 50.9415,
    "In": 114.818,
    "MgO": 40.304,
    "Be": 9.0121831,
    "LiF": 25.939403,
    "Pb": 207.2,
    "Sn": 118.710,
    "Mg": 24.305,
    "CsBr": 212.809451,
    "Ca": 40.078,
    "Tl": 204.38,
    "NaCl": 58.44276928,
    "LiI": 133.84547,
    "LiBr": 86.845,
    "NaBr": 102.89376928,
    "NaI": 149.89423928,
    "KF": 58.096703,
    "RbF": 104.466203,
    "LiCl": 42.394,
    "Li": 6.94,
    "Na": 22.98976928,
    "KI": 166.00277,
    "RbI": 212.37277,
    "RbBr": 165.3718,
    "K": 39.0983,
    "Rb": 85.4678,
    "NaF": 41.988172443,
    "RbCl": 120.9208,
    "Nd": 144.242,
}
MODEL = {
    "sms4": ("SunMorse4", "sun_morse_4", "SMS4", "preferred benchmark family"),
    "sms3": ("SunMorse3", "sun_morse_3", "SMS3", "same-data model sensitivity"),
    "mrs3": ("Morse3", "morse_3", "MRS3", "same-data model sensitivity"),
}


def slug(value: str) -> str:
    value = value.lower().replace("n-h2", "solid_h2")
    return re.sub(r"[^a-z0-9]+", "_", value).strip("_")


def source_rows() -> list[dict[str, str]]:
    resource = resources.files("peritheos").joinpath(
        "data", "datasets", "sun-2010-table2-morse-eos-parameters.csv"
    )
    with resource.open(encoding="utf-8", newline="") as handle:
        return [row for row in csv.DictReader(handle) if row["solid"] != "Brass"]


def selected_curves(limit: int) -> dict[str, list[str]]:
    rows = source_rows()
    available = len(rows) * len(MODEL)
    if not len(rows) <= limit <= available:
        raise ValueError(f"record limit must be in [{len(rows)}, {available}]")
    ordered = [(row["solid"], prefix) for prefix in MODEL for row in rows]
    selected: dict[str, list[str]] = {row["solid"]: [] for row in rows}
    for solid, prefix in ordered[:limit]:
        selected[solid].append(prefix)
    return selected


def reference() -> dict[str, object]:
    return {
        "authors": ["Sun", "Wu", "Guo", "Cai"],
        "year": 2010,
        "title": "Two Universal Equations of State for Solids",
        "source": "Zeitschrift fuer Naturforschung A",
        "volume": "65",
        "locator": "34-44",
        "doi": DOI,
    }


def record(row: dict[str, str], prefix: str) -> dict[str, object]:
    formula = row["formula"]
    material_slug = slug(row["solid"])
    eos_type, model, source_name, role = MODEL[prefix]
    pressure_max = float(row["pressure_max_gpa"])
    return {
        "identifier": f"{material_slug}_sun_2010_{prefix}_1",
        "label": f"Sun et al. (2010), {row['solid']} {source_name} legacy benchmark",
        "record_kind": "published",
        "equation_kind": "isothermal",
        "reference": reference(),
        "volume_basis": {
            "kind": "formula_units",
            "formula_units": 1,
            "molar_mass_g_mol": MOLAR_MASS[formula],
        },
        "eos": {
            "type": eos_type,
            "model": model,
            "parameters": {
                "V0": float(row["v0_cm3_mol"]) * AVOGADRO_CM3_TO_A3,
                "K0": float(row[f"{prefix}_k0_gpa"]),
                "K0_prime": float(row[f"{prefix}_k0_prime"]),
            },
        },
        "parameter_errors": {"V0": None, "K0": None, "K0_prime": None},
        "parameter_error_confidence": None,
        "parameter_covariance": None,
        "fixed_parameters": [],
        "experimental_pressure_range_gpa": [0.0, pressure_max],
        "pressure_range_status": "reported_exactly",
        "validity": {
            "pressure_gpa": [0.0, pressure_max],
            "notes": [
                "This is the exact source fit range, not a verified single-phase stability range.",
                "The source does not republish row-level volumes, weights, or pressure calibrations.",
            ],
        },
        "notes": (
            f"Table 2 {source_name} {role}. The source fits a legacy normalized-"
            "compression compilation and reports no parameter uncertainties or covariance. "
            "The material card is phase-unresolved because the article does not preserve "
            "phase-resolved observation lineage."
        ),
        "pressure_calibration": {
            "status": "unresolved",
            "methods": [],
            "recalculation": {
                "status": "not_possible",
                "notes": "Exact row-wise legacy observations and their pressure calibrations are not republished.",
            },
            "audit_date": "2026-09-05",
            "notes": "No pressure scale is inferred from the compiled material label or publication era.",
        },
        "parameter_provenance": {
            "source_solid": row["solid"],
            "table_prefix": prefix,
            "V0": "Table 2, experimental molar volume converted with exact Avogadro constant",
            "K0": f"Table 2, {source_name}",
            "K0_prime": f"Table 2, {source_name}",
            "fit_range": "Table 1",
            "mean_pressure_error_gpa": float(row[f"{prefix}_mean_pressure_error_gpa"]),
        },
        "scientific_validation": {
            "status": "primary_source_validated",
            "note": "Equation, coefficients, molar-volume basis, fit range, and reported mean pressure error were checked directly against the primary article.",
            "audit_date": "2026-09-05",
            "verified_fields": [
                "equation",
                "parameters",
                "units",
                "reference_state",
                "validity",
            ],
            "primary_source_check": {
                "access_url": f"https://doi.org/{DOI}",
                "locations": ["Equations (3)-(15)", "Tables 1-2"],
                "finding": "The paper explicitly defines all three models and prints V0, K0, K0-prime, the pressure range, and mean pressure error for every fitted solid.",
                "doi": DOI,
            },
            "primary_data_check": {
                "status": "parameterization_only",
                "audit_date": "2026-09-05",
                "source_locations": [
                    "Tables 1-2; cited legacy compression compilations"
                ],
                "finding": "Source coefficients are complete and executable, but the fitted row-level compression observations and regression weights are not republished.",
            },
        },
        "source_lineage": [
            {
                "role": "scientific equation, fit, and coefficients",
                "citation": f"Sun et al. (2010), equations (6)-(9), Tables 1-2 ({source_name})",
                "doi": DOI,
            },
            {
                "role": "legacy normalized-compression observations",
                "citation": "Kennedy and Keeler (1972), except the source's separately cited H2, W, and NaCl datasets",
            },
        ],
    }


def document(row: dict[str, str], prefixes: list[str]) -> dict[str, object]:
    formula = row["formula"]
    material_slug = slug(row["solid"])
    return {
        "format": "peritheos.material",
        "format_version": 3,
        "identifier": f"{material_slug}_sun_2010_legacy",
        "name": f"{row['solid']} legacy compression aggregate (Sun 2010)",
        "formula": formula,
        "phase": "legacy multi-source compression aggregate; phase unresolved",
        "aliases": [row["solid"]],
        "cell_contents": "One chemical formula unit; no crystallographic conventional cell is inferred.",
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        },
        "formula_units_per_cell": 1,
        "atom_sites": [],
        "peaks": [],
        "source": {
            "article_url": f"https://doi.org/{DOI}",
            "article_sha256": PDF_SHA256,
            "parameter_table": "sun-2010-table2-morse-eos-parameters.csv",
            "parameter_table_sha256": TABLE_SHA256,
        },
        "notes": "Separate phase-unresolved card for historical universal-EOS benchmark fits; it must not be merged with a specific crystallographic phase or used as a modern pressure standard.",
        "eos_records": [record(row, prefix) for prefix in prefixes],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--record-limit", type=int, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("peritheos/data/materials"),
    )
    args = parser.parse_args()
    selected = selected_curves(args.record_limit)
    args.output.mkdir(parents=True, exist_ok=True)
    written = 0
    for row in source_rows():
        prefixes = selected[row["solid"]]
        path = args.output / f"{slug(row['solid'])}_sun_2010_legacy.eosmat"
        path.write_text(json.dumps(document(row, prefixes), indent=1) + "\n")
        written += len(prefixes)
    print(f"wrote {written} records across {len(selected)} material documents")


if __name__ == "__main__":
    main()
