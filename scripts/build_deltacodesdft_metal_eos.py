#!/usr/bin/env python3
"""Import 50 elemental-metal BM3 records from the Delta-project archive.

The official archive stores per-atom BM3 parameters and P1-expanded CIFs.  This
helper preserves those inputs in a compact JSON transcription, converts V0 to
each material document's public cell convention, and either extends an exact
phase match or creates a source-structure material card.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"
DATASETS = ROOT / "peritheos" / "data" / "datasets"
SOURCE_JSON = DATASETS / "deltacodesdft-metal-eos-source.json"
PARAMETER_CSV = DATASETS / "deltacodesdft-metal-eos-parameters.csv"

METALS = (
    "Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo "
    "Tc Ru Rh Pd Ag Cd In Sn Cs Ba Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po"
).split()
FLEUR_DUPLICATES = ("Al", "Fe")

EXISTING_PHASES = {
    "Mg": "magnesium_hcp",
    "Al": "aluminum",
    "Ti": "titanium_alpha",
    "Cr": "chromium",
    "Fe": "fe",
    "Co": "cobalt_hcp",
    "Ni": "nickel",
    "Cu": "copper",
    "Zr": "zirconium_alpha",
    "Nb": "niobium",
    "Mo": "molybdenum",
    "Ru": "ruthenium",
    "Rh": "rhodium",
    "Pd": "palladium",
    "Ag": "silver",
    "Ta": "tantalum",
    "W": "tungsten",
    "Re": "rhenium",
    "Os": "osmium",
    "Ir": "iridium",
    "Pt": "platinum",
    "Au": "gold",
    "Pb": "lead_fcc",
}

NEW_PHASES = {
    "Li": ("lithium_9r", "Lithium (9R)", "9R lithium, source P1-expanded cell"),
    "Be": (
        "beryllium_hcp",
        "Beryllium (hcp)",
        "hcp beryllium, source P1-expanded cell",
    ),
    "Na": ("sodium_9r", "Sodium (9R)", "9R sodium, source P1-expanded cell"),
    "K": ("potassium_bcc", "Potassium (bcc)", "bcc potassium, source P1-expanded cell"),
    "Ca": ("calcium_fcc", "Calcium (fcc)", "fcc calcium, source P1-expanded cell"),
    "Sc": ("scandium_hcp", "Scandium (hcp)", "hcp scandium, source P1-expanded cell"),
    "V": ("vanadium_bcc", "Vanadium (bcc)", "bcc vanadium, source P1-expanded cell"),
    "Mn": (
        "manganese_fcc_afm",
        "Manganese (fcc AFM model)",
        "antiferromagnetic fcc manganese model, source P1-expanded cell",
    ),
    "Zn": ("zinc_hcp", "Zinc (hcp)", "hcp zinc, source P1-expanded cell"),
    "Ga": (
        "gallium_alpha",
        "Gallium (alpha)",
        "alpha gallium, source P1-expanded cell",
    ),
    "Rb": ("rubidium_bcc", "Rubidium (bcc)", "bcc rubidium, source P1-expanded cell"),
    "Sr": (
        "strontium_fcc",
        "Strontium (fcc)",
        "fcc strontium, source P1-expanded cell",
    ),
    "Y": ("yttrium_hcp", "Yttrium (hcp)", "hcp yttrium, source P1-expanded cell"),
    "Tc": (
        "technetium_hcp",
        "Technetium (hcp)",
        "hcp technetium, source P1-expanded cell",
    ),
    "Cd": ("cadmium_hcp", "Cadmium (hcp)", "hcp cadmium, source P1-expanded cell"),
    "In": (
        "indium_bct",
        "Indium (body-centered tetragonal)",
        "body-centered tetragonal indium, source P1-expanded cell",
    ),
    "Sn": (
        "tin_alpha",
        "Tin (alpha)",
        "alpha-tin diamond structure, source P1-expanded cell",
    ),
    "Cs": ("cesium_bcc", "Cesium (bcc)", "bcc cesium, source P1-expanded cell"),
    "Ba": ("barium_bcc", "Barium (bcc)", "bcc barium, source P1-expanded cell"),
    "Lu": ("lutetium_hcp", "Lutetium (hcp)", "hcp lutetium, source P1-expanded cell"),
    "Hf": ("hafnium_hcp", "Hafnium (hcp)", "hcp hafnium, source P1-expanded cell"),
    "Hg": (
        "mercury_bct_dft",
        "Mercury (PBE bct model)",
        "PBE body-centered tetragonal mercury model, source P1-expanded cell",
    ),
    "Tl": ("thallium_hcp", "Thallium (hcp)", "hcp thallium, source P1-expanded cell"),
    "Bi": ("bismuth_a7", "Bismuth (A7)", "A7 bismuth, source P1-expanded cell"),
    "Po": (
        "polonium_simple_cubic",
        "Polonium (simple cubic)",
        "simple-cubic polonium, source P1-expanded cell",
    ),
}

MOLAR_MASS = {
    "Li": 6.94,
    "Be": 9.0121831,
    "Na": 22.98976928,
    "Mg": 24.305,
    "Al": 26.9815385,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955908,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938044,
    "Fe": 55.845,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90584,
    "Zr": 91.224,
    "Nb": 92.90637,
    "Mo": 95.95,
    "Tc": 98.0,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.710,
    "Cs": 132.90545196,
    "Ba": 137.327,
    "Lu": 174.9668,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.96657,
    "Hg": 200.592,
    "Tl": 204.38,
    "Pb": 207.2,
    "Bi": 208.98040,
    "Po": 209.0,
}

ARCHIVE_DOI = "10.24435/materialscloud:5e-mv"
METHOD_DOI = "10.1080/10408436.2013.772503"
CAMPAIGN_DOI = "10.1126/science.aad3000"
ARCHIVE_URL = "https://archive.materialscloud.org/record/2023.133"
AUTHOR_PDF = "https://backoffice.biblio.ugent.be/download/4188194/4188200"

CODE_METHOD = {
    "WIEN2k": "WIEN2k 13.1 all-electron APW+lo, scalar-relativistic PBE",
    "FLEUR": "FLEUR 0.26 all-electron LAPW (+lo), scalar-relativistic PBE",
}
MAGNETIC_STATE = {
    "Cr": "antiferromagnetic spin polarization",
    "Mn": "antiferromagnetic spin polarization",
    "Fe": "ferromagnetic spin polarization",
    "Co": "ferromagnetic spin polarization",
    "Ni": "ferromagnetic spin polarization",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_parameter_file(path: Path) -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 4 or fields[0].startswith("#"):
            continue
        try:
            rows[fields[0]] = {
                "V0_a3_atom": float(fields[1]),
                "K0_gpa": float(fields[2]),
                "K0_prime": float(fields[3]),
            }
        except ValueError:
            continue
    return rows


def cif_number(text: str, key: str) -> float:
    match = re.search(rf"^{re.escape(key)}\s+([+-]?[0-9.]+)", text, re.MULTILINE)
    if match is None:
        raise ValueError(f"missing CIF field {key}")
    return float(match.group(1))


def parse_cif(path: Path, element: str) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    lattice = {
        "a": cif_number(text, "_cell_length_a"),
        "b": cif_number(text, "_cell_length_b"),
        "c": cif_number(text, "_cell_length_c"),
        "alpha": cif_number(text, "_cell_angle_alpha"),
        "beta": cif_number(text, "_cell_angle_beta"),
        "gamma": cif_number(text, "_cell_angle_gamma"),
    }
    sites = []
    for match in re.finditer(
        rf"^\s*{element}\d+\s+([0-9.]+)\s+([+-]?[0-9.]+)\s+"
        rf"([+-]?[0-9.]+)\s+([+-]?[0-9.]+)",
        text,
        re.MULTILINE,
    ):
        occupancy, x, y, z = map(float, match.groups())
        sites.append(
            {
                "element": element,
                "wyckoff": "1a",
                "x": x,
                "y": y,
                "z": z,
                "occupancy": occupancy,
            }
        )
    if not sites:
        raise ValueError(f"no sites found in {path}")
    return {
        "file": f"CIFs.tar.gz/{element}.cif",
        "sha256": sha256(path),
        "lattice": lattice,
        "sites": sites,
    }


def cell_volume(lattice: dict[str, float]) -> float:
    alpha, beta, gamma = (
        math.radians(lattice[name]) for name in ("alpha", "beta", "gamma")
    )
    factor = math.sqrt(
        1
        + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
        - math.cos(alpha) ** 2
        - math.cos(beta) ** 2
        - math.cos(gamma) ** 2
    )
    return lattice["a"] * lattice["b"] * lattice["c"] * factor


def bm3_pressure_ratio(ratio: float, k0: float, k0_prime: float) -> float:
    eta = ratio ** (-1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def transcribe_source(source_root: Path) -> dict[str, Any]:
    wien_path = source_root / "WIEN2k.txt"
    fleur_path = source_root / "history" / "history" / "FLEUR.txt"
    readme_path = source_root / "README.txt"
    eosfit_path = source_root / "eosfit.py"
    wien = parse_parameter_file(wien_path)
    fleur = parse_parameter_file(fleur_path)
    structures = {
        element: parse_cif(source_root / "cifs" / f"{element}.cif", element)
        for element in METALS
    }
    records = [
        {"code": "WIEN2k", "element": element, **wien[element]} for element in METALS
    ]
    records.extend(
        {"code": "FLEUR", "element": element, **fleur[element]}
        for element in FLEUR_DUPLICATES
    )
    if len(records) != 50 or len(structures) != 48:
        raise AssertionError("the curated Delta-project selection must remain 50/48")
    return {
        "format": "peritheos.deltacodesdft-metal-eos-source",
        "format_version": 1,
        "archive": {
            "title": "Delta project web site archive",
            "url": ARCHIVE_URL,
            "doi": ARCHIVE_DOI,
            "related_publication_doi": CAMPAIGN_DOI,
            "license": "CC BY 4.0",
            "files": {
                "README.txt": sha256(readme_path),
                "WIEN2k.txt": sha256(wien_path),
                "history/history/FLEUR.txt": sha256(fleur_path),
                "eosfit.py": sha256(eosfit_path),
            },
        },
        "selection": {
            "description": "48 WIEN2k-PBE metal EOS plus FLEUR-PBE cross-code parameterizations for Al and Fe",
            "record_count": 50,
            "element_count": 48,
            "code_methods": CODE_METHOD,
        },
        "records": records,
        "structures": structures,
    }


def write_source_data(source: dict[str, Any]) -> None:
    DATASETS.mkdir(parents=True, exist_ok=True)
    SOURCE_JSON.write_text(
        json.dumps(source, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    fields = (
        "code",
        "element",
        "material_identifier",
        "source_v0_a3_atom",
        "formula_units_per_cell",
        "executable_v0_a3_cell",
        "k0_gpa",
        "k0_prime",
        "source_structure_file",
        "reference_doi",
    )
    with PARAMETER_CSV.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in source["records"]:
            element = row["element"]
            material_id = (
                EXISTING_PHASES[element]
                if element in EXISTING_PHASES
                else NEW_PHASES[element][0]
            )
            structure = source["structures"][element]
            z = len(structure["sites"])
            if material_id in EXISTING_PHASES.values():
                document = json.loads((MATERIALS / f"{material_id}.eosmat").read_text())
                z = document["formula_units_per_cell"]
            writer.writerow(
                {
                    "code": row["code"],
                    "element": element,
                    "material_identifier": material_id,
                    "source_v0_a3_atom": row["V0_a3_atom"],
                    "formula_units_per_cell": z,
                    "executable_v0_a3_cell": round(row["V0_a3_atom"] * z, 8),
                    "k0_gpa": row["K0_gpa"],
                    "k0_prime": row["K0_prime"],
                    "source_structure_file": structure["file"],
                    "reference_doi": CAMPAIGN_DOI,
                }
            )


def reference() -> dict[str, Any]:
    return {
        "authors": ["Lejaeghere", "et al."],
        "year": 2016,
        "title": "Reproducibility in density functional theory calculations of solids",
        "source": "Science",
        "volume": "351",
        "locator": "aad3000",
        "doi": CAMPAIGN_DOI,
    }


def record(
    source: dict[str, Any], row: dict[str, Any], material_id: str, z: float
) -> dict[str, Any]:
    code = row["code"]
    element = row["element"]
    year = 2016
    source_v0 = float(row["V0_a3_atom"])
    structure = source["structures"][element]
    structure_v_atom = cell_volume(structure["lattice"]) / len(structure["sites"])
    ratio_center = structure_v_atom / source_v0
    volume_ratio = [0.94 * ratio_center, 1.06 * ratio_center]
    pressure_range = sorted(
        bm3_pressure_ratio(ratio, row["K0_gpa"], row["K0_prime"])
        for ratio in volume_ratio
    )
    suffix = code.lower()
    magnetic_state = MAGNETIC_STATE.get(element)
    notes = (
        f"Static 0 K {CODE_METHOD[code]} BM3 fit to seven frozen-geometry E(V) "
        "calculations. "
        "The official workflow scales the archived equilibrium CIF to 94-106% "
        "of its volume; spin-orbit coupling is omitted."
    )
    if magnetic_state is not None:
        notes += f" The source requires {magnetic_state} for {element}."
    if element == "Mn":
        notes += (
            " The source deliberately uses a simpler antiferromagnetic fcc Mn "
            "model instead of complex alpha-Mn; K0-prime=-0.21 is source-reported "
            "and must not be extrapolated beyond the fitted volume window."
        )
    locations = [
        f"official archive {code if code == 'WIEN2k' else 'history/history/FLEUR'}.txt",
        "official archive README.txt, eosfit.py, and CIFs.tar.gz",
    ]
    return {
        "identifier": f"{material_id}_lejaeghere_{year}_{suffix}_pbe_bm3",
        "label": f"Lejaeghere et al. ({year}), {code} PBE BM3",
        "record_kind": "published",
        "equation_kind": "isothermal",
        "reference": reference(),
        "sample_composition": (
            f"Elemental {element}; archived Delta-project reference structure"
            + (f" with {magnetic_state}." if magnetic_state is not None else ".")
        ),
        "volume_basis": {
            "kind": "formula_units",
            "formula_units": z,
            "molar_mass_g_mol": MOLAR_MASS[element],
        },
        "eos": {
            "type": "BM3",
            "model": "birch_murnaghan_3",
            "parameters": {
                "V0": round(source_v0 * z, 8),
                "K0": float(row["K0_gpa"]),
                "K0_prime": float(row["K0_prime"]),
            },
        },
        "parameter_errors": {"V0": None, "K0": None, "K0_prime": None},
        "parameter_error_confidence": None,
        "parameter_covariance": None,
        "fixed_parameters": [],
        "temperature_ref": 0.0,
        "pressure_range_status": "theoretical",
        "validity": {
            "pressure_gpa": pressure_range,
            "temperature_k": [0.0, 0.0],
            "volume_ratio": volume_ratio,
            "notes": [
                "Static PBE fitting envelope, not an experimental phase-stability field.",
                "Pressure bounds are derived by evaluating BM3 at the archived volume-envelope endpoints.",
                "Use outside the archived seven-volume fitting envelope is extrapolation.",
            ],
        },
        "parameter_provenance": {
            "equation": "Official eosfit.py version 3.1, standard third-order Birch-Murnaghan E(V) fit.",
            "V0": f"{code} parameter table: {source_v0} A^3/atom; multiplied by Z={z:g} for the executable cell basis.",
            "K0": f"{code} parameter table, GPa.",
            "K0_prime": f"{code} parameter table, dimensionless BP column.",
            "uncertainties": "The official archive publishes no coefficient uncertainties or covariance.",
            "method": (
                f"Static {CODE_METHOD[code]} total energies at seven fixed cell "
                "volumes; no spin-orbit coupling"
                + (f"; {magnetic_state}." if magnetic_state is not None else ".")
            ),
        },
        "pressure_calibration": {
            "status": "not_applicable",
            "methods": [
                {
                    "kind": "ab_initio",
                    "source_location": "Official README.txt and eosfit.py",
                    "scope": "static total-energy calculations",
                }
            ],
            "recalculation": {
                "status": "not_applicable",
                "notes": "No experimental pressure scale applies.",
            },
            "audit_date": "2026-09-06",
        },
        "scientific_validation": {
            "status": "primary_source_validated",
            "note": "Equation, coefficients, units, structures, and seven-volume workflow were checked against the primary paper and official archived supplement.",
            "audit_date": "2026-09-06",
            "verified_fields": [
                "canonical_doi",
                "composition",
                "phase",
                "equation",
                "parameters",
                "units",
                "reference_state",
                "theoretical_method",
                "volume_range",
                "numerical_reproduction",
            ],
            "primary_source_check": {
                "access_url": ARCHIVE_URL,
                "doi": CAMPAIGN_DOI,
                "locations": locations,
                "finding": "The archived table explicitly reports per-atom V0, B0 in GPa, and dimensionless B1/BP from the documented seven-point BM3 workflow.",
            },
            "primary_data_check": {
                "status": "parameterization_only",
                "audit_date": "2026-09-06",
                "source_locations": [
                    f"official archive {code if code == 'WIEN2k' else 'history/history/FLEUR'}.txt",
                    f"official archive {structure['file']}",
                ],
                "finding": "The official archive preserves the complete EOS coefficients and frozen structure, but not the seven row-level E(V) values for this code; direct coefficient refitting is unavailable.",
            },
        },
        "source_lineage": [
            {
                "role": "scientific equation, workflow, and reference implementation",
                "citation": "Lejaeghere et al. (2014), Section IV.A and official eosfit.py",
                "doi": METHOD_DOI,
                "url": AUTHOR_PDF,
            },
            {
                "role": "multi-code reproducibility campaign",
                "citation": "Lejaeghere et al. (2016)",
                "doi": CAMPAIGN_DOI,
            },
            {
                "role": f"{code} coefficients and frozen elemental structure",
                "citation": f"Delta project web site archive, {ARCHIVE_DOI}",
                "url": ARCHIVE_URL,
                "license": "CC BY 4.0",
                "files": [
                    f"{code if code == 'WIEN2k' else 'history/history/FLEUR'}.txt",
                    structure["file"],
                ],
            },
        ],
        "notes": notes,
    }


def new_document(source: dict[str, Any], element: str) -> dict[str, Any]:
    identifier, name, phase = NEW_PHASES[element]
    structure = source["structures"][element]
    z = len(structure["sites"])
    return {
        "format": "peritheos.material",
        "format_version": 3,
        "identifier": identifier,
        "name": name,
        "formula": element,
        "phase": phase,
        "cell_contents": f"{z} {element} atom{'s' if z != 1 else ''} in the archived P1-expanded source cell",
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        },
        "symmetry": "TRICLINIC",
        "lattice": structure["lattice"],
        "formula_units_per_cell": z,
        "space_group": "P1",
        "space_group_number": 1,
        "atom_sites": structure["sites"],
        "peaks": [],
        "source": {
            "structure": structure["file"],
            "structure_sha256": structure["sha256"],
            "archive_url": ARCHIVE_URL,
            "archive_doi": ARCHIVE_DOI,
            "license": "CC BY 4.0",
        },
        "notes": (
            "The exact explicit P1 cell from the official Delta-project CIF is retained "
            "for diffraction and volume-basis reproducibility; the phase field records "
            "the crystallographic prototype. The structural cell and fitted V0 may "
            "differ because the archive combines code-specific EOS fits with a common "
            "reference geometry."
        ),
        "eos_records": [],
    }


def apply_records(source: dict[str, Any]) -> None:
    documents: dict[str, dict[str, Any]] = {}
    for row in source["records"]:
        element = row["element"]
        material_id = (
            EXISTING_PHASES[element]
            if element in EXISTING_PHASES
            else NEW_PHASES[element][0]
        )
        if material_id not in documents:
            path = MATERIALS / f"{material_id}.eosmat"
            documents[material_id] = (
                json.loads(path.read_text(encoding="utf-8"))
                if path.exists()
                else new_document(source, element)
            )
        document = documents[material_id]
        z = float(document["formula_units_per_cell"])
        new_record = record(source, row, material_id, z)
        records = [
            item
            for item in document["eos_records"]
            if item["identifier"]
            not in {
                new_record["identifier"],
                f"{material_id}_lejaeghere_2014_{row['code'].lower()}_pbe_bm3",
            }
        ]
        records.append(new_record)
        document["eos_records"] = records
    for material_id, document in documents.items():
        (MATERIALS / f"{material_id}.eosmat").write_text(
            json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    print(f"wrote 50 records across {len(documents)} material documents")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        type=Path,
        help="extracted official archive directory containing README.txt and WIEN2k.txt",
    )
    parser.add_argument(
        "--apply", action="store_true", help="update bundled material documents"
    )
    args = parser.parse_args()
    if args.source_root is not None:
        source = transcribe_source(args.source_root)
        write_source_data(source)
    else:
        source = json.loads(SOURCE_JSON.read_text(encoding="utf-8"))
    if args.apply:
        apply_records(source)


if __name__ == "__main__":
    main()
