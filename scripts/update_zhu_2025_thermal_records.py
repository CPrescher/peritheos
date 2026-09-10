#!/usr/bin/env python3
"""Split the Ye (2017) 300 K records from the Zhu et al. P-V-T records."""

from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"

YE_REFERENCE = {
    "authors": ["Ye", "Prakapenka", "Meng", "Shim"],
    "year": 2017,
    "title": "Intercomparison of the gold, platinum, and MgO pressure scales up to 140 GPa and 2500 K",
    "source": "Journal of Geophysical Research: Solid Earth",
    "volume": "122",
    "locator": "3450-3464",
    "doi": "10.1002/2016JB013811",
}

ZHU_REFERENCE = {
    "authors": ["Zhu", "Ye", "Sun", "Katsura", "Caracas"],
    "year": 2025,
    "title": "Internally consistent pressure-volume-temperature equations of state among platinum, gold, and MgO: Implications for phase transition boundaries in the mantle",
    "source": "ESSOAr preprint",
    "locator": "version 1",
    "doi": "10.22541/essoar.176236186.65259830/v1",
}

RELEASE = {
    "doi": "10.17632/6kxnhc2g73.3",
    "url": "https://data.mendeley.com/datasets/6kxnhc2g73/3",
    "archive_sha256": "77fbd3e3d02c2e5d27bf66b1ea69a88aed46266e1207b32b461fed33f48e4a0a",
}

CONFIG = {
    "gold": {
        "old_id": "gold_zhu_2025_vinet_300k",
        "ye_id": "gold_ye_2017_vinet_300k",
        "zhu_id": "gold_zhu_2025_pvt",
        "symbol": "Au",
        "static": {"V0": 67.85, "K0": 167.0, "K0_prime": 5.9},
        "thermal": {
            "Tr": 300.0,
            "theta0": 180.0,
            "gamma0": 2.93,
            "a": 0.75,
            "b": 2.7,
            "n": 1.0,
            "beta0": 0.0,
            "m": 1.0,
        },
        "thermal_errors": {"gamma0": 0.18, "a": 0.02, "b": 0.7},
        "fit_datasets": [
            "zhu_2025_au_shock",
            "zhu_2025_au_zero_pressure_thermal_expansion",
        ],
        "pressure_range": [0.0, 197.039],
        "temperature_range": [293.0, 1324.0],
        "volume_ratio": [46.7931 / 67.85, 61.1706 / 67.85],
        "source_files": "Au/Au_shock.dat and Au/Au_high_T.dat",
        "optimizer_note": "The record follows the v3 optimizer's theta0=180 K and its downstream rounded gamma coefficients. The separate v3 pressure calculator instead hard-codes theta0=170 K.",
    },
    "platinum": {
        "old_id": "platinum_zhu_2025_vinet_300k",
        "ye_id": "platinum_ye_2017_vinet_300k",
        "zhu_id": "platinum_zhu_2025_pvt",
        "symbol": "Pt",
        "static": {"V0": 60.38, "K0": 277.3, "K0_prime": 5.23},
        "thermal": {
            "Tr": 300.0,
            "theta0": 240.0,
            "gamma0": 2.75,
            "a": 0.39,
            "b": 5.1,
            "n": 1.0,
            "beta0": 0.002145,
            "m": 0.65,
        },
        "thermal_errors": {},
        "fit_datasets": [
            "zhu_2025_pt_shock",
            "zhu_2025_pt_zero_pressure_thermal_expansion",
        ],
        "pressure_range": [0.0, 682.402],
        "temperature_range": [293.0, 1900.0],
        "volume_ratio": [35.9002 / 60.38, 56.0169 / 60.38],
        "source_files": "Pt/Pt_shock_700.dat and Pt/Pt_high_T.dat",
        "optimizer_note": "The record follows the v3 optimizer's theta0=240 K and its downstream rounded gamma coefficients. The separate v3 pressure calculator instead hard-codes theta0=230 K.",
    },
    "mgo": {
        "old_id": "mgo_zhu_2025_vinet_300k",
        "ye_id": "mgo_ye_2017_vinet_300k",
        "zhu_id": "mgo_zhu_2025_pvt",
        "symbol": "MgO",
        "static": {"V0": 74.71, "K0": 160.3, "K0_prime": 4.182},
        "thermal": {
            "Tr": 300.0,
            "theta0": 761.0,
            "gamma0": 1.53,
            "a": 1.0,
            "b": 1.43,
            "n": 2.0,
            "beta0": -0.0008061,
            "m": 4.8,
        },
        "thermal_errors": {},
        "fit_datasets": ["zhu_2025_mgo_pvt"],
        "pressure_range": [0.0001, 138.104],
        "temperature_range": [298.0, 2986.0],
        "volume_ratio": [51.48 / 74.71, 85.8479 / 74.71],
        "source_files": "MgO/MgO_pvt_40_new.dat",
        "optimizer_note": "The record follows the v3 optimizer and property script (gamma0=1.53 after rounding). The separate v3 pressure calculator instead hard-codes gamma0=1.52.",
    },
}

DATASETS = {
    "gold": [
        {
            "identifier": "zhu_2025_au_shock",
            "kind": "shock_wave",
            "description": "Twelve Au shock states used by the Zhu et al. v3 thermal fit.",
            "source_location": "Mendeley v3: Au/Au_shock.dat",
            "columns": [
                (
                    "volume_a3_conventional_cell",
                    "conventional_unit_cell_volume",
                    "angstrom^3/conventional_unit_cell",
                    "value",
                ),
                ("shock_pressure_gpa", "pressure", "GPa", "value"),
                ("particle_velocity_m_s", "particle_velocity", "m/s", "value"),
            ],
            "path": "datasets/zhu-2025-au-shock.csv",
            "sha256": "ef2fb3cbfc1636d6a513c40783e6c860d5d3ced3a9e48e0e939ef36f8f8ece01",
            "official_sha256": "fab474e9615ba8708d72d187672d93ad55ae6c15a79bc0c31e8f554beed4c2f5",
        },
        {
            "identifier": "zhu_2025_au_zero_pressure_thermal_expansion",
            "kind": "temperature_volume",
            "description": "Ten zero-pressure Au temperature-volume states used by the Zhu et al. v3 thermal fit.",
            "source_location": "Mendeley v3: Au/Au_high_T.dat",
            "columns": [
                ("temperature_k", "temperature", "K", "value"),
                (
                    "volume_a3_conventional_cell",
                    "conventional_unit_cell_volume",
                    "angstrom^3/conventional_unit_cell",
                    "value",
                ),
            ],
            "path": "datasets/zhu-2025-au-zero-pressure-thermal-expansion.csv",
            "sha256": "3baac523dbf58873d01e05499f73934bbebdfc5bbd03de876c9e4f06dfeee611",
            "official_sha256": "a803b6f3805dec09e08b2312af466c80250b91239288924d061d029e35e56d5d",
        },
    ],
    "platinum": [
        {
            "identifier": "zhu_2025_pt_shock",
            "kind": "shock_wave",
            "description": "Sixty-eight Pt shock states used by the Zhu et al. v3 thermal fit.",
            "source_location": "Mendeley v3: Pt/Pt_shock_700.dat",
            "columns": [
                (
                    "volume_a3_conventional_cell",
                    "conventional_unit_cell_volume",
                    "angstrom^3/conventional_unit_cell",
                    "value",
                ),
                ("shock_pressure_gpa", "pressure", "GPa", "value"),
                ("particle_velocity_m_s", "particle_velocity", "m/s", "value"),
            ],
            "path": "datasets/zhu-2025-pt-shock.csv",
            "sha256": "025914a80db03e52fe631c8f18646b722a558026d616079edea3d47c7509b924",
            "official_sha256": "5591f904a2c0e7200f5e9a6f0a3151607247b187d8e5473ee900990ea6526238",
        },
        {
            "identifier": "zhu_2025_pt_zero_pressure_thermal_expansion",
            "kind": "temperature_volume",
            "description": "Seventeen zero-pressure Pt temperature-volume states used by the Zhu et al. v3 thermal fit.",
            "source_location": "Mendeley v3: Pt/Pt_high_T.dat",
            "columns": [
                ("temperature_k", "temperature", "K", "value"),
                (
                    "volume_a3_conventional_cell",
                    "conventional_unit_cell_volume",
                    "angstrom^3/conventional_unit_cell",
                    "value",
                ),
            ],
            "path": "datasets/zhu-2025-pt-zero-pressure-thermal-expansion.csv",
            "sha256": "c338f3db484be4927880b90f7fcf88ce6b29c3c880e726755d3b29cf778cda70",
            "official_sha256": "986ce36b663c065f639fce79ffcd6f3f896d76f2740bad10a0e1d5eeb1647aa8",
        },
    ],
    "mgo": [
        {
            "identifier": "zhu_2025_mgo_pvt",
            "kind": "pressure_volume_temperature",
            "description": "All 213 MgO P-V-T states used by the Zhu et al. v3 thermal fit.",
            "source_location": "Mendeley v3: MgO/MgO_pvt_40_new.dat",
            "columns": [
                ("pressure_gpa", "pressure", "GPa", "value"),
                ("temperature_k", "temperature", "K", "value"),
                (
                    "volume_a3_conventional_cell",
                    "conventional_unit_cell_volume",
                    "angstrom^3/conventional_unit_cell",
                    "value",
                ),
            ],
            "path": "datasets/zhu-2025-mgo-pvt.csv",
            "sha256": "839acb4435c741aaa7266c956d9eba196756e5782ca63d47fb0b1a930fe4285c",
            "official_sha256": "7fce0f2c195d0d3a2bb7f59e22064c283a6da64c259823466bb8616207ae34a0",
        }
    ],
}


def _ye_record(source: dict, config: dict) -> dict:
    record = copy.deepcopy(source)
    record["identifier"] = config["ye_id"]
    record["aliases"] = [config["old_id"]]
    record["label"] = (
        f"Ye et al. (2017), {config['symbol']} Vinet 300 K co-compression fit"
    )
    record["reference"] = YE_REFERENCE
    record["source_lineage"] = [
        line
        for line in record["source_lineage"]
        if line.get("doi") not in {ZHU_REFERENCE["doi"], RELEASE["doi"]}
    ]
    record["validity"]["notes"] = [
        "This is Ye et al.'s fitted 300 K co-compression isotherm.",
        "Zhu et al. later fixed this static branch inside a separate thermal P-V-T model.",
    ]
    record["parameter_provenance"]["equation"] = (
        "Ye et al. (2017), Section 3.3 and Table 1: third-order Vinet EOS."
    )
    record["parameter_provenance"]["K0_prime"] = record["parameter_provenance"][
        "K0_prime"
    ].replace("; Zhu (2025) adopts rather than refits it", "")
    record["parameter_provenance"]["scope"] = (
        "Complete Ye et al. 300 K Vinet fit; no Zhu thermal parameters are part of this record."
    )
    record["scientific_validation"]["note"] = (
        "The equation, parameters, uncertainty, reference state, and corrected co-compression inputs were checked against Ye et al. (2017)."
    )
    record["scientific_validation"]["primary_source_check"].pop(
        "adopting_source_doi", None
    )
    record["scientific_validation"]["primary_source_check"]["finding"] = (
        record["scientific_validation"]["primary_source_check"]["finding"].split(
            "; Zhu"
        )[0]
        + "."
    )
    return record


def _zhu_record(source: dict, config: dict) -> dict:
    thermal = config["thermal"]
    record = copy.deepcopy(source)
    record.pop("aliases", None)
    record["identifier"] = config["zhu_id"]
    record["label"] = (
        f"Zhu et al. (2025), {config['symbol']} internally consistent P-V-T EOS"
    )
    record["reference"] = ZHU_REFERENCE
    record["source_lineage"] = [
        {
            "role": "fixed 300 K Vinet branch",
            "citation": "Ye et al. (2017), Section 3.3 and Table 1",
            "doi": YE_REFERENCE["doi"],
        },
        {
            "role": "thermal model equations and version-1 parameter table",
            "citation": "Zhu et al. (2025), equations (2)-(8) and Table 1",
            "doi": ZHU_REFERENCE["doi"],
        },
        {
            "role": "version-3 fit inputs, revised optimizer, and downstream property coefficients",
            "citation": "Zhu et al. Mendeley Data version 3",
            "doi": RELEASE["doi"],
        },
    ]
    record["equation_kind"] = "thermal"
    record["eos"]["parameters"] = config["static"]
    record["parameter_errors"] = {name: None for name in config["static"]}
    record["fixed_parameters"] = list(config["static"])
    record["fit_datasets"] = config["fit_datasets"]
    record["thermal"] = {
        "type": "AsymptoticPowerLawMieGruneisenDebyeExcess",
        "model": "asymptotic_power_law_mie_gruneisen_debye_excess",
        "parameters": thermal,
        "parameter_errors": {
            name: config["thermal_errors"].get(name) for name in thermal
        },
        "fixed_parameters": ["Tr", "theta0", "a", "n", "beta0", "m"],
    }
    record["experimental_pressure_range_gpa"] = config["pressure_range"]
    record["pressure_range_status"] = "reference_parameterization"
    record["experimental_temperature_range_k"] = config["temperature_range"]
    record["validity"] = {
        "pressure_gpa": config["pressure_range"],
        "temperature_k": config["temperature_range"],
        "volume_ratio": config["volume_ratio"],
        "notes": [
            "Bounds conservatively report the extrema of the released v3 fit-input tables, not a rectangular phase-stability field.",
            config["optimizer_note"],
        ],
    }
    record["parameter_provenance"] = {
        "reference_isotherm": {
            name: "Fixed Ye et al. (2017) 300 K Vinet coefficient used by the released v3 optimizer."
            for name in config["static"]
        },
        "thermal_correction": {
            name: (
                f"Mendeley Data v3 optimizer and downstream property implementation for {config['symbol']}; "
                "beta0 is stored in J mol^-1 K^-2 after conversion from the "
                "source's J g^-1 K^-2 coefficient."
            )
            for name in thermal
        },
        "additional": {
            "thermal_model": "Zhu et al. (2025), equations (2)-(8): asymptotic-power-law Gruneisen parameter, Debye thermal energy, and a volume-dependent T-squared excess Helmholtz term, all referenced to 300 K.",
            "version_note": config["optimizer_note"],
        },
    }
    record["pressure_calibration"] = {
        "status": "not_applicable",
        "methods": [
            {
                "kind": "self_consistent",
                "reference": ZHU_REFERENCE,
                "source_location": "Equations (2)-(8), released v3 fit inputs, and optimizer scripts",
                "scope": "Joint static, shock, and thermal constraints on the internally consistent pressure scale.",
            }
        ],
        "recalculation": {
            "status": "ready",
            "notes": "scripts/reproduce_zhu_2025_pressure_standards.py independently translates the released v3 iterative robust-bisquare optimizer and refits gamma0 and b from every bundled observation.",
        },
        "audit_date": "2026-09-08",
    }
    record["scientific_validation"] = {
        "status": "primary_source_validated",
        "note": "The thermal equations were checked against the preprint and the executable coefficients and fit inputs against the CC BY 4.0 Mendeley v3 release.",
        "audit_date": "2026-09-08",
        "verified_fields": [
            "equation",
            "parameters",
            "units",
            "reference_state",
            "phase",
            "validity",
            "primary_data",
        ],
        "primary_source_check": {
            "access_url": "https://d197for5662m48.cloudfront.net/documents/publicationstatus/288995/preprint_pdf/85b400a506980803732aafbd5816f7fb.pdf",
            "doi": ZHU_REFERENCE["doi"],
            "locations": [
                "Equations (2)-(8)",
                "Table 1",
                f"Mendeley v3 {config['symbol']} optimizer and property scripts",
            ],
            "finding": "The record implements the complete 300 K-referenced Debye plus T-squared-excess pressure form and the v3 optimizer/downstream-property coefficient set.",
        },
        "primary_data_check": {
            "status": "bundled",
            "audit_date": "2026-09-08",
            "dataset_identifiers": config["fit_datasets"],
            "finding": f"Every active numeric row from the released v3 {config['source_files']} inputs is bundled as a stable CSV normalization under CC BY 4.0.",
        },
        "reproduction": {
            "method": "Independent translation of the released Mendeley v3 iterative robust-bisquare optimizer using all bundled shock, thermal-expansion, or P-V-T rows.",
            "note": "This establishes thermal-refit parity against the optimizer/property parameter set. The inconsistent standalone v3 pressure calculators remain documented as a separate implementation diagnostic.",
        },
    }
    record["notes"] = (
        "This Zhu record is the broader temperature-dependent P-V-T model. The separately identified Ye record remains the source of the fitted 300 K co-compression isotherm."
    )
    return record


def _dataset(spec: dict, used_by: str) -> dict:
    return {
        "identifier": spec["identifier"],
        "kind": spec["kind"],
        "description": spec["description"],
        "reference": ZHU_REFERENCE,
        "source_location": spec["source_location"],
        "source_url": RELEASE["url"],
        "license": "CC-BY-4.0",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "license_scope": "Zhu et al. Mendeley Data version 3 and this normalized redistribution.",
        "columns": [
            {"name": name, "quantity": quantity, "unit": unit, "role": role}
            for name, quantity, unit, role in spec["columns"]
        ],
        "resource": {
            "path": spec["path"],
            "sha256": spec["sha256"],
            "media_type": "text/csv",
        },
        "provenance": {
            "type": "lossless_csv_normalization",
            "release_doi": RELEASE["doi"],
            "release_archive_sha256": RELEASE["archive_sha256"],
            "official_filename": spec["source_location"].split(": ", 1)[1],
            "official_sha256": spec["official_sha256"],
            "normalization": "Removed blank and percent-comment lines, replaced whitespace separation with CSV commas, and supplied explicit quantity names; active numeric tokens and row order are unchanged.",
        },
        "used_by_eos_records": [used_by],
    }


def update() -> None:
    for material, config in CONFIG.items():
        path = MATERIALS / f"{material}.eosmat"
        document = json.loads(path.read_text(encoding="utf-8"))
        records = document["eos_records"]
        records[:] = [
            record for record in records if record["identifier"] != config["zhu_id"]
        ]
        matches = [
            index
            for index, record in enumerate(records)
            if record["identifier"] in {config["old_id"], config["ye_id"]}
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one source record for {material}, got {matches}"
            )
        index = matches[0]
        source = records[index]
        records[index : index + 1] = [
            _ye_record(source, config),
            _zhu_record(source, config),
        ]

        for dataset in document.get("datasets", []):
            dataset["used_by_eos_records"] = [
                config["ye_id"] if identifier == config["old_id"] else identifier
                for identifier in dataset.get("used_by_eos_records", [])
            ]
        new_ids = {spec["identifier"] for spec in DATASETS[material]}
        document["datasets"] = [
            dataset
            for dataset in document.get("datasets", [])
            if dataset["identifier"] not in new_ids
        ]
        document["datasets"].extend(
            _dataset(spec, config["zhu_id"]) for spec in DATASETS[material]
        )
        path.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    update()
