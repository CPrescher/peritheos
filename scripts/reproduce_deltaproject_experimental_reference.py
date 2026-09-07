#!/usr/bin/env python3
"""Verify the Delta archived experimental-reference parameterizations.

This is deliberately a coefficient/provenance reproduction, not an
experimental refit: the Delta archive does not contain a unified P-V dataset
for these reference constructions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "peritheos" / "data" / "datasets"
MATERIALS = ROOT / "peritheos" / "data" / "materials"
SOURCE = DATASETS / "deltaproject-experimental-reference-source.json"
PARAMETERS = DATASETS / "deltaproject-experimental-reference-parameters.csv"
KNITTLE_RECONSTRUCTION = (
    DATASETS / "deltaproject-knittle-1995-b0prime-reconstruction.csv"
)
OSMIUM_DATA = DATASETS / "osmium-takemura-2004-table1-compression.csv"
ELEMENTS = {
    "Li",
    "Na",
    "Mg",
    "Al",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Cs",
    "Ba",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Tl",
    "Pb",
    "Bi",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_exp(path: Path) -> dict[str, tuple[float, float, float]]:
    parsed: dict[str, tuple[float, float, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 4 or fields[0].startswith("#"):
            continue
        parsed[fields[0]] = tuple(map(float, fields[1:]))
    return parsed


def _independent_bm3_pressure(
    volume: float, v0: float, k0: float, k0_prime: float
) -> float:
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (
        1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0)
    )


def _reproduce_knittle_b0_prime() -> tuple[int, list[str], list[str]]:
    reconstructed: list[str] = []
    unresolved: list[str] = []
    with KNITTLE_RECONSTRUCTION.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 32:
        raise AssertionError(f"expected 32 Knittle rows, found {len(rows)}")
    for row in rows:
        status = row["status"]
        if status == "unresolved":
            unresolved.append(row["element"])
            continue
        values = [float(value) for value in row["knittle_component_values"].split(";")]
        if row["reduction"] in {"direct", "direct_selected"}:
            result = values[0]
        elif row["reduction"] == "arithmetic_mean":
            result = sum(values) / len(values)
        elif row["reduction"] == "arithmetic_mean_rounded_3dp":
            result = round(sum(values) / len(values), 3)
        else:
            raise AssertionError(f"unknown reduction {row['reduction']!r}")
        if abs(result - float(row["delta_b0_prime"])) > 1.0e-12:
            raise AssertionError(f"Knittle reconstruction failed for {row['element']}")
        reconstructed.append(row["element"])
    return len(rows), reconstructed, unresolved


def _reproduce_osmium_fixed_v0_fit() -> tuple[float, float]:
    with OSMIUM_DATA.open(encoding="utf-8", newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if row["run"] != "2.0"]
    pressure = np.asarray([float(row["pressure_gpa"]) for row in rows])
    volume_ratio = np.asarray([float(row["volume_ratio"]) for row in rows])

    def residual(parameters: np.ndarray) -> np.ndarray:
        k0, k0_prime = parameters
        return np.asarray(
            [
                _independent_bm3_pressure(volume, 1.0, k0, k0_prime)
                for volume in volume_ratio
            ]
        ) - pressure

    fit = least_squares(residual, np.asarray([395.0, 4.5]))
    if not fit.success:
        raise AssertionError(f"osmium fixed-V0 fit failed: {fit.message}")
    return float(fit.x[0]), float(fit.x[1])


def reproduce(source_exp: Path | None = None) -> dict[str, Any]:
    metadata = json.loads(SOURCE.read_text(encoding="utf-8"))
    with PARAMETERS.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 42:
        raise AssertionError(f"expected 42 compiled-reference rows, found {len(rows)}")

    archived: dict[str, tuple[float, float, float]] | None = None
    if source_exp is not None:
        expected = metadata["archive"]["files"][
            "Delta_v3-1_0.zip/history.tar.gz/history/exp.txt"
        ]["sha256"]
        if _sha256(source_exp) != expected:
            raise AssertionError("exp.txt SHA-256 does not match the archived source")
        archived = _parse_exp(source_exp)
        if set(archived).intersection(ELEMENTS) != ELEMENTS:
            raise AssertionError("the archived exp.txt is missing a selected metal row")

    maximum_parameter_error = 0.0
    maximum_pressure_error = 0.0
    identifiers: set[str] = set()
    for row in rows:
        material = json.loads(
            MATERIALS.joinpath(f"{row['material_identifier']}.eosmat").read_text(
                encoding="utf-8"
            )
        )
        record = next(
            item
            for item in material["eos_records"]
            if item["identifier"] == row["record_identifier"]
        )
        identifiers.add(record["identifier"])
        if record["record_kind"] != "derived":
            raise AssertionError(f"{record['identifier']} is not classified as derived")
        if record["pressure_range_status"] != "reference_parameterization":
            raise AssertionError(
                f"{record['identifier']} claims a non-reference pressure range"
            )
        if "experimental_pressure_range_gpa" in record:
            raise AssertionError(
                f"{record['identifier']} incorrectly claims experimental coverage"
            )
        primary_data_status = record["scientific_validation"]["primary_data_check"][
            "status"
        ]
        if primary_data_status != "parameterization_only":
            raise AssertionError(f"{record['identifier']} incorrectly claims raw data")

        parameters = record["eos"]["parameters"]
        expected = (float(row["V0"]), float(row["K0"]), float(row["K0_prime"]))
        actual = (parameters["V0"], parameters["K0"], parameters["K0_prime"])
        maximum_parameter_error = max(
            maximum_parameter_error,
            *(abs(left - right) for left, right in zip(expected, actual)),
        )

        if archived is not None:
            element = record["label"].split(",")[1].strip().split()[0]
            v0_atom, k0, k0_prime = archived[element]
            z = float(material["formula_units_per_cell"])
            source_expected = (round(v0_atom * z, 10), k0, k0_prime)
            maximum_parameter_error = max(
                maximum_parameter_error,
                *(
                    abs(left - right)
                    for left, right in zip(source_expected, actual)
                ),
            )

        eos = BM3(*actual)
        for ratio in (0.94, 1.0, 1.06):
            volume = actual[0] * ratio
            maximum_pressure_error = max(
                maximum_pressure_error,
                abs(
                    eos.pressure(volume)
                    - _independent_bm3_pressure(volume, *actual)
                ),
            )

    conflicts = {item["element"]: item for item in metadata["known_conflicts"]}
    if set(conflicts) != {"Sr", "Tl"}:
        raise AssertionError("unexpected or unrecorded source conflict")
    if conflicts["Tl"]["archive_exp_txt"] != 3.0:
        raise AssertionError("unexpected Tl archive value")
    if conflicts["Sr"]["anderson_1990_table_1"] != 2.41:
        raise AssertionError("unexpected Sr original-source value")

    knittle_rows, knittle_reconstructed, knittle_unresolved = (
        _reproduce_knittle_b0_prime()
    )
    if len(knittle_reconstructed) != 30 or knittle_unresolved != ["Fe", "Au"]:
        raise AssertionError("unexpected Knittle recovery coverage")

    sn_k0_gpa = (0.690 + 2.0 * 0.293) / 3.0 * 100.0
    if round(sn_k0_gpa, 1) != 42.5:
        raise AssertionError("alpha-Sn elastic-constant derivation failed")

    os_k0, os_k0_prime = _reproduce_osmium_fixed_v0_fit()
    if abs(os_k0 - 395.0) > 1.0 or abs(os_k0_prime - 4.5) > 0.02:
        raise AssertionError("rounded Takemura Table I does not reproduce the Os fit")

    return {
        "records": len(rows),
        "unique_identifiers": len(identifiers),
        "archive_transcription_checked": archived is not None,
        "maximum_parameter_error": maximum_parameter_error,
        "maximum_independent_pressure_error_gpa": maximum_pressure_error,
        "raw_observation_refit": "not possible for a unified Delta triplet; upstream availability is reported separately",
        "knittle_b0_prime_rows": knittle_rows,
        "knittle_b0_prime_reconstructed": len(knittle_reconstructed),
        "knittle_b0_prime_unresolved": knittle_unresolved,
        "alpha_sn_k0_from_elastic_constants_gpa": sn_k0_gpa,
        "osmium_fixed_v0_refit": {"K0": os_k0, "K0_prime": os_k0_prime},
        "known_source_conflicts": len(conflicts),
        "known_phase_mismatches": len(metadata["known_phase_mismatches"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-exp",
        type=Path,
        help="optional extracted Materials Cloud Delta_v3-1_0 history/history/exp.txt",
    )
    args = parser.parse_args()
    print(json.dumps(reproduce(args.source_exp), indent=2))


if __name__ == "__main__":
    main()
