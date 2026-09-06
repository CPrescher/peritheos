#!/usr/bin/env python3
"""Reproduce the 50 archived Delta-project metal BM3 parameterizations."""

from __future__ import annotations

import csv
import json
import math
from importlib import resources
from typing import Any

import numpy as np

from peritheos import Material, get_material_document

TABLE = "deltacodesdft-metal-eos-parameters.csv"
SOURCE = "deltacodesdft-metal-eos-source.json"


def bm3_pressure(volume: float, v0: float, k0: float, k0_prime: float) -> float:
    """Independent standard third-order Birch-Murnaghan pressure expression."""
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def cell_volume(lattice: dict[str, float]) -> float:
    alpha, beta, gamma = (
        math.radians(lattice[name]) for name in ("alpha", "beta", "gamma")
    )
    metric = (
        1
        + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
        - math.cos(alpha) ** 2
        - math.cos(beta) ** 2
        - math.cos(gamma) ** 2
    )
    return lattice["a"] * lattice["b"] * lattice["c"] * math.sqrt(metric)


def load_inputs() -> tuple[list[dict[str, str]], dict[str, Any]]:
    root = resources.files("peritheos").joinpath("data", "datasets")
    with root.joinpath(TABLE).open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    source = json.loads(root.joinpath(SOURCE).read_text(encoding="utf-8"))
    return rows, source


def reproduce() -> dict[str, float | int]:
    rows, source = load_inputs()
    source_rows = {(row["code"], row["element"]): row for row in source["records"]}
    largest_parameter_error = 0.0
    largest_pressure_error = 0.0
    largest_p0 = 0.0
    for row in rows:
        code = row["code"]
        element = row["element"]
        source_row = source_rows[(code, element)]
        material_id = row["material_identifier"]
        identifier = f"{material_id}_lejaeghere_2016_{code.lower()}_pbe_bm3"
        document = get_material_document(material_id)
        stored = next(
            record
            for record in document["eos_records"]
            if record["identifier"] == identifier
        )
        parameters = stored["eos"]["parameters"]
        expected = {
            "V0": float(row["executable_v0_a3_cell"]),
            "K0": float(source_row["K0_gpa"]),
            "K0_prime": float(source_row["K0_prime"]),
        }
        largest_parameter_error = max(
            largest_parameter_error,
            *(abs(float(parameters[name]) - value) for name, value in expected.items()),
        )
        executable = (
            Material.from_eosmat(document, record_identifiers=[identifier])
            .eos_records[0]
            .eos
        )
        v0 = expected["V0"]
        largest_p0 = max(largest_p0, abs(float(executable.pressure(v0))))
        for ratio in (0.94, 1.06):
            expected_pressure = bm3_pressure(
                ratio * v0, v0, expected["K0"], expected["K0_prime"]
            )
            largest_pressure_error = max(
                largest_pressure_error,
                abs(float(executable.pressure(ratio * v0)) - expected_pressure),
            )

        structure = source["structures"][element]
        structure_v_atom = cell_volume(structure["lattice"]) / len(structure["sites"])
        center = structure_v_atom / float(source_row["V0_a3_atom"])
        np.testing.assert_allclose(
            stored["validity"]["volume_ratio"],
            [0.94 * center, 1.06 * center],
            rtol=0.0,
            atol=1.0e-14,
        )
    return {
        "records": len(rows),
        "elements": len({row["element"] for row in rows}),
        "wien2k_records": sum(row["code"] == "WIEN2k" for row in rows),
        "fleur_records": sum(row["code"] == "FLEUR" for row in rows),
        "largest_parameter_error": largest_parameter_error,
        "largest_p0_gpa": largest_p0,
        "largest_independent_pressure_error_gpa": largest_pressure_error,
    }


def main() -> None:
    for key, value in reproduce().items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
