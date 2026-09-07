#!/usr/bin/env python3
"""Reproduce the Dorfman et al. (2012) simultaneous Vinet fits.

Supply the publisher file ``jgrb17272-sup-0002-txts01.pdf`` with
``--source-pdf``. The parser verifies its SHA-256 checksum before extracting
Tables S1-S6 with Poppler's ``pdftotext -tsv`` output. It can also regenerate
the bundled factual CSV transcription without copying the source PDF itself.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

import numpy as np
from scipy.optimize import least_squares

DOI = "10.1029/2012JB009292"
CORRECTION_DOI = "10.1029/2012JB009800"
ARTICLE_URL = (
    "https://duffy.princeton.edu/sites/g/files/toruqf616/files/dorfman_jgr_12.pdf"
)
ARTICLE_SHA256 = "1c9fe0f884fdb61f9d52d60308276d63dd434f27e59bcea2c74471566d0ec5e0"
SOURCE_FILENAME = "jgrb17272-sup-0002-txts01.pdf"
SOURCE_URL = (
    "https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?"
    "doi=10.1029%2F2012JB009292&file=jgrb17272-sup-0002-txts01.pdf"
)
SOURCE_SHA256 = "d7f7214cae4dabb3aae8a8585df4af7ada7990cb36b83a42ed7ed97e9eef7f2a"
DATASET_FILENAME = "dorfman-2012-tables-s1-s6-cocompression.csv"
DATASET_LICENSE = "CC0-1.0"
DATASET_LICENSE_URL = "https://creativecommons.org/publicdomain/zero/1.0/"

MATERIAL_ORDER = ("Au", "MgO", "Mo", "NaCl B2", "Ne", "Pt")
REFERENCE_PARAMETERS = {
    "Au": {"V0": 67.85, "K0": 167.0, "K0_prime": 5.88},
    "MgO": {"V0": 74.698, "K0": 160.6, "K0_prime": 4.37},
    "Mo": {"V0": 31.17, "K0": 261.0, "K0_prime": 4.19},
    "NaCl B2": {"V0": 41.35, "K0": 24.2, "K0_prime": 5.76},
    "Ne": {"V0": 88.967, "K0": 1.04, "K0_prime": 8.48},
    "Pt": {"V0": 60.38, "K0": 277.0, "K0_prime": 5.43},
}
PUBLISHED_FREE_PARAMETERS = {
    "Au": {"V0": 67.85, "K0": 167.0, "K0_prime": 5.84},
    "MgO": {"V0": 74.698, "K0": 160.6, "K0_prime": 4.37},
    "Mo": {"V0": 31.17, "K0": 271.0, "K0_prime": 3.89},
    "NaCl B2": {"V0": 41.35, "K0": 26.4, "K0_prime": 5.49},
    "Ne": {"V0": 88.967, "K0": 1.02, "K0_prime": 8.50},
    "Pt": {"V0": 60.38, "K0": 280.0, "K0_prime": 5.29},
}

# Table layout after ``pdftotext -tsv``.  The run boundaries are the printed
# vertical positions of the subheadings in the seven-page auxiliary PDF.
PAGE_LAYOUTS = {
    1: ((202.15, 303.62), ((480.47, "MMH6"), (math.inf, "MMH12"))),
    2: ((151.42, 252.88, 354.35), ((math.inf, "ANN11"),)),
    3: ((151.42, 252.88, 354.35), ((454.08, "MNN13"), (math.inf, "MNN11"))),
    4: ((151.42, 252.88, 354.35), ((358.44, "PMNN7"), (math.inf, "PMN3"))),
    5: ((202.15, 303.62), ((374.07, "PMH6"), (math.inf, "PNH7"))),
}
RUN_MATERIALS = {
    "MMH6": ("Mo", "MgO"),
    "MMH12": ("Mo", "MgO"),
    "ANN11": ("Au", "NaCl B2", "Ne"),
    "MNN13": ("MgO", "NaCl B2", "Ne"),
    "MNN11": ("MgO", "NaCl B2", "Ne"),
    "PMNN7": ("Pt", "MgO", "Ne"),
    "PMN3": ("Pt", "MgO", "Ne"),
    "PMH6": ("Pt", "MgO"),
    "PNH7": ("Pt", "NaCl B2"),
}
RUN_TABLES = {
    "MMH6": "S1",
    "MMH12": "S1",
    "ANN11": "S2",
    "MNN13": "S3",
    "MNN11": "S3",
    "PMNN7": "S4",
    "PMN3": "S4",
    "PMH6": "S5",
    "PNH7": "S6",
}
VOLUME_TOKEN = re.compile(r"^(?P<value>\d+\.\d+)(?:\((?P<error>\d+)\)|(?P<star>\*))$")
PRESSURE_TOKEN = re.compile(r"^\d+\.\d+$")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vinet_pressure(volume: float, parameters: dict[str, float]) -> float:
    x = (volume / parameters["V0"]) ** (1.0 / 3.0)
    eta = 1.5 * (parameters["K0_prime"] - 1.0)
    return 3.0 * parameters["K0"] * (1.0 - x) * math.exp(eta * (1.0 - x)) / x**2


def _run_for_position(page: int, top: float) -> str:
    for boundary, run in PAGE_LAYOUTS[page][1]:
        if top < boundary:
            return run
    raise AssertionError("unreachable")


def _parse_volume(token: str) -> tuple[float, Optional[float], bool]:
    match = VOLUME_TOKEN.fullmatch(token)
    if match is None:
        raise ValueError(f"not a source-table volume: {token!r}")
    value_text = match.group("value")
    error_text = match.group("error")
    error = None
    if error_text is not None:
        error = int(error_text) * 10.0 ** (-len(value_text.partition(".")[2]))
    return float(value_text), error, match.group("star") is not None


def extract_observations(source_pdf: Path) -> list[dict[str, Any]]:
    actual_sha256 = sha256(source_pdf)
    if actual_sha256 != SOURCE_SHA256:
        raise ValueError(
            f"unexpected supporting-PDF SHA-256 {actual_sha256}; expected {SOURCE_SHA256}"
        )
    completed = subprocess.run(
        ["pdftotext", "-f", "1", "-l", "5", "-tsv", str(source_pdf), "-"],
        check=True,
        capture_output=True,
        text=True,
    )
    words = list(csv.DictReader(io.StringIO(completed.stdout), delimiter="\t"))
    by_line: dict[tuple[int, float], list[dict[str, str]]] = defaultdict(list)
    for word in words:
        if word["level"] != "5":
            continue
        page = int(word["page_num"])
        if page not in PAGE_LAYOUTS:
            continue
        by_line[(page, round(float(word["top"]), 2))].append(word)

    observations: list[dict[str, Any]] = []
    for (page, top), line_words in sorted(by_line.items()):
        volume_words = [
            word for word in line_words if VOLUME_TOKEN.fullmatch(word["text"])
        ]
        if not volume_words:
            continue
        run = _run_for_position(page, top)
        materials = RUN_MATERIALS[run]
        volume_columns = PAGE_LAYOUTS[page][0]
        values: dict[str, dict[str, Any]] = {}
        for word in volume_words:
            left = float(word["left"])
            column = min(
                range(len(volume_columns)), key=lambda i: abs(left - volume_columns[i])
            )
            if abs(left - volume_columns[column]) > 1.0:
                raise ValueError(
                    f"unrecognized table column at page {page}, y={top}: {left}"
                )
            material = materials[column]
            volume, volume_sigma, single_peak = _parse_volume(word["text"])
            next_column = (
                volume_columns[column + 1]
                if column + 1 < len(volume_columns)
                else math.inf
            )
            pressure_candidates = [
                item
                for item in line_words
                if PRESSURE_TOKEN.fullmatch(item["text"])
                and left + 35.0 < float(item["left"]) < next_column
            ]
            if len(pressure_candidates) != 1:
                raise ValueError(
                    f"expected one printed pressure at page {page}, y={top}, {material}; "
                    f"found {pressure_candidates}"
                )
            values[material] = {
                "source_volume_token": word["text"],
                "volume_a3_conventional_cell": volume,
                "volume_sigma_a3": volume_sigma,
                "single_peak_mgo": single_peak,
                "source_pressure_token": pressure_candidates[0]["text"],
                "printed_pressure_gpa": float(pressure_candidates[0]["text"]),
            }
        if len(values) < 2:
            raise ValueError(
                f"source row at page {page}, y={top} has fewer than two volumes"
            )
        observations.append({"run": run, "source_page": page, "values": values})
    return observations


def write_dataset(observations: list[dict[str, Any]], output: Path) -> None:
    """Write a lossless long-form transcription of the source table cells."""
    fields = [
        "observation_index",
        "source_table",
        "source_pdf_page",
        "run",
        "row_in_run",
        "material",
        "volume_source_token",
        "volume_a3_conventional_cell",
        "volume_standard_error_a3",
        "reported_pressure_source_token",
        "reported_pressure_gpa",
        "single_peak_mgo",
    ]
    run_rows: dict[str, int] = defaultdict(int)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for observation_index, observation in enumerate(observations, 1):
            run = observation["run"]
            run_rows[run] += 1
            for material, value in observation["values"].items():
                writer.writerow(
                    {
                        "observation_index": observation_index,
                        "source_table": RUN_TABLES[run],
                        "source_pdf_page": observation["source_page"],
                        "run": run,
                        "row_in_run": run_rows[run],
                        "material": material,
                        "volume_source_token": value["source_volume_token"],
                        "volume_a3_conventional_cell": value[
                            "volume_a3_conventional_cell"
                        ],
                        "volume_standard_error_a3": value["volume_sigma_a3"],
                        "reported_pressure_source_token": value[
                            "source_pressure_token"
                        ],
                        "reported_pressure_gpa": value["printed_pressure_gpa"],
                        "single_peak_mgo": int(value["single_peak_mgo"]),
                    }
                )


def load_dataset(path: Path) -> list[dict[str, Any]]:
    """Load the bundled long-form CSV back into aligned source observations."""
    grouped: dict[int, dict[str, Any]] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            observation_index = int(row["observation_index"])
            observation = grouped.setdefault(
                observation_index,
                {
                    "run": row["run"],
                    "source_page": int(row["source_pdf_page"]),
                    "values": {},
                },
            )
            if observation["run"] != row["run"]:
                raise ValueError(f"mixed runs in observation {observation_index}")
            material = row["material"]
            observation["values"][material] = {
                "source_volume_token": row["volume_source_token"],
                "volume_a3_conventional_cell": float(
                    row["volume_a3_conventional_cell"]
                ),
                "volume_sigma_a3": (
                    float(row["volume_standard_error_a3"])
                    if row["volume_standard_error_a3"]
                    else None
                ),
                "single_peak_mgo": bool(int(row["single_peak_mgo"])),
                "source_pressure_token": row["reported_pressure_source_token"],
                "printed_pressure_gpa": float(row["reported_pressure_gpa"]),
            }
    expected = list(range(1, len(grouped) + 1))
    if sorted(grouped) != expected:
        raise ValueError("observation indices are not contiguous")
    return [grouped[index] for index in expected]


def _pairs(observations: list[dict[str, Any]]) -> list[tuple[str, float, str, float]]:
    pairs = []
    for observation in observations:
        row_values = observation["values"]
        for left, right in itertools.combinations(row_values, 2):
            pairs.append(
                (
                    left,
                    row_values[left]["volume_a3_conventional_cell"],
                    right,
                    row_values[right]["volume_a3_conventional_cell"],
                )
            )
    return pairs


def _fit_parameters(mode: str, vector: np.ndarray) -> dict[str, dict[str, float]]:
    result = {name: dict(values) for name, values in REFERENCE_PARAMETERS.items()}
    index = 0
    for material in MATERIAL_ORDER:
        if material == "MgO":
            continue
        if mode == "free" or material in ("NaCl B2", "Ne"):
            result[material]["K0"] = float(vector[index])
            index += 1
        result[material]["K0_prime"] = float(vector[index])
        index += 1
    if index != len(vector):
        raise AssertionError((index, len(vector)))
    return result


def fit(observations: list[dict[str, Any]], mode: str) -> dict[str, Any]:
    if mode not in ("fixed", "free"):
        raise ValueError(mode)
    initial_parameters = (
        REFERENCE_PARAMETERS if mode == "fixed" else PUBLISHED_FREE_PARAMETERS
    )
    initial = []
    parameter_names = []
    for material in MATERIAL_ORDER:
        if material == "MgO":
            continue
        if mode == "free" or material in ("NaCl B2", "Ne"):
            initial.append(initial_parameters[material]["K0"])
            parameter_names.append(f"{material}.K0")
        initial.append(initial_parameters[material]["K0_prime"])
        parameter_names.append(f"{material}.K0_prime")
    pairs = _pairs(observations)

    def residuals_for(parameters: dict[str, dict[str, float]]) -> np.ndarray:
        result = []
        for left, left_volume, right, right_volume in pairs:
            left_pressure = vinet_pressure(left_volume, parameters[left])
            right_pressure = vinet_pressure(right_volume, parameters[right])
            mean_pressure = (left_pressure + right_pressure) / 2.0
            result.append((left_pressure - right_pressure) / math.sqrt(mean_pressure))
        return np.asarray(result)

    def residuals(vector: np.ndarray) -> np.ndarray:
        return residuals_for(_fit_parameters(mode, vector))

    solution = least_squares(residuals, np.asarray(initial), method="lm")
    parameters = _fit_parameters(mode, solution.x)
    published_parameters = (
        REFERENCE_PARAMETERS if mode == "fixed" else PUBLISHED_FREE_PARAMETERS
    )
    comparisons = {}
    for material in MATERIAL_ORDER:
        if material == "MgO":
            continue
        comparisons[material] = {}
        for name in ("K0", "K0_prime"):
            if name == "K0" and mode == "fixed" and material not in ("NaCl B2", "Ne"):
                continue
            fitted = parameters[material][name]
            published = published_parameters[material][name]
            comparisons[material][name] = {
                "published": published,
                "refit": fitted,
                "difference": fitted - published,
                "relative_difference": (fitted - published) / published,
            }
    return {
        "success": bool(solution.success),
        "message": solution.message,
        "objective": float(np.dot(solution.fun, solution.fun)),
        "pair_count": len(pairs),
        "free_parameter_count": len(initial),
        "degrees_of_freedom": len(pairs) - len(initial),
        "parameter_order": parameter_names,
        "parameters": parameters,
        "published_objective": float(
            np.dot(
                residuals_for(published_parameters), residuals_for(published_parameters)
            )
        ),
        "published_parameter_comparison": comparisons,
    }


def reproduce_observations(observations: list[dict[str, Any]]) -> dict[str, Any]:
    run_counts: dict[str, int] = defaultdict(int)
    printed_pressure_ranges: dict[str, list[float]] = {}
    starred_mgo_rows = 0
    source_pressure_errors = []
    for observation in observations:
        run_counts[observation["run"]] += 1
        starred_mgo_rows += any(
            value["single_peak_mgo"] for value in observation["values"].values()
        )
        for material, value in observation["values"].items():
            pressure = value["printed_pressure_gpa"]
            if material not in printed_pressure_ranges:
                printed_pressure_ranges[material] = [pressure, pressure]
            else:
                printed_pressure_ranges[material][0] = min(
                    printed_pressure_ranges[material][0], pressure
                )
                printed_pressure_ranges[material][1] = max(
                    printed_pressure_ranges[material][1], pressure
                )
            source_pressure_errors.append(
                vinet_pressure(
                    value["volume_a3_conventional_cell"],
                    REFERENCE_PARAMETERS[material],
                )
                - value["printed_pressure_gpa"]
            )
    return {
        "source": {
            "doi": DOI,
            "correction_doi": CORRECTION_DOI,
            "article_url": ARTICLE_URL,
            "article_sha256": ARTICLE_SHA256,
            "supporting_filename": SOURCE_FILENAME,
            "supporting_url": SOURCE_URL,
            "supporting_sha256": SOURCE_SHA256,
            "source_locations": [
                "Equation (2)",
                "Equation (3)",
                "Table 2",
                "Tables S1-S6",
            ],
            "license": (
                "No explicit reusable license was identified for the publisher PDF; "
                "it is not relicensed by the CSV CC0 dedication."
            ),
        },
        "reproduction": {
            "status": "source_rows_refitted_but_published_coefficients_not_reproduced",
            "extractor": "Poppler pdftotext -tsv, pages 1-5",
            "optimizer": "scipy.optimize.least_squares(method='lm')",
            "interpretation": (
                "All unordered material pairs present in each source row are fitted "
                "simultaneously. Printed pressures are used only for the independent "
                "transcription check, never as fit observations."
            ),
            "known_limitations": [
                "Run AN012 is listed in article Table 1 but absent from Tables S1-S6.",
                "The article supplies neither fit code nor a covariance matrix.",
                "The published coefficients have a higher Equation (3) objective than "
                "the refit for both reported constraint modes.",
            ],
        },
        "normalization": {
            "volume": "angstrom^3 per conventional unit cell",
            "formula_units_per_cell": {
                "Au": 4,
                "MgO": 4,
                "Mo": 2,
                "NaCl B2": 1,
                "Ne": 4,
                "Pt": 4,
            },
            "temperature_k": 300.0,
        },
        "selection": {
            "observation_rows": len(observations),
            "volume_values": sum(len(row["values"]) for row in observations),
            "pair_count": len(_pairs(observations)),
            "run_counts": dict(sorted(run_counts.items())),
            "printed_pressure_range_gpa_by_material": dict(
                sorted(printed_pressure_ranges.items())
            ),
            "single_peak_mgo_rows_included": starred_mgo_rows,
            "rule": (
                "Use every nonblank paired volume printed in Tables S1-S6; include "
                "the three MgO values marked as determined from the (200) peak alone."
            ),
            "weighting": (
                "Equation (3): each pair contributes (P_a-P_b)^2/mean(P_a,P_b); "
                "the parenthetical volume precision is not an additional fit weight."
            ),
            "missing_run": (
                "Run AN012 is listed in article Table 1 (Au + NaCl, no medium, "
                "1-164 GPa) but is not present in auxiliary Tables S1-S6."
            ),
        },
        "published_pressure_transcription_check": {
            "parameter_set": "Table 2, K0 fixed for metals",
            "prediction_count": len(source_pressure_errors),
            "max_abs_difference_gpa": float(np.max(np.abs(source_pressure_errors))),
            "rms_difference_gpa": float(
                np.sqrt(np.mean(np.square(source_pressure_errors)))
            ),
        },
        "fits": {mode: fit(observations, mode) for mode in ("fixed", "free")},
    }


def reproduce(source_pdf: Path) -> dict[str, Any]:
    return reproduce_observations(extract_observations(source_pdf))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dataset-output", type=Path)
    args = parser.parse_args()
    observations = extract_observations(args.source_pdf)
    if args.dataset_output is not None:
        write_dataset(observations, args.dataset_output)
    result = reproduce_observations(observations)
    if args.dataset_output is not None:
        result["dataset_release"] = {
            "filename": DATASET_FILENAME,
            "sha256": sha256(args.dataset_output),
            "license": DATASET_LICENSE,
            "license_url": DATASET_LICENSE_URL,
            "scope": (
                "Peritheos-created factual CSV transcription, normalization, column "
                "naming, and arrangement, solely to the extent contributors hold rights; "
                "the article, publisher PDF, and third-party rights are excluded."
            ),
        }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
