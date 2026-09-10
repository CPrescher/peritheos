#!/usr/bin/env python3
"""Reconstruct the observable core of the Sokolova et al. (2013) calibration.

The source's complete thermodynamic optimization is not reproducible: it does
not publish all rows, weights, or optimizer details.  This script therefore
performs one explicitly narrower calculation.  It combines all eleven bundled
room-temperature marker series in one objective, converts their published ruby
pressures back to the common measured R1 shift, evaluates the published Table 4
Holzapfel curves, and fits only the two shared Equation (20) ruby coefficients.
It never refits eleven independent P-V curves.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATASET_ROOT = ROOT / "peritheos" / "data" / "datasets"
MANIFEST_PATH = DATASET_ROOT / "sokolova-2013-global-calibration-manifest.json"
OUTPUT_PATH = ROOT / "docs" / "data" / "sokolova-2013-global-calibration.json"
ANGSTROM3_PER_CM3_MOL = 1.0e24 / 6.02214076e23


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _power_law_shift(pressure: float, *, a_gpa: float, b: float) -> float:
    return math.pow(1.0 + b * pressure / a_gpa, 1.0 / b) - 1.0


def _holzapfel_2005_shift(pressure: float) -> float:
    def residual(ratio: float) -> float:
        exponent = (14.7 + 7.5) / 7.5 * (1.0 - ratio**-7.5)
        return 1845.0 / (14.7 + 7.5) * math.expm1(exponent) - pressure

    return brentq(residual, 1.0, 2.0) - 1.0


def _source_shift(calibration: str, value: float) -> float:
    if calibration == "direct_relative_shift":
        return value
    if calibration == "ruby_mao_1986":
        return _power_law_shift(value, a_gpa=1904.0, b=7.665)
    if calibration == "ruby_dewaele_2004":
        return _power_law_shift(value, a_gpa=1904.0, b=9.5)
    if calibration == "ruby_holzapfel_2005":
        return _holzapfel_2005_shift(value)
    raise ValueError(f"unsupported source calibration {calibration!r}")


def _volume_ratio(marker: dict[str, Any], row: dict[str, str]) -> float:
    material = marker["material"]
    value = float(row[marker["volume_column"]])
    if material == "niobium":
        return value
    atomic_reference = marker["V0_cm3_mol"] * ANGSTROM3_PER_CM3_MOL
    if material == "diamond":
        return value**3 / 8.0 / atomic_reference
    if material == "mgo":
        return value / (4.0 * atomic_reference)
    return value / atomic_reference


def _holzapfel_pressure(
    marker: dict[str, Any], volume_ratio: float, *, table: str
) -> float:
    """Evaluate a published Sokolova Holzapfel isotherm directly."""
    k0 = marker[f"{table}_K0_gpa"]
    k0_prime = marker[f"{table}_K0_prime"]
    n = marker["n"]
    z = marker["Z"]
    v0 = marker["V0_cm3_mol"]
    x = volume_ratio ** (1.0 / 3.0)
    fermi_pressure = 1003.6 * (z * n / v0) ** (5.0 / 3.0)
    c0 = -math.log(3.0 * k0 / fermi_pressure)
    c2 = 1.5 * (k0_prime - 3.0) - c0
    return (
        3.0
        * k0
        * math.exp(c0 * (1.0 - x))
        * (x**-5 - x**-4)
        * (1.0 + c2 * x - c2 * x * x)
    )


def _vinet_pressure(marker: dict[str, Any], volume_ratio: float) -> float:
    """Evaluate the Table 3 Vinet isotherm used in the ruby calibration."""
    k0 = marker["table3_K0_gpa"]
    k0_prime = marker["table3_K0_prime"]
    x = volume_ratio ** (1.0 / 3.0)
    eta = 1.5 * k0_prime - 1.5
    return 3.0 * k0 * (1.0 - x) * x**-2 * math.exp(eta * (1.0 - x))


def load_observations(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    observations = []
    for marker in manifest["markers"]:
        with (DATASET_ROOT / marker["resource"]).open(
            newline="", encoding="utf-8"
        ) as handle:
            rows = csv.DictReader(handle)
            for source_row, row in enumerate(rows, start=1):
                try:
                    if marker["material"] == "diamond" and (
                        row["pressure_gauge"] != "ruby"
                        or float(row["temperature_k"]) != 298.0
                    ):
                        continue
                    pressure_or_shift = row[marker["pressure_column"]]
                    if pressure_or_shift == "":
                        continue
                    shift = _source_shift(
                        marker["source_calibration"], float(pressure_or_shift)
                    )
                    volume_ratio = _volume_ratio(marker, row)
                    table2_pressure = _holzapfel_pressure(
                        marker, volume_ratio, table="table2"
                    )
                    table3_pressure = _vinet_pressure(marker, volume_ratio)
                    table4_pressure = _holzapfel_pressure(
                        marker, volume_ratio, table="table4"
                    )
                except (KeyError, TypeError, ValueError):
                    continue
                if (
                    shift < 0.0
                    or min(table2_pressure, table3_pressure, table4_pressure) < 0.0
                ):
                    continue
                observations.append(
                    {
                        "material": marker["material"],
                        "record_identifier": marker["record_identifier"],
                        "dataset_identifier": marker["dataset_identifier"],
                        "source_row": source_row,
                        "relative_r1_shift": shift,
                        "volume_ratio": volume_ratio,
                        "table2_holzapfel_pressure_gpa": table2_pressure,
                        "table3_vinet_pressure_gpa": table3_pressure,
                        "calibration_target_pressure_gpa": 0.5
                        * (table2_pressure + table3_pressure),
                        "table4_pressure_gpa": table4_pressure,
                    }
                )
    return observations


def reconstruct() -> dict[str, Any]:
    manifest = _load_json(MANIFEST_PATH)
    observations = load_observations(manifest)
    counts = Counter(row["material"] for row in observations)
    expected = {row["material"] for row in manifest["markers"]}
    if set(counts) != expected:
        raise ValueError(
            f"missing marker observations: {sorted(expected - set(counts))}"
        )

    def raw_residuals(parameters: tuple[float, float] | np.ndarray) -> np.ndarray:
        a_gpa, m = parameters
        return np.asarray(
            [
                a_gpa * row["relative_r1_shift"] * (1.0 + m * row["relative_r1_shift"])
                - row["calibration_target_pressure_gpa"]
                for row in observations
            ]
        )

    def marker_equal_residuals(
        parameters: tuple[float, float] | np.ndarray,
    ) -> np.ndarray:
        raw = raw_residuals(parameters)
        return np.asarray(
            [
                value / math.sqrt(counts[row["material"]])
                for value, row in zip(raw, observations)
            ]
        )

    published = np.asarray([1870.0, 6.0])
    fit = least_squares(
        marker_equal_residuals,
        x0=published,
        bounds=([1500.0, 0.0], [2200.0, 20.0]),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    fitted = np.asarray(fit.x)

    marker_results = []
    for marker in manifest["markers"]:
        material = marker["material"]
        indices = [
            index
            for index, observation in enumerate(observations)
            if observation["material"] == material
        ]
        published_residual = raw_residuals(published)[indices]
        fitted_residual = raw_residuals(fitted)[indices]
        marker_results.append(
            {
                "material": material,
                "record_identifier": marker["record_identifier"],
                "dataset_identifier": marker["dataset_identifier"],
                "observations": len(indices),
                "relative_r1_shift_range": [
                    min(observations[index]["relative_r1_shift"] for index in indices),
                    max(observations[index]["relative_r1_shift"] for index in indices),
                ],
                "table4_pressure_range_gpa": [
                    min(
                        observations[index]["table4_pressure_gpa"] for index in indices
                    ),
                    max(
                        observations[index]["table4_pressure_gpa"] for index in indices
                    ),
                ],
                "published_calibration_rmse_gpa": float(
                    np.sqrt(np.mean(np.square(published_residual)))
                ),
                "published_calibration_mean_residual_gpa": float(
                    np.mean(published_residual)
                ),
                "fitted_calibration_rmse_gpa": float(
                    np.sqrt(np.mean(np.square(fitted_residual)))
                ),
                "table4_closure_rmse_gpa": float(
                    np.sqrt(
                        np.mean(
                            np.square(
                                [
                                    1870.0
                                    * observations[index]["relative_r1_shift"]
                                    * (
                                        1.0
                                        + 6.0 * observations[index]["relative_r1_shift"]
                                    )
                                    - observations[index]["table4_pressure_gpa"]
                                    for index in indices
                                ]
                            )
                        )
                    )
                ),
            }
        )

    published_equal = marker_equal_residuals(published)
    fitted_equal = marker_equal_residuals(fitted)
    return {
        "format": "peritheos.sokolova-2013-global-calibration-reconstruction",
        "format_version": 1,
        "generated_with": "scripts/reproduce_sokolova_2013_global_calibration.py",
        "manifest": "peritheos/data/datasets/sokolova-2013-global-calibration-manifest.json",
        "classification": "source_constrained_cross_calibration_reconstruction",
        "independent_eos_refit": False,
        "qualification": (
            "All eleven machine-readable comparison series enter one shared "
            "Equation (20) objective against the mean Table 2 Holzapfel and Table 3 "
            "Vinet isotherm. Table 4 is retained as a post-calibration closure test. "
            "The calculation cannot recover the eleven EOS "
            "parameter sets independently because the source omits complete "
            "thermochemical/ultrasonic rows, weights, and covariance."
        ),
        "objective": manifest["reconstructable_cross_calibration"]["objective"],
        "observations": len(observations),
        "markers": len(counts),
        "published_parameters": {"A_gpa": 1870.0, "m": 6.0},
        "sensitivity_fit": {
            "parameters": {"A_gpa": float(fitted[0]), "m": float(fitted[1])},
            "bounds": manifest["reconstructable_cross_calibration"]["constraints"],
            "success": bool(fit.success),
            "message": str(fit.message),
            "A_relative_difference": float((fitted[0] - published[0]) / published[0]),
            "m_relative_difference": float((fitted[1] - published[1]) / published[1]),
        },
        "published_calibration_marker_equal_rmse_gpa": float(
            np.sqrt(np.sum(np.square(published_equal)) / len(counts))
        ),
        "fitted_calibration_marker_equal_rmse_gpa": float(
            np.sqrt(np.sum(np.square(fitted_equal)) / len(counts))
        ),
        "table4_closure_marker_equal_rmse_gpa": float(
            np.sqrt(
                np.mean(
                    [item["table4_closure_rmse_gpa"] ** 2 for item in marker_results]
                )
            )
        ),
        "marker_results": marker_results,
        "full_optimization_blockers": manifest["published_objective"]["not_published"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="compare with committed JSON"
    )
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    result = reconstruct()
    rendered = json.dumps(result, indent=2, sort_keys=False) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale reconstruction: {args.output}")
    else:
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
