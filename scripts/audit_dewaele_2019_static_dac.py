#!/usr/bin/env python3
"""Audit independently refittable Dewaele (2019) static-DAC metal rows.

The paper reports fifteen metals reduced on two ruby pressure scales.  This
script deliberately separates coefficient verification from row-level
reproduction.  It refits only source rows that are currently bundled and
reports the source gaps that prevent the complete 30-record family from being
called reproducible.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import NamedTuple

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"
MATERIALS = ROOT / "peritheos" / "data" / "materials"


class SourceSpec(NamedTuple):
    material_file: str
    record_prefix: str
    dataset: str
    source_scale: str
    pressure_column: str
    volume_column: str
    volume_divisor: float = 1.0
    lattice_111_column: str | None = None
    scope: str = "complete"
    include_column: str | None = None
    include_value: str = "true"


SOURCES = (
    SourceSpec(
        "gold.eosmat",
        "gold_dewaele_2019",
        "gold-takemura-2008-table3-compression.csv",
        "dor",
        "ruby_pressure_dorogokupets_gpa",
        "",
        lattice_111_column="lattice_111_angstrom",
    ),
    SourceSpec(
        "copper.eosmat",
        "copper_dewaele_2019",
        "copper-dewaele-2004-table1-compression.csv",
        "mao",
        "ruby_pressure_classical_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "platinum.eosmat",
        "platinum_dewaele_2019",
        "platinum-dewaele-2004-table1-compression.csv",
        "mao",
        "ruby_pressure_classical_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "tantalum.eosmat",
        "tantalum_dewaele_2019",
        "tantalum-dewaele-2004-table1-compression.csv",
        "mao",
        "ruby_pressure_classical_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "aluminum.eosmat",
        "aluminum_dewaele_2019",
        "aluminum-dewaele-2004-table1-compression.csv",
        "mao",
        "ruby_pressure_classical_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "tungsten.eosmat",
        "tungsten_dewaele_2019",
        "tungsten-dewaele-2004-table1-compression.csv",
        "mao",
        "ruby_pressure_classical_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "silver.eosmat",
        "silver_dewaele_2019",
        "silver-dewaele-2008-table2-compression.csv",
        "mao",
        "ruby_pressure_mao_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "cobalt_hcp.eosmat",
        "cobalt_hcp_dewaele_2019",
        "cobalt-dewaele-2008-table2-compression.csv",
        "mao",
        "ruby_pressure_mao_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "molybdenum.eosmat",
        "molybdenum_dewaele_2019",
        "molybdenum-dewaele-2008-table2-compression.csv",
        "mao",
        "ruby_pressure_mao_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "nickel.eosmat",
        "nickel_dewaele_2019",
        "nickel-dewaele-2008-table2-compression.csv",
        "mao",
        "ruby_pressure_mao_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "zinc_hcp.eosmat",
        "zinc_hcp_dewaele_2019",
        "zinc-dewaele-2008-table2-compression.csv",
        "mao",
        "ruby_pressure_mao_gpa",
        "atomic_volume_a3",
    ),
    SourceSpec(
        "rhenium.eosmat",
        "rhenium_dewaele_2019",
        "rhenium-anzellini-2014-table3-compression.csv",
        "dor",
        "pressure_gpa",
        "volume_a3_conventional_cell",
        volume_divisor=2.0,
    ),
    SourceSpec(
        "iron.eosmat",
        "iron_dewaele_2019",
        "iron-dewaele-2006-epaps-compression.csv",
        "dor",
        "reported_pressure_gpa",
        "atomic_volume_a3_per_atom",
        scope=(
            "complete 53-row hcp selection; reported pressure is the source's "
            "ruby-linked effective pressure for both ruby and W-gauge rows"
        ),
        include_column="dewaele_2019_vinet_fit_included",
    ),
    SourceSpec(
        "lead_hcp.eosmat",
        "lead_hcp_dewaele_2019",
        "lead-dewaele-2019-table2-compression.csv",
        "dor",
        "pressure_dorogokupets_gpa",
        "atomic_volume_a3",
        scope=(
            "current_run_only; a separately bundled Figure 3 digitization of "
            "Kuznetsov et al. (2002) is not treated as exact fit input"
        ),
    ),
)

SOURCE_GAPS = {
    "beryllium_hcp": (
        "Lazicki et al. (2012) Table I prints an internally inconsistent a=2.262 A "
        "at 24.6 GPa; the unrounded/corrected author row is unavailable."
    ),
    "lead_hcp": (
        "Dewaele (2019) fits its Table 2 run together with Kuznetsov et al. "
        "(2002). Twenty-one hcp square markers are now bundled as an explicit "
        "Figure 3 digitization, but the exact numerical rows, selected subset, "
        "weights, and Brown-NaCl pressure-scale treatment remain unavailable."
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _vinet_pressure(parameters: np.ndarray, volume: np.ndarray | float):
    v0, k0, k0_prime = parameters
    x = (np.asarray(volume) / v0) ** (1.0 / 3.0)
    return 3.0 * k0 * (1.0 - x) / x**2 * np.exp(1.5 * (k0_prime - 1.0) * (1.0 - x))


def _mao_to_dor(pressure: np.ndarray) -> np.ndarray:
    ratio = (1.0 + 7.665 * pressure / 1904.0) ** (1.0 / 7.665)
    shift = ratio - 1.0
    return 1884.0 * shift * (1.0 + 5.5 * shift)


def _dor_to_mao(pressure: np.ndarray) -> np.ndarray:
    shift = (-1.0 + np.sqrt(1.0 + 22.0 * pressure / 1884.0)) / 11.0
    ratio = 1.0 + shift
    return 1904.0 / 7.665 * (ratio**7.665 - 1.0)


def _target_pressure(source: np.ndarray, source_scale: str, target_scale: str):
    if source_scale == target_scale:
        return source.copy()
    if source_scale == "mao" and target_scale == "dor":
        return _mao_to_dor(source)
    if source_scale == "dor" and target_scale == "mao":
        return _dor_to_mao(source)
    raise ValueError(f"unsupported scale conversion: {source_scale} -> {target_scale}")


def _load_rows(spec: SourceSpec) -> tuple[np.ndarray, np.ndarray, str]:
    path = DATA / spec.dataset
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if spec.include_column:
        rows = [
            row
            for row in rows
            if row[spec.include_column].strip().lower() == spec.include_value.lower()
        ]
    pressure = np.array([float(row[spec.pressure_column]) for row in rows])
    if spec.lattice_111_column:
        lattice = np.array([float(row[spec.lattice_111_column]) for row in rows])
        volume = lattice**3 / 4.0
    else:
        volume = np.array([float(row[spec.volume_column]) for row in rows])
        volume /= spec.volume_divisor
    return pressure, volume, _sha256(path)


def _published(spec: SourceSpec, scale: str):
    path = MATERIALS / spec.material_file
    document = json.loads(path.read_text(encoding="utf-8"))
    identifier = f"{spec.record_prefix}_{scale}_vinet"
    record = next(
        row for row in document["eos_records"] if row["identifier"] == identifier
    )
    z = float(document["formula_units_per_cell"])
    parameters = record["eos"]["parameters"]
    errors = record["parameter_errors"]
    return (
        identifier,
        np.array(
            [
                float(parameters["V0"]) / z,
                float(parameters["K0"]),
                float(parameters["K0_prime"]),
            ]
        ),
        np.array(
            [
                float(errors["V0"]) / z,
                float(errors["K0"]),
                float(errors["K0_prime"]),
            ]
        ),
    )


def _fit_pressure(volume: np.ndarray, pressure: np.ndarray, initial: np.ndarray):
    result = least_squares(
        lambda parameters: _vinet_pressure(parameters, volume) - pressure,
        initial,
        bounds=([0.8 * initial[0], 1.0, 0.5], [1.3 * initial[0], 1000.0, 15.0]),
    )
    residual = _vinet_pressure(result.x, volume) - pressure
    return result.x, float(np.sqrt(np.mean(residual**2))), float(np.max(abs(residual)))


def _fit_volume(volume: np.ndarray, pressure: np.ndarray, initial: np.ndarray):
    def predicted_volumes(parameters: np.ndarray) -> np.ndarray:
        return np.array(
            [
                brentq(
                    lambda trial: float(_vinet_pressure(parameters, trial) - target),
                    0.25 * parameters[0],
                    parameters[0],
                )
                for target in pressure
            ]
        )

    result = least_squares(
        lambda parameters: predicted_volumes(parameters) - volume,
        initial,
        bounds=([0.8 * initial[0], 1.0, 0.5], [1.3 * initial[0], 1000.0, 15.0]),
    )
    residual = predicted_volumes(result.x) - volume
    return result.x, float(np.sqrt(np.mean(residual**2))), float(np.max(abs(residual)))


def _original_fe_2006_reproduction() -> dict[str, object]:
    path = DATA / "iron-dewaele-2006-epaps-compression.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [
            row
            for row in csv.DictReader(stream)
            if row["original_2006_helium_vinet_fit_included"] == "true"
        ]
    pressure = np.array([float(row["reported_pressure_gpa"]) for row in rows])
    volume = np.array([float(row["atomic_volume_a3_per_atom"]) for row in rows])
    published = np.array([11.214, 163.4, 5.38])
    errors = np.array([0.049, 7.9, 0.16])
    parameters, rmse, maximum = _fit_pressure(volume, pressure, published)
    return {
        "dataset": path.name,
        "dataset_sha256": _sha256(path),
        "rows": len(rows),
        "pressure_range_gpa": [float(pressure.min()), float(pressure.max())],
        "published_atomic_parameters": published.tolist(),
        "reported_95pct_errors": errors.tolist(),
        "unweighted_pressure_residual_fit": {
            "parameters": parameters.tolist(),
            "rmse_gpa": rmse,
            "max_abs_gpa": maximum,
            "absolute_parameter_difference": abs(parameters - published).tolist(),
            "within_reported_95pct_errors": bool(
                np.all(abs(parameters - published) <= errors)
            ),
        },
    }


def reproduce() -> dict[str, object]:
    fits: dict[str, object] = {}
    for spec in SOURCES:
        source_pressure, volume, checksum = _load_rows(spec)
        for scale in ("mao", "dor"):
            identifier, published, errors = _published(spec, scale)
            pressure = _target_pressure(source_pressure, spec.source_scale, scale)
            p_parameters, p_rmse, p_max = _fit_pressure(volume, pressure, published)
            v_parameters, v_rmse, v_max = _fit_volume(volume, pressure, published)
            fits[identifier] = {
                "dataset": spec.dataset,
                "dataset_sha256": checksum,
                "rows": len(volume),
                "scope": spec.scope,
                "pressure_range_gpa": [float(pressure.min()), float(pressure.max())],
                "published_atomic_parameters": published.tolist(),
                "reported_95pct_errors": errors.tolist(),
                "unweighted_pressure_residual_fit": {
                    "parameters": p_parameters.tolist(),
                    "rmse_gpa": p_rmse,
                    "max_abs_gpa": p_max,
                    "absolute_parameter_difference": abs(
                        p_parameters - published
                    ).tolist(),
                    "within_reported_95pct_errors": bool(
                        np.all(abs(p_parameters - published) <= errors)
                    ),
                },
                "unweighted_volume_residual_fit": {
                    "parameters": v_parameters.tolist(),
                    "rmse_a3_per_atom": v_rmse,
                    "max_abs_a3_per_atom": v_max,
                    "absolute_parameter_difference": abs(
                        v_parameters - published
                    ).tolist(),
                    "within_reported_95pct_errors": bool(
                        np.all(abs(v_parameters - published) <= errors)
                    ),
                },
            }

    identifiers = []
    for path in MATERIALS.glob("*.eosmat"):
        document = json.loads(path.read_text(encoding="utf-8"))
        identifiers.extend(
            row["identifier"]
            for row in document.get("eos_records", [])
            if "dewaele_2019" in row.get("identifier", "")
            and row["identifier"].endswith(("_mao_vinet", "_dor_vinet"))
        )
    return {
        "catalog_record_count": len(identifiers),
        "catalog_identifiers": sorted(identifiers),
        "row_level_refits": fits,
        "source_publication_checks": {
            "iron_dewaele_2006_helium_vinet": _original_fe_2006_reproduction()
        },
        "source_gaps": SOURCE_GAPS,
        "complete_family_reproducible": False,
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
