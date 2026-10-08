"""Normalize bundled Dewaele KCl pressures and refit only the reference isotherm.

Run as ``python -m scripts.audit_kcl_pressure_normalization``. Derived fits are
diagnostics, never replacements for published EOS coefficients. No raw ruby
wavelengths or measurement uncertainty model are reconstructed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from peritheos import (
    get_material_document,
    get_pressure_calibration,
    recalculate_ruby_pressure,
    resolve_dataset_pressure,
)
from peritheos.eos.rt import Vinet
from peritheos.fitting import fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/kcl-pressure-normalization.json"
SOURCE = "ruby_dorogokupets_oganov_2007"
CASES = (
    ("kcl", "kcl_b2_dewaele_2012_vinet_3", "kcl_dewaele_2012_table1_compression"),
    ("kcl_b1", "kcl_b1_dewaele_2012_vinet_1", "kcl_b1_dewaele_2012_table2_compression"),
)


def audit() -> dict:
    results = []
    for material, record_id, dataset_id in CASES:
        document = get_material_document(material)
        record = next(
            r for r in document["eos_records"] if r["identifier"] == record_id
        )
        metadata = next(
            d for d in document["datasets"] if d["identifier"] == dataset_id
        )
        resolved = resolve_dataset_pressure(document, record_id, dataset_id)
        # These tables declare formula-unit volumes; the material EOS uses cells.
        volume = (
            resolved.dataset["volume_a3_per_formula_unit"]
            * document["formula_units_per_cell"]
        )
        initial = record["eos"]["parameters"]
        fixed = {key: initial[key] for key in record["fixed_parameters"]}
        reductions = {}
        for target in (SOURCE, "ruby_mao_1986"):
            pressure = (
                resolved.pressure_gpa
                if target == SOURCE
                else recalculate_ruby_pressure(resolved.pressure_gpa, SOURCE, target)
            )
            fit = fit_rt_eos(
                Vinet,
                volume,
                pressure,
                {key: value for key, value in initial.items() if key not in fixed},
                fixed=fixed,
            )
            if not fit.success:
                raise RuntimeError(f"KCl diagnostic fit failed: {record_id}, {target}")
            limits = get_pressure_calibration(target).validity["pressure_gpa"]
            outside = (pressure < limits[0]) | (pressure > limits[1])
            reductions[target] = {
                "pressure_gpa": [round(float(p), 10) for p in pressure],
                "pressure_extrapolated_row_indices": np.flatnonzero(outside).tolist(),
                "conditional_unweighted_vinet_refit": {
                    "parameters_conventional_cell": {
                        key: round(value, 10) for key, value in fit.parameters.items()
                    },
                    "fixed_parameters": fixed,
                    "rmse_gpa": round(float(np.sqrt(np.mean(fit.residuals**2))), 10),
                    "observation_temperature_k": 298.0,
                    "uncertainty_treatment": "not_propagated",
                },
            }
        results.append(
            {
                "record_identifier": record_id,
                "dataset_identifier": dataset_id,
                "resource_sha256": metadata["resource"]["sha256"],
                "observations": len(volume),
                "row_indices": resolved.row_indices.tolist(),
                "source_calibration_record": SOURCE,
                "normalization_convention": "same_corrected_ruby_r1_ratio",
                "original_published_parameters": initial,
                "reductions": reductions,
            }
        )
    return {
        "generated_with": "scripts/audit_kcl_pressure_normalization.py",
        "qualification": (
            "Rounded reported source-scale pressures, all bundled rows, conditional "
            "unweighted 298 K Vinet fits with source-fixed parameters. Not original "
            "raw-wavelength re-reduction or exact weighted-fit reproduction. Row-wise "
            "measurement uncertainties, source fit weights and covariance are not "
            "available. B2's MD thermal term is not refitted or transferred to a new "
            "scale. Published EOS coefficients and CSV observations are unchanged. "
            "Mao rows above 80 GPa are retained with explicit extrapolation indices."
        ),
        "records": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = json.dumps(audit(), indent=2, allow_nan=False) + "\n"
    if args.check:
        if OUTPUT.read_text() != expected:
            raise SystemExit(f"stale audit: {OUTPUT}")
    else:
        OUTPUT.write_text(expected)


if __name__ == "__main__":
    main()
