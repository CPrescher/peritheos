#!/usr/bin/env python3
"""Reproduce the executable pressure-volume EOS branches in Fu et al. (2023)."""

from __future__ import annotations

import csv
import json
from importlib import resources
from math import isclose, sqrt

from peritheos import Material, get_material_document

BRIDGMANITE = "mg088fe010al014si090o3_bridgmanite"
BRIDGMANITE_RECORDS = (
    "mg088fe010al014si090o3_bridgmanite_fu_2023_bm2",
    "mg088fe010al014si090o3_bridgmanite_fu_2023_bm3",
)
CA_RECORD = "ca_perovskite_fu_2023_bm3_mgd_refit"


def reproduce() -> dict[str, object]:
    document = get_material_document(BRIDGMANITE)
    dataset = next(
        row
        for row in document["datasets"]
        if row["identifier"] == "mg088fe010al014si090o3_bridgmanite_fu_2023_table_s2_pv"
    )
    resource = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with resource.open("r", encoding="utf-8", newline="") as stream:
        observations = list(csv.DictReader(stream))
    assert len(observations) == 43

    results: dict[str, object] = {}
    for identifier in BRIDGMANITE_RECORDS:
        model = Material.from_eosmat(
            document, record_identifiers=[identifier]
        ).get_eos_record(identifier)
        residuals = [
            float(model.pressure(float(row["volume_a3"]))) - float(row["pressure_gpa"])
            for row in observations
        ]
        rmse = sqrt(sum(value * value for value in residuals) / len(residuals))
        assert rmse < 0.7
        test_volume = 0.86 * 163.75
        pressure = float(model.pressure(test_volume))
        assert isclose(float(model.volume(pressure)), test_volume, rel_tol=2e-10)
        results[identifier] = {
            "rows": len(observations),
            "pressure_rmse_gpa": rmse,
            "max_abs_pressure_residual_gpa": max(map(abs, residuals)),
        }

    ca_document = get_material_document("ca_perovskite")
    ca_model = Material.from_eosmat(
        ca_document, record_identifiers=[CA_RECORD]
    ).get_eos_record(CA_RECORD)
    checkpoints = {}
    for temperature in (300.0, 1200.0, 2200.0):
        pressure = float(ca_model.pressure(38.0, temperature))
        recovered = float(ca_model.volume(pressure, temperature))
        assert isclose(recovered, 38.0, rel_tol=2e-10)
        checkpoints[str(int(temperature))] = pressure
    assert checkpoints["2200"] > checkpoints["300"]
    results[CA_RECORD] = {"pressure_at_38_a3_gpa": checkpoints}
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
