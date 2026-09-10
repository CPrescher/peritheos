#!/usr/bin/env python3
"""Reproduce the executable pressure-volume EOS branches in Fu et al. (2023)."""

from __future__ import annotations

import csv
import json
from importlib import resources
from math import isclose, sqrt
from pathlib import Path

from peritheos import Material, get_material_document

BRIDGMANITE = "mg088fe010al014si090o3_bridgmanite"
BRIDGMANITE_RECORDS = (
    "mg088fe010al014si090o3_bridgmanite_fu_2023_bm2",
    "mg088fe010al014si090o3_bridgmanite_fu_2023_bm3",
)
CA_RECORD = "ca_perovskite_fu_2023_bm3_mgd_refit"
CA_REFIT_RECORD = "ca_perovskite_fu_2023_candidate_data_unweighted_bm3_mgd_refit"
AUDIT_PATH = (
    Path(__file__).resolve().parents[1] / "docs/data/fu-2023-casio3-refit-audit.json"
)


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
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    assert audit["fit_observation_count"] == 174
    assert audit["fit_output_count"] == 242
    assert all(item["success"] for item in audit["sensitivity_fits"])
    assert all(
        not source.get("redistributed", False)
        for source in audit["sources"].values()
        if source.get("role", "").startswith("candidate fit observations")
    )
    results[CA_RECORD] = {
        "analytical_checkpoints": {
            "kind": "published_parameterization_round_trip",
            "fit_observations": False,
            "pressure_at_38_a3_gpa": checkpoints,
        },
        "coefficient_audit": {
            "status": "parity_not_achieved",
            "candidate_fit_observations": audit["fit_observation_count"],
            "fit_outputs": audit["fit_output_count"],
            "weighting_disclosure": audit["weighting_disclosure"],
            "sensitivity_fits": audit["sensitivity_fits"],
        },
    }
    assert audit["registered_refit"]["record_identifier"] == CA_REFIT_RECORD
    diagnostic = next(
        item
        for item in audit["sensitivity_fits"]
        if item["objective"] == audit["registered_refit"]["objective"]
    )
    refit_source = next(
        item
        for item in ca_document["eos_records"]
        if item["identifier"] == CA_REFIT_RECORD
    )
    expected = diagnostic["parameters"]
    assert isclose(refit_source["eos"]["parameters"]["V0"], expected["V0"])
    assert isclose(refit_source["eos"]["parameters"]["K0"], expected["K0"])
    assert isclose(refit_source["thermal"]["parameters"]["gamma0"], expected["gamma0"])
    assert isclose(refit_source["thermal"]["parameters"]["q"], expected["q"])
    refit_model = Material.from_eosmat(
        ca_document, record_identifiers=[CA_REFIT_RECORD]
    ).get_eos_record(CA_REFIT_RECORD)
    refit_checkpoints = {}
    for temperature in (700.0, 1200.0, 2200.0):
        pressure = float(refit_model.pressure(38.0, temperature))
        recovered = float(refit_model.volume(pressure, temperature))
        assert isclose(recovered, 38.0, rel_tol=2e-10)
        refit_checkpoints[str(int(temperature))] = pressure
    results[CA_REFIT_RECORD] = {
        "record_kind": refit_source["record_kind"],
        "objective": audit["registered_refit"]["objective"],
        "fit_observations": audit["fit_observation_count"],
        "fit_outputs": audit["fit_output_count"],
        "pressure_at_38_a3_gpa": refit_checkpoints,
    }
    return {"records": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
