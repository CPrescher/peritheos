"""Keep FeS source discrepancies visible and verify the conditional evaluator."""

import hashlib
import json

import numpy as np
import pytest
from jsonschema import validate

from peritheos import (
    Material,
    eosmat_schema,
    get_material_document,
    list_eos_records,
    validate_eosmat_document,
)
from scripts.audit_morard_2026_fes import (
    CELL_PER_MOLAR,
    DATA,
    DATASET_ID,
    OUTPUT,
    RECORD_ID,
    ledger_outcome,
    observations,
    pressure,
    reproduce,
)


def candidate():
    return Material.from_eosmat(
        get_material_document("fes_vi"), require_primary_validation=False
    ).get_eos_record(RECORD_ID)


def test_source_record_is_preserved_but_deferred_from_normal_catalog():
    document = get_material_document("fes_vi")
    validate_eosmat_document(document)
    validate(document, eosmat_schema())
    record = document["eos_records"][0]
    assert record["scientific_validation"]["status"] == "deferred"
    assert record["default"] is False
    assert record["eos"]["parameters"] == {
        "V0": pytest.approx(102.28920653790897),
        "K0": 115.5,
        "K0_prime": 4.99,
    }
    assert record["parameter_errors"]["V0"] == pytest.approx(0.45 * CELL_PER_MOLAR)
    assert record["thermal"]["parameters"]["gamma0"] == 2.42
    assert record["thermal"]["parameters"]["n"] == 2
    assert RECORD_ID not in {r.identifier for r in list_eos_records()}
    with pytest.raises(ValueError, match="deferred"):
        Material.from_eosmat(document)
    material = Material.from_eosmat(document, require_primary_validation=False)
    exported = material.to_eosmat()
    assert (
        exported["eos_records"][0]["scientific_validation"]
        == record["scientific_validation"]
    )
    assert exported["datasets"][0] == document["datasets"][0]
    assert (
        Material.from_eosmat(exported, require_primary_validation=False).to_eosmat()
        == exported
    )


def test_complete_source_table_and_its_actual_uncertainty_conventions():
    d = observations()
    document = get_material_document("fes_vi")
    resource = document["datasets"][0]["resource"]
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == resource["sha256"]
    assert document["datasets"][0]["identifier"] == DATASET_ID
    assert len(d["source_row"]) == 146
    np.testing.assert_array_equal(d["source_row"], np.arange(3, 149))
    assert d["pressure_gpa"][[0, -1]].tolist() == [29.59643816464321, 149.9480701694589]
    np.testing.assert_array_equal(d["sigma_temperature_k"], np.full(146, 150))
    a, b, c = (d[f"{axis}_angstrom"] for axis in "abc")
    np.testing.assert_array_equal(d["volume_angstrom3"], a * b * c)
    linear_error = d["volume_angstrom3"] * sum(
        d[f"sigma_{axis}_angstrom"] / d[f"{axis}_angstrom"] for axis in "abc"
    )
    np.testing.assert_allclose(d["sigma_volume_angstrom3"], linear_error, rtol=1e-14)
    np.testing.assert_allclose(d["order_parameter_c"], c / (np.sqrt(3) * b), rtol=1e-14)
    assert not np.allclose(d["order_parameter_c"], np.sqrt(3) * b / c)


def test_independent_quadrature_native_pressure_inversion_and_validity():
    record = candidate()
    for volume, temperature in [(15.4, 300), (12, 300), (12, 1500), (10.1, 2500)]:
        expected = float(pressure(np.array(volume), np.array(temperature)))
        actual = record.pressure(
            volume * CELL_PER_MOLAR, temperature, check_validity=False
        )
        assert actual == pytest.approx(expected, abs=2e-7)
        assert record.volume(
            actual, temperature, check_validity=False
        ) == pytest.approx(volume * CELL_PER_MOLAR, rel=1e-8)
    for p, t in [(40, 1300), (100, 2000), (145, 1900)]:
        v = record.volume(p, t)
        assert record.pressure(v, t) == pytest.approx(p, abs=1e-7)
    with pytest.raises(ValueError, match="outside the published calibration"):
        record.volume(200, 3500, check_validity=True)


def test_replay_does_not_hide_failed_published_residual_claim():
    report = reproduce()
    assert (
        json.loads(OUTPUT.read_text(encoding="utf-8"))["source_rows"]
        == report["source_rows"]
    )
    integrated = report["published_coefficient_replays"]["integrated_gruneisen"]
    assert integrated["rmse_gpa"] == pytest.approx(1.9110674244, abs=1e-8)
    assert integrated["max_abs_gpa"] == pytest.approx(6.1493197517, abs=1e-8)
    assert len(integrated["rows_exceeding_3_gpa"]) == 15
    assert integrated["over_3_gpa_at_pressure_ge_40"] == 13
    diagnostic = report["conditional_gamma_only_fits"][
        "integrated_gruneisen:unweighted_pressure"
    ]
    assert diagnostic["gamma0"] == pytest.approx(2.28197939, abs=1e-7)
    assert diagnostic["conditional_standard_error"] > 0
    assert diagnostic["converged"]
    assert (
        report["published_coefficient_replays"]["variable_exponent"]["max_abs_gpa"] > 6
    )


def test_all_final_pressures_reproduce_from_paired_gsas_markers():
    report = reproduce()["calibration_replay"]
    assert report["matched_rows"] == 146
    assert report["max_abs_pressure_difference_gpa"] < 2e-12
    material = Material.from_eosmat(
        get_material_document("fes_vi"), require_primary_validation=False
    )
    paired = material.get_dataset("fes_morard_2026_paired_kcl")
    assert len(paired) == 146
    np.testing.assert_allclose(
        paired.values("replayed_pressure_gpa"),
        observations()["pressure_gpa"],
        atol=2e-12,
        rtol=0,
    )
    assert (
        material.eos_records[0].pressure_calibration["recalculation"]["status"]
        == "ready"
    )


def test_experimental_structure_and_source_files_are_independent():
    document = get_material_document("fes_vi")
    assert document["formula_units_per_cell"] == 4
    assert document["space_group"] == "Pnma"
    assert document["space_group_number"] == 62
    assert [
        (s["element"], s["multiplicity"], s["wyckoff"]) for s in document["atom_sites"]
    ] == [("Fe", 4, "4c"), ("S", 4, "4c")]
    assert (
        document["source"]["structure_reference"]["doi"] == "10.1016/j.epsl.2008.05.017"
    )
    lattice = document["lattice"]
    assert lattice["a"] * lattice["b"] * lattice["c"] == pytest.approx(78.54, abs=0.015)
    source_dir = DATA.parent / "morard_2026_sources"
    manifest = json.loads((source_dir / "manifest.json").read_text(encoding="utf-8"))
    for source in manifest["sources"]:
        if source["bundled"]:
            assert (
                hashlib.sha256(
                    (source_dir / source["filename"]).read_bytes()
                ).hexdigest()
                == source["sha256"]
            )
            assert source["license"] == "CC BY 4.0"


def test_staged_reproduction_is_similar_with_separate_residual_finding():
    record = get_material_document("fes_vi")["eos_records"][0]
    outcome = ledger_outcome(record)
    assert record["scientific_validation"]["reproduction_status"] == "similar"
    assert outcome["status"] == outcome["reproduction_status"] == "similar"
    assert outcome["residual_claim_reproduction_status"] == "not_reproduced"
    assert outcome["stage_observations"] == {"cold": 21, "thermal": 146, "quenched": 0}
    assert outcome["observations"] == 167
    assert all(
        p["similar"] and p["reported_intervals_overlap"] for p in outcome["parameters"]
    )
    assert all(p["within_combined_2sigma"] is None for p in outcome["parameters"])
    assert outcome["parameters"][-1]["refit"] == pytest.approx(2.46418)
