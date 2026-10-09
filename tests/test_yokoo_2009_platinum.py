"""Checks distinguish Yokoo's Vinet fit, thermal checkpoints and observations."""

import hashlib
import json

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_yokoo_2009_platinum import (
    DATA,
    OUTPUT,
    RECORD,
    SOURCE,
    V0,
    check_saved,
    gamma_theta,
    ledger_outcome,
    phonon_pressure,
    reproduce,
    source_vinet,
    table_rows,
)


def test_vinet_source_branch_reconstructs_and_roundtrips():
    record = get_eos_record(RECORD)
    volume = V0 * np.array([1, 0.98, 0.9, 0.8, 0.7, 0.6])
    expected = source_vinet(volume)
    assert record.pressure(volume) == pytest.approx(expected, abs=1e-9)
    assert record.volume(expected) == pytest.approx(volume, rel=1e-9)
    restored = Material.from_eosmat(
        get_material_document("platinum"), record_identifiers=[RECORD]
    ).eos_records[0]
    assert restored.pressure(volume) == pytest.approx(expected, abs=1e-9)
    raw = next(
        r
        for r in get_material_document("platinum")["eos_records"]
        if r["identifier"] == RECORD
    )
    assert raw["equation_kind"] == "isothermal"
    assert "thermal" not in raw
    assert raw["parameter_errors"] == dict(V0=None, K0=None, K0_prime=None)
    assert raw["parameter_error_confidence"] is None
    # The cold BM3 error width must not migrate to the derived Vinet fit.
    assert (
        json.loads((DATA / SOURCE).read_text(encoding="utf-8"))["table4"][
            "B0_prime_0k_error"
        ]
        == 0.10
    )


def test_table5_is_derived_output_and_retains_phase_annotations():
    rows = table_rows()
    assert len(rows) == 168
    assert len({r["volume_ratio"] for r in rows}) == 21
    liquid = [r for r in rows if r["source_phase_annotation"] == "first_liquid_state"]
    assert [(r["volume_ratio"], r["temperature_k"]) for r in liquid] == [
        ("1.00", "3000"),
        ("0.98", "3000"),
    ]
    audit = reproduce()
    assert audit["max_abs_vinet_minus_full_model_300k_table_gpa"] == pytest.approx(
        6.00925735085, abs=1e-8
    )
    assert audit["cold_checkpoint_diagnostic"]["max_abs_residual_gpa"] < 0.006
    assert (
        audit["cold_checkpoint_diagnostic"]["kind"]
        == "derived_output_only_not_independent_fit"
    )
    assert ledger_outcome({"identifier": RECORD})["status"] == "not_refittable"
    assert ledger_outcome({"identifier": RECORD})["observations"] == 0


def test_phonon_law_has_integrated_gamma_and_requires_additional_term():
    for ratio in [1, 0.8, 0.6]:
        delta = 1e-5
        derivative = (
            np.log(gamma_theta(ratio * np.exp(delta))[1])
            - np.log(gamma_theta(ratio * np.exp(-delta))[1])
        ) / (2 * delta)
        assert derivative == pytest.approx(-gamma_theta(ratio)[0], rel=1e-8)
    assert phonon_pressure(1, 0) == 0
    audit = reproduce()
    endpoint = audit["table5_thermal_term_diagnostic"][-1]
    assert float(endpoint["temperature_k"]) == 3000
    assert endpoint["table_total_minus_0k_minus_phonon_gpa"] > 1.9
    assert audit["independent_equation_max_difference_gpa"] < 1e-9
    assert all(
        row["rounding_intervals_overlap"]
        for row in audit["upstream_shock_momentum_checks"]
    )


def test_saved_audit_and_inputs_remain_reproducible():
    audit = reproduce()
    saved = json.loads(OUTPUT.read_text(encoding="utf-8"))
    check_saved(saved, audit)
    assert saved["input_sha256"] == audit["input_sha256"]
    for filename, digest in audit["input_sha256"].items():
        assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == digest
    for key in [
        "independent_equation_max_difference_gpa",
        "max_abs_vinet_minus_full_model_300k_table_gpa",
        "density_derived_V0_cell_a3",
    ]:
        assert audit[key] == pytest.approx(saved[key], abs=1e-9)
    for current, stored in zip(
        audit["table5_thermal_term_diagnostic"],
        saved["table5_thermal_term_diagnostic"],
    ):
        assert current["table_total_minus_0k_minus_phonon_gpa"] == pytest.approx(
            stored["table_total_minus_0k_minus_phonon_gpa"], abs=1e-9
        )


def test_sakai_2018_links_the_vinet_branch_without_claiming_paired_data():
    re_record = next(
        r
        for r in get_material_document("rhenium")["eos_records"]
        if r["identifier"] == "rhenium_sakai_2018_yokoo_pt_vinet"
    )
    calibration = re_record["pressure_calibration"]
    assert calibration["methods"][0]["reference_eos_record"] == RECORD
    assert calibration["status"] == "resolved"
    assert calibration["recalculation"]["status"] == "missing_calibrant_observations"
    assert get_eos_record(RECORD).pressure(V0) == pytest.approx(0)


def test_paired_volumes_retain_original_pressure_calibration():
    pairs = reproduce()["dewaele_simultaneous_volume_pairs"]
    assert len(pairs) == 36
    assert pairs[0]["gold_atomic_volume_a3"] == 16.716
    assert pairs[0]["platinum_atomic_volume_a3"] == 14.979
    assert all(p["role"] == "comparison_only_not_yokoo_fit_target" for p in pairs)
    assert any(
        abs(
            p["derived_yokoo_vinet_pt_pressure_gpa"]
            - p["original_revised_ruby_pressure_gpa"]
        )
        > 1
        for p in pairs
    )
