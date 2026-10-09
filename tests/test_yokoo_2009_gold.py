"""Verify Au reference branches without treating model grids as observations."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_yokoo_2009_gold import (
    DATA,
    OUTPUT,
    RECORDS,
    SOURCE,
    V0,
    check_saved,
    gamma_theta,
    ledger_outcome,
    phonon_pressure,
    reproduce,
    shock_relation_checkpoints,
    source_bm3,
    source_vinet,
    table_rows,
)


@pytest.mark.parametrize(
    "identifier,equation,kp",
    [
        ("gold_yokoo_2009_bm3_300k", source_bm3, 5.79),
        ("gold_yokoo_2009_vinet_300k", source_vinet, 5.94),
    ],
)
def test_300k_branches_match_source_equations_and_roundtrip(identifier, equation, kp):
    record = get_eos_record(identifier)
    ratios = np.array([1, 0.98, 0.9, 0.8, 0.7, 0.6])
    expected = equation(ratios)
    assert record.pressure(V0 * ratios) == pytest.approx(expected, abs=1e-9)
    assert record.volume(expected) == pytest.approx(V0 * ratios, rel=1e-9)
    document = get_material_document("gold")
    restored = Material.from_eosmat(document, record_identifiers=[identifier])
    assert restored.eos_records[0].pressure(V0 * ratios) == pytest.approx(
        expected, abs=1e-9
    )
    raw = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    assert raw["eos"]["parameters"] == dict(V0=67.72, K0=167.5, K0_prime=kp)
    assert raw["equation_kind"] == "isothermal"
    assert "thermal" not in raw
    assert raw["parameter_errors"] == dict(V0=None, K0=None, K0_prime=None)
    assert raw["parameter_error_confidence"] is None
    assert raw["parameter_covariance"] is None
    assert raw["validity"]["pressure_gpa"][1] == pytest.approx(float(equation(0.6)))


def test_density_conversion_and_error_conventions_are_explicit():
    source = json.loads((DATA / SOURCE).read_text())
    assert source["volume_basis"]["V0_cell_a3_derived_from_printed_density"] == (
        pytest.approx(67.71649780791556)
    )
    assert source["volume_basis"]["V0_cell_a3_rounded"] == V0
    assert source["table2"]["B0_prime_0k_error"] == 0.10
    assert source["table2"]["eos_parameter_error_confidence"] is None
    assert "2 sigma" in source["table2"]["shock_regression_error_confidence"]


def test_table3_preserves_liquid_markers_and_does_not_fill_blank_cells():
    rows = table_rows()
    assert len(rows) == 162
    liquid = [r for r in rows if r["source_phase_annotation"] == "first_liquid_state"]
    assert [(r["volume_ratio"], r["temperature_k"]) for r in liquid] == [
        ("1.00", "2000"),
        ("0.98", "2000"),
        ("0.96", "2500"),
        ("0.94", "2500"),
        ("0.92", "3000"),
        ("0.90", "3000"),
    ]
    source = json.loads((DATA / SOURCE).read_text())
    states = {(float(r["volume_ratio"]), float(r["temperature_k"])) for r in rows}
    assert len(states) == len(rows)
    assert source["table3"]["observations"] is False
    assert len(source["table3"]["omitted_cells"]) == 6
    assert all(
        (cell["volume_ratio"], cell["temperature_k"]) not in states
        for cell in source["table3"]["omitted_cells"]
    )
    document = get_material_document("gold")
    dataset = next(
        d
        for d in document["datasets"]
        if d["identifier"] == "gold_yokoo_2009_table3_isochores"
    )
    assert dataset["kind"] == "source_derived_eos_checkpoints"
    assert dataset["used_by_eos_records"] == [
        *RECORDS,
        "gold_yokoo_2009_pvt_reconstruction",
    ]


def test_fitted_300k_branches_remain_distinct_from_full_model_output():
    audit = reproduce()
    branches = audit["branches"]
    assert branches["gold_yokoo_2009_bm3_300k"][
        "max_abs_fit_minus_full_model_table3_300k_gpa"
    ] == pytest.approx(2.95987378691, abs=1e-9)
    assert branches["gold_yokoo_2009_vinet_300k"][
        "max_abs_fit_minus_full_model_table3_300k_gpa"
    ] == pytest.approx(13.100754896054, abs=1e-9)
    assert audit["cold_checkpoint_diagnostic"]["max_abs_residual_gpa"] < 0.01
    assert audit["cold_checkpoint_diagnostic"]["kind"] == (
        "derived_output_only_not_independent_fit"
    )
    for identifier in RECORDS:
        outcome = ledger_outcome({"identifier": identifier})
        assert outcome["status"] == "not_refittable"
        assert outcome["observations"] == 0
    ledger = json.loads(
        (Path(__file__).parents[1] / "docs/data/primary-eos-refits.json").read_text()
    )
    outcomes = {r["record_identifier"]: r for r in ledger["records"]}
    for identifier in RECORDS:
        assert outcomes[identifier]["status"] == "not_refittable"
        assert outcomes[identifier]["observations"] == 0


def test_au_phonon_law_is_integrated_but_not_a_complete_thermal_eos():
    assert gamma_theta(1) == pytest.approx((2.96, 170.0))
    for ratio in [1, 0.8, 0.6]:
        delta = 1e-5
        derivative = (
            np.log(gamma_theta(ratio * np.exp(delta))[1])
            - np.log(gamma_theta(ratio * np.exp(-delta))[1])
        ) / (2 * delta)
        assert derivative == pytest.approx(-gamma_theta(ratio)[0], rel=1e-8)
    assert phonon_pressure(1, 0) == 0
    endpoint = reproduce()["table3_thermal_term_diagnostic"][-1]
    assert endpoint["volume_ratio"] == "0.60"
    assert endpoint["temperature_k"] == "3000"
    # A rounded-table residual cannot become an exact electronic fit coefficient.
    assert endpoint["table_total_minus_0k_minus_phonon_gpa"] < -0.04


def test_shock_regression_and_simultaneous_volumes_keep_their_source_roles():
    shock = shock_relation_checkpoints()
    assert all(
        s["kind"] == "derived_shock_regression_output_not_observation" for s in shock
    )
    assert shock[0]["pressure_increment_gpa"] == 0
    assert shock[-1]["particle_velocity_km_s"] == 3.5
    assert shock[-1]["shock_velocity_km_s"] == pytest.approx(8.62125)
    pairs = reproduce()["dewaele_simultaneous_volume_pairs"]
    assert len(pairs) == 36
    assert all(p["role"] == "comparison_only_not_yokoo_fit_target" for p in pairs)
    assert pairs[0]["original_revised_ruby_pressure_gpa"] == 2.27
    assert pairs[0]["gold_atomic_volume_a3"] == 16.716
    assert pairs[-1]["gold_csv_row_one_based"] == 37


def test_saved_audit_and_source_hashes_are_reproducible():
    audit = reproduce()
    check_saved(json.loads(OUTPUT.read_text()), audit)
    for filename, digest in audit["input_sha256"].items():
        assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == digest


def test_electronic_source_preserves_original_and_erratum_values():
    audit = reproduce()
    nodes = audit["electronic_source_nodes"]
    assert len(nodes) == 51
    changed = [n for n in nodes if n["gold_erratum_applied"] == "1"]
    assert [int(n["temperature_k"]) for n in changed] == list(range(3100, 4000, 100))
    node = next(n for n in nodes if n["temperature_k"] == "3100")
    assert float(node["gold_original_electronic_pressure_gpa"]) == 0.19
    assert float(node["gold_corrected_electronic_pressure_gpa"]) == 0.13
    assert all(
        float(row["temperature_k"]) <= 3000
        for row in audit["table3_thermal_term_diagnostic"]
    )
    endpoint = audit["table3_thermal_term_diagnostic"][-1]
    assert endpoint["upstream_electronic_pressure_at_source_node_gpa"] == 0.12
    assert audit["max_abs_table_minus_0k_phonon_electronic_gpa"] < 0.22
    source = json.loads((DATA / SOURCE).read_text())
    assert (
        source["electronic_source"]["sha256"]
        == audit["input_sha256"][source["electronic_source"]["resource"]]
    )
