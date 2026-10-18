"""Primary observations, composition-safe indexing and bundled EOS acceptance."""

import copy
import hashlib
import json
import math
import subprocess
import sys

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_material, get_material_document
from peritheos.errors import MaterialError
from scripts.eosmat_io import read_material
from scripts.reproduce_phase_e_crichton_ross_2000 import (
    CARD,
    ERRORS,
    PARAMETERS,
    PUBLISHED,
    RECORD,
    ROOT,
    SOURCE,
    TABLE,
    diagnostic_fit,
    load_observations,
    parse_reported,
    reproduce,
    source_dp_dv,
    source_pressure,
)


def document():
    return read_material(CARD)


def test_bundled_identity_published_coefficients_and_reference_state():
    doc = document()
    Draft202012Validator(eosmat_schema()).validate(doc)
    assert doc["identifier"] == "phase_e_fe0072"
    assert doc["formula"] == "Mg1.96Fe0.072Si1.04H3.7O6"
    assert doc["sample_composition"]["reported_errors"] == {
        "Mg": 0.07,
        "Fe": 0.005,
        "Si": 0.05,
        "H": 0.8,
        "O": None,
    }
    assert doc["formula_units_per_cell"] == 1
    assert doc["space_group_number"] == 166
    assert doc["crystal_system"] == "trigonal"
    assert doc["symmetry"] == "HEXAGONAL"  # Conventional axes, not crystal system.
    assert len(doc["eos_records"]) == 1  # Axial fits are not physical volume EOSs.
    record = doc["eos_records"][0]
    assert record["eos"]["parameters"] == dict(zip(PARAMETERS, PUBLISHED))
    assert record["parameter_errors"] == dict(zip(PARAMETERS, ERRORS))
    assert record["fixed_parameters"] == []
    assert record["parameter_covariance"] is None
    assert record["record_kind"] == "published"
    assert record["temperature_ref"] == 298.0
    assert "experimental_temperature_range_k" not in record
    assert record["temperature_ref_provenance"]["kind"] == "assumed_room_temperature"
    material = Material.from_eosmat(doc)
    assert get_material_document("phase_e_fe0072") == doc
    assert get_material("phase_e_fe0072").get_eos_record(RECORD).is_default
    eos = material.get_eos_record(RECORD)
    assert eos.pressure(105.627) == pytest.approx(0, abs=1e-12)
    assert eos.eos.bulk_modulus(105.627) == pytest.approx(92.9, abs=1e-10)
    normalized = material.to_eosmat()
    restored = Material.from_eosmat(normalized).to_eosmat()
    assert restored == normalized
    for key in ["sample_composition", "datasets", "source", "peaks"]:
        assert restored[key] == doc[key]


def test_complete_primary_table_null_errors_and_independent_csv():
    rows = load_observations()
    assert rows.shape == (11, 8)
    assert np.isnan(rows[0, 1])
    assert rows[0, [0, 2, 3, 4, 5, 6, 7]].tolist() == [
        0.0001,
        2.9655,
        0.0002,
        13.8703,
        0.0005,
        105.633,
        0.016,
    ]
    np.testing.assert_allclose(
        rows[-1], [6.671, 0.008, 2.9090, 0.0003, 13.5972, 0.0006, 99.648, 0.018]
    )
    assert rows[:, 0].tolist() == [
        0.0001,
        0.296,
        1.128,
        1.987,
        2.973,
        3.944,
        4.443,
        5.301,
        5.539,
        6.139,
        6.671,
    ]
    assert rows[:, 6].tolist() == [
        105.633,
        105.294,
        104.401,
        103.529,
        102.623,
        101.794,
        101.367,
        100.665,
        100.492,
        100.025,
        99.648,
    ]
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert (
        hashlib.sha256(TABLE.read_bytes()).hexdigest()
        == source["transcription"]["sha256"]
    )
    doc = document()
    material = Material.from_eosmat(doc)
    dataset = doc["datasets"][0]
    assert dataset["rows"][0][1] is None
    assert (
        material.get_dataset(dataset["identifier"]).values("volume_a3").tolist()
        == rows[:, 6].tolist()
    )
    assert dataset["used_by_eos_records"] == [RECORD]
    assert parse_reported("105.633(16)") == [105.633, 0.016]
    assert parse_reported("0.0001") == [0.0001, None]


def test_source_equation_residual_and_chi_checkpoints():
    result = reproduce()
    check = result["source_equation_validation"]
    assert check["max_abs_pressure_residual_gpa"] == pytest.approx(0.0236970945293)
    assert check["matches_source_max_residual_at_printed_precision"]
    # The publication independently reports maximum |dP|=0.02 and chi_w²=0.4.
    assert check["max_abs_pressure_residual_gpa"] == pytest.approx(0.02, abs=0.005)
    assert check["conditional_effective_variance_chi2_per_dof"] == pytest.approx(
        0.4, abs=0.05
    )
    eos = Material.from_eosmat(document()).get_eos_record(RECORD)
    for volume, pressure in [(104.0, 1.52587797937949), (100.0, 6.19649221371672)]:
        assert source_pressure(volume) == pytest.approx(pressure, abs=1e-12)
        assert eos.pressure(volume) == pytest.approx(pressure, abs=1e-11)
    volumes = np.linspace(99.648, 105.627, 20)
    finite_difference = (
        source_pressure(volumes + 1e-4) - source_pressure(volumes - 1e-4)
    ) / 2e-4
    np.testing.assert_allclose(
        source_dp_dv(volumes, PUBLISHED), finite_difference, rtol=1e-8
    )
    pressures = np.array([0.0, 0.296, 2.973, 6.671])
    np.testing.assert_allclose(
        eos.pressure(eos.volume(pressures)), pressures, atol=1e-10
    )
    eos.volume(6.671, check_validity=True)
    with pytest.raises(MaterialError):
        eos.volume(6.8, check_validity=True)


def test_conditional_refits_do_not_replace_published_parameters():
    observations = load_observations()
    before = copy.deepcopy(document())
    for mode in [
        "unweighted_pressure",
        "iterated_effective_variance",
        "direct_effective_variance",
    ]:
        fit = diagnostic_fit(observations, mode)
        assert fit["all_within_published_errors"]
        assert fit["fixed_parameters"] == []
        assert fit["covariance"] is None
    fit = diagnostic_fit(observations, "iterated_effective_variance")
    np.testing.assert_allclose(
        list(fit["parameters"].values()), [105.6275239, 92.81451, 7.277217], rtol=2e-6
    )
    assert fit["effective_variance_chi2_per_dof"] == pytest.approx(0.33974793, abs=1e-7)
    without = diagnostic_fit(
        observations, "iterated_effective_variance", include_ambient=False
    )
    assert without["observations"] == 10
    assert not without["all_within_published_errors"]
    assert document() == before
    parity = reproduce()["fit_recovery"]
    assert parity["status"] == "parity_within_reported_uncertainty"
    assert parity["max_coefficient_shift_over_published_error"] < 0.13
    assert parity["max_curve_difference_gpa_on_2000_point_volume_grid"] < 0.01
    assert parity["all_curve_differences_within_combined_observation_errors"]
    assert (
        document()["eos_records"][0]["scientific_validation"]["reproduction_scope"][
            "fit_recovery"
        ]
        == parity["status"]
    )
    workflow = reproduce()["exact_author_workflow_reproduction"]
    assert workflow["provenance_only"]
    assert workflow["does_not_negate_fit_parity"]
    assert (
        reproduce()["exact_author_workflow_reproduction"]["status"] == "not_reproduced"
    )


def test_diffraction_fallback_uses_target_cell_and_separate_provenance():
    doc = document()
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert doc["atom_sites"] == []
    assert len(doc["peaks"]) == 8
    assert [peak[:3] for peak in doc["peaks"]] == source["diffraction_fallback"]["hkl"]
    a, c = doc["lattice"]["a"], doc["lattice"]["c"]
    for h, k, ell, d, marker in doc["peaks"]:
        assert (-h + k + ell) % 3 == 0
        assert d == pytest.approx(
            1
            / math.sqrt(
                4 * (h * h + h * k + k * k) / (3 * a * a) + ell * ell / (c * c)
            ),
            abs=1e-10,
        )
        assert marker == 1.0
    assert doc["peaks"][0][3] == pytest.approx(c / 3)
    assert doc["peaks"][5][3] == pytest.approx(a / 2)
    structure = doc["source"]["structure"]
    assert structure["reference"]["doi"] == "10.1007/BF00202972"
    assert structure["fallback_reference"]["doi"] == "10.2138/am-2000-5-616"
    assert structure["cell_reference"]["doi"] == "10.1180/002646100549427"
    assert reproduce()["cell_consistency"]["all_within_printed_rounding_bounds"]
    # Original structure inconsistencies are retained, not repaired by occupancy edits.
    checks = reproduce()["structure_site_audit"]
    assert checks[0]["table1_site_totals"]["Mg"] == pytest.approx(2.031)
    assert checks[1]["table1_site_totals"]["Mg"] == pytest.approx(2.169)
    assert (
        source["structure_samples"]["samples"][0]["reported_formula"]
        == "Mg2.08Si1.16H3.20O6"
    )


def test_calibration_axial_and_thermal_scope():
    doc = document()
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    calibration = doc["eos_records"][0]["pressure_calibration"]
    assert calibration["methods"][0]["reference"]["doi"] == "10.1107/S0021889897000861"
    assert (
        calibration["methods"][0]["reference_eos_record"]
        == "alpha_quartz_angel_1997_bm3_1"
    )
    assert len(calibration["methods"]) == 1
    assert calibration["auxiliary_checks"][0]["kind"] == "ruby_fluorescence"
    assert calibration["auxiliary_checks"][0]["status"] == "unresolved"
    assert calibration["recalculation"]["status"] == "missing_calibrant_observations"
    assert source["axial_parameterizations"]["not_volume_eos"]
    assert source["axial_parameterizations"]["a"]["K_effective_gpa"] == 94.8
    assert source["thermal_assessment"]["status"] == "deferred_distinct_composition"
    assert "thermal" not in doc["eos_records"][0]
    assert source["rights"]["source_license"] == "unspecified"


def test_author_thesis_version_is_preserved_without_merging_tables():
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    thesis = source["author_thesis_crosscheck"]
    assert len(thesis["raw_rows"]) == 11
    assert len(thesis["differences_from_journal"]) == 4
    assert thesis["raw_rows"][0][1] == "2.9654(2)"
    assert thesis["raw_rows"][1][0] == "0.296(4)"
    assert thesis["raw_rows"][5][1] == "2.9295(2)"
    assert thesis["raw_rows"][-1][1] == "2.9091(3)"
    journal = document()["datasets"][0]["rows"]
    assert journal[0][2] == 2.9655
    assert journal[1][1] == 0.003
    assert journal[5][2] == 2.9296
    assert journal[-1][2] == 2.9090
    assert thesis["bm2_comparison"]["fixed_parameters"] == ["K0_prime"]
    assert thesis["bm2_comparison"]["chi_w_squared"] == 7.3


def test_reproduction_report_is_current():
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/reproduce_phase_e_crichton_ross_2000.py"),
            "--check",
        ],
        cwd=ROOT,
        check=True,
    )
