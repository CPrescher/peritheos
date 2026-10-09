"""Separate Matsui's verified lattice equation, observations and model output."""

import hashlib
import json

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import DebyeTabulatedThermalPressure
from peritheos.eosmat import validate_eosmat_document
from peritheos.errors import EosmatError, EosValidationError, UnsupportedOperationError
from scripts.audit_matsui_2009_platinum import (
    DATA,
    OUTPUT,
    RECORD,
    THERMAL_RECORD,
    V0,
    audit,
    electronic_increment,
    rows,
    source_components,
)


def test_independent_lattice_equation_and_cell_normalization():
    result = audit()
    assert result["independent_cold_max_difference_gpa"] < 1e-9
    assert result["independent_phonon_max_difference_gpa"] < 1e-9
    cold, phonon, theta, gamma = source_components(V0, 3000)
    assert cold == pytest.approx(0, abs=1e-12)
    assert theta == pytest.approx(230)
    assert gamma == pytest.approx(2.7)
    assert phonon == pytest.approx(19.9448402173, abs=1e-8)
    assert cold + phonon + 1.60 == pytest.approx(21.54, abs=0.005)
    assert source_components(np.array([1, 0.8]) * V0, 300)[1] == pytest.approx([0, 0])
    v, step = 0.8 * V0, 1e-5
    theta_plus = source_components(v * np.exp(step), 1000)[2]
    theta_minus = source_components(v * np.exp(-step), 1000)[2]
    derivative = -(np.log(theta_plus) - np.log(theta_minus)) / (2 * step)
    assert derivative == pytest.approx(source_components(v, 1000)[3], rel=1e-9)


def test_rounded_outputs_constrain_electronic_pressure_without_fitting_function():
    result = audit()
    for item in result["known_electronic_table3_checks"].values():
        assert item["max_abs_pressure_difference_gpa"] < 0.008
    for item in result["table3"]["inferred_electronic"]:
        assert item["volume_spread_gpa"] < 0.01
        assert item["common_rounding_interval_gpa"] is not None
    inferred = result["table3"]["inferred_electronic"]
    assert inferred[3]["common_rounding_interval_gpa"][0] > 0.80
    assert inferred[4]["common_rounding_interval_gpa"][0] > 1.59
    examples = result["electronic_examples_reference_subtracted_gpa"]
    # A single quadratic coefficient cannot explain all examples even allowing
    # their +/-0.005 GPa last-digit rounding uncertainty.
    bounds = [
        ((p - 0.005) / (int(t) ** 2 - 300**2), (p + 0.005) / (int(t) ** 2 - 300**2))
        for t, p in examples.items()
    ]
    assert max(lo for lo, hi in bounds) > min(hi for lo, hi in bounds)
    assert (
        result["continuous_thermal_record_status"]
        == "registered_with_explicit_linear_interpolation"
    )
    assert result["author_fit_reproduction_status"] == "not_reproduced"


def test_static_and_shock_observations_keep_their_own_pressure_provenance():
    observed = rows("table1")
    assert len(observed) == 12
    for row in observed:
        assert float(row["observed_pressure_gpa"]) - float(
            row["calculated_pressure_gpa"]
        ) == pytest.approx(
            float(row["printed_observed_minus_calculated_gpa"]), abs=1e-12
        )
        assert float(row["observed_pressure_error_gpa"]) > 0
        assert float(row["volume_ratio_error"]) > 0
    result = audit()
    assert result["table1_printed_residual_rmse_gpa"] == pytest.approx(0.092285788)
    assert {r["shot"] for r in result["holmes_subset"]} == {"TaPt5", "TaPt6", "TaPt8"}
    for row in result["holmes_subset"]:
        assert row["momentum_pressure_gpa"] == pytest.approx(
            row["published_pressure_gpa"], abs=0.1
        )
        assert row["mass_conservation_volume_ratio"] == pytest.approx(
            row["published_density_volume_ratio"], abs=1e-4
        )
    assert len(rows("expansion")) == 32
    assert (
        result["arblaster_expansion"][14]["volume_ratio_from_rounded_lattice_parameter"]
        == 1
    )


def test_catalog_retains_cold_branch_and_evidence_is_not_claimed_as_fit_input():
    document = get_material_document("platinum")
    raw = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert "thermal" not in raw
    assert "fit_datasets" not in raw
    assert raw["parameter_errors"] == {"V0": None, "K0": None, "K0_prime": None}
    record = get_eos_record(RECORD)
    restored = Material.from_eosmat(document, record_identifiers=[RECORD]).eos_records[
        0
    ]
    volumes = V0 * np.array([1, 0.9, 0.7])
    expected = source_components(volumes, 300)[0]
    assert restored.pressure(volumes) == pytest.approx(expected, abs=1e-9)
    assert record.volume(expected) == pytest.approx(volumes, rel=1e-8)
    evidence = {d["identifier"]: d for d in document["datasets"]}
    table1 = evidence["platinum_matsui_2009_table1_pvt"]
    assert "MgO" in table1["pressure_scale"]["observed_pressure_gpa"]
    for identifier in (
        "platinum_matsui_2009_table1_pvt",
        "platinum_matsui_2009_table3_grid",
        "platinum_arblaster_1997_table2_expansion",
    ):
        dataset = evidence[identifier]
        path = DATA.parent / dataset["resource"]["path"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )


def test_saved_audit_matches_independent_replay():
    result, saved = audit(), json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert result["input_sha256"] == saved["input_sha256"]
    assert result["parameters"] == saved["parameters"]
    assert np.array(result["table3"]["phonon_increment_gpa"]) == pytest.approx(
        np.array(saved["table3"]["phonon_increment_gpa"]), abs=1e-9
    )
    for filename, expected in saved["input_sha256"].items():
        assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == expected


def test_recovered_pressure_nodes_and_au_only_erratum():
    source = rows("electronic")
    assert len(source) == 51
    assert [int(r["temperature_k"]) for r in source] == list(range(0, 5001, 100))
    corrected = [r for r in source if r["gold_erratum_applied"] == "1"]
    assert [int(r["temperature_k"]) for r in corrected] == list(range(3100, 3901, 100))
    assert [float(r["gold_corrected_electronic_pressure_gpa"]) for r in corrected] == [
        0.13,
        0.14,
        0.15,
        0.17,
        0.18,
        0.19,
        0.21,
        0.22,
        0.24,
    ]
    for r in source:
        if r["gold_erratum_applied"] == "0":
            assert (
                r["gold_original_electronic_pressure_gpa"]
                == r["gold_corrected_electronic_pressure_gpa"]
            )
    assert float(source[21]["platinum_electronic_pressure_gpa"]) == 0.95
    eos = get_eos_record(THERMAL_RECORD).eos
    t = np.arange(100, 5001, 100)
    assert eos.electronic_pressure_increment(t) == pytest.approx(
        electronic_increment(t)
    )
    assert eos.electronic_pressure_increment([300, 1000, 3000, 5000]) == pytest.approx(
        [0, 0.21, 1.6, 3.78]
    )
    # Sub-grid pressure is qualified interpolation, not a published source node.
    assert eos.electronic_pressure_increment(2050) == pytest.approx(0.86)
    assert "electronic_pressure_gpa" not in eos.parameter_values()
    assert (
        eos.with_parameters(q=1.2).configuration_values() == eos.configuration_values()
    )


def test_complete_pressure_equation_grid_and_bounded_inversions():
    record = get_eos_record(THERMAL_RECORD)
    t = np.array([300, 500, 1000, 2000, 3000])
    v = V0 * np.array([1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7])[:, None]
    cold, phonon, _, _ = source_components(v, t)
    oracle = cold + phonon + electronic_increment(t)
    assert record.pressure(v, t) == pytest.approx(oracle, abs=1e-9)
    result = audit()
    assert result["independent_full_max_difference_gpa"] < 1e-9
    assert result["table3"]["max_abs_pressure_difference_gpa"] < 0.0071
    assert result["table1_max_abs_calculated_reconstruction_difference_gpa"] < 0.0086
    assert record.volume(oracle, t) == pytest.approx(
        np.broadcast_to(v, oracle.shape), rel=1e-8
    )
    eos = record.eos
    assert eos.temperature(oracle, v * record.volume_scale) == pytest.approx(
        np.broadcast_to(t, oracle.shape), abs=1e-6
    )
    for temperature in [300, 2050, 5000]:
        p = record.pressure(0.85 * V0, temperature)
        assert eos.temperature(p, 0.85 * V0 * record.volume_scale) == pytest.approx(
            temperature, abs=1e-6
        )
    for temperature in [0, -1, 5001, np.nan]:
        with pytest.raises(EosValidationError):
            eos.pressure(V0 * record.volume_scale, temperature)
    with pytest.raises(EosValidationError, match="root"):
        eos.temperature(record.pressure(V0, 5000) + 1, V0 * record.volume_scale)
    with pytest.raises(UnsupportedOperationError):
        eos.heat_capacity_v(V0 * record.volume_scale, 3000)
    va = 0.85 * V0 * record.volume_scale
    vh = eos.volume_with_dac_confinement(eos.rt_eos.pressure(va), 2050, f_dac=0.3)
    assert eos.temperature_from_volumes(va, vh, f_dac=0.3) == pytest.approx(
        2050, abs=1e-5
    )


def test_tabulated_pressure_round_trip_preserves_configuration_and_provenance():
    document = get_material_document("platinum")
    material = Material.from_eosmat(document, record_identifiers=[THERMAL_RECORD])
    serialized = json.loads(json.dumps(material.to_eosmat()))
    validate_eosmat_document(serialized)
    restored = Material.from_eosmat(serialized).eos_records[0]
    original = next(
        r for r in document["eos_records"] if r["identifier"] == THERMAL_RECORD
    )
    output = serialized["eos_records"][0]
    assert output["thermal"]["configuration"] == original["thermal"]["configuration"]
    assert output["determination_method"] == "hybrid"
    assert "fit_datasets" not in output
    assert output["parameter_errors"] == original["parameter_errors"]
    assert restored.pressure(V0 * 0.9, 2050) == pytest.approx(
        material.eos_records[0].pressure(V0 * 0.9, 2050), abs=1e-9
    )
    output["thermal"]["configuration"]["electronic_temperature_k"][1] = 0
    with pytest.raises(EosmatError):
        validate_eosmat_document(serialized)


@pytest.mark.parametrize(
    "temperatures, pressures",
    [
        ([0, 100, 100], [0, 0.01, 0.02]),
        ([0, 100], [0]),
        ([0, 100], [0, np.nan]),
        ([0, 100], [0, -0.01]),
    ],
)
def test_invalid_electronic_configuration_is_rejected(temperatures, pressures):
    with pytest.raises(EosValidationError):
        DebyeTabulatedThermalPressure(
            Vinet(1, 273, 5.2),
            50,
            230,
            2.7,
            1.1,
            1,
            electronic_temperature_k=temperatures,
            electronic_pressure_gpa=pressures,
            interpolation="linear",
        )
