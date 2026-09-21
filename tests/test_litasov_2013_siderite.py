"""Primary-table transcription, source benchmark, and acceptance boundaries."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material, get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.fitting import fit_joint_eos
from scripts.reproduce_litasov_2013_siderite import (
    DATASET,
    RECORD,
    THERMAL_RECORD,
    bm3,
    ledger_outcome,
    observations,
    reproduce,
    thermal_pressure,
)

ROOT = Path(__file__).resolve().parents[1]


def test_complete_table_and_printed_uncertainties():
    document = get_material_document("siderite_fe095mn005")
    dataset = document["datasets"][0]
    resource = ROOT / "peritheos/data" / dataset["resource"]["path"]
    assert (
        hashlib.sha256(resource.read_bytes()).hexdigest()
        == dataset["resource"]["sha256"]
    )
    data = observations()
    assert len(data) == 111
    assert Counter(data["run"]) == {1: 66, 2: 12, 3: 6, 4: 27}
    assert sum(data["temperature_k"] == 300) == 27
    assert data["row"].tolist() == list(range(1, 112))
    assert data["pressure_gpa"].min() == 0
    assert data["pressure_gpa"].max() == 33.01
    assert data["temperature_k"].max() == 1673
    # Independent anchors from each source run, including its uncertainties.
    for row, p, v, sv in [
        (1, 0, 293.49, 0.03),
        (67, 28.09, 246.71, 0.03),
        (79, 33.01, 241.54, 0.04),
        (111, 0, 293.11, 0.02),
    ]:
        value = data[row - 1]
        assert value["pressure_gpa"] == p
        assert value["volume_a3"] == v
        assert value["volume_sigma_a3"] == sv
    assert data[0]["a_sigma_a"] == 0.0002
    assert data[0]["c_sigma_a"] == 0.0017
    assert data[78]["gold_volume_a3"] == 59.34
    assert set(data["source_pdf_page"]) == {14, 15, 16}
    assert (
        max(abs(np.sqrt(3) / 2 * data["a_a"] ** 2 * data["c_a"] - data["volume_a3"]))
        < 0.006
    )
    assert not any(
        c["name"] in {"pressure_sigma_gpa", "temperature_sigma_k"}
        for c in dataset["columns"]
    )
    loaded = get_material("siderite_fe095mn005").get_dataset(DATASET)
    assert loaded is not None


def test_composition_structure_and_source_parameters():
    document = get_material_document("siderite_fe095mn005")
    assert document["formula"] == "Fe0.95Mn0.05CO3"
    assert document["space_group_number"] == 167
    assert document["formula_units_per_cell"] == 6
    contents = {
        s["element"]: int(s["wyckoff"][:-1]) * s["occupancy"]
        for s in document["atom_sites"]
    }
    assert contents == pytest.approx({"Fe": 5.7, "Mn": 0.3, "C": 6, "O": 18})
    lattice = document["lattice"]
    assert np.sqrt(3) / 2 * lattice["a"] ** 2 * lattice["c"] == pytest.approx(
        293.49, abs=0.005
    )
    assert len(document["eos_records"]) == 2
    record = document["eos_records"][0]
    assert record["eos"]["parameters"] == {"V0": 293.4, "K0": 120, "K0_prime": 3.57}
    assert record["parameter_errors"] == {"V0": 0.1, "K0": 1, "K0_prime": 0.09}
    assert record["fixed_parameters"] == ["V0"]
    assert "thermal" not in record
    assert record["experimental_temperature_range_k"] == [300, 300]
    assert record["pressure_calibration"]["status"] == "partially_resolved"
    assert record["pressure_calibration"]["recalculation"]["status"] == "not_possible"
    assert all(
        "reference_eos_record" not in method
        for method in record["pressure_calibration"]["methods"]
    )


def test_highest_pressure_primary_benchmark_and_supported_domain():
    record = get_eos_record(RECORD)
    assert abs(record.pressure(241.54, 300) - 33.01) < 0.1
    volumes = np.linspace(241.54, 293.4, 30)
    pressures = record.pressure(volumes, 300)
    np.testing.assert_allclose(pressures, bm3(volumes, 293.4, 120, 3.57), atol=1e-11)
    np.testing.assert_allclose(record.volume(pressures, 300), volumes, atol=1e-8)
    with pytest.raises(ValueError, match="isothermal 300 K"):
        record.volume(10, 1000, check_validity=True)
    with pytest.raises(ValueError, match="outside the published"):
        record.volume(45, 300, check_validity=True)


def test_refit_and_unresolved_source_differences():
    result = reproduce()
    rt = result["rt_refit"]
    assert rt["parameters"] == pytest.approx(
        {"K0": 120.224213, "K0_prime": 3.572214}, rel=1e-7
    )
    assert rt["rmse_gpa"] == pytest.approx(0.22730866)
    assert result["thermal_joint_refit"]["solver_success"]
    assert result["thermal_joint_refit"]["observations"] == 111
    assert result["thermal_joint_refit"]["parameters"]["alpha0"] == pytest.approx(
        3.58784e-5, rel=1e-5
    )
    assert result["thermal_staged_published_cold_sensitivity"]["observations"] == 84
    assert result["thermal_published_table"]["rmse_gpa"] == pytest.approx(0.232563185)
    assert result["thermal_published_abstract"]["rmse_gpa"] == pytest.approx(
        0.299370193
    )
    assert result["gold_calibration_diagnostic"]["max_abs_residual_gpa"] > 0.7
    assert not result["gold_calibration_diagnostic"]["recalibration_applied"]
    document = get_material_document("siderite_fe095mn005")
    assert ledger_outcome(document["eos_records"][0])["status"] == "similar"
    audit = json.loads(
        (ROOT / "docs/data/litasov-2013-siderite-source-audit.json").read_text()
    )
    assert audit["accepted_record_identifiers"] == [RECORD, THERMAL_RECORD]
    assert (
        audit["candidates"][1]["disposition"] == "withheld_source_coefficient_conflict"
    )
    assert audit["candidates"][1]["alpha0_abstract"] == 3.57e-5
    assert audit["candidates"][1]["alpha0_table1"] == 3.77e-5


def test_thermal_source_variants_are_evidence_only_and_refit_is_opt_in():
    material = get_material("siderite_fe095mn005")
    document = material.to_eosmat()
    variants = document["source"]["thermal_parameterizations"]
    assert [v["parameters"]["alpha0"] for v in variants] == [3.57e-5, 3.77e-5]
    for variant in variants:
        assert variant["executable"] is False
        with pytest.raises(KeyError):
            get_eos_record(variant["identifier"])
    assert [r.identifier for r in material.eos_records if r.is_default] == [RECORD]
    raw = next(r for r in document["eos_records"] if r["identifier"] == THERMAL_RECORD)
    assert raw["record_kind"] == "refit"
    assert raw["derived_from_record"] == RECORD
    assert raw["fit_provenance"]["selection"]["included_rows"] == 111
    assert raw["fit_provenance"]["selection"]["excluded_source_rows_1_based"] == []
    assert raw["fixed_parameters"] == ["V0"]
    assert raw["thermal"]["fixed_parameters"] == ["Tr"]
    assert raw["parameter_errors"]["V0"] is None
    assert raw["eos"]["parameters"]["K0"] == pytest.approx(121.90288)
    assert raw["eos"]["parameters"]["K0_prime"] == pytest.approx(3.427642)
    assert raw["thermal"]["reference_volume_law"] == "integrated_expansivity"
    assert raw["thermal"]["thermal_expansion_law"] == "linear_temperature"
    assert ledger_outcome(raw)["status"] == "parity"
    reloaded = Material.from_eosmat(document).get_eos_record(THERMAL_RECORD)
    assert reloaded.pressure(247.82, 1673) == pytest.approx(32.3574785491)


def test_thermal_mapping_covariance_and_inversion_on_primary_states():
    record = get_eos_record(THERMAL_RECORD)
    data = observations()
    result = reproduce()["thermal_joint_refit"]
    stored = get_material_document("siderite_fe095mn005")
    raw = next(r for r in stored["eos_records"] if r["identifier"] == THERMAL_RECORD)
    params = {**raw["eos"]["parameters"], **raw["thermal"]["parameters"]}
    scaled = np.array(
        [params[n] for n in ["K0", "K0_prime", "alpha0", "alpha1", "dK_dT"]]
    ) / np.array([1, 1, 1e-5, 1e-8, 1])
    calculated = record.pressure(data["volume_a3"], data["temperature_k"])
    reference = thermal_pressure(scaled, data["volume_a3"], data["temperature_k"])
    # Equation mapping uses identical stored coefficients. A fresh independent
    # optimizer run is checked separately at bounded pressure precision.
    np.testing.assert_allclose(calculated, reference, atol=2e-11, rtol=0)
    refitted_scaled = np.array(
        [
            result["parameters"][n]
            for n in ["K0", "K0_prime", "alpha0", "alpha1", "dK_dT"]
        ]
    ) / np.array([1, 1, 1e-5, 1e-8, 1])
    np.testing.assert_allclose(
        calculated,
        thermal_pressure(refitted_scaled, data["volume_a3"], data["temperature_k"]),
        atol=5e-6,
        rtol=0,
    )
    assert np.sqrt(np.mean((calculated - data["pressure_gpa"]) ** 2)) == pytest.approx(
        0.221211666
    )
    np.testing.assert_allclose(
        record.volume(calculated, data["temperature_k"]), data["volume_a3"], atol=1e-8
    )
    assert record.pressure(293.4, 300) == pytest.approx(0, abs=1e-12)
    # High-T source state (row 84) compared with the conditional fit error and
    # the printed pressure bound, not claimed as independent holdout validation.
    prediction = record.pressure_with_uncertainty(247.82, 1673)
    assert prediction.standard_error == pytest.approx(0.30756958, rel=1e-5)
    assert abs(record.pressure(247.82, 1673) - 32.10) < prediction.standard_error + 0.1
    covariance = np.array(record.parameter_covariance)
    np.testing.assert_allclose(
        covariance,
        result["estimated_covariance_free_parameters"],
        rtol=2e-4,
        atol=1e-20,
    )
    assert record.covariance_parameters == (
        "rt_eos.K0",
        "rt_eos.K0_prime",
        "alpha0",
        "alpha1",
        "dK_dT",
    )
    np.testing.assert_allclose(
        np.sqrt(covariance.diagonal()),
        result["estimated_standard_errors"],
        rtol=2e-4,
        atol=1e-20,
    )
    record.volume(32.1, 1673, check_validity=True)
    for pressure, temperature in [(45, 300), (20, 1800)]:
        with pytest.raises(ValueError, match="outside the published"):
            record.volume(pressure, temperature, check_validity=True)


def test_refit_parity_requires_pressure_agreement_as_well_as_coefficients(monkeypatch):
    document = get_material_document("siderite_fe095mn005")
    record = next(
        r for r in document["eos_records"] if r["identifier"] == THERMAL_RECORD
    )
    report = reproduce()
    stored = {**record["eos"]["parameters"], **record["thermal"]["parameters"]}
    parameters = report["thermal_joint_refit"]["parameters"]
    for name in parameters:
        parameters[name] = stored[name]
    # This fits inside alpha1's coefficient tolerance but produces a detectable
    # pressure change, which must prevent a parity classification.
    parameters["alpha1"] += 1.5e-13
    monkeypatch.setattr(
        "scripts.reproduce_litasov_2013_siderite.reproduce", lambda: report
    )
    outcome = ledger_outcome(record)
    assert all(p["similar"] for p in outcome["parameters"])
    assert outcome["status"] == "parity_not_achieved"


def test_native_joint_fit_recovers_the_independent_reference_result():
    data = observations()
    result = fit_joint_eos(
        ThermalReferenceStateEOS,
        BM3,
        data["volume_a3"],
        data["temperature_k"],
        data["pressure_gpa"],
        initial={
            "rt_eos.K0": 120.0,
            "rt_eos.K0_prime": 3.57,
            "alpha0": 3.77e-5,
            "alpha1": 6e-10,
            "dK_dT": -0.015,
        },
        fixed={"rt_eos.V0": 293.4, "Tr": 300.0},
        configuration={
            "thermal_expansion_law": "linear_temperature",
            "reference_volume_law": "integrated_expansivity",
        },
        max_nfev=5000,
    )
    assert result.success
    expected = reproduce()["thermal_joint_refit"]["parameters"]
    for name, value in expected.items():
        key = f"rt_eos.{name}" if name in {"K0", "K0_prime"} else name
        assert result.parameters[key] == pytest.approx(value, rel=5e-4)
