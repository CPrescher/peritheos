import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, get_eos_record, get_material, get_material_document
from peritheos.eosmat import eosmat_schema
from peritheos.errors import EosValidationError, MaterialError
from scripts.diopside_reproduce import (
    DATA,
    bm3,
    metric_differences,
    normalized_csv,
    reproduce,
    thermal_pressure,
)

STATIC = "diopside_thompson_2008_bm3"
THERMAL = "diopside_zhao_1998_thermal_bm3"
ROOT = Path(__file__).resolve().parents[1]


def test_diopside_schema_catalog_and_lossless_round_trip():
    document = get_material_document("diopside")
    Draft202012Validator(eosmat_schema()).validate(document)
    material = get_material("diopside")
    assert material.default_record().identifier == STATIC
    assert get_material("CaMgSi2O6").identifier == "diopside"
    assert get_eos_record(THERMAL).is_thermal
    restored = Material.from_eosmat(material.to_eosmat())
    assert restored.to_eosmat() == material.to_eosmat()
    for identifier in (STATIC, THERMAL):
        a = material.get_eos_record(identifier)
        b = restored.get_eos_record(identifier)
        assert b.pressure(425.0, 300.0) == pytest.approx(a.pressure(425.0, 300.0))
    assert document["eos_records"][0]["parameter_errors"]["K0_prime"] == 0.3
    assert document["eos_records"][1]["parameter_errors"]["V0"] is None
    assert document["eos_records"][1]["thermal"]["parameter_errors"]["dK_dT"] == 0.0004
    assert (
        document["eos_records"][1]["scientific_validation"]["numerical_reproduction"][
            "author_fit_reproduced"
        ]
        is False
    )


def test_diopside_observed_structure_cell_contents_and_provenance():
    document = get_material_document("diopside")
    assert document["formula_units_per_cell"] == 4
    assert document["space_group_number"] == 15
    assert document["space_group"] == "C2/c"
    contents = Counter()
    for site in document["atom_sites"]:
        contents[site["element"]] += int(site["wyckoff"][:-1]) * site["occupancy"]
    assert contents == {"Ca": 4, "Mg": 4, "Si": 8, "O": 24}
    assert document["atom_sites"][0]["y"] == 0.908083
    assert document["atom_sites"][1]["y"] == 0.301478
    assert document["atom_sites"][2]["x"] == 0.286187
    lattice = document["lattice"]
    volume = (
        lattice["a"] * lattice["b"] * lattice["c"] * np.sin(np.deg2rad(lattice["beta"]))
    )
    assert volume == pytest.approx(438.63, abs=0.003)
    source = document["source"]
    assert source["structure_reference"]["doi"] == "10.2138/am.2008.2684"
    resource = ROOT / "peritheos/data" / source["structure_dataset"]
    assert (
        hashlib.sha256(resource.read_bytes()).hexdigest()
        == source["structure_dataset_sha256"]
    )
    with resource.open(newline="") as stream:
        sites = list(csv.DictReader(stream))
    assert len(sites) == 6
    assert sites[0]["y_error"] == "0.000097"
    assert sites[2]["z_error"] == "0.000099"


def test_diopside_original_observations_errors_and_inconsistencies_are_retained():
    document = get_material_document("diopside")
    expected_counts = {
        "diopside_thompson_2008_table1": 14,
        "diopside_zhao_1998_table1": 32,
        "diopside_zhao_1998_table1_corrected_pvt": 32,
    }
    assert len(document["datasets"]) == len(expected_counts)
    for dataset in document["datasets"]:
        count = expected_counts[dataset["identifier"]]
        resource = ROOT / "peritheos/data" / dataset["resource"]["path"]
        assert (
            hashlib.sha256(resource.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
        if dataset["identifier"] != "diopside_zhao_1998_table1_corrected_pvt":
            assert normalized_csv(resource.stem).encode() == resource.read_bytes()
        rows = list(csv.DictReader(resource.read_text().splitlines()))
        assert len(rows) == count
        assert list(rows[0]) == [c["name"] for c in dataset["columns"]]
    rows = list(
        csv.DictReader(
            (DATA / "diopside-zhao-1998-table1.csv").read_text().splitlines()
        )
    )
    assert rows[0]["cell_volume_a3"] == rows[-1]["cell_volume_a3"] == "438.67"
    assert rows[14]["pressure_gpa"] == "6.95"
    assert rows[14]["temperature_k"] == "1072"
    assert rows[3]["lattice_a_error_angstrom"] == "0.007"
    assert all(r["pressure_error_gpa"] == r["temperature_error_k"] == "" for r in rows)
    assert max(float(r["temperature_k"]) for r in rows) == 1280
    assert max(float(r["pressure_gpa"]) for r in rows) == 8.18
    assert metric_differences(rows)[0] == pytest.approx(4.905287185142356)
    assert metric_differences(rows)[15] == pytest.approx(0.6645623390308515)


def test_diopside_source_checkpoints_and_independent_equations():
    static = get_eos_record(STATIC)
    thermal = get_eos_record(THERMAL)
    # Off-reference observations, tolerances based on printed source P precision.
    assert static.pressure(407.16, 300.0) == pytest.approx(10.16, abs=0.04)
    assert thermal.pressure(421.43, 1272.0) == pytest.approx(8.18, abs=0.12)
    assert thermal.pressure(426.65, 1280.0) == pytest.approx(6.88, abs=0.12)
    # Page 27 independently reports mean alpha over 300-1300 K as 3.82(8)e-5.
    alpha = thermal.eos.alpha0 + thermal.eos.alpha1 * 800.0
    assert alpha == pytest.approx(3.82e-5, abs=0.08e-5)
    volumes = np.array([416.8, 421.43, 426.65, 444.61])
    temperatures = np.array([308.0, 1272.0, 1280.0, 1272.0])
    assert thermal.pressure(volumes, temperatures) == pytest.approx(
        thermal_pressure(volumes, temperatures), abs=2e-12
    )
    assert static.pressure(volumes, 300.0) == pytest.approx(
        bm3(volumes, 438.66, 118, 3.8), abs=2e-12
    )
    assert thermal.eos._native is not None
    assert static.eos._native is not None


def test_diopside_reference_temperature_law_inversion_and_modulus():
    record = get_eos_record(THERMAL)
    model = record.eos
    assert model.pressure(438.67, 300) == pytest.approx(0, abs=1e-12)
    assert model.bulk_modulus(438.67, 300) == pytest.approx(109.1)
    for temperature in (300.0, 800.0, 1280.0):
        dt = temperature - 300.0
        v0_t = 438.67 * np.exp(2.32e-5 * dt + 0.5 * 1.88e-8 * (temperature**2 - 300**2))
        assert model.pressure(v0_t, temperature) == pytest.approx(0, abs=1e-12)
        assert model.bulk_modulus(v0_t, temperature) == pytest.approx(
            109.1 - 0.0205 * dt
        )
        for pressure in (0.0, 4.0, 8.18):
            volume = record.volume(pressure, temperature, check_validity=True)
            assert record.pressure(
                volume, temperature, check_validity=(pressure == 4.0)
            ) == pytest.approx(pressure, abs=1e-9)
            step = 1e-3
            numerical = (
                -volume
                * (
                    model.pressure(volume + step, temperature)
                    - model.pressure(volume - step, temperature)
                )
                / (2 * step)
            )
            assert model.bulk_modulus(volume, temperature) == pytest.approx(
                numerical, rel=1e-8
            )
    static = get_eos_record(STATIC)
    for pressure in (0.0001, 5.0, 10.16):
        volume = static.volume(pressure, 300.0, check_validity=True)
        assert static.pressure(
            volume, 300.0, check_validity=(pressure == 5.0)
        ) == pytest.approx(pressure, abs=1e-9)
    with pytest.raises((MaterialError, EosValidationError)):
        record.volume(8.3, 800, check_validity=True)
    with pytest.raises((MaterialError, EosValidationError)):
        record.volume(5, 1300, check_validity=True)
    with pytest.raises((MaterialError, EosValidationError)):
        static.volume(5, 1000, check_validity=True)


def test_diopside_diagnostic_refits_and_reproducibility_qualifications():
    result = reproduce()
    stored = json.loads((ROOT / "docs/data/diopside-reproduction.json").read_text())
    assert result["resources"] == stored["resources"]
    zhao = result["zhao_1998"]
    assert zhao["published_curve_all_rows"]["pressure_rmse_gpa"] == pytest.approx(
        0.426861235534794
    )
    assert zhao["published_curve_without_row15_diagnostic_only"][
        "pressure_rmse_gpa"
    ] == pytest.approx(0.06552064301005417)
    assert zhao["author_fit_reproduced"] is False
    assert zhao["author_workflow_parity"] is False
    full, sensitivity = zhao["diagnostic_refits"]
    assert full["observations"] == 32 and full["excluded_source_rows"] == []
    assert sensitivity["observations"] == 31 and sensitivity[
        "excluded_source_rows"
    ] == [15]
    assert full["coefficients"]["K0_prime"] < 1
    for name, published, error in [
        ("K0", 109.1, 4),
        ("K0_prime", 4.84, 0.6),
        ("alpha0", 2.32e-5, 0.05e-5),
        ("alpha1", 1.88e-8, 0.07e-8),
        ("dK_dT", -0.0205, 0.0004),
    ]:
        assert sensitivity["coefficients"][name] == pytest.approx(published, abs=error)
    thompson = result["thompson_2008"]
    assert thompson["published_curve_all_rows"]["pressure_rmse_gpa"] == pytest.approx(
        0.03739721229323705
    )
    for name, published, error in [
        ("V0", 438.66, 0.02),
        ("K0", 118, 1),
        ("K0_prime", 3.8, 0.3),
    ]:
        assert thompson["diagnostic_refit"]["coefficients"][name] == pytest.approx(
            published, abs=error
        )
    assert thompson["author_workflow_parity"] is False


def test_diopside_calibrations_and_rights_scope():
    records = get_material_document("diopside")["eos_records"]
    assert (
        records[0]["pressure_calibration"]["methods"][0]["reference_calibration_record"]
        == "ruby_mao_1978"
    )
    assert (
        records[1]["pressure_calibration"]["methods"][0]["reference"]["doi"]
        == "10.1063/1.1660714"
    )
    for record in records:
        assert (
            record["pressure_calibration"]["recalculation"]["status"]
            == "missing_calibrant_observations"
        )
    license_text = (DATA / "diopside.LICENSE.md").read_text()
    assert "contributor rights alone" in license_text
    assert "Source rights remain unspecified" in license_text


def test_diopside_plot_evidence_supports_table_error_without_rewriting_observations():
    from scripts.diopside_check_plots import inspect_plots

    result = inspect_plots()
    stored = json.loads((ROOT / "docs/data/diopside-plot-check.json").read_text())
    assert result["row15"] == stored["row15"]
    parity = result["parity_assessment"]
    assert parity["parameter_parity_with_documented_table_error"] is True
    assert parity["published_curves_agree_at_graphical_resolution"] is True
    assert parity["source_table_errors_block_published_eos_acceptance"] is False
    assert parity["exact_original_author_workflow_recovered"] is False
    assert all(
        check["within_reported_error"] for check in parity["parameter_checks"].values()
    )
    validation = get_material_document("diopside")["eos_records"][1][
        "scientific_validation"
    ]
    summary = validation["parity_assessment"]
    for key, value in parity.items():
        if key != "parameter_checks":
            assert summary[key] == value
    assert summary["parameters_within_reported_errors"] == list(
        parity["parameter_checks"]
    )
    assert summary["parameter_check_report"] == "docs/data/diopside-plot-check.json"
    point = result["row15"]
    assert point["digitized_pressure_gpa"] == pytest.approx(4.58, abs=0.05)
    assert point["digitized_volume_a3"] == pytest.approx(431.66, abs=0.2)
    assert point["published_eos_pressure_gpa"] == pytest.approx(
        point["digitized_pressure_gpa"], abs=0.05
    )
    assert point["table_pressure_gpa"] == 6.95
    assert result["ambient"]["digitized_volume_a3"] == pytest.approx(438.67, abs=0.2)
    assert (
        max(abs(x) for x in result["figure2_curves"]["digitized_minus_model_a3"]) < 0.3
    )
    diagnostic = result["graphical_interpretation_diagnostic"]
    assert diagnostic["published_curve_pressure_rmse_gpa"] < 0.07
    assert diagnostic["original_author_workflow_parity"] is False
    for name, published, error in [
        ("K0", 109.1, 4),
        ("K0_prime", 4.84, 0.6),
        ("alpha0", 2.32e-5, 0.05e-5),
        ("alpha1", 1.88e-8, 0.07e-8),
        ("dK_dT", -0.0205, 0.0004),
    ]:
        assert diagnostic["refit_coefficients"][name] == pytest.approx(
            published, abs=error
        )
    rows = list(
        csv.DictReader(
            (DATA / "diopside-zhao-1998-table1.csv").read_text().splitlines()
        )
    )
    assert rows[14]["pressure_gpa"] == "6.95"
    assert get_material("diopside").default_record().identifier == STATIC
    assert get_eos_record(THERMAL).eos.rt_eos.K0 == 109.1


def test_diopside_corrected_dataset_is_reusable_and_preserves_correction_provenance():
    from scripts.diopside_check_plots import inspect_plots
    from scripts.diopside_correct_dataset import IDENTIFIER, OUTPUT, corrected_csv

    report = inspect_plots()
    assert (
        OUTPUT.read_bytes()
        == corrected_csv(report["row15"]["digitized_pressure_gpa"]).encode()
    )
    corrected = list(csv.DictReader(OUTPUT.read_text().splitlines()))
    original = list(
        csv.DictReader(
            (DATA / "diopside-zhao-1998-table1.csv").read_text().splitlines()
        )
    )
    assert len(corrected) == len(original) == 32
    for source, analysis in zip(original, corrected):
        for field in (
            "source_row",
            "temperature_k",
            "temperature_error_k",
            "cell_volume_a3",
            "cell_volume_error_a3",
            "pressure_error_gpa",
        ):
            assert analysis[field] == source[field]
        assert analysis["source_pressure_gpa"] == source["pressure_gpa"]
        if source["source_row"] != "15":
            assert analysis["pressure_gpa"] == source["pressure_gpa"]
            assert analysis["pressure_correction_applied"] == "0"
            assert analysis["pressure_graphical_reading_bound_gpa"] == ""
    assert corrected[14]["pressure_gpa"] == "4.58"
    assert corrected[14]["source_pressure_gpa"] == "6.95"
    assert corrected[14]["pressure_correction_applied"] == "1"
    assert corrected[14]["pressure_source"] == "figure_2_digitized"
    assert corrected[14]["pressure_graphical_reading_bound_gpa"] == "0.05"
    assert corrected[14]["pressure_error_gpa"] == ""
    assert corrected[-1]["duplicate_of_source_row"] == "1"
    assert "lattice_a_angstrom" not in corrected[0]
    assert "source_lattice_volume_mismatch" in corrected[0]["quality_notes"]
    pressures = np.array([float(r["pressure_gpa"]) for r in corrected])
    volumes = np.array([float(r["cell_volume_a3"]) for r in corrected])
    temperatures = np.array([float(r["temperature_k"]) for r in corrected])
    residual = thermal_pressure(volumes, temperatures) - pressures
    assert np.sqrt(np.mean(residual**2)) == pytest.approx(
        report["graphical_interpretation_diagnostic"][
            "published_curve_pressure_rmse_gpa"
        ]
    )
    document = get_material_document("diopside")
    metadata = next(d for d in document["datasets"] if d["identifier"] == IDENTIFIER)
    assert metadata["preferred_for_refitting"] is True
    assert metadata["used_by_eos_records"] == []
    assert document["eos_records"][1]["recommended_refit_dataset"] == IDENTIFIER
    assert document["eos_records"][1]["fit_datasets"] == ["diopside_zhao_1998_table1"]
