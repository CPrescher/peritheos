"""Final-source custody, actual fit scope, independent benchmarks and limitations."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import Avogadro

from peritheos import Material, get_material, get_material_document
from scripts.reproduce_ma_2024_kcl import (
    DATA,
    RECORD_ID,
    SOURCE,
    SOURCE_SHA256,
    SUPPLEMENT,
    SUPPLEMENT_SHA256,
    TABLES,
    acoustic_reduction,
    debye_fit,
    extract,
    joint_fit,
    load,
    matsui_nacl_pressure,
    pressure,
    publisher_supplement_check,
    reproduce,
    thermal_pressure,
    walker_input_check,
    workbook_cells,
)

ROOT = Path(__file__).resolve().parents[1]


def test_publisher_supplement_custody_and_printed_table_agreement(tmp_path):
    assert hashlib.sha256(SUPPLEMENT.read_bytes()).hexdigest() == SUPPLEMENT_SHA256
    check = publisher_supplement_check()
    assert check["table_rows_including_headers"] == [22, 13, 48, 50]
    assert check["checked_cells"] > 1400
    assert check["differences_beyond_printed_rounding"] == []
    assert check["embedded_files"] == []
    invalid = tmp_path / "wrong.docx"
    invalid.write_bytes(b"not the publisher supplement")
    with pytest.raises(ValueError, match="recovered publisher supplement"):
        publisher_supplement_check(invalid)


@pytest.mark.parametrize(
    "relative_volume,temperature,source_pressure",
    [
        # Independent Matsui 2012 Table 1 values, pp. 1672-1673. Volumes rounded
        # to four decimals and pressures to two; not values from Ma's reductions.
        (0.7669, 300, 12.04),
        (0.7702, 300, 11.73),
        (0.7742, 300, 11.35),
        (0.7860, 300, 10.31),
        (0.8115, 300, 8.30),
        (0.8291, 300, 7.07),
        (0.8355, 300, 6.66),
        (0.8310, 473, 7.43),
        (0.8358, 473, 7.12),
        (0.8340, 673, 7.81),
        (0.8379, 673, 7.56),
    ],
)
def test_matsui_primary_table_benchmarks(relative_volume, temperature, source_pressure):
    assert matsui_nacl_pressure(
        relative_volume * 179.425, temperature
    ) == pytest.approx(source_pressure, abs=0.02)


def test_final_author_workbook_is_checksummed_and_extractions_are_lossless(
    tmp_path, monkeypatch
):
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    original = workbook_cells(SOURCE)
    assert original["Table S2"]["D3"] == "32.48(9)"
    assert "2 standard deviations" in original["Table S1"]["A23"]
    assert original["Table S1"]["E3"] == pytest.approx(24.321784093513394)
    assert original["Table S1"]["H15"] == pytest.approx(3.1578608426387715)
    # Re-extract to an isolated directory; never rewrite checked-in observations.
    from scripts import reproduce_ma_2024_kcl as audit

    monkeypatch.setattr(audit, "DATA", tmp_path)
    monkeypatch.setattr(audit, "SOURCE", tmp_path / "source" / "tables.xlsx")
    audit.SOURCE.parent.mkdir()
    rows = extract(SOURCE)
    assert {key: len(value) for key, value in rows.items()} == {
        "acoustic": 11,
        "walker": 8,
        "cold_grid": 186,
        "thermal_grid": 371,
    }
    for filename, _ in TABLES.values():
        assert (tmp_path / filename).read_bytes() == (DATA / filename).read_bytes()
    assert json.loads((tmp_path / "source/cached-cells.json").read_text()) == original
    invalid = tmp_path / "wrong-version.xlsx"
    invalid.write_bytes(b"not the final workbook")
    with pytest.raises(ValueError, match="final 7svmv9hvft"):
        workbook_cells(invalid)


def test_published_record_units_configuration_scope_and_existing_default():
    document = get_material_document("kcl")
    records = {r["identifier"]: r for r in document["eos_records"]}
    source = records[RECORD_ID]
    assert source["default"] is False
    assert [r["identifier"] for r in records.values() if r.get("default")] == [
        "kcl_b2_dewaele_2012_vinet_3"
    ]
    assert source["eos"]["parameters"]["V0"] * Avogadro * 1e-24 == pytest.approx(32.48)
    assert source["parameter_errors"]["V0"] * Avogadro * 1e-24 == pytest.approx(0.09)
    assert source["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 251.0,
        "gamma0": 1.92,
        "q": 1.0,
        "n": 2,
    }
    assert source["thermal"]["fixed_parameters"] == ["Tr", "q", "n"]
    assert source["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert source["thermal"]["thermal_pressure_reference"] == "reference_temperature"
    assert source["experimental_temperature_range_k"] == [300, 300]
    assert source["experimental_pressure_range_gpa"][-1] == 85
    assert "parameter_covariance" not in source
    assert source["uncertainty_provenance"]["status"] == "not_deposited"
    assert source["fit_datasets"] == ["kcl_ma_2024_acoustic", "kcl_ma_2024_walker"]
    assert source["pressure_calibration"]["status"] == "partially_resolved"
    assert (
        source["pressure_calibration"]["methods"][1]["reference"]["doi"]
        == "10.2138/am.2012.4136"
    )
    assert "reference_eos_record" not in source["pressure_calibration"]["methods"][1]
    assert (
        source["scientific_validation"]["primary_source_check"]["doi"]
        == "10.1029/2024JB028819"
    )
    assert "10.17632/7svmv9hvft.1" in [r["doi"] for r in source["source_lineage"]]


def test_published_source_grids_are_predictions_and_independent_benchmarks():
    record = get_material("kcl").get_eos_record(RECORD_ID)
    cold, thermal = load("cold_grid"), load("thermal_grid")
    assert min(r["molar_volume_cm3_mol"] for r in cold) == pytest.approx(14.0)
    assert max(r["model_pressure_gpa"] for r in cold) == pytest.approx(
        144.18264431138437
    )
    volumes = np.array([r["molar_volume_cm3_mol"] for r in cold])
    published = np.array([r["model_pressure_gpa"] for r in cold])
    result = record.pressure(volumes / (Avogadro * 1e-24), 300.0)
    # Printed Table 1 rounding, independently diagnosed from deposited S3.
    assert np.max(np.abs(result - published)) < 0.107
    assert np.all(
        np.abs(result - published)
        < np.array([r["model_pressure_error_gpa"] for r in cold])
    )
    assert result == pytest.approx(pressure(volumes), rel=2e-13, abs=2e-13)
    v0 = 32.48 / (Avogadro * 1e-24)
    temperatures = np.array([r["temperature_k"] for r in thermal])
    increments = record.pressure(v0, temperatures)
    assert temperatures[[0, -1]] == pytest.approx([300, 4000])
    assert (
        np.max(
            np.abs(
                increments
                - np.array([r["model_thermal_pressure_gpa"] for r in thermal])
            )
        )
        < 0.014
    )
    assert increments[0] == pytest.approx(0.0, abs=1e-13)
    for volume, temperature in [(15.0, 2000.0), (20.0, 1000.0), (32.48, 4000.0)]:
        expected = float(pressure(volume)) + thermal_pressure(volume, temperature)
        assert record.pressure(
            volume / (Avogadro * 1e-24), temperature
        ) == pytest.approx(expected, rel=2e-12)


def test_joint_fit_uses_reduced_kt_and_recalibrated_pv_with_honest_weights(monkeypatch):
    from scripts import reproduce_ma_2024_kcl as audit

    original = audit.load

    def forbid_model_inputs(key):
        assert key in {"acoustic", "walker"}
        rows = original(key)
        # Corrupt downstream EOS pressures: their values must not affect fits.
        for row in rows:
            row["model_pressure_gpa"] = -999999.0
        return rows

    baseline = joint_fit()
    monkeypatch.setattr(audit, "load", forbid_model_inputs)
    replay = joint_fit()
    assert replay["parameters_molar"] == baseline["parameters_molar"]
    assert replay["observations"] == 19
    assert replay["acoustic_rows"] == 11 and replay["walker_rows"] == 8
    assert replay["success"]
    assert replay["parameters_molar"] == pytest.approx(
        {
            "V0_cm3_mol": 32.47937363348643,
            "K0_gpa": 21.348727595446597,
            "K0_prime": 4.846212436598056,
        },
        rel=2e-7,
    )
    assert abs(replay["difference_from_published"][0]) < 0.09
    assert abs(replay["difference_from_published"][1]) < 0.70
    assert abs(replay["difference_from_published"][2]) < 0.083
    assert joint_fit("unweighted")["parameters_molar"]["K0_gpa"] < 17.0
    assert (
        abs(joint_fit("effective_errors")["parameters_molar"]["K0_gpa"] - 21.33) < 0.70
    )
    assert "not the source covariance" in replay["covariance_provenance"]
    with pytest.raises(ValueError, match="unknown joint-fit"):
        joint_fit("invented_source_weights")


def test_acoustic_debye_reduction_and_thermal_stage_expose_printed_eq9_discrepancy():
    for row in load("acoustic"):
        reduced = acoustic_reduction(row)
        assert abs(reduced["source_theta_difference_k"]) < 4e-6
        assert abs(reduced["source_kt_difference_gpa"]) < 0.021
        assert row["theta_d_k"] / reduced["printed_eq9_theta_d_k"] == pytest.approx(
            2 * 2 ** (1 / 3), rel=1e-8
        )
    fit = debye_fit()
    assert abs(fit["gamma0"] - 1.92) < 0.11
    assert abs(fit["theta0_k"] - 251) < 22


def test_original_walker_pairs_are_retained_without_claiming_upstream_parity():
    audit = walker_input_check()
    assert (
        audit["status"]
        == "deposited_inputs_preserved_upstream_recalculation_not_reproduced"
    )
    rows = audit["rows"]
    assert len(rows) == 8
    assert rows[0]["nacl_file"] == "r35189"
    assert rows[0]["kcl_file"] == "r35193"
    assert rows[0]["difference_gpa"] == pytest.approx(-0.147057334, abs=1e-8)
    assert all(abs(r["difference_gpa"]) > 0.04 for r in rows)
    for original, deposited in zip(rows, load("walker")):
        assert original["kcl_molar_volume_from_original_cell_cm3_mol"] == pytest.approx(
            deposited["molar_volume_cm3_mol"], abs=7e-5
        )
        assert original["nacl_lattice_a_esd_angstrom"] > 0
        assert original["b2_kcl_cell_volume_esd_a3"] > 0
        assert (
            round(original["kcl_molar_volume_from_original_cell_cm3_mol"], 4)
            == deposited["molar_volume_cm3_mol"]
        )
        assert original["source_temperature_k"] == deposited["source_temperature_k"]
        assert original["printed_lattice_half_digit_pressure_bound_gpa"] < 0.0017
        assert (
            abs(original["difference_gpa"])
            > 25 * original["printed_lattice_half_digit_pressure_bound_gpa"]
        )
        assert matsui_nacl_pressure(
            original["implied_nacl_lattice_at_300k_angstrom"] ** 3
        ) == pytest.approx(deposited["matsui_300k_pressure_gpa"], abs=1e-9)


def test_api_native_python_inversion_and_roundtrip():
    material = get_material("kcl")
    record = material.get_eos_record(RECORD_ID)
    restored = Material.from_eosmat(material.to_eosmat()).get_eos_record(RECORD_ID)
    volumes = np.array([15.0, 20.0, 24.0]) / (Avogadro * 1e-24)
    temperatures = np.array([300.0, 1200.0, 3000.0])
    pressures = record.pressure(volumes, temperatures)
    assert record.volume(pressures, temperatures) == pytest.approx(volumes, rel=2e-10)
    assert restored.pressure(volumes, temperatures) == pytest.approx(
        pressures, rel=2e-13
    )
    native = record.eos
    assert native._native is not None
    # Force the independently implemented Python fallback for the same model.
    from peritheos.eos.rt import BM3
    from peritheos.eos.thermal import MieGruneisenDebye

    python = MieGruneisenDebye(BM3(3.248, 21.33, 4.836), 300, 251, 1.92, 1, 2)
    python._native = None
    assert python.pressure(np.array([1.5, 2.0, 2.4]), temperatures) == pytest.approx(
        pressures, rel=2e-13
    )


def test_committed_report_has_source_grid_precision_and_scientific_limits():
    report = json.loads((ROOT / "docs/data/ma-2024-kcl-reproduction.json").read_text())
    replay = reproduce()
    assert replay["source"] == report["source"]
    for name in ["cold_max_abs_error_gpa", "thermal_max_abs_error_gpa"]:
        assert replay["source_grid_validation"][name] == pytest.approx(
            report["source_grid_validation"][name], abs=1e-8
        )
    grid = report["source_grid_validation"]["grid_precision_diagnostic"]
    assert replay["publisher_supplement_check"] == report["publisher_supplement_check"]
    assert (
        replay["cold_grid_covariance_identifiability"][
            "first_order_variance_design_rank"
        ]
        == 5
    )
    assert grid["cold_max_abs_residual_gpa"] < 1e-10
    assert grid["thermal_max_abs_residual_gpa"] < 1e-10
    assert report["source"]["excluded_older_deposit"] == "10.17632/mws6hnp49j.1"
    assert "not a refit of observations" in grid["meaning"]
    for mode in report["fits"]:
        assert replay["fits"][mode]["parameters_molar"] == pytest.approx(
            report["fits"][mode]["parameters_molar"], rel=1e-6
        )


def test_dataset_metadata_classifies_model_outputs_and_uncertainty():
    document = get_material_document("kcl")
    datasets = {r["identifier"]: r for r in document["datasets"]}
    for key, (filename, columns) in TABLES.items():
        dataset = datasets["kcl_ma_2024_" + key]
        assert dataset["license"] == "CC-BY-4.0"
        assert (
            dataset["resource"]["sha256"]
            == hashlib.sha256((DATA / filename).read_bytes()).hexdigest()
        )
        with (DATA / filename).open(newline="") as stream:
            assert next(csv.reader(stream)) == columns
        assert [r["name"] for r in dataset["columns"]] == columns
        for column in dataset["columns"]:
            if column["name"].startswith("model_pressure"):
                assert column["quantity"] == "modeled_pressure"
            if "_error_" in column["name"]:
                assert column["role"] == "uncertainty"
        if "grid" in key:
            assert dataset["kind"] == "model_grid"
