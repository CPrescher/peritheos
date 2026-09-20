"""Source-level regressions for Campbell's four experimental thermal EOS."""

import json

import numpy as np
import pytest
from scipy.constants import Avogadro

from peritheos import (
    Material,
    get_eos_record,
    get_material,
    get_material_document,
    list_eos_records,
    list_materials,
    validate_eosmat_document,
)
from scripts.apply_primary_source_audit import audit_record
from scripts.reconstruct_seagle_2008_pressures import (
    calibration_checks,
    hcp_pressure,
    reconstructed_rows,
)
from scripts.reproduce_campbell_2009_buffers import (
    OUTPUT,
    RECORDS,
    SOURCE,
    buffer_checks,
    current_data,
    pressure,
    read_rows,
    seagle_nominal_data,
    seagle_reconstructed_data,
    volume,
)


@pytest.mark.parametrize("material", ["nickel", "nickel_oxide_b1"])
def test_campbell_native_equation_matches_independent_molar_quadrature(material):
    record = get_eos_record(RECORDS[material])
    parameters = SOURCE[material]
    cell_per_molar = 4e24 / Avogadro
    # Compression and heating both matter: this detects use of the integrated
    # Debye-temperature law or an incorrect atoms/formula-units normalization.
    for ratio, temperature in [(1, 295), (0.8, 295), (0.9, 1500), (0.75, 2400)]:
        molar_volume = parameters[0] * ratio
        expected = float(pressure(molar_volume, temperature, parameters))
        actual = record.pressure(
            molar_volume * cell_per_molar, temperature, check_validity=False
        )
        assert actual == pytest.approx(expected, abs=2e-7)
    expected_volume = volume(10, 1500, parameters) * cell_per_molar
    assert record.volume(10, 1500) == pytest.approx(expected_volume, rel=1e-8)
    with pytest.raises(ValueError, match="outside the published calibration"):
        record.volume(200, 4000, check_validity=True)


@pytest.mark.parametrize("material", SOURCE)
def test_campbell_interchange_preserves_source_semantics(material):
    document = get_material_document(material)
    if material in {"fe_fcc", "feo"}:
        # Retain source parameters losslessly without constructing an executable
        # model from the explicitly deferred Fe/FeO records.
        exported = json.loads(json.dumps(document))
        validate_eosmat_document(exported)
    else:
        model = Material.from_eosmat(document)
        exported = model.to_eosmat()
        assert Material.from_eosmat(exported).to_eosmat() == exported
    original = next(
        r for r in document["eos_records"] if r["identifier"] == RECORDS[material]
    )
    restored = next(
        r for r in exported["eos_records"] if r["identifier"] == RECORDS[material]
    )
    for key, value in original.items():
        if key == "thermal":
            assert all(restored[key][k] == v for k, v in value.items())
        else:
            assert restored[key] == value
    assert document["formula_units_per_cell"] == 4
    assert document["space_group_number"] == 225
    assert original["thermal"]["debye_temperature_law"] == "variable_exponent"
    assert original["parameter_errors"]["V0"] is None
    assert original["temperature_ref"] == 295
    assert original["parameter_error_confidence"] is None
    assert "campbell_2009_tables_s4_s5_buffer_results" not in original["fit_datasets"]
    if material in {"fe_fcc", "feo"}:
        assert "K0_prime" in original["fixed_parameters"]
        assert "q" not in original["thermal"]["fixed_parameters"]
    else:
        assert "K0_prime" not in original["fixed_parameters"]
        assert "q" in original["thermal"]["fixed_parameters"]


@pytest.mark.parametrize("material", ["fe_fcc", "feo"])
def test_campbell_mismatch_records_are_source_only_and_stay_deferred(material):
    identifier = RECORDS[material]
    document = get_material_document(material)
    record = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    validation = record["scientific_validation"]
    assert validation["status"] == "deferred"
    assert "discrepancies" in validation["note"]
    assert record["eos"]["parameters"]["K0"] == SOURCE[material][1]
    assert identifier not in {r.identifier for r in list_eos_records()}
    with pytest.raises(KeyError):
        get_eos_record(identifier)
    with pytest.raises(ValueError, match="deferred"):
        Material.from_eosmat(document, record_identifiers=[identifier])
    if material == "fe_fcc":
        assert material not in {m.identifier for m in list_materials()}
        with pytest.raises(KeyError):
            get_material(material)
    else:
        # One deferred record must not disable the other accepted FeO models.
        assert get_material(material).eos_records
        assert identifier not in {
            r.identifier for r in get_material(material).eos_records
        }
    rebuilt = audit_record(record, f"{material}.eosmat")["scientific_validation"]
    assert rebuilt["status"] == "deferred"
    assert rebuilt["primary_source_check"] == validation["primary_source_check"]
    assert rebuilt["unresolved"] == validation["unresolved"]


def test_campbell_observations_keep_phases_errors_and_temperature_convention():
    rows = read_rows("ni-nio-campbell-2009-table-s3-pvt.csv")
    assert len(rows) == 101
    for row in rows:
        assert float(row["delta_nio_ni_molar_volume_cm3_mol"]) == pytest.approx(
            float(row["nio_molar_volume_cm3_mol"])
            - float(row["nickel_molar_volume_cm3_mol"])
        )
        assert float(row["reported_pressure_uncertainty_gpa"]) > 0
    # This source row uses a 293 K ambient baseline; others use 295 K.
    first = rows[0]
    measured_temperature = 293 + (float(first["sample_temperature_k"]) - 293) / 0.97
    assert measured_temperature == pytest.approx(2524)
    assert float(first["nacl_temperature_k"]) == pytest.approx(
        (3 * measured_temperature + 293) / 4
    )
    # The workbook uses /4; the article's Equation (4) prints /2.
    assert float(first["nacl_temperature_uncertainty_k"]) == pytest.approx(
        (measured_temperature - 293) / 4
    )
    assert any(
        r["nacl_phase"] == "B2" and float(r["reported_pressure_gpa"]) < 26 for r in rows
    )
    assert current_data("fe_fcc").shape == (6, 15)
    assert current_data("feo").shape == (6, 25)
    assert current_data("nickel").shape == (6, 92)
    assert current_data("nickel_oxide_b1", include_quenched=True).shape == (6, 101)


def test_campbell_independent_published_buffer_checkpoints_and_refit_limits():
    checks = buffer_checks()
    assert len(checks) == 5
    assert all(abs(row["difference_cm3_mol"]) < 0.02 for row in checks)
    assert max(abs(row["difference_cm3_mol"]) for row in checks) < 0.0072
    report = json.loads(OUTPUT.read_text())
    for row in report["record_refits"].values():
        assert row["solver_success"]
        assert row["rmse_gpa"] < row["published_rmse_gpa"]
        assert "not an exact recovery" in row["qualification"]
        assert np.linalg.eigvalsh(row["covariance"]).min() > 0
    nio = report["record_refits"][RECORDS["nickel_oxide_b1"]]
    assert nio["parameters"]["rt_eos.K0"] == pytest.approx(209.2349, abs=0.001)
    assert nio["sensitivity"]["pressure_uncertainty"]["parameters"][
        "rt_eos.K0"
    ] == pytest.approx(192.5593, abs=0.001)


def test_seagle_hcp_calibration_reproduces_independent_published_pressures():
    checks = calibration_checks()
    original = checks["seagle_2006_table_3"]
    assert original["observations"] == 41
    assert original["density_text"]["within_reported_pressure_uncertainty"] == 41
    assert original["density_text"]["rmse_gpa"] < 0.25
    # Detect accidental use of the conflicting printed Table 2 reference volume.
    assert original["table_2"]["rmse_gpa"] > 2.4
    later = checks["seagle_2008_table_2_lower_eutectic"]
    assert later["observations"] == 5
    assert later["density_text"]["max_absolute_difference_gpa"] < 0.60
    assert float(hcp_pressure(18.15, 2501.5)) == pytest.approx(82, abs=2)


def test_seagle_reconstruction_preserves_missing_pressures_and_phase_selection():
    rows = reconstructed_rows(seagle_nominal_data("feo")[1])
    assert len(rows) == 81
    hcp = [r for r in rows if r["status"] == "hcp_seagle_2006_density_text"]
    assert len(hcp) == 65
    unresolved = [r for r in rows if r not in hcp]
    assert len(unresolved) == 16
    assert all(r["reconstructed_pressure_gpa"] == "" for r in unresolved)
    assert sum(r["status"] == "fcc_calibration_incomplete" for r in rows) == 14
    assert all(
        r["reconstructed_pressure_gpa"] > r["table_2_reference_pressure_gpa"]
        for r in hcp
    )
    assert seagle_reconstructed_data().shape == (6, 65)
    assert seagle_reconstructed_data(include_nominal=True).shape == (6, 81)
    report = json.loads(OUTPUT.read_text())
    fits = report["record_refits"][RECORDS["feo"]]["sensitivity"]
    assert fits["seagle_hcp_density_text"]["observations"] == 90
    combined = fits["seagle_hcp_density_text_plus_nominal"]
    assert combined["observations"] == 106
    assert combined["active_parameter_bounds"] == ["q"]
    assert combined["seagle_pressure_reconstruction"]["nominal_only_rows"] == 16
