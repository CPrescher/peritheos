import csv
import hashlib
import io
import math
from importlib import resources

import pytest

from peritheos import Material, get_material_document, validate_eosmat_document

MATERIAL_ID = "sio2_stv_andr"
RECORD_ID = "sio2_stv_andr_schoelmerich_2020_bm3_3"
DATASET_ID = "stishovite_schoelmerich_2020_table1_shock"


def _source_record_and_rows():
    document = get_material_document(MATERIAL_ID)
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_ID
    )
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_bytes()
    )
    rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
    return document, source, dataset, payload, rows


def test_schoelmerich_2020_record_is_the_new_shock_reduced_300_k_fit():
    document, source, _, _, _ = _source_record_and_rows()

    validate_eosmat_document(document)
    assert source["reference"]["doi"] == "10.1038/s41598-020-66340-y"
    assert source["equation_kind"] == "isothermal"
    assert source["record_kind"] == "published"
    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 46.5, "K0": 307.0, "K0_prime": 4.66},
    }
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": 4.0,
        "K0_prime": 0.15,
    }
    assert source["fixed_parameters"] == ["V0"]
    assert source["temperature_ref"] == 300.0
    assert source["experimental_pressure_range_gpa"] == [0.0, 336.0]
    assert source["experimental_temperature_range_k"] == [324.0, 4757.0]
    assert source["fit_datasets"] == [DATASET_ID]

    validation = source["scientific_validation"]
    assert validation["status"] == "primary_source_validated"
    assert validation["primary_source_check"]["doi"] == ("10.1038/s41598-020-66340-y")
    assert validation["primary_data_check"]["status"] == "bundled"
    assert source["thermal_reduction"]["reproduction_status"] == (
        "not_directly_refittable"
    )
    assert source["pressure_calibration"]["methods"][0]["kind"] == "shock_wave"
    assert any(
        "not an equilibrium phase-stability field" in note
        for note in source["validity"]["notes"]
    )
    assert {item["location"] for item in validation["reported_inconsistencies"]} == {
        "Article EOS paragraph",
        "Article Table 1",
    }


def test_schoelmerich_2020_bm3_reproduces_the_published_figure_s5_curve():
    document, _, _, _, _ = _source_record_and_rows()
    record = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]

    assert record.reference_volume == pytest.approx(46.5)
    assert record.pressure(46.5, 300.0, check_validity=True) == pytest.approx(0.0)
    assert record.eos.bulk_modulus(46.5) == pytest.approx(307.0)
    assert record.eos.bulk_modulus_derivative(46.5) == pytest.approx(4.66, rel=2e-6)

    # The solid 300 K curve and grey reduced points in Supplementary Figure S5
    # are approximately 114 GPa at V/V0=37.3/46.5 and 324 GPa at 31.0/46.5.
    assert record.pressure(37.3, 300.0, check_validity=True) == pytest.approx(
        114.0, abs=5.0
    )
    assert record.pressure(31.0, 300.0, check_validity=True) == pytest.approx(
        324.0, abs=5.0
    )
    pressure = record.pressure(31.0, 300.0)
    assert record.volume(pressure, 300.0) == pytest.approx(31.0)
    assert record.within_calibration_range(31.0, 300.0)
    assert not record.within_calibration_range(30.9, 300.0)


def test_schoelmerich_2020_table1_transcription_and_shock_pressure_reduction():
    _, _, dataset, payload, rows = _source_record_and_rows()

    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    assert len(rows) == 6
    assert len(rows[0]) == len(dataset["columns"]) == 17
    assert [row["run_id"] for row in rows] == [
        "SACLA-795997",
        "SACLA-796485",
        "SACLA-796491",
        "LCLS-235",
        "LCLS-233",
        "LCLS-239",
    ]

    ambient = rows[0]
    assert ambient["pressure_gpa"] == ""
    assert float(ambient["volume_a3_conventional_cell"]) == pytest.approx(46.5)
    assert float(ambient["volume_uncertainty_a3_conventional_cell"]) == 0.2

    run_233 = rows[4]
    assert run_233["lattice_a_angstrom"] == ""
    assert run_233["lattice_c_angstrom"] == ""
    assert float(run_233["pressure_gpa"]) == 317.0
    assert float(run_233["volume_a3_conventional_cell"]) == 31.3

    rho0 = float(ambient["density_g_cm3"])
    for row in rows[1:]:
        pressure = float(row["pressure_gpa"])
        pressure_uncertainty = float(row["pressure_uncertainty_gpa"])
        calculated = (
            rho0
            * float(row["shock_velocity_km_s"])
            * float(row["particle_velocity_km_s"])
        )
        assert calculated == pytest.approx(pressure, abs=pressure_uncertainty)

        # The printed energy numbers follow the Rankine-Hugoniot specific-energy
        # result in MJ/kg, exposing the source table's inconsistent kJ/mol label.
        density = float(row["density_g_cm3"])
        calculated_specific_energy = 0.5 * pressure * (1.0 / rho0 - 1.0 / density)
        assert calculated_specific_energy == pytest.approx(
            float(row["reported_energy_increment_kj_mol"]),
            abs=float(row["reported_energy_increment_uncertainty_kj_mol"]),
        )

    for row in rows[:4] + rows[5:]:
        a = float(row["lattice_a_angstrom"])
        c = float(row["lattice_c_angstrom"])
        volume = float(row["volume_a3_conventional_cell"])
        volume_uncertainty = float(row["volume_uncertainty_a3_conventional_cell"])
        assert math.prod((a, a, c)) == pytest.approx(volume, abs=volume_uncertainty)
