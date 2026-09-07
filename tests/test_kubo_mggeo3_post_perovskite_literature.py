import csv
import hashlib
import math
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document
from peritheos.errors import MaterialError

MATERIAL_IDENTIFIER = "mggeo3_post_perovskite"
RECORD_IDENTIFIER = "mggeo3_post_perovskite_kubo_2006_bm3_1"
DATASET_IDENTIFIER = "mggeo3_post_perovskite_kubo_2006_table_s1_pv"
SOURCE_TABLE_SHA256 = "d6be700950a3c9421a649701b88c6bd0cd5ce4a84d5608beea4417c3d758684d"


def _source_record_dataset_and_rows():
    document = get_material_document(MATERIAL_IDENTIFIER)
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_IDENTIFIER
    )
    dataset = next(
        dataset
        for dataset in document["datasets"]
        if dataset["identifier"] == DATASET_IDENTIFIER
    )
    resource = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with resource.open("rb") as stream:
        content = stream.read()
    rows = list(csv.DictReader(content.decode("utf-8").splitlines()))
    return document, source, dataset, rows, content


def test_kubo_2006_record_structure_and_table_s1_transcription():
    document, source, dataset, rows, content = _source_record_dataset_and_rows()

    assert document["space_group"] == "Cmcm"
    assert document["formula_units_per_cell"] == 4
    assert [site["wyckoff"] for site in document["atom_sites"]] == [
        "4a",
        "4c",
        "4c",
        "8f",
    ]
    assert math.prod(document["lattice"][axis] for axis in "abc") == pytest.approx(
        139.899086334
    )
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(179.2),
        "K0": pytest.approx(207.0),
        "K0_prime": pytest.approx(4.4),
    }
    assert source["parameter_errors"] == {
        "V0": pytest.approx(0.7),
        "K0": pytest.approx(5.0),
        "K0_prime": None,
    }
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert source["fixed_parameters"] == ["K0_prime"]
    assert source["temperature_ref"] == pytest.approx(300.0)
    assert source["experimental_pressure_range_gpa"] == [47.3, 196.3]
    assert dataset["used_by_eos_records"] == [RECORD_IDENTIFIER]
    assert dataset["source_file_sha256"] == SOURCE_TABLE_SHA256
    assert hashlib.sha256(content).hexdigest() == dataset["resource"]["sha256"]

    assert len(rows) == 25
    assert sum(int(row["used_in_published_fit"]) for row in rows) == 22
    assert rows[0] == {
        "run": "1",
        "lattice_a_angstrom": "2.5982",
        "lattice_a_standard_deviation_angstrom": "0.0011",
        "lattice_b_angstrom": "8.4140",
        "lattice_b_standard_deviation_angstrom": "0.0065",
        "lattice_c_angstrom": "6.4043",
        "lattice_c_standard_deviation_angstrom": "0.0037",
        "unit_cell_volume_a3": "140.00603412",
        "pt_111_d_spacing_angstrom": "2.1248",
        "pt_200_d_spacing_angstrom": "1.8451",
        "pressure_gpa": "87.8",
        "differential_stress_gpa": "1.8",
        "exposure_timing": "just after heating",
        "used_in_published_fit": "1",
    }
    assert [
        float(row["pressure_gpa"])
        for row in rows
        if row["used_in_published_fit"] == "0"
    ] == [45.0, 36.3, 7.4]
    for row in rows:
        volume = math.prod(float(row[f"lattice_{axis}_angstrom"]) for axis in "abc")
        assert float(row["unit_cell_volume_a3"]) == pytest.approx(volume, abs=5.0e-9)


def test_kubo_2006_bm3_reproduces_table_s1_and_round_trips():
    document, _, _, rows, _ = _source_record_dataset_and_rows()
    record = Material.from_eosmat(
        document, record_identifiers=[RECORD_IDENTIFIER]
    ).eos_records[0]
    fitted_rows = [row for row in rows if row["used_in_published_fit"] == "1"]
    volumes = np.array([float(row["unit_cell_volume_a3"]) for row in fitted_rows])
    source_pressures = np.array([float(row["pressure_gpa"]) for row in fitted_rows])

    calculated_pressures = record.pressure(volumes, 300.0, check_validity=False)
    residuals = calculated_pressures - source_pressures

    assert calculated_pressures[0] == pytest.approx(88.2910956, abs=1.0e-8)
    assert np.sqrt(np.mean(residuals**2)) == pytest.approx(1.3945065, abs=1.0e-7)
    assert np.max(np.abs(residuals)) == pytest.approx(3.9378893, abs=1.0e-7)
    assert record.pressure(record.reference_volume, 300.0) == pytest.approx(
        0.0, abs=1.0e-14
    )
    model_volume = record.volume(100.0, 300.0, check_validity=True)
    assert record.pressure(model_volume, 300.0, check_validity=True) == pytest.approx(
        100.0
    )
    with pytest.raises(MaterialError, match="outside the published"):
        record.pressure(160.0, 300.0, check_validity=True)


def test_kubo_2006_unweighted_refit_recovers_published_parameters():
    _, _, _, rows, _ = _source_record_dataset_and_rows()
    fitted_rows = [row for row in rows if row["used_in_published_fit"] == "1"]
    volumes = np.array([float(row["unit_cell_volume_a3"]) for row in fitted_rows])
    pressures = np.array([float(row["pressure_gpa"]) for row in fitted_rows])

    def residuals(parameters):
        volume_zero, bulk_modulus = parameters
        eta = (volume_zero / volumes) ** (1.0 / 3.0)
        calculated = (
            1.5
            * bulk_modulus
            * (eta**7 - eta**5)
            * (1.0 + 0.75 * (4.4 - 4.0) * (eta**2 - 1.0))
        )
        return calculated - pressures

    fit = least_squares(
        residuals,
        [179.2, 207.0],
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    volume_zero, bulk_modulus = fit.x

    assert volume_zero == pytest.approx(179.2005459, abs=1.0e-6)
    assert bulk_modulus == pytest.approx(206.7112286, abs=1.0e-6)
    assert abs(volume_zero - 179.2) < 0.7
    assert abs(bulk_modulus - 207.0) < 5.0


def test_kubo_2006_pt111_observations_reproduce_holmes_pressures():
    _, _, _, rows, _ = _source_record_dataset_and_rows()
    platinum = Material.from_eosmat(
        get_material_document("platinum"),
        record_identifiers=["platinum_holmes_1989_vinet_1"],
    ).eos_records[0]
    d_spacings = np.array([float(row["pt_111_d_spacing_angstrom"]) for row in rows])
    platinum_cell_volumes = (np.sqrt(3.0) * d_spacings) ** 3
    source_pressures = np.array([float(row["pressure_gpa"]) for row in rows])

    recalculated = platinum.pressure(platinum_cell_volumes, 300.0, check_validity=False)

    assert recalculated[0] == pytest.approx(87.97145598, abs=1.0e-8)
    assert np.max(np.abs(recalculated - source_pressures)) < 0.446
