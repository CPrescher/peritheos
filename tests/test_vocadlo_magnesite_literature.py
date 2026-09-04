import csv
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document

RECORD_IDENTIFIER = "magnesite_vocadlo_1999_bm3_2"
DATASET_IDENTIFIER = "magnesite_vocadlo_1999_table1_energy_volume"
GPA_A3_TO_EV = 1.0e-21 / 1.602176634e-19


def _source_record_and_rows():
    document = get_material_document("magnesite")
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
    with resource.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    return document, source, dataset, rows


def test_vocadlo_1999_record_and_table_i_transcription():
    document, source, dataset, rows = _source_record_and_rows()

    assert document["space_group"] == "R-3c"
    assert document["formula_units_per_cell"] == 6
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(287.49),
        "K0": pytest.approx(99.0),
        "K0_prime": pytest.approx(4.28),
    }
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(0.5),
        "K0_prime": pytest.approx(0.01),
    }
    assert source["fixed_parameters"] == []
    assert source["temperature_ref"] == pytest.approx(0.0)
    assert source["experimental_pressure_range_gpa"] == [-2.41596, 207.419]
    assert source["pressure_range_status"] == "theoretical"
    assert dataset["used_by_eos_records"] == [RECORD_IDENTIFIER]
    assert len(rows) == 9
    assert rows[0] == {
        "energy_ev": "-216.32189",
        "volume_a3": "295",
        "pressure_gpa": "-2.41596",
    }
    assert rows[-1] == {
        "energy_ev": "-167.715033",
        "volume_a3": "160",
        "pressure_gpa": "207.419",
    }


def test_vocadlo_1999_bm3_reproduces_all_table_i_pressures():
    document, _, _, rows = _source_record_and_rows()
    record = Material.from_eosmat(
        document, record_identifiers=[RECORD_IDENTIFIER]
    ).eos_records[0]
    volumes = np.array([float(row["volume_a3"]) for row in rows])
    source_pressures = np.array([float(row["pressure_gpa"]) for row in rows])

    calculated_pressures = record.pressure(volumes, check_validity=False)

    np.testing.assert_allclose(calculated_pressures, source_pressures, atol=1.2e-4)
    assert record.pressure(record.reference_volume) == pytest.approx(0.0, abs=1e-14)
    assert record.volume(207.419, check_validity=False) == pytest.approx(
        160.0, abs=5.0e-5
    )


def test_vocadlo_1999_energy_volume_refit_recovers_published_bm3():
    _, _, _, rows = _source_record_and_rows()
    volumes = np.array([float(row["volume_a3"]) for row in rows])
    energies = np.array([float(row["energy_ev"]) for row in rows])

    def residuals(parameters):
        energy_zero, volume_zero, bulk_modulus, bulk_modulus_prime = parameters
        strain_term = (volume_zero / volumes) ** (2.0 / 3.0)
        calculated = energy_zero + (
            9.0
            * volume_zero
            * bulk_modulus
            * GPA_A3_TO_EV
            / 16.0
            * (
                (strain_term - 1.0) ** 3 * bulk_modulus_prime
                + (strain_term - 1.0) ** 2 * (6.0 - 4.0 * strain_term)
            )
        )
        return calculated - energies

    fit = least_squares(
        residuals,
        [-216.38, 287.49, 99.0, 4.28],
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    energy_zero, volume_zero, bulk_modulus, bulk_modulus_prime = fit.x

    assert energy_zero == pytest.approx(-216.381962, abs=1.0e-6)
    assert volume_zero == pytest.approx(287.49198, abs=2.0e-5)
    assert bulk_modulus == pytest.approx(98.9208, abs=1.0e-4)
    assert bulk_modulus_prime == pytest.approx(4.279108, abs=1.0e-6)
    assert np.sqrt(np.mean(residuals(fit.x) ** 2)) < 0.0042
