import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_sun_2016_2022_casio3 import (
    fit_sun_2016,
    fit_sun_2022,
)

ROOT = Path(__file__).resolve().parents[1]


def test_sun_tetragonal_casio3_volume_basis_and_primary_regression():
    document = get_material_document("ca_perovskite_tetragonal")
    source = document["eos_records"][0]
    record = Material.from_eosmat(document).eos_records[0]

    assert document["formula_units_per_cell"] == 4
    assert document["space_group"] == "I4/mcm"
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(4.0 * 45.6),
        "K0": pytest.approx(229.0),
        "K0_prime": pytest.approx(4.0),
    }
    assert source["parameter_errors"] == {
        "V0": pytest.approx(4.0 * 0.2),
        "K0": pytest.approx(4.0),
        "K0_prime": None,
    }
    assert source["fixed_parameters"] == ["K0_prime"]
    assert source["experimental_pressure_range_gpa"] == [21.5, 199.2]
    assert record.validity.pressure_gpa == (30.0, 150.0)
    assert record.pressure(140.6, 300.0) == pytest.approx(100.4396514506)
    with pytest.raises(ValueError, match="outside the published calibration/data"):
        record.volume(21.5, 300.0, check_validity=True)
    assert source["pressure_calibration"]["methods"][0]["material"] == "Pt"
    assert source["pressure_calibration"]["recalculation"]["status"] == (
        "missing_calibrant_observations"
    )
    assert source["fit_datasets"] == ["ca_perovskite_tetragonal_sun_2022_table1_pv"]


def test_sun_cubic_casio3_thermal_parameters_and_primary_regression():
    document = get_material_document("ca_perovskite")
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == "ca_perovskite_sun_2016_bm3_3"
    )
    record = Material.from_eosmat(
        document,
        record_identifiers=["ca_perovskite_sun_2016_bm3_3"],
    ).eos_records[0]

    assert document["formula_units_per_cell"] == 1
    assert record.reference_volume == pytest.approx(45.4)
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(45.4),
        "K0": pytest.approx(249.0),
        "K0_prime": pytest.approx(4.0),
    }
    assert source["thermal"]["parameters"] == {
        "Tr": pytest.approx(300.0),
        "theta0": pytest.approx(1000.0),
        "gamma0": pytest.approx(1.8),
        "q": pytest.approx(1.1),
        "n": pytest.approx(5.0),
    }
    assert record.validity.temperature_k == (1200.0, 2600.0)
    pressure = record.pressure(36.68, 2200.0)
    assert pressure == pytest.approx(94.9028891889)
    assert record.volume(pressure, 2200.0) == pytest.approx(36.68)
    with pytest.raises(ValueError, match="outside the published calibration/data"):
        record.pressure(36.68, 300.0, check_validity=True)
    assert source["pressure_calibration"]["methods"][0]["material"] == "Pt"
    assert source["pressure_calibration"]["recalculation"]["status"] == (
        "missing_calibrant_observations"
    )
    assert source["fit_datasets"] == ["ca_perovskite_sun_2016_table1_pvt"]


@pytest.mark.parametrize(
    ("material_id", "dataset_id", "filename", "row_count", "checksum"),
    [
        (
            "ca_perovskite",
            "ca_perovskite_sun_2016_table1_pvt",
            "ca-perovskite-sun-2016-table1-pvt.csv",
            144,
            "0d1bff5418ef54eb97e8d608121e0154a773871450a32824c45da9b66df1bb2f",
        ),
        (
            "ca_perovskite_tetragonal",
            "ca_perovskite_tetragonal_sun_2022_table1_pv",
            "ca-perovskite-tetragonal-sun-2022-table1-pv.csv",
            23,
            "c05a0bd7e76d8d969ce8cb4e2de0a00b2541b96cd4ebefc211187621ff584028",
        ),
    ],
)
def test_sun_casio3_primary_tables_are_complete_and_checksummed(
    material_id, dataset_id, filename, row_count, checksum
):
    document = get_material_document(material_id)
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == dataset_id
    )
    path = ROOT / "peritheos" / "data" / "datasets" / filename
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == row_count
    assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum
    assert dataset["resource"]["sha256"] == checksum
    assert "not asserted to be openly licensed" in dataset["license"]


def test_sun_2016_unweighted_source_protocol_reproduces_model_1():
    result = fit_sun_2016()
    scale = get_eos_record("ca_perovskite_sun_2016_bm3_3").volume_scale

    assert result.success
    assert [
        result.parameters["rt_eos.V0"] / scale,
        result.parameters["rt_eos.K0"],
        result.parameters["gamma0"],
        result.parameters["q"],
    ] == pytest.approx(
        [45.4631688751, 250.264586487, 1.87264323881, 1.04254373813],
        rel=2.0e-8,
    )
    assert [
        result.standard_errors["rt_eos.V0"] / scale,
        result.standard_errors["rt_eos.K0"],
        result.standard_errors["gamma0"],
        result.standard_errors["q"],
    ] == pytest.approx(
        [0.127774556706, 4.46330160974, 0.161015674523, 0.404541872717],
        rel=2.0e-7,
    )
    assert np.sqrt(np.mean(result.residuals**2)) == pytest.approx(0.936048593565)


def test_sun_2022_unweighted_source_protocol_reproduces_both_table_2_fits():
    fixed = fit_sun_2022(fixed_k0_prime=True)
    free = fit_sun_2022(fixed_k0_prime=False)

    assert fixed.success and free.success
    assert [fixed.parameters["V0"] / 4.0, fixed.parameters["K0"]] == pytest.approx(
        [45.5620119969, 228.983858403], rel=2.0e-8
    )
    assert [
        fixed.standard_errors["V0"] / 4.0,
        fixed.standard_errors["K0"],
    ] == pytest.approx([0.15792732064, 4.15911556756], rel=2.0e-7)
    assert [
        free.parameters["V0"] / 4.0,
        free.parameters["K0"],
        free.parameters["K0_prime"],
    ] == pytest.approx([45.6042497227, 226.687151416, 4.02795199335], rel=2.0e-8)
    assert [
        free.standard_errors["V0"] / 4.0,
        free.standard_errors["K0"],
        free.standard_errors["K0_prime"],
    ] == pytest.approx([0.42277993001, 21.4924952316, 0.257598311607], rel=2.0e-7)


def test_sun_2022_table_subsets_and_lattice_volumes_are_not_conflated():
    path = (
        ROOT
        / "peritheos"
        / "data"
        / "datasets"
        / "ca-perovskite-tetragonal-sun-2022-table1-pv.csv"
    )
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert [row["source_subset"] for row in rows].count("sun_2016_reanalyzed") == 14
    assert [row["source_subset"] for row in rows].count("sun_2022_new") == 9
    for row in rows:
        lattice_volume = float(row["a_angstrom"]) ** 2 * float(row["c_angstrom"])
        assert lattice_volume == pytest.approx(
            float(row["volume_a3_conventional_cell"]), abs=0.15
        )


def test_sun_pressure_scale_parameters_are_executable_but_rows_are_missing():
    platinum = get_eos_record("pt_fcc_fei_2007")
    assert platinum.validity.pressure_gpa == (0.0, 94.0)
    assert platinum.validity.temperature_k == (300.0, 1873.0)
    assert platinum.eos.rt_eos.parameter_values() == {
        "V0": pytest.approx(60.38 * platinum.volume_scale),
        "K0": pytest.approx(277.0),
        "K0_prime": pytest.approx(5.08),
    }
