import pytest

from peritheos import Material, get_material_document


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
