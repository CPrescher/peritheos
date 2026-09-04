import pytest

from peritheos import Material, get_material_document

RECORD_ID = "ca_perovskite_akber_knutson_2002_bm3_4"


def _source_and_record():
    document = get_material_document("ca_perovskite")
    source = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    record = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]
    return document, source, record


def _bm3_pressure(volume):
    v0 = 45.92
    k0 = 229.0
    k0_prime = 4.3
    compression = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5
        * k0
        * (compression**7 - compression**5)
        * (1.0 + 0.75 * (k0_prime - 4.0) * (compression**2 - 1.0))
    )


def test_akber_knutson_cubic_vibc_record_preserves_source_choices():
    document, source, _ = _source_and_record()

    assert document["formula"] == "CaSiO3"
    assert document["space_group"] == "Pm-3m"
    assert document["space_group_number"] == 221
    assert document["formula_units_per_cell"] == 1
    assert source["reference"]["doi"] == "10.1029/2001GL013523"
    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 45.92, "K0": 229.0, "K0_prime": 4.3},
    }
    assert source["parameter_errors"] == {
        "V0": 0.02,
        "K0": 2.0,
        "K0_prime": 0.1,
    }
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert source["fixed_parameters"] == []
    assert source["fit_datasets"] == []
    assert source["temperature_ref"] == 300.0
    assert source["pressure_range_status"] == "theoretical"
    assert source["experimental_pressure_range_gpa"] == [0.0, 150.0]


def test_akber_knutson_bm3_reproduces_figure_1_curve_and_inverts():
    _, source, record = _source_and_record()
    checks = source["scientific_validation"]["numerical_reproduction"][
        "figure_1_curve_checks"
    ]

    for check in checks:
        volume = check["volume_a3_per_formula_unit"]
        pressure = record.pressure(volume)
        assert pressure == pytest.approx(_bm3_pressure(volume), abs=1.0e-12)
        assert pressure == pytest.approx(check["calculated_pressure_gpa"], abs=5.0e-7)
        assert pressure == pytest.approx(
            check["figure_pressure_gpa_approx"],
            abs=check["absolute_tolerance_gpa"],
        )

    assert record.pressure(45.92) == pytest.approx(0.0, abs=1.0e-12)
    assert record.eos.bulk_modulus(45.92) == pytest.approx(229.0)
    for pressure in (10.0, 50.0, 100.0, 149.9):
        volume = record.volume(pressure, check_validity=True)
        assert record.pressure(volume, check_validity=True) == pytest.approx(
            pressure, rel=1.0e-11
        )


def test_akber_knutson_phase_model_and_calibration_boundaries_are_explicit():
    document, source, _ = _source_and_record()
    validation = source["scientific_validation"]
    parameterizations = validation["reported_parameterizations"]

    assert [item["structure"] for item in parameterizations] == [
        "Pnma",
        "Pm3m",
        "Pm3m",
    ]
    assert [item["disposition"] for item in parameterizations] == [
        "excluded_distinct_phase",
        "selected",
        "audit_only",
    ]
    assert parameterizations[0]["V0_a3_per_formula_unit"] == 45.9
    assert parameterizations[1]["model"] == "VIBC"
    assert parameterizations[2]["model"] == "VIB without covalent correction"
    assert validation["calibration_model_boundary"]["material"] == "stishovite SiO2"
    assert (
        validation["primary_data_check"]["status"]
        == "theoretical_parameterization_only"
    )
    assert source["pressure_calibration"]["status"] == "not_applicable"
    assert source["pressure_calibration"]["methods"][0]["kind"] == "other"
    assert (
        "does not distinguish"
        in validation["unresolved_issues"]["figure_1_curve_symmetry"]
    )

    assert not any(
        RECORD_ID in dataset.get("used_by_eos_records", [])
        for dataset in document.get("datasets", [])
    )
    assert not any(
        "akber_knutson_2002" in record["identifier"]
        and record["identifier"] != RECORD_ID
        for record in document["eos_records"]
    )


def test_akber_knutson_figure_values_are_not_reference_state_identities():
    _, source, _ = _source_and_record()
    checks = source["scientific_validation"]["numerical_reproduction"][
        "figure_1_curve_checks"
    ]

    assert all(check["volume_a3_per_formula_unit"] != 45.92 for check in checks)
    assert checks[0]["figure_pressure_gpa_approx"] == 43.0
