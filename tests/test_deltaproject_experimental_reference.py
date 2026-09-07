import csv
from importlib import resources

import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_deltaproject_experimental_reference import reproduce


def _rows():
    resource = resources.files("peritheos").joinpath(
        "data", "datasets", "deltaproject-experimental-reference-parameters.csv"
    )
    with resource.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _record(row):
    document = get_material_document(row["material_identifier"])
    return document, next(
        item
        for item in document["eos_records"]
        if item["identifier"] == row["record_identifier"]
    )


def test_delta_rows_are_derived_reference_parameterizations():
    rows = _rows()
    assert len(rows) == 42
    assert len({row["record_identifier"] for row in rows}) == 42
    for row in rows:
        document, record = _record(row)
        assert record["record_kind"] == "derived"
        assert record["pressure_range_status"] == "reference_parameterization"
        assert "experimental_pressure_range_gpa" not in record
        assert record["validity"]["volume_ratio"] == [0.94, 1.06]
        assert (
            record["scientific_validation"]["primary_data_check"]["status"]
            == "parameterization_only"
        )
        executable = Material.from_eosmat(
            document, record_identifiers=[record["identifier"]]
        ).eos_records[0]
        assert executable.pressure(record["eos"]["parameters"]["V0"], 0.0) == (
            pytest.approx(0.0, abs=1.0e-8)
        )


def test_delta_compiled_reference_reproduction():
    result = reproduce()
    assert result["records"] == 42
    assert result["unique_identifiers"] == 42
    assert result["maximum_parameter_error"] == 0.0
    assert result["maximum_independent_pressure_error_gpa"] < 1.0e-12
    assert result["knittle_b0_prime_rows"] == 32
    assert result["knittle_b0_prime_reconstructed"] == 30
    assert result["knittle_b0_prime_unresolved"] == ["Fe", "Au"]
    assert result["alpha_sn_k0_from_elastic_constants_gpa"] == pytest.approx(
        42.53333333333333
    )
    assert result["osmium_fixed_v0_refit"] == pytest.approx(
        {"K0": 395.5658321642279, "K0_prime": 4.4988725479206835}
    )
    assert result["known_source_conflicts"] == 2
    assert result["known_phase_mismatches"] == 1
    assert result["raw_observation_refit"].startswith("not possible")


def test_delta_records_expose_property_level_upstream_recovery():
    expected_status = {
        "manganese_alpha": ("published_fit_recovered_plot_only", None),
        "tin_alpha": ("derived_elastic_property_recovered", "exact_phase_mismatch"),
        "osmium": ("secondary_table_recovered", "row_level_data_recovered_and_refit"),
        "strontium_fcc": ("secondary_table_recovered", "exact_with_upstream_conflict"),
        "thallium_hcp": (
            "secondary_table_recovered",
            "citation_resolved_conflict_unresolved",
        ),
    }
    rows = {row["material_identifier"]: row for row in _rows()}
    for material_identifier, (k0_status, k0_prime_status) in expected_status.items():
        _, record = _record(rows[material_identifier])
        recovery = record["scientific_validation"]["upstream_data_recovery"]
        assert recovery["K0"]["status"] == k0_status
        if k0_prime_status == "exact_phase_mismatch":
            assert recovery["K0_prime"]["element_reconstruction"]["status"] == (
                k0_prime_status
            )
            assert recovery["phase_mismatch"]["element"] == "Sn"
        elif k0_prime_status == "exact_with_upstream_conflict":
            assert recovery["K0_prime"]["element_reconstruction"]["status"] == (
                k0_prime_status
            )
            assert recovery["source_conflict"]["anderson_1990_table_1"] == 2.41
        elif k0_prime_status is not None:
            assert recovery["K0_prime"]["status"] == k0_prime_status
