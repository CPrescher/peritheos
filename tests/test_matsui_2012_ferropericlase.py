import csv
import hashlib
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.materials import Material, MaterialError
from scripts.reproduce_matsui_2012_ferropericlase import diagnostics, load_data

CASES = {
    "mg083fe017o": {
        "record": "mg083fe017o_matsui_2012_bm3_mgd_1",
        "dataset": "mg083fe017o_matsui_2012_table1_pvt",
        "rows": 34,
        "included": 23,
        "counts": {300.0: 12, 700.0: 11, 1100.0: 11},
        "v0": 75.849,
        "kp": 4.08,
        "gamma0": 1.53,
        "q": 0.7,
        "lattice": 4.2331,
        "sha256": "55e5f1580599a3306784d333176d8b2cf264c9c11ea3d4de86af38fd234fc1f3",
        "observed_rmse": 0.1233027487,
        "observed_max": 0.2163884252,
        "table_rmse": 0.0054742591,
        "table_max": 0.0139179342,
        "phase": "high-spin B1 ferropericlase",
    },
    "mg075fe025o": {
        "record": "mg075fe025o_matsui_2012_bm3_mgd_1",
        "dataset": "mg075fe025o_matsui_2012_table2_pvt",
        "rows": 39,
        "included": 30,
        "counts": {300.0: 14, 700.0: 13, 1100.0: 12},
        "v0": 76.372,
        "kp": 4.22,
        "gamma0": 1.64,
        "q": 0.7,
        "lattice": 4.2427,
        "sha256": "c9202d12b55f1e9892c10d6fea43fa7625a74f9288bc2d9cbb8d28ce512533cd",
        "observed_rmse": 0.1426097852,
        "observed_max": 0.2669266040,
        "table_rmse": 0.0028967125,
        "table_max": 0.0052404910,
        "phase": "B1 ferropericlase",
    },
}


@pytest.mark.parametrize("material_identifier", CASES)
def test_matsui_records_preserve_published_parameters_and_reference_state(
    material_identifier,
):
    expected = CASES[material_identifier]
    document = get_material_document(material_identifier)
    record = document["eos_records"][0]

    assert document["phase"] == expected["phase"]
    assert document["space_group"] == "Fm-3m"
    assert document["space_group_number"] == 225
    assert document["formula_units_per_cell"] == 4
    assert document["lattice"]["a"] == pytest.approx(expected["lattice"])
    assert record["identifier"] == expected["record"]
    assert record["reference"]["doi"] == "10.2138/am.2012.3937"
    assert record["eos"]["parameters"] == {
        "V0": expected["v0"],
        "K0": 160.0,
        "K0_prime": expected["kp"],
    }
    assert record["fixed_parameters"] == ["K0"]
    assert record["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 500.0,
        "gamma0": expected["gamma0"],
        "q": expected["q"],
        "n": 2.0,
    }
    assert record["thermal"]["fixed_parameters"] == ["Tr", "theta0", "n"]
    assert record["thermal"]["debye_temperature_law"] == ("integrated_gruneisen")
    assert record["thermal"]["thermal_pressure_reference"] == ("reference_temperature")
    assert record["validity"]["pressure_gpa"] == [0.0, 47.0]
    assert record["validity"]["temperature_k"] == [300.0, 1100.0]


@pytest.mark.parametrize("material_identifier", CASES)
def test_matsui_complete_tables_and_fit_selection(material_identifier):
    expected = CASES[material_identifier]
    document = get_material_document(material_identifier)
    dataset = document["datasets"][0]
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert dataset["identifier"] == expected["dataset"]
    assert dataset["resource"]["sha256"] == expected["sha256"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected["sha256"]
    assert len(rows) == expected["rows"]
    assert sum(row["fit_included"] == "1" for row in rows) == expected["included"]
    assert {
        temperature: sum(float(row["temperature_k"]) == temperature for row in rows)
        for temperature in (300.0, 700.0, 1100.0)
    } == expected["counts"]
    assert all(
        float(row["pressure_gpa"]) < 47.0 for row in rows if row["fit_included"] == "1"
    )
    assert all(
        float(row["pressure_gpa"]) >= 47.0 and row["exclusion_reason"]
        for row in rows
        if row["fit_included"] == "0"
    )


@pytest.mark.parametrize("material_identifier", CASES)
def test_matsui_reference_subtracted_curve_and_high_spin_validity(
    material_identifier,
):
    expected = CASES[material_identifier]
    document = get_material_document(material_identifier)
    loaded = Material.from_eosmat(
        document,
        record_identifiers=[expected["record"]],
    ).eos_records[0]
    data = load_data(material_identifier)

    assert loaded.thermal_pressure_increment(
        data["volume"],
        300.0,
    ) == pytest.approx(np.zeros(len(data)), abs=1.0e-14)
    loaded.volume(47.0, 300.0, check_validity=True)
    with pytest.raises(MaterialError, match="outside the published"):
        loaded.volume(47.01, 300.0, check_validity=True)


@pytest.mark.parametrize("material_identifier", CASES)
def test_matsui_published_curves_reproduce_printed_tables(material_identifier):
    expected = CASES[material_identifier]
    result = diagnostics(material_identifier)

    assert result.rows == expected["rows"]
    assert result.included_rows == expected["included"]
    assert result.included_observed_rmse_gpa == pytest.approx(
        expected["observed_rmse"], abs=1.0e-10
    )
    assert result.included_observed_max_abs_gpa == pytest.approx(
        expected["observed_max"], abs=1.0e-10
    )
    assert result.table_calculated_rmse_gpa == pytest.approx(
        expected["table_rmse"], abs=1.0e-10
    )
    assert result.table_calculated_max_abs_gpa == pytest.approx(
        expected["table_max"], abs=1.0e-10
    )
