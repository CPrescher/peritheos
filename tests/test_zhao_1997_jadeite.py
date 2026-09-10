import csv
from pathlib import Path

import pytest

from peritheos import Material, get_material_document
from peritheos.eos.thermal import ThermalReferenceStateEOS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/jadeite-zhao-1997-table1-pvt.csv"


def test_zhao_record_is_the_complete_preferred_thermal_parameterization():
    document = get_material_document("naalsi2o6")
    record = document["eos_records"][0]

    assert record["label"] == "Zhao et al. (1997), jadeite thermal BM3"
    assert record["eos"]["parameters"] == {
        "V0": 403.0,
        "K0": 124.5,
        "K0_prime": 5.0,
    }
    assert record["fixed_parameters"] == ["V0", "K0_prime"]
    assert record["parameter_errors"] == {
        "V0": None,
        "K0": 4.0,
        "K0_prime": None,
    }
    assert record["thermal"] == {
        "type": "AlphaKT",
        "model": "thermal_reference_state",
        "thermal_expansion_law": "linear_temperature",
        "reference_volume_law": "integrated_expansivity",
        "parameters": {
            "Tr": 300.0,
            "alpha0": 2.56e-5,
            "alpha1": 2.6e-9,
            "dK_dT": -0.0165,
        },
        "parameter_errors": {
            "Tr": None,
            "alpha0": 2.2e-6,
            "alpha1": 1.8e-9,
            "dK_dT": 0.0049,
        },
        "fixed_parameters": ["Tr"],
    }
    assert record["experimental_pressure_range_gpa"] == [0.0, 8.16]
    assert record["experimental_temperature_range_k"] == [300.0, 1280.0]
    assert record["fit_provenance"]["selection"] == {
        "dataset": "jadeite_zhao_1997_table1_pvt",
        "included_rows": 31,
        "excluded_bundled_rows": 0,
        "source_preselection": (
            "Only observations collected after heating removed the cold-compression "
            "deviatoric stress were admitted to Table 1"
        ),
    }
    assert record["fit_provenance"]["source_objective"] == "not_reported"
    assert record["fit_provenance"]["source_weights"] == "not_reported"
    assert record["fit_provenance"]["source_covariance"] == "not_reported"
    assert record["fit_provenance"]["software"] == {
        "name": "not reported",
        "version": "not reported",
    }
    assert record["fit_provenance"]["refined_parameters"] == [
        "K0",
        "alpha0",
        "alpha1",
        "dK_dT",
    ]
    assert record["fit_provenance"]["fixed_parameters"] == [
        "V0",
        "K0_prime",
        "Tr",
    ]
    assert record["fit_provenance"]["statistics"] == {"observations": 31}


def test_zhao_table_is_complete_and_uses_no_post_bundling_exclusions():
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 31
    assert [float(row["pressure_gpa"]) for row in rows][::6] == pytest.approx(
        [0.0, 6.33, 4.73, 3.17, 1.77, 0.41]
    )
    assert min(float(row["temperature_k"]) for row in rows) == 300.0
    assert max(float(row["temperature_k"]) for row in rows) == 1280.0
    assert min(float(row["unit_cell_volume_a3"]) for row in rows) == 385.41
    assert max(float(row["unit_cell_volume_a3"]) for row in rows) == 406.70


def test_zhao_equation_one_is_executable_over_pressure_and_temperature():
    material = Material.from_eosmat(get_material_document("naalsi2o6"))
    eos = material.eos_records[0].eos

    assert isinstance(eos, ThermalReferenceStateEOS)
    assert eos.pressure(403.0, 300.0) == pytest.approx(0.0, abs=1e-14)
    assert eos.pressure(390.0, 1000.0) == pytest.approx(6.674297363244782)
    assert eos.pressure(388.54, 1280.0) == pytest.approx(8.085208941433727)


def test_zhao_pressure_calibration_is_identified_but_not_recalculable():
    record = get_material_document("naalsi2o6")["eos_records"][0]
    calibration = record["pressure_calibration"]

    assert calibration["status"] == "resolved"
    method = calibration["methods"][0]
    assert method["material"] == "NaCl"
    assert method["reference"]["doi"] == "10.1063/1.1660714"
    assert calibration["recalculation"]["status"] == ("missing_calibrant_observations")
