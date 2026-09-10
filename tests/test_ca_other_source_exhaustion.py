import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from peritheos import Material, validate_eosmat_document

ROOT = Path(__file__).parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
DATA = ROOT / "peritheos/data/datasets"


def _document(name: str):
    document = json.loads((MATERIALS / name).read_text(encoding="utf-8"))
    validate_eosmat_document(document)
    return document


def _source_record(name: str, identifier: str):
    document = _document(name)
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )
    executable = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).get_eos_record(identifier)
    return source, executable


def test_karki_1998_cubic_casio3_source_record():
    source, executable = _source_record(
        "ca_perovskite.eosmat", "ca_perovskite_karki_crain_1998_static_bm3"
    )
    assert source["reference"]["doi"] == "10.1029/98GL51952"
    assert source["eos"]["parameters"] == {
        "V0": 45.34,
        "K0": 241.0,
        "K0_prime": 4.14,
    }
    assert source["temperature_ref"] == 0.0
    assert executable.pressure(executable.volume(140.0)) == pytest.approx(140.0)


def test_shim_2002_tetragonal_casio3_source_record():
    source, executable = _source_record(
        "casio3_perovskite_tetragonal.eosmat",
        "casio3_perovskite_tetragonal_shim_2002_bm3_1",
    )
    assert source["reference"]["doi"] == "10.1029/2002GL016148"
    assert source["fixed_parameters"] == ["V0", "K0_prime"]
    assert source["eos"]["parameters"] == {
        "V0": 45.58,
        "K0": 255.0,
        "K0_prime": 4.0,
    }
    assert source["parameter_errors"]["K0"] == 5.0
    assert executable.pressure(executable.volume(45.8)) == pytest.approx(45.8)


def test_ono_2013_thermal_extension_keeps_source_lineage_separate():
    source, executable = _source_record(
        "ca_perovskite.eosmat", "ca_perovskite_ono_2013_bm3_log_thermal"
    )
    assert source["reference"]["doi"] == "10.3390/e15104300"
    assert source["thermal"]["model"] == "log_volume_thermal_pressure"
    assert source["thermal"]["parameters"] == {
        "Tr": 300.0,
        "alpha_KT_ref": 0.0083,
        "dK_dT_V": -0.0031,
    }
    assert source["fixed_parameters"] == ["V0", "K0", "K0_prime"]
    assert [entry["doi"] for entry in source["source_lineage"]] == [
        "10.3390/e15104300",
        "10.1016/S0031-9201(00)00154-0",
    ]
    assert executable.eos.thermal_pressure(40.0, 2000.0) == pytest.approx(
        13.421792967662018
    )


def test_sun_2019_fesio3_liquid_reference_and_table():
    source, executable = _source_record(
        "fesio3_liquid.eosmat", "fesio3_liquid_sun_2019_2500k_bm4_1"
    )
    assert source["reference"]["doi"] == "10.1029/2018GL081421"
    assert source["eos"]["parameters"] == {
        "V0": 78.46047092,
        "K0": 2.217,
        "K0_prime": 22.913,
        "K0_double_prime": -175.628,
    }
    assert source["fit_datasets"] == ["fesio3_liquid_sun_2019_table1_pvt"]
    classification = source["numerical_reproduction"]["state_classification"]
    assert classification["liquid_circle_states"] == 40
    assert classification["nonliquid_square_states"] == 6
    volume_a3 = 40.72 * 1.0e24 / 6.02214076e23
    assert executable.pressure(volume_a3) == pytest.approx(1.084850164170989)

    path = DATA / "fesio3-liquid-sun-2019-table1-pvt.csv"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == (
        "f0f4ce1eb820fdb91d2c988baeb367c25f808028b99ed7c0040ed1ce75b335f7"
    )
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 46
    assert sum(row["phase_state"] == "liquid" for row in rows) == 40
    assert all(
        (row["phase_state"] == "liquid") == (row["used_in_published_eos_fit"] == "yes")
        for row in rows
    )
    assert {
        (float(row["volume_ratio_to_vx"]), int(row["temperature_k"]))
        for row in rows
        if row["phase_state"] == "nonliquid"
    } == {
        (0.4, 2500),
        (0.4, 3000),
        (0.4, 4000),
        (0.5, 2500),
        (0.5, 3000),
        (0.6, 2500),
    }
    checkpoint = next(
        row
        for row in rows
        if row["volume_ratio_to_vx"] == "1.0" and row["temperature_k"] == "2500"
    )
    assert float(checkpoint["pressure_gpa"]) == 1.07
    assert float(checkpoint["pressure_standard_error_gpa"]) == 0.12


def test_source_exhaustion_reproduction_script():
    path = ROOT / "scripts/reproduce_ca_other_source_exhaustion.py"
    spec = importlib.util.spec_from_file_location("reproduce_ca_other_exhaustion", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["sun_table1_states"] == 46
    assert result["sun_state_classification"] == {
        "liquid_states_used_in_published_eos_fit": 40,
        "nonliquid_states_excluded": 6,
        "nonliquid_volume_ratio_temperature_pairs": [
            [0.4, 2500],
            [0.4, 3000],
            [0.4, 4000],
            [0.5, 2500],
            [0.5, 3000],
            [0.6, 2500],
        ],
    }
    assert result["sun_published_rounded_model_liquid_state_rmse_gpa"] == (
        pytest.approx(0.6409268656)
    )
    assert (
        result["sun_2500k_table1_checkpoint"]["absolute_difference_gpa"]
        < result["sun_2500k_table1_checkpoint"]["published_pressure_standard_error_gpa"]
    )
    assert result["ono_2000k_checkpoint"]["thermal_pressure_gpa"] == pytest.approx(
        result["ono_2000k_checkpoint"]["expected_thermal_pressure_gpa"]
    )
