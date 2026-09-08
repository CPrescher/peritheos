import csv
import importlib.util
from pathlib import Path

import pytest

from peritheos import get_material_document

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "reproduce_duffy_ahrens_1995_mgo_hugoniot.py"
SPEC = importlib.util.spec_from_file_location("duffy_ahrens_1995", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_duffy_ahrens_table3_transcription_and_hugoniot_selection():
    document = get_material_document("mgo")
    record = next(
        row
        for row in document["eos_records"]
        if row["identifier"] == "mgo_b1_duffy_ahrens_1995_hugoniot_5"
    )
    dataset = next(
        row
        for row in document["datasets"]
        if row["identifier"] == "mgo_duffy_ahrens_1995_table3_hugoniot"
    )
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert [row["shot_number"] for row in rows] == ["233", "840", "841", "843"]
    assert [float(row["particle_velocity_km_s"]) for row in rows] == [
        3.367,
        0.893,
        1.742,
        0.513,
    ]
    assert [float(row["shock_velocity_km_s"]) for row in rows] == [
        11.042,
        7.96,
        9.01,
        7.53,
    ]
    assert all(row["used_in_published_hugoniot_fit"] == "1" for row in rows)
    assert dataset["used_by_eos_records"] == [record["identifier"]]
    assert record["equation_kind"] == "hugoniot"
    assert record["branch_kind"] == "untransformed"
    assert "isotherm" in record["notes"]


def test_duffy_ahrens_uncertainty_aware_fit_reproduces_source_rounding():
    result = MODULE.reproduce()
    eiv = result["fits"]["errors_in_variables"]
    wls = result["fits"]["shock_velocity_wls"]

    assert eiv["parameters"]["c0"] == pytest.approx(6.87043987023)
    assert eiv["parameters"]["s"] == pytest.approx(1.23798631265)
    assert eiv["standard_errors"]["c0"] == pytest.approx(0.100349188)
    assert eiv["standard_errors"]["s"] == pytest.approx(0.041203119)
    assert wls["parameters"]["c0"] == pytest.approx(6.87029565579)
    assert wls["parameters"]["s"] == pytest.approx(1.23803732073)
    assert result["published_rounding_reproduced"] == {"c0": True, "s": True}
    assert result["rankine_hugoniot_checks"]["max_abs_pressure_difference_gpa"] < 0.1
    assert result["rankine_hugoniot_checks"]["max_abs_density_difference_g_cm3"] < 0.001
