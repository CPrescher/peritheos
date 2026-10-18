"""Keep published MgO evaluation distinct from historical Au thermal ancestry."""

import csv
from pathlib import Path

import pytest

from peritheos import (
    find_pressure_calibration_path,
    get_cross_calibration_edge,
    get_eos_record,
    get_eos_record_document,
)
from peritheos.eos.thermal import VariableQMieGruneisenDebye
from peritheos.errors import ValidationError

MGO = "mgo_speziale_2001_bm3_2"
THERMAL_AU = (
    "gold_fei_2004_bm3_mgd",
    "gold_fei_2007_vinet_2",
    "gold_hirose_2008_bm3_fit1",
    "gold_hirose_2008_bm3_fit2",
)
DATA = Path(__file__).resolve().parents[1] / "peritheos/data/datasets"


@pytest.mark.parametrize("identifier", THERMAL_AU)
def test_unrecovered_thermal_reducers_do_not_generate_executable_ancestry(identifier):
    calibration = get_eos_record_document(identifier)["pressure_calibration"]
    assert calibration["status"] == "partially_resolved"
    assert calibration["recalculation"]["status"] == "not_possible"
    mgo_methods = [
        method
        for method in calibration["methods"]
        if method.get("reference", {}).get("doi") == "10.1029/2000JB900318"
    ]
    assert len(mgo_methods) == 1
    # These fields generate executable pressure-identity edges regardless of
    # recalculation status. Published-model availability is insufficient here.
    assert "reference_eos_record" not in mgo_methods[0]
    assert "reference_calibration_record" not in mgo_methods[0]


def test_cold_ancestry_remains_executable_and_historical_thermal_edge_is_evidence():
    cold = get_eos_record_document("gold_hirose_2008_bm3_300k")
    assert cold["pressure_calibration"]["status"] == "resolved"
    assert cold["pressure_calibration"]["recalculation"]["status"] == "ready"
    assert find_pressure_calibration_path(cold["identifier"], MGO)
    fei = get_eos_record_document("gold_fei_2007_vinet_2")
    assert (
        fei["pressure_calibration"]["methods"][0]["reference_calibration_record"]
        == "ruby_dewaele_2004"
    )
    thermal = get_cross_calibration_edge(
        "gold_hirose_2008_to_mgo_speziale_2001_thermal"
    )
    assert thermal["executable"] is False
    assert thermal["target_node"] == "mgo_speziale_2001_thermal"
    with pytest.raises(ValidationError, match="No executable"):
        find_pressure_calibration_path(thermal["source_node"], thermal["target_node"])


def test_executable_published_mgo_does_not_reproduce_later_hot_pressure_targets():
    record = get_eos_record(MGO)
    assert isinstance(record.eos, VariableQMieGruneisenDebye)
    with (DATA / "gold-hirose-2006-table1.csv").open() as stream:
        rows = [
            row
            for row in csv.DictReader(stream)
            if row["use_in_fei2007_thermal_replay"] == "true"
        ]
    reported = [float(row["mgo_speziale_pressure_gpa"]) for row in rows]
    predicted = [
        record.pressure(float(row["mgo_volume_a3"]), float(row["temperature_k"]))
        for row in rows
    ]
    assert reported == [106.5, 117.9]
    assert predicted == pytest.approx([105.609013, 115.947872], abs=1e-6)
    assert all(abs(a - b) > 0.8 for a, b in zip(reported, predicted))
