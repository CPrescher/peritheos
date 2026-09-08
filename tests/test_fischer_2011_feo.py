"""Primary-source and refit checks for Fischer et al. (2011) FeO."""

import csv
from importlib import resources

import pytest

from peritheos import get_material_document
from peritheos.materials import Material
from scripts import validate_primary_eos_refits as refits

B1_ID = "feo_fischer_2011_bm3_2"
B8_ID = "feo_b8_2_fischer_2011_bm3_1"
SEAGLE_DATASET_ID = "feo_seagle_2008_supplement_volume_temperature"
OZAWA_DATASET_ID = "feo_ozawa_2010_table2_pvt"
CAMPBELL_DATASET_ID = "feo_campbell_2009_table_s2_pvt"


def _record(material, identifier):
    document = get_material_document(material)
    return document, next(
        row for row in document["eos_records"] if row["identifier"] == identifier
    )


def _dataset(document, identifier):
    return next(row for row in document["datasets"] if row["identifier"] == identifier)


@pytest.mark.parametrize(
    ("material", "identifier", "gamma0", "q"),
    [("feo", B1_ID, 1.41, 0.5), ("feo_b8_2", B8_ID, 1.73, 1.0)],
)
def test_fischer_records_encode_published_mgd(material, identifier, gamma0, q):
    document, record = _record(material, identifier)
    thermal = record["thermal"]
    assert thermal["type"] == "MieGruneisenDebye"
    assert thermal["debye_temperature_law"] == "integrated_gruneisen"
    assert thermal["parameters"] == {
        "Tr": 300.0,
        "theta0": 417.0,
        "gamma0": gamma0,
        "q": q,
        "n": 2.0,
    }
    assert thermal["fixed_parameters"] == ["Tr", "theta0", "q", "n"]
    assert record["pressure_calibration"]["methods"][0]["material"] == "hcp Fe"
    loaded = Material.from_eosmat(document, record_identifiers=[identifier])
    eos = loaded.eos_records[0]
    volume = 55.0 if material == "feo" else 26.0
    assert eos.pressure(volume, 2000.0, check_validity=False) > (
        eos.pressure(volume, 300.0, check_validity=False) + 5.0
    )


def test_fischer_table_s1_phase_rows_and_b8_fallback_are_explicit():
    document, _ = _record("feo_b8_2", B8_ID)
    dataset = document["datasets"][0]
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    b1 = [row for row in rows if row["b1_feo_molar_volume_cm3_mol"]]
    b8 = [row for row in rows if row["b8_feo_molar_volume_cm3_mol"]]
    unconstrained = [
        row
        for row in b8
        if row["b8_observation"] == "measured_two_peaks_uncertainty_unconstrained"
    ]
    assert (len(rows), len(b1), len(b8), len(unconstrained)) == (52, 42, 21, 3)
    assert min(float(row["temperature_k"]) for row in b1) == 1143.9
    assert {
        row["b8_feo_molar_volume_uncertainty_cm3_mol"] for row in unconstrained
    } == {""}
    assert {
        float(row["b8_feo_molar_volume_fit_uncertainty_cm3_mol"])
        for row in unconstrained
    } == {0.1}
    assert all(
        not row["b8_feo_molar_volume_cm3_mol"]
        for row in rows
        if row["b8_observation"] == "detected_one_peak_no_volume"
    )


def test_fischer_refits_use_recovered_numerical_phase_rows():
    outcomes = {}
    for material, identifier in (("feo", B1_ID), ("feo_b8_2", B8_ID)):
        document, record = _record(material, identifier)
        outcomes[identifier] = refits._fit_record(
            document, record, document["datasets"][0]
        )
    assert outcomes[B1_ID]["status"] == "parity"
    assert outcomes[B1_ID]["observations"] == 42
    assert [p["refit"] for p in outcomes[B1_ID]["parameters"]] == pytest.approx(
        [144.5079729305, 3.7156235852, 1.6301955006]
    )
    assert outcomes[B8_ID]["status"] == "parity"
    assert outcomes[B8_ID]["observations"] == 29
    assert [p["refit"] for p in outcomes[B8_ID]["parameters"]] == pytest.approx(
        [39.8745039519, 137.8628936218, 1.6661584486]
    )


def test_seagle_supplement_recovers_all_subsolidus_calibrant_rows():
    document, _ = _record("feo", B1_ID)
    dataset = _dataset(document, SEAGLE_DATASET_ID)
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 81
    assert sum(row["iron_phase"] == "hcp" for row in rows) == 65
    assert sum(row["iron_phase"] == "fcc" for row in rows) == 14
    assert sum(not row["iron_phase"] for row in rows) == 2
    assert sum(
        row["pressure_recalculation_status"]
        == "ready_seagle_table1_series_pressure"
        for row in rows
    ) == 2
    assert rows[48] == {
        "source_row": "49",
        "feo_unit_cell_volume_a3": "60.676",
        "feo_unit_cell_volume_uncertainty_a3": "0.048",
        "iron_unit_cell_volume_a3": "17.729",
        "iron_unit_cell_volume_uncertainty_a3": "0.133",
        "iron_phase": "hcp",
        "temperature_k": "2210",
        "source_series_pressure_gpa": "",
        "source_series_pressure_uncertainty_gpa": "",
        "pressure_recalculation_status": "ready_dewaele_2006_hcp_fe",
    }
    assert {
        (
            row["source_row"],
            row["feo_unit_cell_volume_a3"],
            row["temperature_k"],
            row["source_series_pressure_gpa"],
            row["source_series_pressure_uncertainty_gpa"],
        )
        for row in rows
        if row["pressure_recalculation_status"]
        == "ready_seagle_table1_series_pressure"
    } == {
        ("61", "68.480", "2580", "50", "4"),
        ("62", "68.404", "2590", "50", "4"),
    }


def test_ozawa_printed_table_supplies_b1_and_b8_fit_candidates():
    b1_document, _ = _record("feo", B1_ID)
    b8_document, _ = _record("feo_b8_2", B8_ID)
    b1_dataset = _dataset(b1_document, OZAWA_DATASET_ID)
    b8_dataset = _dataset(b8_document, OZAWA_DATASET_ID)
    assert b1_dataset["resource"] == b8_dataset["resource"]

    path = resources.files("peritheos.data").joinpath(
        b1_dataset["resource"]["path"]
    )
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 20
    assert sum(row["phase"] == "B1" for row in rows) == 12
    assert sum(row["phase"] == "B8" for row in rows) == 8

    iron = Material.from_eosmat(
        get_material_document("iron"),
        record_identifiers=["iron_dewaele_2006_vinet_thermal"],
    ).eos_records[0]
    pressure_residuals = [
        iron.pressure(
            float(row["iron_hcp_unit_cell_volume_a3"]),
            float(row["temperature_k"]),
        )
        - float(row["pressure_gpa"])
        for row in rows
    ]
    assert max(abs(value) for value in pressure_residuals) < 3.0


def test_campbell_supplement_recovers_complete_fe_feo_table():
    document, _ = _record("feo", B1_ID)
    dataset = _dataset(document, CAMPBELL_DATASET_ID)
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 25
    assert sum(row["iron_phase"] == "fcc" for row in rows) == 15
    assert sum(row["iron_phase"] == "hcp" for row in rows) == 10
    assert all(row["nacl_phase"] in {"B1", "B2"} for row in rows)
    assert rows[0]["run_id"] == "FeFeO_030"
    assert float(rows[0]["feo_molar_volume_cm3_mol"]) == pytest.approx(
        11.4433538834759
    )
    assert float(rows[-1]["reported_pressure_gpa"]) == pytest.approx(
        55.1248168632253
    )
