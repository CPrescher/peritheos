"""Primary-source and refit checks for Fischer et al. (2011) FeO."""

import csv
from importlib import resources

import pytest

from peritheos import get_material_document
from peritheos.materials import Material
from scripts import validate_primary_eos_refits as refits

B1_ID = "feo_fischer_2011_bm3_2"
B8_ID = "feo_b8_2_fischer_2011_bm3_1"


def _record(material, identifier):
    document = get_material_document(material)
    return document, next(
        row for row in document["eos_records"] if row["identifier"] == identifier
    )


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


def test_fischer_conditional_refits_use_all_numerical_phase_rows():
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
    assert outcomes[B8_ID]["status"] == "parity_not_achieved"
    assert outcomes[B8_ID]["observations"] == 21
    assert [p["refit"] for p in outcomes[B8_ID]["parameters"]] == pytest.approx(
        [37.2814998117, 183.8222355571, 1.9961519777]
    )
