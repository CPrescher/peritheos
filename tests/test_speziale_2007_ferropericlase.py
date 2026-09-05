import csv
import hashlib
from importlib import resources

import pytest

from peritheos import get_material_document
from peritheos.materials import Material, MaterialError
from scripts.reproduce_speziale_2007_ferropericlase import reproduce

DATASET_ID = "mg080fe020o_speziale_2007_table1_pv"
RECORD_ID = "mg080fe020o_speziale_2007_high_spin_bm3_1"
CHECKSUM = "0b7c788f9d5ade52431338aec1da6745c7f528c16e551470de953ed9823f7a42"


def _load():
    document = get_material_document("mg080fe020o")
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_ID
    )
    loaded = next(
        record
        for record in Material.from_eosmat(document).eos_records
        if record.identifier == RECORD_ID
    )
    dataset = next(
        dataset
        for dataset in document["datasets"]
        if dataset["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, source, loaded, dataset, path, rows


def test_speziale_high_spin_record_preserves_source_parameters():
    document, source, _, _, _, _ = _load()

    assert document["formula"] == "Mg0.80Fe0.20O"
    assert document["phase"] == "high-spin B1 ferropericlase"
    assert document["space_group"] == "Fm-3m"
    assert document["formula_units_per_cell"] == 4
    assert source["reference"]["doi"] == "10.1029/2006JB004730"
    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 76.03, "K0": 158.0, "K0_prime": 4.4},
    }
    assert source["parameter_errors"] == {
        "V0": 0.09,
        "K0": 3.0,
        "K0_prime": 0.2,
    }
    assert source["fixed_parameters"] == []
    assert source["fit_datasets"] == [DATASET_ID]


def test_speziale_table1_transcription_is_complete_and_checksummed():
    _, _, _, dataset, path, rows = _load()

    assert hashlib.sha256(path.read_bytes()).hexdigest() == CHECKSUM
    assert dataset["resource"]["sha256"] == CHECKSUM
    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert len(rows) == 30
    assert [int(row["source_row"]) for row in rows] == list(range(1, 31))
    assert sum(row["fit_included"] == "1" for row in rows) == 18
    assert sum(row["fit_included"] == "0" for row in rows) == 12
    assert (
        max(float(row["pressure_gpa"]) for row in rows if row["fit_included"] == "1")
        == 39.5
    )
    assert (
        min(float(row["pressure_gpa"]) for row in rows if row["fit_included"] == "0")
        == 42.1
    )


def test_speziale_published_curve_matches_table_within_uncertainty():
    result = reproduce()

    assert result["observations_in_source_table"] == 30
    assert result["observations_in_high_spin_fit"] == 18
    assert result["excluded_spin_transition_observations"] == 12
    assert result["published_curve_pressure_rmse_gpa"] == pytest.approx(
        1.2608936698, abs=2.0e-10
    )
    assert result["published_curve_effective_sigma_rms"] == pytest.approx(
        0.7034358017, abs=2.0e-10
    )
    assert result["published_curve_max_abs_effective_sigma"] < 1.65


def test_speziale_record_is_executable_only_on_high_spin_branch():
    _, source, loaded, _, _, _ = _load()

    pressure = loaded.pressure(64.08, 300.0, check_validity=True)
    assert pressure == pytest.approx(39.4326808828, abs=2.0e-9)
    assert loaded.volume(pressure, 300.0) == pytest.approx(64.08)
    with pytest.raises(MaterialError, match="outside the published"):
        loaded.volume(42.0, 300.0, check_validity=True)

    notes = " ".join(source["validity"]["notes"])
    assert "continuous-spin-crossover" in notes
    assert source["pressure_calibration"]["status"] == "unresolved"
