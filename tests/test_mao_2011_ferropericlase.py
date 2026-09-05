import csv
import hashlib
from collections import Counter
from importlib import resources

import pytest

from peritheos import get_material_document
from peritheos.materials import Material, MaterialError
from scripts.reproduce_mao_2011_ferropericlase import reproduce_branch

DATASET_ID = "mg075fe025o_mao_2011_figure1_300k_digitized"
HIGH_SPIN_ID = "mg075fe025o_mao_2011_high_spin_bm3_1"
LOW_SPIN_ID = "mg075fe025o_mao_2011_low_spin_bm3_1"
CHECKSUM = "dd8a928c47ceb80c7c595e92fedbe00b223bfcd5eeb76daf17fd4e7100455eac"


def _load():
    document = get_material_document("mg075fe025o")
    sources = {record["identifier"]: record for record in document["eos_records"]}
    loaded = {
        record.identifier: record
        for record in Material.from_eosmat(document).eos_records
    }
    dataset = next(
        dataset
        for dataset in document["datasets"]
        if dataset["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, sources, loaded, dataset, path, rows


def test_mao_spin_branches_preserve_corrected_source_parameters():
    document, sources, _, _, _, _ = _load()
    high_spin = sources[HIGH_SPIN_ID]
    low_spin = sources[LOW_SPIN_ID]

    assert document["phase"] == "B1 ferropericlase"
    assert document["space_group"] == "Fm-3m"
    assert document["space_group_number"] == 225
    assert document["formula_units_per_cell"] == 4
    assert high_spin["reference"]["doi"] == "10.1029/2011GL049915"
    assert high_spin["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 76.34, "K0": 162.0, "K0_prime": 4.0},
    }
    assert high_spin["parameter_errors"] == {
        "V0": 0.01,
        "K0": 1.0,
        "K0_prime": None,
    }
    assert high_spin["fixed_parameters"] == ["V0", "K0_prime"]
    assert low_spin["eos"]["parameters"] == {
        "V0": 74.4,
        "K0": 166.0,
        "K0_prime": 4.0,
    }
    assert low_spin["parameter_errors"] == {
        "V0": 0.6,
        "K0": 7.0,
        "K0_prime": None,
    }
    assert low_spin["fixed_parameters"] == ["K0_prime"]
    assert high_spin["fit_datasets"] == [DATASET_ID]
    assert low_spin["fit_datasets"] == [DATASET_ID]


def test_mao_figure1_digitization_is_complete_and_checksummed():
    _, _, _, dataset, path, rows = _load()

    assert hashlib.sha256(path.read_bytes()).hexdigest() == CHECKSUM
    assert dataset["resource"]["sha256"] == CHECKSUM
    assert dataset["used_by_eos_records"] == [HIGH_SPIN_ID, LOW_SPIN_ID]
    assert len(rows) == 42
    assert [int(row["source_order"]) for row in rows] == list(range(1, 43))
    assert Counter(row["spin_regime"] for row in rows) == {
        "high_spin": 15,
        "mixed_spin": 9,
        "low_spin": 18,
    }
    assert sum(row["used_in_high_spin_fit"] == "1" for row in rows) == 15
    assert sum(row["used_in_low_spin_fit"] == "1" for row in rows) == 18
    assert all(
        row["used_in_high_spin_fit"] == "0" and row["used_in_low_spin_fit"] == "0"
        for row in rows
        if row["spin_regime"] == "mixed_spin"
    )


def test_mao_digitized_refits_recover_both_branch_parameterizations():
    high_spin = reproduce_branch("high_spin")
    low_spin = reproduce_branch("low_spin")

    assert high_spin["observations"] == 15
    assert high_spin["refit_parameters"] == pytest.approx(
        {"V0": 76.34, "K0": 163.0460834463, "K0_prime": 4.0}, abs=2.0e-9
    )
    assert high_spin["published_curve_pressure_rmse_gpa"] == pytest.approx(
        0.2491706696, abs=2.0e-10
    )
    assert abs(high_spin["refit_parameters"]["K0"] - 162.0) < 1.1

    assert low_spin["observations"] == 18
    assert low_spin["refit_parameters"] == pytest.approx(
        {"V0": 74.47146131, "K0": 165.6618754184, "K0_prime": 4.0},
        abs=2.0e-9,
    )
    assert low_spin["published_curve_pressure_rmse_gpa"] == pytest.approx(
        0.4942818332, abs=2.0e-10
    )
    assert abs(low_spin["refit_parameters"]["V0"] - 74.4) < 0.6
    assert abs(low_spin["refit_parameters"]["K0"] - 166.0) < 7.0


def test_mao_branch_validity_excludes_the_mixed_spin_interval():
    _, sources, loaded, _, _, _ = _load()
    high_spin = loaded[HIGH_SPIN_ID]
    low_spin = loaded[LOW_SPIN_ID]

    high_pressure = high_spin.pressure(64.12, 300.0, check_validity=True)
    assert high_pressure == pytest.approx(40.0797376217, abs=2.0e-9)
    assert high_spin.volume(high_pressure, 300.0) == pytest.approx(64.12)
    with pytest.raises(MaterialError, match="outside the published"):
        high_spin.volume(50.01, 300.0, check_validity=True)

    low_pressure = low_spin.pressure(54.25, 300.0, check_validity=True)
    assert low_pressure == pytest.approx(98.7964769732, abs=2.0e-9)
    assert low_spin.volume(low_pressure, 300.0) == pytest.approx(54.25)
    with pytest.raises(MaterialError, match="outside the published"):
        low_spin.volume(76.0, 300.0, check_validity=True)

    for identifier in (HIGH_SPIN_ID, LOW_SPIN_ID):
        notes = " ".join(sources[identifier]["validity"]["notes"])
        assert "continuous spin-crossover model" in notes
