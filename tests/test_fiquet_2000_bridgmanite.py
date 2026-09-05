import csv
import hashlib
from importlib import resources

import pytest

from peritheos import get_material_document
from peritheos.materials import Material, MaterialError
from scripts.reproduce_fiquet_2000_bridgmanite import reproduce

RECORD_ID = "bridgmanite_fiquet_2000_bm3_1"
DATASET_ID = "bridgmanite_fiquet_2000_table1_298k_pv"
CHECKSUM = "92e3a39df52f6f5a30a60667d446e8415bb5965f6266c6e2cc51ac37063929e1"


def _load():
    document = get_material_document("bridgmanite")
    source = next(r for r in document["eos_records"] if r["identifier"] == RECORD_ID)
    loaded = next(
        r
        for r in Material.from_eosmat(document).eos_records
        if r.identifier == RECORD_ID
    )
    dataset = next(d for d in document["datasets"] if d["identifier"] == DATASET_ID)
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return source, loaded, dataset, path, rows


def test_fiquet_record_preserves_source_bm3():
    source, _, _, _, _ = _load()
    assert source["reference"]["doi"] == "10.1029/1999GL008397"
    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 162.27, "K0": 253.0, "K0_prime": 3.9},
    }
    assert source["parameter_errors"] == {
        "V0": 0.01,
        "K0": 9.0,
        "K0_prime": 0.2,
    }
    assert source["fixed_parameters"] == []
    assert source["fit_datasets"] == [DATASET_ID]


def test_fiquet_298k_table_subset_is_exact_and_checksummed():
    _, _, dataset, path, rows = _load()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == CHECKSUM
    assert dataset["resource"]["sha256"] == CHECKSUM
    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert len(rows) == 25
    assert {row["temperature_k"] for row in rows} == {"298"}
    assert all(row["fit_included"] == "1" for row in rows)
    assert {row["table_block"] for row in rows} == {"left", "right"}


def test_fiquet_refit_recovers_reported_moduli():
    result = reproduce()
    assert result["observations"] == 25
    assert result["published_curve_pressure_rmse_gpa"] == pytest.approx(
        0.7070408617, abs=2.0e-10
    )
    assert result["fixed_v0_refit"] == pytest.approx(
        {
            "V0": 162.27,
            "K0": 253.3244173577,
            "K0_prime": 3.9161902310,
            "pressure_rmse_gpa": 0.6837616952,
        },
        abs=2.0e-9,
    )
    assert abs(result["fixed_v0_refit"]["K0"] - 253.0) < 9.0
    assert abs(result["fixed_v0_refit"]["K0_prime"] - 3.9) < 0.2


def test_fiquet_isotherm_is_executable_only_in_source_range():
    source, loaded, _, _, _ = _load()
    pressure = loaded.pressure(140.827, 298.0, check_validity=True)
    assert pressure == pytest.approx(47.2723242731, abs=2.0e-9)
    assert loaded.volume(pressure, 298.0) == pytest.approx(140.827)
    with pytest.raises(MaterialError, match="outside the published"):
        loaded.volume(20.0, 298.0, check_validity=True)
    assert "not either composite thermal inversion" in " ".join(
        source["validity"]["notes"]
    )
