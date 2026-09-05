import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MATERIAL = (
    ROOT / "peritheos" / "data" / "materials" / "mg05fe05al05si05o3_bridgmanite.eosmat"
)
DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "fa50-bridgmanite-koemets-2023-supplement-table3.csv"
)
SCRIPT = ROOT / "scripts" / "reproduce_koemets_2023_fa50_bridgmanite.py"
RECORD_ID = "mg05fe05al05si05o3_bridgmanite_koemets_2023_bm2_4"
DATASET_ID = "mg05fe05al05si05o3_bridgmanite_koemets_2023_table_s3_pv"


def _load_reproduction_module():
    spec = importlib.util.spec_from_file_location("koemets_reproduction", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_koemets_record_and_dataset_metadata_are_linked():
    document = json.loads(MATERIAL.read_text(encoding="utf-8"))
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )

    assert record["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 171.2, "K0": 221.0},
    }
    assert record["fit_datasets"] == [DATASET_ID]
    assert record["experimental_pressure_range_gpa"] == [10.8, 60.0]
    assert record["scientific_validation"]["independent_refit"]["result"] == ("parity")
    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert (
        dataset["resource"]["sha256"]
        == hashlib.sha256(DATASET.read_bytes()).hexdigest()
    )


def test_koemets_table_s3_selection_and_reproduction():
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    selected = [row for row in rows if row["used_in_bm2_fit"] == "1"]
    excluded = [row for row in rows if row["used_in_bm2_fit"] == "0"]
    assert len(rows) == 18
    assert len(selected) == 16
    assert [float(row["pressure_gpa"]) for row in excluded] == [7.5, 8.6]

    result = _load_reproduction_module().reproduce()
    assert result["selected_pressure_range_gpa"] == [10.8, 60.0]
    assert result["published_curve_pressure_rmse_gpa"] == pytest.approx(0.6776309819)
    assert result["errors_in_variables_refit"]["parameters"] == pytest.approx(
        {"V0": 171.3069453112, "K0": 220.0083064600}
    )
