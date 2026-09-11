from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

import pytest

from peritheos import get_material_document

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/audit_schoelmerich_2020_stishovite.py"
AUDIT = ROOT / "docs/data/schoelmerich-2020-stishovite-audit.json"


def _load_script():
    spec = importlib.util.spec_from_file_location("schoelmerich_2020_audit", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_table1_is_preserved_without_a_production_eos():
    document = get_material_document("sio2_stv_andr")
    records = {item["identifier"] for item in document["eos_records"]}
    assert "stishovite_schoelmerich_2020_shock_300k_bm3" not in records

    datasets = {item["identifier"]: item for item in document["datasets"]}
    assert "stishovite_schoelmerich_2020_figure3_corrected" not in datasets
    assert "stishovite_schoelmerich_2020_table1_shock" not in datasets
    path = (
        ROOT / "peritheos/data/datasets/stishovite-schoelmerich-2020-table1-shock.csv"
    )
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 6
    assert hashlib.sha256(path.read_bytes()).hexdigest() == (
        "72a73c8812b01c11ffb4dea7958459f8b4e782f5d7f83acf7793c3198be19726"
    )
    assert Counter(row["state_origin"] for row in rows) == {
        "ambient_xrd": 1,
        "xrd_refined_shock_state": 4,
        "visar_rankine_hugoniot_only": 1,
    }


def test_deterministic_audit_records_rejection_and_source_checks(assert_audit_close):
    result = _load_script().audit()
    assert_audit_close(
        result, json.loads(AUDIT.read_text(encoding="utf-8")), rel=1e-12, abs=1e-12
    )

    table = result["reported_table"]
    assert table["rankine_hugoniot_pressure_check"][
        "max_abs_difference_gpa"
    ] == pytest.approx(0.8314)
    energy = table["internal_energy_unit_audit"]
    assert energy["rmse_treating_reported_values_as_mj_kg"] < 0.1
    assert energy["rmse_treating_reported_values_as_kj_mol"] > 500.0

    decision = result["catalog_decision"]
    assert decision["outcome"] == "direct_refit_unavailable"
    assert decision["production_eos_record"] is False
    assert (
        "Figure digitization cannot provide an independent reproduction"
        in (decision["reason"])
    )
