import json
from collections import Counter
from pathlib import Path

import pytest

from peritheos import list_eos_record_documents

ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "docs" / "data" / "primary-eos-refits.json"
MARKDOWN_PATH = ROOT / "docs" / "primary-eos-refits.md"


def load_ledger():
    return json.loads(LEDGER_PATH.read_text(encoding="utf-8"))


def test_primary_refit_ledger_covers_every_bundled_record_once():
    ledger = load_ledger()
    identifiers = [item["record_identifier"] for item in ledger["records"]]

    assert ledger["format"] == "peritheos.primary-eos-refit-validation"
    assert ledger["format_version"] == 1
    assert len(identifiers) == len(set(identifiers)) == 160
    assert set(identifiers) == set(list_eos_record_documents())


def test_primary_refit_summary_and_results_are_internally_consistent():
    ledger = load_ledger()
    statuses = Counter(item["status"] for item in ledger["records"])

    assert ledger["summary"] == {"total": 160, **dict(sorted(statuses.items()))}
    assert statuses == {
        "parity": 76,
        "similar": 32,
        "parity_not_achieved": 13,
        "not_refittable": 39,
    }
    assert all(
        item.get("reason")
        for item in ledger["records"]
        if item["status"] == "not_refittable"
    )


def test_primary_refit_regression_examples_and_documentation_coverage():
    ledger = load_ledger()
    by_identifier = {item["record_identifier"]: item for item in ledger["records"]}
    markdown = MARKDOWN_PATH.read_text(encoding="utf-8")

    assert by_identifier["aragonite_martinez_1996_bm2_2"]["status"] == "parity"
    assert by_identifier["kcl_campbell_1991_bm2_1"]["status"] == "parity"
    assert by_identifier["b4c_somayazulu_2023_bm3_1"]["status"] == (
        "parity_not_achieved"
    )
    assert by_identifier["b4c_somayazulu_2023_berman_refit"]["status"] == ("parity")
    assert by_identifier["gold_shen_2026_vinet_3"]["status"] == "not_refittable"
    neon_bm3 = by_identifier["neon_fcc_fei_2007_bm3_1"]
    neon_vinet = by_identifier["neon_fcc_fei_2007_vinet_2"]
    assert neon_bm3["status"] == "parity"
    assert neon_vinet["status"] == "parity"
    assert neon_bm3["observations"] == 34
    assert neon_vinet["observations"] == 34
    assert neon_bm3["dataset_identifiers"] == [
        "neon_fei_2007_figure5_digitized",
        "neon_hemley_1989_table1_fei_recalculated",
    ]
    assert neon_bm3["observed_pressure_range_gpa"] == pytest.approx(
        [10.047704716, 115.715945233]
    )
    assert [item["refit"] for item in neon_bm3["parameters"]] == pytest.approx(
        [1.4775172078, 7.8503378205]
    )
    assert [item["refit"] for item in neon_vinet["parameters"]] == pytest.approx(
        [1.1439324121, 8.2580743609]
    )
    assert neon_bm3["free_parameters"] == ["K0", "K0_prime"]
    assert neon_vinet["free_parameters"] == ["rt_eos.K0", "rt_eos.K0_prime"]
    assert "Conditional partial reproduction" in neon_bm3["qualification"]
    assert "Conditional partial reference-isotherm" in neon_vinet["qualification"]
    assert "Finger's low-pressure rows remain unavailable" in (
        neon_bm3["qualification"]
    )
    hemley = by_identifier["neon_fcc_hemley_1989_bm3_refit"]
    assert hemley["status"] == "parity"
    assert hemley["observations"] == 21
    assert hemley["free_parameters"] == ["K0", "K0_prime"]
    assert markdown.count("### `") == 45
    assert all(identifier in markdown for identifier in by_identifier)
