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
    assert len(identifiers) == len(set(identifiers)) == 163
    assert set(identifiers) == set(list_eos_record_documents())


def test_primary_refit_summary_and_results_are_internally_consistent():
    ledger = load_ledger()
    statuses = Counter(item["status"] for item in ledger["records"])

    assert ledger["summary"] == {"total": 163, **dict(sorted(statuses.items()))}
    assert statuses == {
        "parity": 81,
        "similar": 32,
        "parity_not_achieved": 10,
        "not_refittable": 40,
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
    walker = by_identifier["kcl_walker_2002_bm3_2"]
    assert walker["status"] == "similar"
    assert walker["free_parameters"] == [
        "rt_eos.K0",
        "rt_eos.K0_prime",
        "alpha_KT",
    ]
    assert [item["refit"] for item in walker["parameters"]] == pytest.approx(
        [23.7753978405, 4.4158737620, 0.002766245688]
    )
    assert walker["stages"][0]["observations"] == 8
    assert walker["stages"][1]["observations"] == 39
    assert "not a simultaneous four-parameter fit" in walker["qualification"]
    simultaneous = walker["simultaneous_fit_diagnostic"]
    simultaneous_parameters = simultaneous["peritheos_four_parameter_refit"]
    assert simultaneous_parameters["rt_eos.V0"] == pytest.approx(55.39232501)
    assert simultaneous_parameters["rt_eos.K0"] == pytest.approx(14.188725516)
    assert simultaneous_parameters["rt_eos.K0_prime"] == pytest.approx(7.156559906)
    assert simultaneous_parameters["alpha_KT"] == pytest.approx(0.002696920764)
    tateno = by_identifier["kcl_b2_tateno_2019_vinet_4"]
    assert tateno["status"] == "parity"
    assert tateno["dataset_identifiers"] == ["kcl_tateno_2019_table_s1_pvt"]
    assert [item["refit"] for item in tateno["parameters"]] == pytest.approx(
        [18.344616371, 5.600955151, 2.295185193, 0.824899722]
    )
    assert "Corrected final-publication reproduction" in tateno["qualification"]
    chidester = by_identifier["kcl_b2_chidester_2021_bm3_5"]
    assert chidester["status"] == "parity"
    assert chidester["observations"] == 278
    assert chidester["dataset_identifiers"] == [
        "kcl_dewaele_2012_table1_compression",
        "kcl_chidester_2021_supplemental_pvt",
    ]
    assert [item["refit"] for item in chidester["parameters"]] == pytest.approx(
        [53.203552571, 23.972141965, 4.557980658, 2.917137053, 0.965242582]
    )
    assert chidester["high_temperature_refit_rmse_gpa"] == pytest.approx(
        1.581924568
    )
    assert chidester["source_reported_pressure_rmse_gpa"] == 1.6
    assert "155 new high-temperature rows" in chidester["qualification"]
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
    assert (
        "Finger's low-pressure rows remain unavailable" in (neon_bm3["qualification"])
    )
    hemley = by_identifier["neon_fcc_hemley_1989_bm3_refit"]
    assert hemley["status"] == "parity"
    assert hemley["observations"] == 21
    assert hemley["free_parameters"] == ["K0", "K0_prime"]
    rbcl = by_identifier["rbcl_b2_campbell_1994_bm3_1"]
    assert rbcl["status"] == "parity"
    assert rbcl["observations"] == 24
    assert rbcl["dataset_identifiers"] == [
        "rbcl_campbell_1994_table1_compression"
    ]
    assert [item["refit"] for item in rbcl["parameters"]] == pytest.approx(
        [17.8808436600, 5.2381519592]
    )
    cscl = by_identifier["cscl_campbell_1994_bm3_1"]
    assert cscl["status"] == "parity"
    assert cscl["observations"] == 22
    assert cscl["dataset_identifiers"] == [
        "cscl_campbell_1994_table1_compression",
        "cscl_yagi_1978_table1_compression",
    ]
    assert [item["refit"] for item in cscl["parameters"]] == pytest.approx(
        [17.6967334851, 5.2026008383]
    )
    assert "Complete source-data reproduction" in cscl["qualification"]
    explained = [
        item
        for item in ledger["records"]
        if item["status"] in {"similar", "parity_not_achieved", "refit_failed"}
    ]
    assert markdown.count("### `") == len(explained) == 42
    assert all(identifier in markdown for identifier in by_identifier)
    failed = [
        item
        for item in ledger["records"]
        if item["status"] in {"parity_not_achieved", "refit_failed"}
    ]
    for item in failed:
        anchor = f"investigation-{item['record_identifier']}"
        assert markdown.count(f'<a id="{anchor}"></a>') == 1
        assert markdown.count(f"](#{anchor})") >= 2
