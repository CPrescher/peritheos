"""Source preservation and conditional reconstruction limits for supplied Fe3C data."""

import json

import numpy as np
import pytest

from peritheos import get_material, get_material_document
from scripts.reproduce_prescher_2012_fe3c import ROOT, observations, reproduce


@pytest.fixture(scope="module")
def result():
    return reproduce()


def test_supplement_rows_errors_and_atomic_tokens_survive_public_loading():
    raw = json.loads(
        (
            ROOT
            / "peritheos/data/datasets/cementite-prescher-2012-table-s1-tokens.json"
        ).read_text()
    )
    dataset = get_material("cementite").get_dataset("cementite_prescher_2012_table_s1")
    assert len(dataset) == len(raw["rows"]) == 19
    assert dataset.checksum_verified
    np.testing.assert_array_equal(
        dataset["pressure_gpa"],
        [
            1.8,
            4.7,
            6.7,
            9.1,
            11.4,
            13.4,
            16.2,
            18.4,
            20.7,
            22.6,
            25.3,
            27.5,
            29.7,
            32.1,
            34.7,
            37,
            39.6,
            41.9,
            47.4,
        ],
    )
    assert dataset["volume_a3"][0] == 152.98
    assert dataset["volume_error_a3"][0] == 0.07
    assert dataset["volume_a3"][-1] == 129.50
    assert dataset["volume_error_a3"][-1] == 0.05
    assert dataset["r1_percent"].max() == 8.2
    assert sum(len(row["atomic_tokens_xyzu"]) for row in raw["rows"]) == 57
    assert raw["rows"][-1]["atomic_tokens_xyzu"]["Fe1"] == [
        "0.3428(5)",
        "0.1857(5)",
        "0.0629(6)",
        "0.0107(9)",
    ]
    doc = get_material_document("cementite", format_version=4)
    support = next(
        d
        for d in doc["datasets"]
        if d["identifier"] == "cementite_prescher_2012_table_s1"
    )
    assert set(support["used_by_eos_records"]) == {
        "cementite_prescher_2012_full_range_bm3",
        "cementite_prescher_2012_below_22_gpa_bm3",
        "cementite_prescher_2012_22_47_gpa_bm3",
    }
    errors = [c for c in support["columns"] if c.get("of")]
    assert errors and all(c["role"] == "uncertainty" for c in errors)


def test_high_pressure_structure_is_independent_of_eos_intercepts():
    d = get_material_document("cementite_nonmagnetic", format_version=4)
    assert d["space_group"] == "Pnma"
    assert d["space_group_number"] == 62 and d["formula_units_per_cell"] == 4
    lattice = d["lattice"]
    assert [lattice[k] for k in ("a", "b", "c")] == [4.820, 6.357, 4.226]
    assert lattice["a"] * lattice["b"] * lattice["c"] == pytest.approx(129.50, abs=0.05)
    sites = d["atom_sites"]
    assert [(s["element"], s["wyckoff"]) for s in sites] == [
        ("Fe", "8d"),
        ("Fe", "4c"),
        ("C", "4c"),
    ]
    assert [sites[0][k] for k in ("x", "y", "z")] == [0.1857, 0.0629, 0.3428]
    assert sites[1]["y"] == sites[2]["y"] == 0.25
    assert all(s["occupancy"] == 1 for s in sites)
    provenance = d["structure_provenance"]
    assert provenance["pressure_gpa"] == 47.4 and provenance["temperature_k"] == 293
    assert "not refined occupancies" in d["notes"]
    assert [r["parameters"]["V0"]["value"] for r in d["eos_records"]] == [148.8, 147.3]


def test_conditional_coefficient_agreement_does_not_invent_source_v0(result):
    assert result["observations"] == 19
    for name, count in [("full_range", 19), ("below_22_gpa", 9), ("above_22_gpa", 10)]:
        row = result["source_subsets"][name]
        assert row["published_parameters"]["V0"] is None
        assert row["published_parameter_errors"]["V0"] is None
        assert row["published_curve_replay"] == "unavailable_without_source_V0"
        assert row["author_workflow_reproduced"] is False
        fit = row["fits"]["reported_coordinate_widths"]
        assert fit["observations"] == count
        assert fit["K0_and_K0_prime_within_source_widths"] is (name == "full_range")
    fit = result["source_subsets"]["full_range"]["fits"]["reported_coordinate_widths"]
    assert fit["parameters"]["K0"] == pytest.approx(160.53548, rel=1e-6)
    assert fit["parameters"]["K0_prime"] == pytest.approx(6.031741, rel=1e-6)
    assert observations()["pressure_gpa"].min() > 0
    # This historical table-only diagnostic remains unchanged after explicitly
    # adopted graph-derived V0 records were integrated separately.
    records = [
        r
        for r in get_material_document("cementite", format_version=4)["eos_records"]
        if r["reference"].get("doi") == "10.1103/PhysRevB.85.140402"
    ]
    assert len(records) == 3
    assert all(
        r["evidence"]["parameter_provenance"]["V0"]["source_value_published"] is False
        for r in records
    )


def test_archive_working_versions_do_not_replace_final_source_observations():
    archive = json.loads(
        (
            ROOT / "peritheos/data/research/cementite/prescher-2012-archive.json"
        ).read_text()
    )
    assert archive["published_fit_constraints"] is None
    assert archive["published_V0_status"].startswith("Unresolved:")
    book = archive["workbooks"][0]
    cells = {cell["address"]: cell["value"] for cell in book["sheets"][0]["cells"]}
    assert cells["B2"] == 0 and cells["I2"] == pytest.approx(163.40230309104)
    assert cells["B3"] == 1.8 and cells["I3"] == pytest.approx(153.6255042696)
    origin = archive["origin_projects"][0]
    assert len(origin["finite_observation_rows"]) == 10
    assert origin["graphs"] == origin["functions"] == origin["notes"] == 0
    delta = [
        r["archive_minus_published_volume_a3"]
        for r in origin["comparison_to_final_table_s1"]
    ]
    assert max(abs(x) for x in delta) == pytest.approx(0.04)
    assert observations()["volume_a3"][0] == 152.98
    assert len(observations()["volume_a3"]) == 19
