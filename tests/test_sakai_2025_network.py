"""Recover measured pairs and distinguish conditional from sequential calibration."""

import csv
import json
from collections import Counter

import numpy as np
import pytest

from scripts.reproduce_sakai_2025_gold import pressure
from scripts.reproduce_sakai_2025_network import (
    ARCHIVE,
    PAIRS,
    historical_nacl_pairs,
    reproduce,
    workbook_pairs,
)


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_every_workbook_pair_uses_its_source_cells_and_keeps_missing_errors():
    rows = workbook_pairs()
    assert Counter(r["source_table"] for r in rows) == {
        "S1": 139,
        "S2": 33,
        "S3": 54,
        "S4": 19,
        "S5": 48,
        "S6": 19,
        "S7": 81,
    }
    source = json.loads(ARCHIVE.read_text(encoding="utf-8"))["sheets"]
    for r in rows:
        cells = source[r["source_sheet"]]["cells"]
        for role in ("target", "calibrant"):
            cell = cells[r[f"{role}_volume_source_cell"]]
            assert r[f"{role}_volume_a3"] == cell["value"]
            assert r[f"{role}_volume_source_token"] == cell["source_token"]
    iron = next(r for r in rows if r["source_table"] == "S6")
    assert iron["calibrant_volume_source_cell"] == "J5"
    assert iron["calibrant_volume_a3"] == 29.082625
    assert iron["target_error_a3"] is None
    with PAIRS.open(newline="", encoding="utf-8") as stream:
        bundled = list(csv.DictReader(stream))
    assert len(bundled) == 393
    for raw, stored in zip(rows, bundled):
        assert stored == {k: "" if v is None else str(v) for k, v in raw.items()}


def test_sequential_stages_actually_use_fitted_upstream_coefficients(report):
    stages = report["modes"]["sequential_measured_subsets"]
    by_target = {s["target"]: s for s in stages}
    for target, calibrant in (
        ("platinum", "rhenium"),
        ("mgo", "platinum"),
        ("gold", "tungsten"),
        ("molybdenum", "tungsten"),
        ("iron", "tungsten"),
    ):
        assert (
            by_target[target]["upstream_parameters"][calibrant]
            == by_target[calibrant]["fitted_parameters"]
        )
    assert by_target["rhenium"]["fitted_parameters"][2] == pytest.approx(
        4.4109690561, abs=1e-8
    )
    assert by_target["gold"]["fitted_parameters"][2] == pytest.approx(
        5.715517182, abs=1e-8
    )
    for stages in report["modes"].values():
        for stage in stages:
            assert stage["optimizer_success"]
            assert stage["jacobian_rank"] == len(stage["free_parameters"])
            assert not any(stage["active_bounds"])
            assert stage["native_max_difference_gpa"] < 1e-8
            rows = stage["reductions"]
            residual = np.array(
                [
                    float(pressure(r["target_volume_a3"], stage["fitted_parameters"]))
                    - r["calculated_pressure_gpa"]
                    for r in rows
                ]
            )
            if stage["anchor"]:
                a = stage["anchor"]
                residual = np.append(
                    residual,
                    pressure(a["volume_a3"], stage["fitted_parameters"])
                    - a["pressure_gpa"],
                )
            assert np.sqrt(np.mean(residual**2)) == pytest.approx(
                stage["refit_rms_gpa"]
            )
            assert not stage["author_workflow_reproduced"]
            assert not stage["source_uncertainty_reproduced"]
    assert not report["full_network_parity_established"]


def test_nacl_uses_only_named_existing_pairs_and_keeps_reference_point_separate(report):
    rows = historical_nacl_pairs()
    assert len(rows) == 103
    assert all(r["calibrant"] in ("platinum", "gold", "mgo") for r in rows)
    assert all("S9" not in r["source"] for r in rows)
    stage = report["modes"]["published_upstream"][-1]
    assert stage["observations"] == stage["measured_pairs"] + 1
    assert stage["anchor"]["volume_a3"] == 27.25
    assert stage["anchor"]["pressure_gpa"] == 32.3


def test_new_measured_tables_are_available_through_public_dataset_api():
    from peritheos import Material, get_material_document
    from scripts.reproduce_sakai_2025_network import BRANCHES

    for table, target, _, _, _ in BRANCHES:
        if target == "gold":
            continue  # Its existing source-specific Au dataset is preserved.
        document = get_material_document(target)
        record = next(
            r for r in document["eos_records"] if "_sakai_2025_" in r["identifier"]
        )
        material = Material.from_eosmat(
            document, record_identifiers=[record["identifier"]]
        )
        dataset = material.get_dataset(
            f"{target}_sakai_2025_table_s{table}_measured_pairs"
        )
        assert dataset.checksum_verified
        assert len(dataset) == sum(
            r["source_table"] == f"S{table}" for r in workbook_pairs()
        )
        assert set(dataset["target"].tolist()) == {target}


def test_recovered_history_is_reduced_in_dependency_order_without_new_markers(report):
    from scripts.sakai_2025_historical_inputs import historical_rows

    stages = report["modes"]["sequential_with_recovered_history"]
    by_target = {s["target"]: s for s in stages}
    assert [s["target"] for s in stages] == [
        "tungsten",
        "rhenium",
        "platinum",
        "mgo",
        "gold",
        "molybdenum",
        "nacl_b2",
        "iron",
    ]
    assert (
        by_target["rhenium"]["upstream_parameters"]["tungsten"]
        == by_target["tungsten"]["fitted_parameters"]
    )
    assert (
        by_target["iron"]["upstream_parameters"]["nacl_b2"]
        == by_target["nacl_b2"]["fitted_parameters"]
    )
    assert by_target["gold"]["observations"] == 210
    assert by_target["tungsten"]["observations"] == 96
    assert by_target["iron"]["observations"] == 72
    assert by_target["gold"]["fitted_parameters"][2] == pytest.approx(
        5.78075052, abs=1e-7
    )
    # Real W coordinates are used for Re run 3; helium run 2 is omitted.
    re = historical_rows("rhenium")
    assert len(re) == 48
    assert sum(r["calibrant"] == "tungsten" for r in re) == 27
    ruby = [r for r in re if r["calibrant"] == "ruby2020"]
    assert len(ruby) == 21
    assert all(r["calibrant_volume_a3"] is None for r in ruby)
    assert all("rounded" in r["qualification"] for r in ruby)
    assert all(
        r["source"].endswith("300K") or r["source"].startswith("Dorfman")
        for r in historical_rows("mgo")
    )
