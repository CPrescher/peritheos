import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material
from scripts.reproduce_solomatova_2016_ferropericlase import (
    TABLE_7,
    bm3_pressure,
    reproduce,
)

ROOT = Path(__file__).parents[1]
DOI = "10.2138/am-2016-5510"
MATERIALS = (
    "mg090fe010o",
    "mg083fe017o",
    "mg075fe025o",
    "mg065fe035o",
    "mg061fe039o",
    "mg0490fe0483ti0027o",
    "mgfe60o",
)


def _documents() -> list[dict]:
    return [
        json.loads(
            (ROOT / f"peritheos/data/materials/{name}.eosmat").read_text(
                encoding="utf-8"
            )
        )
        for name in MATERIALS
    ]


def _records() -> list[tuple[dict, dict]]:
    return [
        (document, record)
        for document in _documents()
        for record in document["eos_records"]
        if record["reference"]["doi"].lower() == DOI
    ]


def test_solomatova_preserves_all_sixteen_table7_reference_branches():
    records = _records()
    assert len(records) == 16
    observed = sorted(
        tuple(record["eos"]["parameters"].values()) for _document, record in records
    )
    assert observed == sorted(TABLE_7.values())
    assert sum(record["fixed_parameters"] == ["K0_prime"] for _, record in records) == 8
    assert all(
        record["spin_crossover_context"]["mixer_implemented"] is False
        for _, record in records
    )


def test_solomatova_all_branches_match_independent_bm3_equation():
    for document, record in _records():
        parameters = record["eos"]["parameters"]
        volumes = parameters["V0"] * np.array([0.95, 0.85, 0.75])
        expected = bm3_pressure(
            volumes,
            parameters["V0"],
            parameters["K0"],
            parameters["K0_prime"],
        )
        actual = (
            Material.from_eosmat(document)
            .get_eos_record(record["identifier"])
            .pressure(volumes)
        )
        assert np.asarray(actual) == pytest.approx(expected)


def test_solomatova_fp48_primary_tables_are_complete_and_unchanged():
    resources = {
        "mg0490fe0483ti0027o-solomatova-2016-table1-pv.csv": (
            45,
            "8cffb2ec0c7889a4e9d23f8565ba23c32588ddc15f8351d08d4e3c3513ce8af6",
        ),
        "mg0490fe0483ti0027o-solomatova-2016-table3-crossover-grid.csv": (
            36,
            "2ca2cafeff2be1be7eb165876021cb9bf4e2807ff6db904a78ac879dddfb3e3b",
        ),
    }
    for name, (count, checksum) in resources.items():
        path = ROOT / "peritheos/data/datasets" / name
        with path.open(newline="", encoding="utf-8") as stream:
            assert len(list(csv.DictReader(stream))) == count
        assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum

    diagnostics = reproduce()
    assert (
        diagnostics["fp48_endmember_diagnostics"]["high_spin_pressure_rmse_gpa"] < 0.52
    )
    assert (
        diagnostics["fp48_endmember_diagnostics"]["low_spin_pressure_rmse_gpa"] < 0.31
    )


def test_solomatova_audit_disposes_all_sixteen_rows_once():
    audit = (
        ROOT / "docs/literature-reproductions/solomatova-2016-ferropericlase.md"
    ).read_text(encoding="utf-8")
    candidates = (
        "litcurate_d84461e22abf5080",
        "litcurate_1524ace141e7f7f0",
        "litcurate_99bd28480585bdaa",
        "litcurate_b0f8a936ede6bc5e",
        "litcurate_ae5f29f994c6d8f3",
        "litcurate_c7874f561714da3d",
        "litcurate_eb3e9425325e0c18",
        "litcurate_7f3cefd535845c56",
        "litcurate_28e6643f0cf0c399",
        "litcurate_85a318658060d6c0",
        "litcurate_67984b8ee2408ca4",
        "litcurate_2fbe372e10a362f1",
        "litcurate_4fe52dbe32e1dfdb",
        "litcurate_72b31a925a02dab0",
        "litcurate_1b209d91f037a754",
        "litcurate_112ff70f49853ccc",
    )
    assert all(audit.count(candidate) == 1 for candidate in candidates)
    assert audit.count("| ACCEPT |") == 16
