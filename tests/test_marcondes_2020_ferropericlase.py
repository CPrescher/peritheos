import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material
from scripts.reproduce_marcondes_2020_ferropericlase import (
    TABLE_I,
    bm3_pressure,
    reproduce,
)

ROOT = Path(__file__).parents[1]
DOI = "10.1103/physrevb.102.104112"


def _documents() -> list[dict]:
    return [
        json.loads(
            (ROOT / f"peritheos/data/materials/{name}.eosmat").read_text(
                encoding="utf-8"
            )
        )
        for name in ("mg096875fe003125o", "mg09375fe00625o")
    ]


def test_marcondes_keeps_all_eight_distinct_configuration_spin_branches():
    records = [
        record
        for document in _documents()
        for record in document["eos_records"]
        if record["reference"]["doi"].lower() == DOI
    ]
    assert len(records) == 8
    assert {
        record["identifier"].split("marcondes_2020_")[1].rsplit("_bm3_", 1)[0]
        for record in records
    } == {
        "hs",
        "ls",
        "11nn_hs",
        "11nn_ls",
        "11nn_ms",
        "2nn_hs",
        "2nn_ls",
        "2nn_ms",
    }


def test_marcondes_source_parameters_match_executable_curves():
    by_identifier = {
        record["identifier"]: (document, record)
        for document in _documents()
        for record in document["eos_records"]
    }
    for key, parameters in TABLE_I.items():
        identifier_key = key.split("_", 1)[1]
        matches = [
            value
            for identifier, value in by_identifier.items()
            if f"_{identifier_key}_bm3_" in identifier
            and (key.startswith("3fp_") == identifier.startswith("mg096875"))
        ]
        assert len(matches) == 1
        document, record = matches[0]
        assert tuple(record["eos"]["parameters"].values()) == parameters
        v0, k0, k0_prime = parameters
        volumes = v0 * np.array([0.95, 0.85, 0.75])
        executable = Material.from_eosmat(document).get_eos_record(record["identifier"])
        assert np.asarray(executable.pressure(volumes)) == pytest.approx(
            bm3_pressure(volumes, v0, k0, k0_prime)
        )


def test_marcondes_reproduction_and_audit_dispose_all_rows():
    assert reproduce()["accepted_records"] == 8
    audit = (
        ROOT / "docs/literature-reproductions/marcondes-2020-ferropericlase.md"
    ).read_text(encoding="utf-8")
    candidates = (
        "litcurate_8b61325b3373dcc9",
        "litcurate_b5a171ea0852a37e",
        "litcurate_77c20002800f4789",
        "litcurate_ec75255e7d6eae38",
        "litcurate_2116ff41406a0b93",
        "litcurate_cd07018138aa391c",
        "litcurate_c871ad45b4ff997f",
        "litcurate_442cdb585530380d",
    )
    assert all(audit.count(candidate) == 1 for candidate in candidates)
    assert audit.count("| ACCEPT |") == 8
