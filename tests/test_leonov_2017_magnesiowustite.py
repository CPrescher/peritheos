import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material
from scripts.reproduce_leonov_2017_magnesiowustite import (
    TABLE_I,
    bm3_pressure,
    conventional_cell_volume,
    reproduce,
)

ROOT = Path(__file__).parents[1]
DOI = "10.1103/physrevb.96.075136"
MATERIALS = {
    0.000: "feo",
    0.125: "mg0125fe0875o",
    0.250: "mg025fe075o",
    0.375: "mg0375fe0625o",
    0.500: "mg05fe05o",
    0.625: "mg0625fe0375o",
    0.750: "mg075fe025o",
    0.875: "mg0875fe0125o",
}


def _document(identifier: str) -> dict:
    path = ROOT / "peritheos" / "data" / "materials" / f"{identifier}.eosmat"
    return json.loads(path.read_text(encoding="utf-8"))


def _records() -> dict[float, tuple[dict, dict]]:
    result = {}
    for fraction, material in MATERIALS.items():
        document = _document(material)
        matches = [
            row
            for row in document["eos_records"]
            if row["reference"]["doi"].lower() == DOI
        ]
        assert len(matches) == 1
        result[fraction] = (document, matches[0])
    return result


def test_leonov_preserves_all_complete_hs_table_rows_and_correct_volume_basis():
    for fraction, (document, record) in _records().items():
        v0_bohr3, hs_k0, _ls_k0, _transition = TABLE_I[fraction]
        assert record["eos"]["parameters"] == pytest.approx(
            {
                "V0": conventional_cell_volume(v0_bohr3),
                "K0": hs_k0,
                "K0_prime": 4.1,
            },
            abs=1e-9,
        )
        assert record["fixed_parameters"] == ["K0_prime"]
        assert record["temperature_ref"] == 1160
        assert document["formula_units_per_cell"] == 4


def test_leonov_records_match_independent_bm3_equation():
    for _fraction, (document, record) in _records().items():
        executable = Material.from_eosmat(document).get_eos_record(record["identifier"])
        parameters = record["eos"]["parameters"]
        volumes = parameters["V0"] * np.array([0.95, 0.85, 0.75])
        expected = bm3_pressure(
            volumes,
            parameters["V0"],
            parameters["K0"],
            parameters["K0_prime"],
        )
        assert np.asarray(executable.pressure(volumes)) == pytest.approx(expected)


def test_leonov_reproduction_and_audit_dispose_all_rows_without_inventing_ls_v0():
    result = reproduce()
    assert result["accepted_high_spin_records"] == 8
    assert result["rejected_incomplete_low_spin_rows"] == 8
    audit = (
        ROOT / "docs/literature-reproductions/leonov-2017-magnesiowustite.md"
    ).read_text()
    candidates = [
        "litcurate_5b0bfea83d12c0fb",
        "litcurate_e92cab581af0e56b",
        "litcurate_47f21997eaa8b634",
        "litcurate_2f0698ae65c320b2",
        "litcurate_5fa3c1713ba7d985",
        "litcurate_9bd8d5b53f395ddf",
        "litcurate_6776cec20630c908",
        "litcurate_b9d8bb0bcbf290ea",
        "litcurate_00fe1f0d9c6653b4",
        "litcurate_183858ab1b656cf9",
        "litcurate_e2fee3875c64e323",
        "litcurate_edb317e72af3b863",
        "litcurate_68b3a06dda40c4fe",
        "litcurate_edeba2970d95ba8d",
        "litcurate_2c3aa7dfbb7a0b84",
        "litcurate_19730defa86587a4",
    ]
    assert all(audit.count(candidate) == 1 for candidate in candidates)
    assert audit.count("| ACCEPT |") == 8
    assert audit.count("| REJECT |") == 8
