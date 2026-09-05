import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_mookherjee_2015_365a_phase import (
    finite_strain_pressure,
    reproduce,
)

ROOT = Path(__file__).parents[1]
MATERIAL = "mgsioh6_365a_phase"
DOI = "10.2138/am-2015-5312"
RESOURCE = "mgsioh6-365a-mookherjee-2015-supplement-table1-pv.csv"
CHECKSUM = "85484dad4825baeea039b488f42f17ab42b3a5de341f4406c97b46ea24dd748e"


def _rows():
    path = ROOT / "peritheos" / "data" / "datasets" / RESOURCE
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_mookherjee_material_is_exact_mgsioh6_p21_z2():
    document = get_material_document(MATERIAL)
    assert document["formula"] == "MgSiH6O6"
    assert document["phase"].startswith("monoclinic hydroxide perovskite")
    assert document["space_group"] == "P21"
    assert document["space_group_number"] == 4
    assert document["formula_units_per_cell"] == 2
    lattice = document["lattice"]
    volume = (
        lattice["a"]
        * lattice["b"]
        * lattice["c"]
        * math.sin(math.radians(lattice["beta"]))
    )
    assert volume == pytest.approx(194.52, abs=0.01)

    contents = Counter()
    for site in document["atom_sites"]:
        contents[site["element"]] += site["site_multiplicity"] * site["occupancy"]
    assert contents == {"Mg": 2.0, "Si": 2.0, "H": 12.0, "O": 12.0}
    assert sum(contents.values()) == 28.0
    assert document["source"]["topology_proxy"] is True


def test_mookherjee_preserves_three_distinct_source_parameterizations():
    records = get_material_document(MATERIAL)["eos_records"]
    assert len(records) == 3
    assert records[0]["default_for"] == "equilibrium"
    assert "default_for" not in records[1]
    assert "default_for" not in records[2]
    assert records[0]["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 194.52, "K0": 83.0, "K0_prime": 4.9},
    }
    assert records[1]["eos"] == {
        "type": "BM4",
        "model": "birch_murnaghan_4",
        "parameters": {
            "V0": 194.52,
            "K0": 77.0,
            "K0_prime": 7.9,
            "K0_double_prime": -0.7,
        },
    }
    assert records[2]["eos"] == {
        "type": "BM4",
        "model": "birch_murnaghan_4",
        "parameters": {
            "V0": 202.02,
            "K0": 80.0,
            "K0_prime": 3.4,
            "K0_double_prime": -0.05,
        },
    }
    assert records[1]["fixed_parameters"] == ["V0"]
    assert records[2]["fixed_parameters"] == []
    assert "sensitivity" in records[1]["label"]
    assert "GGA model-crystal" in records[2]["label"]
    for record in records:
        assert record["reference"]["doi"].lower() == DOI
        assert record["record_kind"] == "published"
        assert record["scientific_validation"]["status"] == ("primary_source_validated")


def test_mookherjee_official_supplement_is_complete_and_bundled():
    document = get_material_document(MATERIAL)
    dataset = document["datasets"][0]
    path = ROOT / "peritheos" / "data" / "datasets" / RESOURCE
    assert dataset["resource"]["sha256"] == CHECKSUM
    assert hashlib.sha256(path.read_bytes()).hexdigest() == CHECKSUM
    rows = _rows()
    assert len(rows) == 95
    assert Counter(row["experiment"] for row in rows) == {
        "ambient": 1,
        "I": 50,
        "II": 21,
        "III": 23,
    }
    assert Counter(row["path"] for row in rows) == {
        "ambient": 1,
        "compression": 80,
        "decompression": 14,
    }
    assert min(float(row["pressure_gpa"]) for row in rows) == 0.0
    assert max(float(row["pressure_gpa"]) for row in rows) == 41.0


def test_primary_source_finite_strain_equations_match_executable_records():
    document = get_material_document(MATERIAL)
    executable = Material.from_eosmat(document)
    volumes = np.array([155.0, 170.0, 185.0])
    for stored in document["eos_records"]:
        parameters = stored["eos"]["parameters"]
        expected = finite_strain_pressure(
            volumes,
            parameters["V0"],
            parameters["K0"],
            parameters["K0_prime"],
            parameters.get("K0_double_prime"),
        )
        actual = executable.get_eos_record(stored["identifier"]).pressure(volumes)
        assert np.asarray(actual) == pytest.approx(expected, abs=1.0e-12)


def test_mookherjee_reproduction_recovers_curve_and_coefficient_identity():
    result = reproduce()
    assert result["official_supplement"]["observations"] == 95
    assert (
        result["official_supplement"]["compression_observations_with_uncertainties"]
        == 80
    )
    curves = result["published_curve_residuals_all_rows"]
    assert curves["experimental_bm3"]["pressure_rmse_gpa"] < 0.43
    assert curves["experimental_bm4_sensitivity"]["pressure_rmse_gpa"] < 0.41

    refits = result["compression_odr_diagnostic"]
    bm3 = refits["experimental_bm3"]["parameters"]
    assert bm3["K0"] == pytest.approx(83.0, abs=1.0)
    assert bm3["K0_prime"] == pytest.approx(4.9, abs=0.2)
    bm4 = refits["experimental_bm4_sensitivity"]["parameters"]
    assert bm4["K0"] == pytest.approx(77.0, abs=2.0)
    assert bm4["K0_prime"] == pytest.approx(7.9, abs=0.8)
    assert bm4["K0_double_prime"] == pytest.approx(-0.7, abs=0.2)


def test_mookherjee_audit_disposes_all_seven_litcurate_rows_once():
    audit = (
        ROOT / "docs" / "literature-reproductions" / "mookherjee-2015-365a-phase.md"
    ).read_text(encoding="utf-8")
    for candidate in (
        "litcurate_b98adada3a08d56f",
        "litcurate_292be0594883e2b5",
        "litcurate_e091013c8442ddc8",
        "litcurate_efccbe0ac02295bc",
        "litcurate_d87b5d28e5c194af",
        "litcurate_f7574e211806046d",
        "litcurate_7a13e5afc23a7fe5",
    ):
        assert audit.count(candidate) == 1
