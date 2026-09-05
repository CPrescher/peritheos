import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
DOI = "10.1029/2011GL049519"
FORMULA = "Ca0.988Mg0.918Fe0.078Mn0.016C2O6"
DOLOMITE = "ca0988mg0918fe0078mn0016c2o6_dolomite"
DOLOMITE_III = "ca0988mg0918fe0078mn0016c2o6_dolomite_iii"
DOLOMITE_RESOURCE = (
    "ca0988mg0918fe0078mn0016c2o6-mao-2011-figure3-fe-dolomite-digitized.csv"
)
DOLOMITE_III_RESOURCE = (
    "ca0988mg0918fe0078mn0016c2o6-mao-2011-figure3-dolomite-iii-digitized.csv"
)


def _rows(resource):
    path = ROOT / "peritheos" / "data" / "datasets" / resource
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_mao_dolomite_material_identity_and_unresolved_phase_iii_structure():
    dolomite = get_material_document(DOLOMITE)
    phase_iii = get_material_document(DOLOMITE_III)
    assert dolomite["formula"] == phase_iii["formula"] == FORMULA
    assert dolomite["space_group"] == "R-3"
    assert dolomite["space_group_number"] == 148
    assert dolomite["formula_units_per_cell"] == 3
    lattice = dolomite["lattice"]
    assert (math.sqrt(3.0) / 2.0 * lattice["a"] ** 2 * lattice["c"]) == pytest.approx(
        321.75827586, abs=1.0e-8
    )

    assert phase_iii["symmetry"] == "MONOCLINIC"
    assert phase_iii["source"]["indexing_only"] is True
    assert "space_group" not in phase_iii
    assert "lattice" not in phase_iii
    assert phase_iii["atom_sites"] == []
    assert "alternative indexings" in phase_iii["source"]["structure_location"]


def test_mao_dolomite_preserves_three_published_bm3_parameterizations():
    low = get_material_document(DOLOMITE)["eos_records"]
    high_pressure = get_material_document(DOLOMITE_III)["eos_records"]
    assert len(low) == 1
    assert len(high_pressure) == 2
    assert low[0]["eos"]["parameters"] == {
        "V0": 321.77,
        "K0": 94.1,
        "K0_prime": 4.0,
    }
    assert high_pressure[0]["eos"]["parameters"] == {
        "V0": 239.2,
        "K0": 164.0,
        "K0_prime": 4.0,
    }
    assert high_pressure[1]["eos"]["parameters"] == {
        "V0": 231.8,
        "K0": 184.0,
        "K0_prime": 4.0,
    }
    for record in [*low, *high_pressure]:
        assert record["eos"]["type"] == "BM3"
        assert record["fixed_parameters"] == ["K0_prime"]
        assert record["reference"]["doi"].lower() == DOI.lower()
        assert record["scientific_validation"]["status"] == "primary_source_validated"
        assert record["pressure_calibration"]["methods"][0]["material"] == "Pt"
        material = Material.from_eosmat(
            get_material_document(DOLOMITE if record in low else DOLOMITE_III),
            record_identifiers=[record["identifier"]],
        )
        assert material.get_eos_record(record["identifier"]).pressure(
            record["eos"]["parameters"]["V0"]
        ) == pytest.approx(0.0)


@pytest.mark.parametrize(
    ("material", "resource", "checksum", "count"),
    [
        (
            DOLOMITE,
            DOLOMITE_RESOURCE,
            "cd7e22cce57c621402ecba2f577ca87aea697b3e31698d46ce20871cc8c0c4d5",
            7,
        ),
        (
            DOLOMITE_III,
            DOLOMITE_III_RESOURCE,
            "13977fd72c669eb5f02eb0c0633cb238ca8a602d2f90b82fc55a5ba04701fe5d",
            25,
        ),
    ],
)
def test_mao_figure3_digitizations_are_bundled(material, resource, checksum, count):
    dataset = get_material_document(material)["datasets"][0]
    path = ROOT / "peritheos" / "data" / "datasets" / resource
    assert dataset["resource"]["sha256"] == checksum
    assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum
    assert len(_rows(resource)) == count
    assert dataset["provenance"]["source_figure_sha256"] == (
        "cb20211a4551d2ec74a1ece62ec9c65754e63f1331cac5c57498d43bb6b388fe"
    )
    assert dataset["provenance"]["type"] == "digitized_from_figure"


@pytest.mark.parametrize(
    ("material", "record_index", "resource", "flag", "expected_rmse"),
    [
        (DOLOMITE, 0, DOLOMITE_RESOURCE, "fit_included", 0.3991598272),
        (
            DOLOMITE_III,
            0,
            DOLOMITE_III_RESOURCE,
            "high_spin_fit_included",
            2.9550573964,
        ),
        (
            DOLOMITE_III,
            1,
            DOLOMITE_III_RESOURCE,
            "low_spin_fit_included",
            2.9152095407,
        ),
    ],
)
def test_mao_published_curves_reproduce_digitized_figure3(
    material, record_index, resource, flag, expected_rmse
):
    document = get_material_document(material)
    stored = document["eos_records"][record_index]
    executable = Material.from_eosmat(
        document, record_identifiers=[stored["identifier"]]
    ).get_eos_record(stored["identifier"])
    rows = [row for row in _rows(resource) if row[flag] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    residual = np.asarray(executable.pressure(volume)) - pressure
    assert math.sqrt(float(np.mean(residual**2))) == pytest.approx(
        expected_rmse, abs=1.0e-9
    )


def test_mao_audit_disposes_all_three_litcurate_candidates_once():
    audit = (
        ROOT / "docs" / "literature-reproductions" / "mao-2011-dolomite-iii.md"
    ).read_text(encoding="utf-8")
    for candidate in (
        "litcurate_ca7038ed5d3c92bc",
        "litcurate_883e64429604ae2d",
        "litcurate_9fea9bd5e6a1abd0",
    ):
        assert audit.count(candidate) == 1
