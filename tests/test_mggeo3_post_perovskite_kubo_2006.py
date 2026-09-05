import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.eosmat import list_eos_record_documents
from scripts.reproduce_kubo_2006_mggeo3_post_perovskite import (
    PREFERRED,
    SENSITIVITY,
    curve_metrics,
    fit_with_fixed_derivative,
    platinum_reduction_metrics,
)

ROOT = Path(__file__).parents[1]
DOI = "10.1029/2006GL025686"
DATA_RESOURCE = "mggeo3-post-perovskite-kubo-2006-table-s1-pv.csv"
DATA_SHA256 = "2dadf9072636404baf49fa03c04e07d65c093e7d85c6cb1e2508d985056773ce"
PREFERRED_IDENTIFIER = "mggeo3_post_perovskite_kubo_2006_bm3_1"
SENSITIVITY_IDENTIFIER = "mggeo3_post_perovskite_kubo_2006_bm2_sensitivity_2"
LITCURATE_IDENTIFIERS = {
    "litcurate_257ec8acbb2a2877",
    "litcurate_973395417d1d42d7",
    "litcurate_24cd91aaacc9abd7",
    "litcurate_49869576e372eb74",
    "litcurate_21f4bb278eb404db",
    "litcurate_0a9cb64914d6f525",
    "litcurate_89524810ffc6d316",
    "litcurate_7000ad6e95c15aea",
}


def _rows():
    path = ROOT / "peritheos" / "data" / "datasets" / DATA_RESOURCE
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_kubo_mggeo3_phase_structure_and_record_selection():
    document = get_material_document("mggeo3_post_perovskite")
    assert document["formula"] == "MgGeO3"
    assert document["space_group"] == "Cmcm"
    assert document["space_group_number"] == 63
    assert document["formula_units_per_cell"] == 4
    assert math.prod(document["lattice"][axis] for axis in ("a", "b", "c")) == (
        pytest.approx(140.006034, abs=5.0e-7)
    )

    cell_contents = Counter()
    for site in document["atom_sites"]:
        cell_contents[site["element"]] += site["site_multiplicity"] * site["occupancy"]
    assert cell_contents == {"Mg": 4.0, "Ge": 4.0, "O": 12.0}

    records = {record["identifier"]: record for record in document["eos_records"]}
    assert set(records) == {PREFERRED_IDENTIFIER, SENSITIVITY_IDENTIFIER}
    assert records[PREFERRED_IDENTIFIER]["default"] is True
    assert records[PREFERRED_IDENTIFIER]["default_for"] == "equilibrium"
    assert records[SENSITIVITY_IDENTIFIER]["default"] is False


def test_kubo_mggeo3_published_parameters_and_fixed_status():
    records = {
        record["identifier"]: record
        for record in get_material_document("mggeo3_post_perovskite")["eos_records"]
    }
    preferred = records[PREFERRED_IDENTIFIER]
    assert preferred["reference"]["doi"].lower() == DOI.lower()
    assert preferred["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": PREFERRED,
    }
    assert preferred["parameter_errors"] == {
        "V0": 0.7,
        "K0": 5.0,
        "K0_prime": None,
    }
    assert preferred["fixed_parameters"] == ["K0_prime"]
    assert preferred["experimental_pressure_range_gpa"] == [47.3, 196.3]
    assert preferred["temperature_ref"] == 300.0

    sensitivity = records[SENSITIVITY_IDENTIFIER]
    assert sensitivity["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 175.9, "K0": 245.0},
    }
    assert sensitivity["parameter_errors"] == {"V0": 0.6, "K0": 5.0}
    assert sensitivity["fixed_parameters"] == []


def test_kubo_mggeo3_table_s1_is_complete_and_checksum_locked():
    document = get_material_document("mggeo3_post_perovskite")
    dataset = document["datasets"][0]
    assert dataset["identifier"] == "mggeo3_post_perovskite_kubo_2006_table_s1_pv"
    assert dataset["used_by_eos_records"] == [
        PREFERRED_IDENTIFIER,
        SENSITIVITY_IDENTIFIER,
    ]

    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == DATA_SHA256
    assert dataset["resource"]["sha256"] == DATA_SHA256

    rows = _rows()
    assert len(rows) == 25
    assert len(rows[0]) == 17
    included = [row for row in rows if row["fit_included"] == "1"]
    excluded = [row for row in rows if row["fit_included"] == "0"]
    assert len(included) == 22
    assert [float(row["pressure_gpa"]) for row in excluded] == [45.0, 36.3, 7.4]
    assert {int(row["run"]) for row in rows} == {1, 2, 3}

    for row in rows:
        a = float(row["a_angstrom"])
        b = float(row["b_angstrom"])
        c = float(row["c_angstrom"])
        sa = float(row["a_sigma_angstrom"])
        sb = float(row["b_sigma_angstrom"])
        sc = float(row["c_sigma_angstrom"])
        volume = a * b * c
        volume_sigma = volume * math.sqrt((sa / a) ** 2 + (sb / b) ** 2 + (sc / c) ** 2)
        assert float(row["volume_a3"]) == pytest.approx(volume, abs=5.0e-7)
        assert float(row["volume_sigma_a3"]) == pytest.approx(volume_sigma, abs=5.0e-7)


def test_kubo_mggeo3_source_curves_and_refits_reproduce():
    preferred_rmse, preferred_maximum = curve_metrics(PREFERRED)
    assert preferred_rmse == pytest.approx(1.3945069526, abs=5.0e-10)
    assert preferred_maximum == pytest.approx(3.9378912089, abs=5.0e-10)
    assert fit_with_fixed_derivative(4.4) == pytest.approx(
        (179.20054626, 206.71122613), abs=5.0e-7
    )

    sensitivity_rmse, sensitivity_maximum = curve_metrics(SENSITIVITY)
    assert sensitivity_rmse == pytest.approx(1.3561285588, abs=5.0e-10)
    assert sensitivity_maximum == pytest.approx(3.9722613183, abs=5.0e-10)
    assert fit_with_fixed_derivative(4.0) == pytest.approx(
        (175.92371146, 244.61773418), abs=5.0e-7
    )

    published = np.array([PREFERRED["V0"], PREFERRED["K0"]])
    published_errors = np.array([0.7, 5.0])
    assert np.all(np.abs(fit_with_fixed_derivative(4.4) - published) < published_errors)


def test_kubo_mggeo3_peritheos_and_pt_pressure_scale_are_executable():
    document = get_material_document("mggeo3_post_perovskite")
    material = Material.from_eosmat(
        document,
        record_identifiers=[PREFERRED_IDENTIFIER, SENSITIVITY_IDENTIFIER],
    )
    preferred = material.get_eos_record(PREFERRED_IDENTIFIER)
    sensitivity = material.get_eos_record(SENSITIVITY_IDENTIFIER)
    assert preferred.pressure(140.0) == pytest.approx(88.3150121125, abs=5.0e-10)
    assert sensitivity.pressure(140.0) == pytest.approx(88.3722020596, abs=5.0e-10)

    pressure_grid = np.linspace(47.3, 196.3, 11)
    for record in (preferred, sensitivity):
        volumes = record.volume(pressure_grid)
        assert record.pressure(volumes) == pytest.approx(pressure_grid, rel=1.0e-10)

    pt_rmse, pt_maximum = platinum_reduction_metrics()
    assert pt_rmse == pytest.approx(0.2535081479, abs=5.0e-10)
    assert pt_maximum == pytest.approx(0.4452120107, abs=5.0e-10)
    available_records = set(list_eos_record_documents())
    for record in document["eos_records"]:
        method = record["pressure_calibration"]["methods"][0]
        assert method["reference_eos_record"] == "platinum_holmes_1989_vinet_1"
        assert method["reference_eos_record"] in available_records


def test_kubo_mggeo3_reproduction_doc_disposes_every_same_doi_candidate():
    path = (
        ROOT
        / "docs"
        / "literature-reproductions"
        / "kubo-2006-mggeo3-post-perovskite.md"
    )
    content = path.read_text(encoding="utf-8")
    assert all(identifier in content for identifier in LITCURATE_IDENTIFIERS)
    assert "accepted, preferred/default" in content
    assert "accepted, nonpreferred sensitivity" in content
    assert "held / not refittable" in content
    assert content.count("rejected as a Kubo-source record") == 4
