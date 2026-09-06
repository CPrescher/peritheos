import csv
import json
from importlib import resources

import numpy as np
import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from scripts.reproduce_deltacodesdft_metal_eos import reproduce


def source_inputs():
    root = resources.files("peritheos").joinpath("data", "datasets")
    source = json.loads(
        root.joinpath("deltacodesdft-metal-eos-source.json").read_text(encoding="utf-8")
    )
    with root.joinpath("deltacodesdft-metal-eos-parameters.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    return source, rows


def test_delta_project_metal_selection_and_archive_transcription():
    source, rows = source_inputs()

    assert len(rows) == 50
    assert len({row["element"] for row in rows}) == 48
    assert sum(row["code"] == "WIEN2k" for row in rows) == 48
    assert sum(row["code"] == "FLEUR" for row in rows) == 2
    assert {row["element"] for row in rows if row["code"] == "FLEUR"} == {
        "Al",
        "Fe",
    }
    assert source["archive"]["doi"] == "10.24435/materialscloud:5e-mv"
    assert source["archive"]["related_publication_doi"] == ("10.1126/science.aad3000")
    assert source["archive"]["license"] == "CC BY 4.0"
    assert source["archive"]["files"] == {
        "README.txt": "4bcf23218baba567f74997a50033c5b6d08dc9661f882b0ac7f8555505f74eb4",
        "WIEN2k.txt": "d8fb392d8e8a6325175d74a06565d8cfe0a04869deb3bbae7c22cfb5c6a8b977",
        "history/history/FLEUR.txt": "dd924c02bdb31f4a355c55d62f8cc2ef96a0d710e27c0d45d7cb5e9aed09a172",
        "eosfit.py": "c8bddeb1493cbf2a925382f80aa13b29c7fa1913a6ff68b8e862ccb55e4b8b3c",
    }


@pytest.mark.parametrize("ratio", [0.94, 1.0, 1.06])
def test_delta_project_records_are_exact_and_executable(ratio):
    source, rows = source_inputs()
    source_rows = {(row["code"], row["element"]): row for row in source["records"]}

    for row in rows:
        material_id = row["material_identifier"]
        code = row["code"]
        identifier = f"{material_id}_lejaeghere_2016_{code.lower()}_pbe_bm3"
        document = get_material_document(material_id)
        validate_eosmat_document(document)
        stored = next(
            item for item in document["eos_records"] if item["identifier"] == identifier
        )
        expected = source_rows[(code, row["element"])]
        assert stored["eos"]["parameters"] == {
            "V0": float(row["executable_v0_a3_cell"]),
            "K0": expected["K0_gpa"],
            "K0_prime": expected["K0_prime"],
        }
        assert stored["volume_basis"]["formula_units"] == float(
            row["formula_units_per_cell"]
        )
        assert stored["validity"]["pressure_gpa"][0] <= 0.0
        assert stored["validity"]["pressure_gpa"][1] >= 0.0
        loaded = Material.from_eosmat(
            document, record_identifiers=[identifier]
        ).eos_records[0]
        pressure = loaded.pressure(ratio * loaded.reference_volume, 0.0)
        assert np.isfinite(pressure)
        if ratio == 1.0:
            assert pressure == pytest.approx(0.0, abs=1.0e-13)


def test_delta_project_special_phase_assignments_are_not_conflated():
    lithium = get_material_document("lithium_9r")
    sodium = get_material_document("sodium_9r")
    manganese = get_material_document("manganese_fcc_afm")

    assert lithium["formula_units_per_cell"] == 3
    assert sodium["formula_units_per_cell"] == 3
    assert "9R" in lithium["phase"]
    assert "9R" in sodium["phase"]
    assert "antiferromagnetic fcc" in manganese["phase"]
    record = manganese["eos_records"][0]
    assert record["eos"]["parameters"]["K0_prime"] == -0.21
    assert "must not be extrapolated" in record["notes"]
    assert "alpha" not in manganese["identifier"]


def test_delta_project_code_and_magnetic_provenance_is_explicit():
    for material_id, magnetic_state in {
        "chromium": "antiferromagnetic",
        "manganese_fcc_afm": "antiferromagnetic",
        "fe": "ferromagnetic",
        "cobalt_hcp": "ferromagnetic",
        "nickel": "ferromagnetic",
    }.items():
        document = get_material_document(material_id)
        records = [
            record
            for record in document["eos_records"]
            if "_lejaeghere_2016_" in record["identifier"]
        ]
        assert records
        for record in records:
            assert record["reference"]["doi"] == "10.1126/science.aad3000"
            assert magnetic_state in record["sample_composition"]
            assert magnetic_state in record["parameter_provenance"]["method"]
            assert "scalar-relativistic" in record["parameter_provenance"]["method"]

    source, rows = source_inputs()
    assert source["selection"]["code_methods"] == {
        "WIEN2k": "WIEN2k 13.1 all-electron APW+lo, scalar-relativistic PBE",
        "FLEUR": "FLEUR 0.26 all-electron LAPW (+lo), scalar-relativistic PBE",
    }
    assert {row["reference_doi"] for row in rows} == {"10.1126/science.aad3000"}


def test_delta_project_independent_bm3_reproduction():
    result = reproduce()

    assert result["records"] == 50
    assert result["elements"] == 48
    assert result["largest_parameter_error"] == 0.0
    assert result["largest_p0_gpa"] == 0.0
    assert result["largest_independent_pressure_error_gpa"] < 1.0e-10
