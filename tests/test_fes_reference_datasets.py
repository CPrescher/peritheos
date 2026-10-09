"""Reference data are accessible without conflating author/source versions."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from jsonschema import validate

from peritheos import Material, eosmat_schema, get_material_document
from scripts.audit_sata_2010_fes_vi import observations
from scripts.bundle_fes_reference_datasets import AUTHOR, DATA, PREFIX, SOURCE


@pytest.fixture(scope="module")
def material():
    document = get_material_document("fes_vi")
    validate(document, eosmat_schema())
    return Material.from_eosmat(document, require_primary_validation=False)


def test_original_files_packaged_byte_for_byte_and_checksums_verified(material):
    manifest = json.loads((SOURCE / "manifest.json").read_text(encoding="utf-8"))
    for source in manifest["sources"]:
        original = (
            Path("docs/data/fes-sata-2010-table1.csv")
            if source["filename"].endswith(".csv")
            else AUTHOR / source["filename"]
        )
        packaged = SOURCE / source["filename"]
        assert packaged.read_bytes() == original.read_bytes()
        assert hashlib.sha256(packaged.read_bytes()).hexdigest() == source["sha256"]
    assert len(manifest["datasets"]) == 3
    for identifier in manifest["datasets"]:
        dataset = material.get_dataset(identifier)
        assert dataset.checksum_verified
        assert dataset.used_by_eos_records == ()


def test_selected_input_and_cold_subset_exclude_quenched_points(material):
    short = material.get_dataset(PREFIX + "without_additional_cold_pvt")
    literature = material.get_dataset(PREFIX + "literature_cold_pv")
    assert len(short) == 167 and len(literature) == 21
    ids = short["observation_id"].tolist()
    positions = [ids.index(x) for x in literature["observation_id"]]
    for column in [
        "pressure_gpa",
        "temperature_k",
        "molar_volume_cm3_mol",
        "pressure_uncertainty_gpa",
        "temperature_uncertainty_k",
        "molar_volume_uncertainty_cm3_mol",
    ]:
        np.testing.assert_array_equal(short[column][positions], literature[column])
    assert np.sum(short["source_group"] == "morard_thermal") == 146
    assert "additional_cold" not in short["source_group"]
    assert not np.any(
        (short["temperature_k"] == 300) & (short["pressure_uncertainty_gpa"] == 0.0056)
    )
    for dataset in (short, literature):
        assert dataset.metadata["provenance"]["final_author_fit_selection_confirmed"]
        assert (
            "personal communication"
            in dataset.metadata["provenance"]["selection_provenance"]
        )
    identifiers = {x["identifier"] for x in material.datasets}
    assert PREFIX + "combined_pvt" not in identifiers
    assert PREFIX + "additional_cold_pv" not in identifiers
    for path in (
        DATA / "fes-morard-2026-author-combined.csv",
        DATA / "fes-morard-2026-author-additional-cold.csv",
        SOURCE / "FeS6EOSfittot-SataOhfuji.dat",
    ):
        assert not path.exists()
    pv = literature.as_pressure_volume()
    assert pv.volume_unit == "cm^3/mol"
    assert pv.pressure_uncertainty_role == "uncertainty"
    assert pv.pressure_sigma is None and pv.volume_sigma is None
    np.testing.assert_allclose(pv.pressure_uncertainty, 0.01 * pv.pressure, atol=1e-12)
    np.testing.assert_allclose(pv.volume_uncertainty, 0.005 * pv.volume, atol=6e-10)


def test_sata_original_errors_and_pressure_are_not_replaced_by_author_input(material):
    original = material.get_dataset("fes_sata_2010_table1_vi_pv")
    reference = [x for x in observations() if x["included"]]
    assert len(original) == 13
    for column, field in [
        ("pressure_gpa", "pressure_gpa"),
        ("pressure_uncertainty_gpa", "error_pressure_gpa"),
        ("volume_angstrom3", "volume_angstrom3"),
        ("volume_uncertainty_angstrom3", "error_volume_angstrom3"),
    ]:
        np.testing.assert_array_equal(original[column], [x[field] for x in reference])
    assert original["pressure_raw"].tolist() == [x["pressure_raw"] for x in reference]
    assert original["run"].tolist() == [x["run"] for x in reference]
    assert not original.find_columns(role="uncertainty", quantity="temperature")
    combined = material.get_dataset(PREFIX + "literature_cold_pv")
    matched = np.isfinite(combined["matched_sata2010_row"])
    assert np.sum(matched) == 13 and np.sum(~matched) == 8
    row = combined["matched_sata2010_row"] == 3
    assert combined["pressure_gpa"][row].tolist() == [101.2]
    assert original["pressure_gpa"][2] == 101.1


def test_eos_source_parameters_remain_published_and_deferred(material):
    document = material.to_eosmat()
    source = get_material_document("fes_vi")
    assert document["datasets"] == source["datasets"]
    assert (
        document["eos_records"][0]["author_input_followup"]
        == source["eos_records"][0]["author_input_followup"]
    )
    assert (
        Material.from_eosmat(document, require_primary_validation=False).to_eosmat()
        == document
    )
    record = document["eos_records"][0]
    assert record["scientific_validation"]["status"] == "deferred"
    assert record["eos"]["parameters"]["K0"] == 115.5
    assert record["thermal"]["parameters"]["gamma0"] == 2.42
    assert record["thermal"]["parameters"]["Tr"] == 300
    saved = record["author_input_followup"]["saved_author_parameters"]
    assert saved["gamma0"] == 2.41671 and saved["Tr_k"] == 298
    assert saved["q_compromise"] is False
    assert record["fit_datasets"] == ["fes_morard_2026_table_s1_pvt"]
    assert len(document["datasets"]) == 5
    for metadata in document["datasets"]:
        assert (DATA.parent / metadata["resource"]["path"]).is_file()
