import importlib.util
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import eosmat_schema

ROOT = Path(__file__).parents[1]

ZHANG = {
    "bridgmanite_zhang_wentzcovitch_2022_phq_lda_300k_bm3",
    "bridgmanite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3",
    "mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_lda_300k_bm3",
    "mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3",
}
DRIVER = {
    "alpha_quartz_driver_2010_qmc_300k_vinet",
    "sio2_stv_andr_driver_2010_qmc_300k_vinet",
    "seifertite_driver_2010_qmc_300k_vinet",
}
HOLLAND = {
    "bridgmanite_holland_2013_mpv_modified_tait",
    "fesio3_bridgmanite_holland_2013_fpv_modified_tait",
    "al2o3_perovskite_holland_2013_apv_modified_tait",
    "ca_perovskite_holland_2013_cpv_modified_tait",
    "mgo_holland_2013_per_modified_tait",
    "feo_holland_2013_fper_modified_tait",
    "sio2_stv_andr_holland_2013_stv_modified_tait",
}
MATERIALS = {
    "bridgmanite",
    "mgsio3_post_perovskite",
    "alpha_quartz",
    "sio2_stv_andr",
    "seifertite",
    "fesio3_bridgmanite",
    "al2o3_perovskite",
    "ca_perovskite",
    "mgo",
    "feo",
}


def _documents():
    for material in MATERIALS:
        path = ROOT / "peritheos/data/materials" / f"{material}.eosmat"
        yield material, json.loads(path.read_text(encoding="utf-8"))


def _load_script(filename):
    path = ROOT / "scripts" / filename
    spec = importlib.util.spec_from_file_location(filename, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_exact_second_tranche_production_set_and_doi_counts():
    by_doi = {}
    for _, document in _documents():
        for record in document["eos_records"]:
            doi = record.get("reference", {}).get("doi", "").lower()
            if doi in {
                "10.1103/physrevb.106.054103",
                "10.1073/pnas.0912130107",
                "10.1093/petrology/egt035",
                "10.1002/2016gl069332",
            }:
                by_doi.setdefault(doi, set()).add(record["identifier"])

    assert by_doi["10.1103/physrevb.106.054103"] == ZHANG
    assert by_doi["10.1073/pnas.0912130107"] == DRIVER
    assert by_doi["10.1093/petrology/egt035"] == HOLLAND
    assert "10.1002/2016gl069332" not in by_doi
    assert len(ZHANG | DRIVER | HOLLAND) == 14


def test_new_volume_bases_have_explicit_positive_molar_mass():
    expected = ZHANG | DRIVER | HOLLAND
    found = set()
    for _, document in _documents():
        for record in document["eos_records"]:
            if record["identifier"] not in expected:
                continue
            found.add(record["identifier"])
            assert record["volume_basis"]["formula_units"] > 0
            assert record["volume_basis"]["molar_mass_g_mol"] > 0
    assert found == expected


def test_new_records_satisfy_normative_json_schema():
    schema = eosmat_schema()
    record_validator = Draft202012Validator(
        {
            "$schema": schema["$schema"],
            "$ref": "#/$defs/eos_record",
            "$defs": schema["$defs"],
        }
    )
    expected = ZHANG | DRIVER | HOLLAND
    found = set()
    for material, document in _documents():
        for record in document["eos_records"]:
            if record["identifier"] not in expected:
                continue
            found.add(record["identifier"])
            errors = list(record_validator.iter_errors(record))
            assert not errors, (material, [error.message for error in errors])
    assert found == expected


def test_zhang_driver_reproduction_roundtrips_all_seven_records():
    module = _load_script("reproduce_zhang_driver_thermoelastic_eos.py")
    result = module.reproduce()
    assert set(result["records"]) == ZHANG | DRIVER
    assert (
        result["records"]["alpha_quartz_driver_2010_qmc_300k_vinet"][
            "checkpoint_pressure_gpa"
        ]
        == 20.0
    )


def test_holland_reproduction_converts_and_executes_all_seven_records():
    module = _load_script("reproduce_holland_2013_modified_tait.py")
    result = module.reproduce()
    assert set(result["records"]) == {"mpv", "fpv", "apv", "cpv", "per", "fper", "stv"}
    assert result["records"]["fpv"]["source_parameters"]["V0_j_bar_mol"] == 2.548
    assert result["records"]["stv"]["roundtrip_pressure_gpa"] == pytest.approx(30.0)


@pytest.mark.parametrize(
    ("filename", "candidate_ids"),
    [
        (
            "zhang-wentzcovitch-2022-anharmonic-pv-ppv.md",
            [
                "litcurate_de4c3c7e3aa017f4",
                "litcurate_720fd80b9e9294bd",
                "litcurate_081bf45b6736075c",
                "litcurate_faa7c12c60bab5df",
                "litcurate_e591fafb49a5023c",
                "litcurate_5c4acb917edcc312",
            ],
        ),
        (
            "driver-2010-qmc-silica.md",
            [
                "litcurate_fc269c4a99bde326",
                "litcurate_c4d862798e2cba40",
                "litcurate_3bebcc8155727eba",
                "litcurate_3181bf31d9a5246c",
            ],
        ),
        (
            "holland-2013-ncfmas-modified-tait.md",
            [
                "litcurate_2d2e9920973cac04",
                "litcurate_30a0eba3590777a5",
                "litcurate_6f93835769d546c6",
                "litcurate_aee72f71b40279a3",
                "litcurate_ec9908a767c7323f",
                "litcurate_2d3f926a8e10b0e0",
                "litcurate_555f952c38b721f1",
                "litcurate_d9edde4e235c2ee1",
            ],
        ),
        (
            "shukla-2016-ferric-al-bridgmanite.md",
            [
                "litcurate_ae6fe425cad16fdd",
                "litcurate_07a276bc1acd6da1",
                "litcurate_fc233ec735390ab4",
                "litcurate_7b772f4f43012b67",
                "litcurate_477fce3ed82ea93a",
                "litcurate_668d496ea6b360eb",
            ],
        ),
    ],
)
def test_every_litcurate_candidate_is_explicitly_documented(filename, candidate_ids):
    text = (ROOT / "docs/literature-reproductions" / filename).read_text(
        encoding="utf-8"
    )
    for candidate_id in candidate_ids:
        assert text.count(candidate_id) == 1
