import copy
import csv
import hashlib
import io
import json
import math
import re
from importlib import resources
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    EOSMAT_FORMAT,
    EOSMAT_FORMAT_VERSION,
    Material,
    eosmat_schema,
    get_eos_record_document,
    get_material_document,
    list_eos_record_documents,
    list_material_documents,
    list_materials,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
    validate_pressure_calibration_references,
)
from scripts.import_dioptas_eos_database import migrate_document

NINETY_FIVE_PERCENT_INTERVAL_RECORDS = (
    ("alumina", "alumina_dewaele_2013_vinet_1"),
    ("aluminum", "aluminum_dewaele_2004_vinet_1"),
    ("copper", "copper_dewaele_2004_vinet_1"),
    ("diamond", "diamond_dewaele_2008_vinet_2"),
    ("lif_b1", "lif_b1_dewaele_2019_vinet_1"),
    ("nickel", "nickel_dewaele_2008_vinet_1"),
    ("rhenium", "rhenium_anzellini_2014_vinet_1"),
    ("silicon_v", "silicon_v_anzellini_2019_vinet_1"),
    ("silicon_vii", "silicon_vii_anzellini_2019_vinet_1"),
    ("silicon_x", "silicon_x_anzellini_2019_vinet_1"),
    ("silver", "silver_dewaele_2008_vinet_1"),
    ("titanium_alpha", "titanium_alpha_dewaele_2015_vinet_1"),
    ("titanium_omega", "titanium_omega_dewaele_2015_vinet_1"),
    ("tungsten", "tungsten_dewaele_2004_vinet_2"),
)

SHEN_SMITH_2026_RECORDS = (
    ("fe", "fe_shen_2026_vinet_1", 23.506, 116.15919903947005, 2.5, 0.5),
    ("gold", "gold_shen_2026_vinet_3", 67.792, 128.41042483384334, 0.3, 0.01),
    ("iron", "iron_shen_2026_vinet_2", 22.352, 123.15954992657934, 0.6, 0.03),
    ("mgo", "mgo_shen_2026_vinet_3", 74.636, 95.1759313304533, 0.7, 0.03),
    (
        "molybdenum",
        "molybdenum_shen_2026_vinet_1",
        31.14,
        152.7396582018993,
        0.5,
        0.02,
    ),
    ("nacl_b1", "nacl_b1_shen_2026_vinet_1", 179.41, 16.415871613452587, 0.3, 0.06),
    ("nacl_b2", "nacl_b2_shen_2026_vinet_2", 41.35, 18.484248202987555, 0.11, 0.02),
    (
        "platinum",
        "platinum_shen_2026_vinet_2",
        60.364,
        197.26236099266487,
        0.5,
        0.02,
    ),
    (
        "tantalum",
        "tantalum_shen_2026_vinet_2",
        36.14,
        107.10106281373046,
        0.4,
        0.02,
    ),
    (
        "tungsten",
        "tungsten_shen_2026_vinet_3",
        31.704,
        178.85867602216723,
        0.8,
        0.03,
    ),
)


def test_complete_migrated_dioptas_library_is_bundled_and_valid():
    identifiers = list_material_documents()
    documents = [get_material_document(identifier) for identifier in identifiers]

    assert len(identifiers) == 116
    assert len(set(identifiers)) == 116
    assert sum(len(document["eos_records"]) for document in documents) == 163
    assert all(document["eos_records"] for document in documents)
    assert all(document["format"] == EOSMAT_FORMAT for document in documents)
    assert all(
        document["format_version"] == EOSMAT_FORMAT_VERSION for document in documents
    )
    assert all(document["identifier"] in identifiers for document in documents)
    assert all(
        document["units"]
        == {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        }
        for document in documents
    )
    json.dumps(documents, allow_nan=False)


def test_migrated_records_have_completed_primary_source_audit():
    records = [
        record
        for identifier in list_material_documents()
        for record in get_material_document(identifier)["eos_records"]
    ]

    assert len({record["identifier"] for record in records}) == 163
    statuses = [record["scientific_validation"]["status"] for record in records]
    assert set(statuses) == {"primary_source_validated"}
    assert statuses.count("primary_source_validated") == 163
    audit_dates = {
        record["identifier"]: record["scientific_validation"]["audit_date"]
        for record in records
    }
    current_audit_identifiers = {
        "b4c_somayazulu_2023_berman_2",
        "b4c_somayazulu_2023_berman_refit",
        "b4c_somayazulu_2023_bm3_1",
        "ca_perovskite_caracas_2005_bm3_3",
        "ca_perovskite_kawai_2014_vinet_mgd_3",
        "diamond_benedict_2014_dewaele_anchored",
        "diamond_correa_2008_dewaele_anchored",
        "gold_dewaele_2004_vinet_5",
        "gold_fratanduono_2021_vinet_7",
        "gold_takemura_2008_vinet_6",
        "goethite_gleason_2008_bm3_1",
        "kcl_b2_chidester_2021_bm3_5",
        "kcl_b2_tateno_2019_vinet_4",
        "mgo_b1_tange_2009_vinet",
        "mgo_b1_luo_2023_vinet_thermal_5",
        "molybenum_carbide_mo2c_haines_2001_bm3_refit",
        "platinum_dorogokupets_oganov_2007_vinet_4",
        "neon_fcc_hemley_1989_bm3_refit",
        "rbcl_b2_campbell_1994_bm3_1",
        "mgo_dewaele_2000_bm3_mgd_5",
    }
    latest_audit_identifiers = {
        "ca_perovskite_caracas_2005_bm3_3",
        "ca_perovskite_kawai_2014_vinet_mgd_3",
        "molybenum_carbide_mo2c_haines_2001_bm3_refit",
        "neon_fcc_hemley_1989_bm3_refit",
        "kcl_b2_tateno_2019_vinet_4",
        "kcl_b2_chidester_2021_bm3_5",
        "goethite_gleason_2008_bm3_1",
        "rbcl_b2_campbell_1994_bm3_1",
        "mgo_b1_luo_2023_vinet_thermal_5",
        "mgo_dewaele_2000_bm3_mgd_5",
    }
    assert {
        audit_dates[identifier] for identifier in latest_audit_identifiers
    } == {"2026-09-04"}
    goethite = next(
        record
        for record in records
        if record["identifier"] == "goethite_gleason_2008_bm3_1"
    )
    assert goethite["scientific_validation"]["usage_recommendation"] == (
        "not_recommended_for_quantitative_use"
    )
    assert {
        audit_dates[identifier]
        for identifier in current_audit_identifiers
        if identifier not in latest_audit_identifiers
    } == {"2026-09-03"}
    assert {
        date
        for identifier, date in audit_dates.items()
        if identifier not in current_audit_identifiers
    } == {"2026-09-01"}
    assert all(
        record["scientific_validation"]["primary_source_check"] for record in records
    )
    assert all(
        len({item["path"] for item in record.get("audit_corrections", [])})
        == len(record.get("audit_corrections", []))
        for record in records
    )
    migrated_records = [
        record
        for record in records
        if "migration_source" in record["scientific_validation"]
    ]
    native_records = [
        record
        for record in records
        if "migration_source" not in record["scientific_validation"]
    ]
    assert {record["identifier"] for record in native_records} == {
        "aragonite_martinez_1996_bm2_2",
        "b4c_somayazulu_2023_berman_2",
        "b4c_somayazulu_2023_berman_refit",
        "ca_perovskite_caracas_2005_bm3_3",
        "ca_perovskite_kawai_2014_vinet_mgd_3",
        "diamond_benedict_2014_double_debye_4",
        "diamond_benedict_2014_dewaele_anchored",
        "diamond_correa_2008_double_debye_log_moment_5",
        "diamond_correa_2008_dewaele_anchored",
        "kcl_b2_dewaele_2012_vinet_3",
        "gold_dewaele_2004_vinet_5",
        "gold_fratanduono_2021_vinet_7",
        "gold_takemura_2008_vinet_6",
        "kcl_b2_chidester_2021_bm3_5",
        "mgo_b1_luo_2023_vinet_thermal_5",
        "mgo_b1_tange_2009_vinet",
        "mgo_dewaele_2000_bm3_mgd_5",
        "molybenum_carbide_mo2c_haines_2001_bm3_refit",
        "kcl_b2_tateno_2019_vinet_4",
        "platinum_dorogokupets_oganov_2007_vinet_4",
        "neon_fcc_hemley_1989_bm3_refit",
        "rbcl_b2_campbell_1994_bm3_1",
    }
    assert {
        record["scientific_validation"]["migration_source"]["version"]
        for record in migrated_records
    } == {"0.10.0"}
    assert {
        record["scientific_validation"]["migration_source"]["commit"]
        for record in migrated_records
    } == {"5a8bfd81d10bfab3499039603380aae34576d60a"}
    assert all(record["eos"]["model"] for record in records)
    assert all(
        not record.get("thermal") or record["thermal"]["model"] for record in records
    )


def test_primary_source_audit_report_covers_every_migrated_record():
    report = json.loads(
        resources.files("peritheos.data")
        .joinpath("primary-source-audit.json")
        .read_text(encoding="utf-8")
    )
    bundled_ids = {
        record["identifier"]
        for identifier in list_material_documents()
        for record in get_material_document(identifier)["eos_records"]
    }

    assert report["summary"] == {
        "records": 163,
        "primary_source_validated": 163,
    }
    assert report["audit_date"] == "2026-09-04"
    assert {entry["record"] for entry in report["records"]} == bundled_ids
    assert len(report["records"]) == len(bundled_ids)


def test_bundled_primary_dataset_resources_match_metadata():
    data_root = resources.files("peritheos.data")
    checked = 0

    for material_identifier in list_material_documents():
        document = get_material_document(material_identifier)
        for dataset in document.get("datasets", []):
            resource = dataset.get("resource")
            if resource is None:
                continue

            payload = data_root.joinpath(resource["path"]).read_bytes()
            assert hashlib.sha256(payload).hexdigest() == resource["sha256"]

            rows = list(csv.reader(io.StringIO(payload.decode("utf-8"))))
            assert len(rows[0]) == len(dataset["columns"])
            assert len(rows) > 1
            assert all(len(row) == len(dataset["columns"]) for row in rows[1:])
            checked += 1

    assert checked >= 39


def test_ono_cubic_sno2_primary_data_transcription_is_complete():
    documents = [
        get_material_document("sno2_cubic_27gpa"),
        get_material_document("sno2_pa_3_at_48gpa"),
    ]
    datasets = [
        next(
            item
            for item in document["datasets"]
            if item["identifier"] == "sno2_ono_2000_table2_pvt"
        )
        for document in documents
    ]

    assert datasets[0]["resource"] == datasets[1]["resource"]
    assert datasets[0]["used_by_eos_records"] == ["sno2_cubic_27gpa_ono_2000_bm3_1"]
    assert datasets[1]["used_by_eos_records"] == ["sno2_pa_3_at_48gpa_ono_2000_bm3_1"]

    payload = (
        resources.files("peritheos.data")
        .joinpath(datasets[0]["resource"]["path"])
        .read_text(encoding="utf-8")
    )
    rows = list(csv.DictReader(io.StringIO(payload)))
    assert len(rows) == 40
    assert sum(row["experiment"] == "S145" for row in rows) == 34
    assert sum(row["experiment"] == "S192" for row in rows) == 6
    assert min(float(row["pressure_gpa"]) for row in rows) == pytest.approx(16.09)
    assert max(float(row["pressure_gpa"]) for row in rows) == pytest.approx(28.85)
    representative = next(
        row
        for row in rows
        if row["temperature_k"] == "1400" and row["pressure_gpa"] == "28.09"
    )
    assert float(representative["unit_cell_volume_a3"]) == pytest.approx(121.16)
    assert float(representative["pressure_uncertainty_gpa"]) == pytest.approx(0.03)


def test_shen_smith_simultaneous_volume_dataset_is_complete():
    record_by_material = {
        "fe": "fe_shen_2026_vinet_1",
        "gold": "gold_shen_2026_vinet_3",
        "iron": "iron_shen_2026_vinet_2",
        "mgo": "mgo_shen_2026_vinet_3",
        "molybdenum": "molybdenum_shen_2026_vinet_1",
        "nacl_b1": "nacl_b1_shen_2026_vinet_1",
        "nacl_b2": "nacl_b2_shen_2026_vinet_2",
        "platinum": "platinum_shen_2026_vinet_2",
        "tantalum": "tantalum_shen_2026_vinet_2",
        "tungsten": "tungsten_shen_2026_vinet_3",
    }
    datasets = []
    for material, record in record_by_material.items():
        document = get_material_document(material)
        dataset = next(
            item
            for item in document["datasets"]
            if item["identifier"] == "shen_smith_2026_table_s1_simultaneous_volumes"
        )
        assert dataset["used_by_eos_records"] == [record]
        datasets.append(dataset)

    assert all(dataset["resource"] == datasets[0]["resource"] for dataset in datasets)

    payload = (
        resources.files("peritheos.data")
        .joinpath(datasets[0]["resource"]["path"])
        .read_text(encoding="utf-8")
    )
    rows = list(csv.DictReader(io.StringIO(payload)))
    assert len(rows) == 3511
    assert {
        int(row["run_number"]) for row in rows if row["experiment"] == "DAC-1"
    } == set(range(1, 229))
    assert {
        int(row["run_number"]) for row in rows if row["experiment"] == "DAC-2"
    } == set(range(1, 195))
    assert sum(row["phase"] == "Cu" for row in rows) == 422
    assert sum(row["phase"] == "NaCl-B1" for row in rows) == 67
    assert sum(row["phase"] == "NaCl-B2" for row in rows) == 131
    assert (
        sum(
            row["phase"] == "Pt" and row["measurement_sequence"] == "first"
            for row in rows
        )
        == 421
    )
    assert (
        sum(
            row["phase"] == "Pt" and row["measurement_sequence"] == "last"
            for row in rows
        )
        == 422
    )

    first = rows[0]
    assert first == {
        "experiment": "DAC-1",
        "run_number": "1",
        "phase": "Pt",
        "measurement_sequence": "first",
        "unit_cell_volume_a3": "59.54115",
        "unit_cell_volume_standard_error_a3": "0.01721",
    }
    last = rows[-1]
    assert last == {
        "experiment": "DAC-2",
        "run_number": "194",
        "phase": "Pt",
        "measurement_sequence": "last",
        "unit_cell_volume_a3": "48.23100",
        "unit_cell_volume_standard_error_a3": "0.03486",
    }


def test_campbell_cscl_and_rbcl_table_blocks_are_separate_and_complete():
    expected = {
        "cscl": {
            "dataset": "cscl_campbell_1994_table1_compression",
            "record": "cscl_campbell_1994_bm3_1",
            "rows": 13,
            "first": ("6.97", "3.8405"),
            "last": ("28.7", "3.5308"),
        },
        "rbcl": {
            "dataset": "rbcl_campbell_1994_table1_compression",
            "record": "rbcl_b2_campbell_1994_bm3_1",
            "rows": 24,
            "first": ("1.11", "3.8799"),
            "last": ("32.3", "3.3301"),
        },
    }
    resource_paths = set()

    for material_identifier, values in expected.items():
        document = get_material_document(material_identifier)
        dataset = next(
            item
            for item in document["datasets"]
            if item["identifier"] == values["dataset"]
        )
        resource_paths.add(dataset["resource"]["path"])
        payload = (
            resources.files("peritheos.data")
            .joinpath(dataset["resource"]["path"])
            .read_text(encoding="utf-8")
        )
        rows = list(csv.DictReader(io.StringIO(payload)))

        assert len(rows) == values["rows"]
        assert (rows[0]["pressure_gpa"], rows[0]["lattice_parameter_a"]) == values[
            "first"
        ]
        assert (rows[-1]["pressure_gpa"], rows[-1]["lattice_parameter_a"]) == values[
            "last"
        ]
        assert dataset["used_by_eos_records"] == [values["record"]]
        assert {
            column["role"]
            for column in dataset["columns"]
            if column["name"].endswith("uncertainty_gpa")
            or column["name"].endswith("_uncertainty")
        } == {"standard_deviation"}

    assert len(resource_paths) == 2

    cscl = get_material_document("cscl")
    yagi = next(
        item
        for item in cscl["datasets"]
        if item["identifier"] == "cscl_yagi_1978_table1_compression"
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(yagi["resource"]["path"])
        .read_text(encoding="utf-8")
    )
    rows = list(csv.DictReader(io.StringIO(payload)))
    correction = (4.118 / 4.123) ** 3
    assert len(rows) == 9
    assert [float(row["pressure_kbar"]) for row in rows] == list(
        np.arange(10.0, 100.0, 10.0)
    )
    assert [float(row["volume_ratio_reported"]) for row in rows] == pytest.approx(
        [0.9526, 0.9159, 0.8860, 0.8609, 0.8391, 0.8201, 0.8031, 0.7878, 0.7739]
    )
    assert [float(row["volume_ratio_corrected"]) for row in rows] == pytest.approx(
        [float(row["volume_ratio_reported"]) * correction for row in rows],
        abs=5.0e-13,
    )
    assert cscl["eos_records"][0]["fit_datasets"] == [
        "cscl_campbell_1994_table1_compression",
        "cscl_yagi_1978_table1_compression",
    ]


def test_chidester_kcl_dataset_and_new_pressure_markers_match_primary_values():
    document = get_material_document("kcl")
    dataset = next(
        item
        for item in document["datasets"]
        if item["identifier"] == "kcl_chidester_2021_supplemental_pvt"
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_text(encoding="utf-8")
    )
    rows = list(csv.DictReader(io.StringIO(payload)))
    assert len(rows) == 155
    assert rows[0]["Sample/pattern number"] == "B1_578"
    assert rows[-1]["Sample/pattern number"] == "B86_382"
    assert dataset["notes"].endswith("Pt-derived rather than as a raw Pt-volume pair.")

    records = {record["identifier"]: record for record in document["eos_records"]}
    tateno = records["kcl_b2_tateno_2019_vinet_4"]
    assert tateno["eos"]["parameters"] == {
        "V0": 54.5,
        "K0": 18.3,
        "K0_prime": 5.6,
    }
    assert tateno["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 235.0,
        "gamma0": 2.3,
        "q": 0.8,
        "n": 2,
    }
    assert tateno["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert tateno["thermal"]["parameter_errors"]["gamma0"] == pytest.approx(0.2)
    tateno_dataset = next(
        item
        for item in document["datasets"]
        if item["identifier"] == "kcl_tateno_2019_table_s1_pvt"
    )
    tateno_payload = (
        resources.files("peritheos.data")
        .joinpath(tateno_dataset["resource"]["path"])
        .read_text(encoding="utf-8")
    )
    tateno_rows = list(csv.DictReader(io.StringIO(tateno_payload)))
    assert len(tateno_rows) == 39
    six_gpa = next(
        row for row in tateno_rows if row["run"] == "3" and row["pressure_gpa"] == "6.0"
    )
    assert six_gpa["kcl_unit_cell_volume_a3"] == "45.0463"
    assert "official MSA XLSX deposit" in tateno_dataset["notes"]
    assert (
        tateno["pressure_calibration"]["methods"][0]["reference_eos_record"]
        == "platinum_sokolova_2013_holzapfel_3"
    )

    chidester = records["kcl_b2_chidester_2021_bm3_5"]
    assert chidester["eos"]["parameters"] == {
        "V0": pytest.approx(53.137250149563094),
        "K0": 24.0,
        "K0_prime": 4.56,
    }
    assert chidester["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 235.0,
        "gamma0": 2.9,
        "q": 1.0,
        "n": 2,
    }
    assert chidester["thermal"]["debye_temperature_law"] == (
        "integrated_gruneisen"
    )
    assert chidester["fit_datasets"] == [
        "kcl_dewaele_2012_table1_compression",
        "kcl_chidester_2021_supplemental_pvt",
    ]
    assert chidester["scientific_validation"]["primary_data_check"][
        "dataset_identifiers"
    ] == chidester["fit_datasets"]
    dewaele_dataset = next(
        item
        for item in document["datasets"]
        if item["identifier"] == "kcl_dewaele_2012_table1_compression"
    )
    assert chidester["identifier"] in dewaele_dataset["used_by_eos_records"]
    assert chidester["experimental_configuration"] == {
        "geometry": "laser-heated diamond-anvil cell",
        "sample_and_laser_absorber": "approximately 5 micrometer Pt foil",
        "pressure_reference": "simultaneously diffracting Pt foil",
        "pressure_medium": "KCl",
        "thermal_insulator": "KCl",
        "arrangement": (
            "Pt foil secured between two approximately 10 micrometer layers of "
            "dried KCl"
        ),
        "temperature_convention": (
            "KCl diffraction samples a thermal gradient between the hot Pt surface "
            "and the approximately 295 K diamond anvils; the EOS uses the authors' "
            "estimated average KCl temperature"
        ),
    }

    gold = get_eos_record_document("gold_fratanduono_2021_vinet_7")
    assert gold["eos"]["parameters"] == {
        "V0": pytest.approx(4.0 * 16.929),
        "K0": 170.9,
        "K0_prime": 5.88,
    }


def test_pressure_calibration_audit_covers_every_eos_record_and_links_resolve():
    records = [
        record
        for material_identifier in list_material_documents()
        for record in get_material_document(material_identifier)["eos_records"]
    ]

    assert len(records) == 163
    assert set(list_eos_record_documents()) == {
        record["identifier"] for record in records
    }
    assert all("pressure_calibration" in record for record in records)
    assert {record["pressure_calibration"]["audit_date"] for record in records} == {
        "2026-09-03",
        "2026-09-04",
    }
    manifest = json.loads(
        resources.files("peritheos.data.materials")
        .joinpath("manifest.json")
        .read_text(encoding="utf-8")
    )
    summary = manifest["pressure_calibration"]
    assert summary["status_counts"] == {
        status: sum(
            record["pressure_calibration"]["status"] == status for record in records
        )
        for status in ("resolved", "partially_resolved", "not_applicable", "unresolved")
    }
    assert summary["resolvable_reference_eos_uses"] == sum(
        "reference_eos_record" in method
        for record in records
        for method in record["pressure_calibration"]["methods"]
    )
    assert summary["resolvable_reference_calibration_uses"] == sum(
        "reference_calibration_record" in method
        for record in records
        for method in record["pressure_calibration"]["methods"]
    )
    validate_pressure_calibration_references()

    b4c = get_eos_record_document("b4c_somayazulu_2023_bm3_1")
    method = b4c["pressure_calibration"]["methods"][0]
    assert method["reference_eos_record"] == "mgo_b1_tange_2009_vinet"
    reference = get_eos_record_document(method["reference_eos_record"])
    assert reference["eos"]["model"] == "vinet"
    assert reference["thermal"]["model"] == ("asymptotic_power_law_mie_gruneisen_debye")


def test_pressure_calibration_round_trips_through_material_api():
    document = get_material_document("forsterite")
    source = document["eos_records"][0]["pressure_calibration"]
    material = Material.from_eosmat(document)

    assert material.eos_records[0].pressure_calibration == source
    assert material.to_eosmat()["eos_records"][0]["pressure_calibration"] == source


def test_every_primary_validated_migrated_record_is_executable():
    failures = []
    checked = 0
    for identifier in list_material_documents():
        document = get_material_document(identifier)
        for record in document["eos_records"]:
            if record["scientific_validation"]["status"] != "primary_source_validated":
                continue
            checked += 1
            try:
                Material.from_eosmat(
                    document, record_identifiers=[record["identifier"]]
                )
            except (TypeError, ValueError) as error:
                failures.append(f"{record['identifier']}: {error}")

    assert checked == 163
    assert failures == []


def test_benedict_diamond_eosmat_record_loads_and_reproduces_150_gpa():
    document = get_material_document("diamond")
    identifier = "diamond_benedict_2014_double_debye_4"
    material = Material.from_eosmat(document, record_identifiers=[identifier])
    record = material.get_eos_record(identifier)
    source = next(
        item for item in document["eos_records"] if item["identifier"] == identifier
    )

    assert source["thermal"]["model"] == "double_debye_helmholtz"
    assert record.reference_volume == pytest.approx(45.6272)
    assert record.pressure(8.0 * 4.654270411587497, 3000.0) == pytest.approx(
        150.0, rel=5.0e-11
    )
    assert record.volume(150.0, 3000.0) == pytest.approx(
        8.0 * 4.654270411587497, rel=1.0e-11
    )


def test_correa_diamond_eosmat_record_loads_and_reproduces_library_regression():
    document = get_material_document("diamond")
    identifier = "diamond_correa_2008_double_debye_log_moment_5"
    material = Material.from_eosmat(document, record_identifiers=[identifier])
    record = material.get_eos_record(identifier)
    source = next(
        item for item in document["eos_records"] if item["identifier"] == identifier
    )

    assert source["thermal"]["model"] == "double_debye_log_moment_helmholtz"
    assert record.reference_volume == pytest.approx(46.28)
    assert record.pressure(8.0 * 4.43, 5000.0) == pytest.approx(
        202.62811519774186, rel=5.0e-11
    )
    assert record.volume(202.62811519774186, 5000.0) == pytest.approx(
        8.0 * 4.43, rel=1.0e-11
    )


@pytest.mark.parametrize(
    ("absolute_identifier", "anchored_identifier"),
    [
        (
            "diamond_benedict_2014_double_debye_4",
            "diamond_benedict_2014_dewaele_anchored",
        ),
        (
            "diamond_correa_2008_double_debye_log_moment_5",
            "diamond_correa_2008_dewaele_anchored",
        ),
    ],
)
def test_diamond_eosmat_reference_temperature_controls_isotherm_anchoring(
    absolute_identifier, anchored_identifier
):
    document = get_material_document("diamond")
    material = Material.from_eosmat(
        document,
        record_identifiers=[
            "diamond_dewaele_2008_vinet_2",
            absolute_identifier,
            anchored_identifier,
        ],
    )
    source_by_identifier = {
        record["identifier"]: record for record in document["eos_records"]
    }
    reference = material.get_eos_record("diamond_dewaele_2008_vinet_2")
    absolute = material.get_eos_record(absolute_identifier)
    anchored = material.get_eos_record(anchored_identifier)
    volume = 40.0
    temperature = 3000.0

    assert (
        source_by_identifier[absolute_identifier]["thermal"]["parameters"]["Tr"] is None
    )
    assert (
        source_by_identifier[anchored_identifier]["thermal"]["parameters"]["Tr"]
        == 298.0
    )
    assert anchored.pressure(volume, 298.0) == pytest.approx(
        reference.pressure(volume, 298.0)
    )
    assert anchored.pressure(volume, temperature) == pytest.approx(
        reference.pressure(volume, 298.0)
        + absolute.pressure(volume, temperature)
        - absolute.pressure(volume, 298.0)
    )


def test_all_sokolova_eosmat_records_explain_fit_and_workbook_lineage():
    sources = [
        record
        for identifier in list_material_documents()
        for record in get_material_document(identifier)["eos_records"]
        if (record.get("thermal") or {}).get("type") == "Sokolova2016"
    ]

    assert len(sources) == 11
    for source in sources:
        assert source["reference"]["year"] == 2013
        assert source["reference"]["doi"] == "10.1016/j.rgg.2013.01.005"
        assert source["scientific_validation"]["primary_source_check"]["doi"] == (
            "10.1016/j.rgg.2013.01.005"
        )
        lineage = {item["role"]: item["doi"] for item in source["source_lineage"]}
        assert lineage["reference volume, composition, and shock inputs"] == (
            "10.1016/j.rgg.2013.01.005"
        )
        assert (
            lineage["final cross-calibrated Holzapfel and thermal coefficients"]
            == "10.1016/j.rgg.2013.01.005"
        )
        assert lineage["spreadsheet implementation, conventions, and corrections"] == (
            "10.1016/j.cageo.2016.06.002"
        )
        assert "not derived from one experimental dataset" in source["notes"]

    mgo = next(source for source in sources if source["identifier"].startswith("mgo_"))
    assert any(
        item["role"] == "implemented MgO anharmonic-coefficient correction"
        for item in mgo["source_lineage"]
    )


def test_primary_audit_records_corrections_and_known_source_limitations():
    graphite = get_material_document("graphite")["eos_records"][0]
    b4c = get_eos_record_document("b4c_somayazulu_2023_bm3_1")
    zircon = get_material_document("zircon")["eos_records"][0]
    shen = next(
        record
        for record in get_material_document("gold")["eos_records"]
        if record["identifier"] == "gold_shen_2026_vinet_3"
    )

    assert graphite["parameter_errors"]["V0"] == pytest.approx(0.02)
    assert graphite["audit_corrections"][0]["primary_reference"]["location"] == (
        "page 12599, Murnaghan-fit paragraph"
    )
    assert b4c["scientific_validation"]["reported_inconsistencies"][0] == {
        "field": "eos.order",
        "abstract": "third-order Birch-Murnaghan",
        "figure_1_caption": "second-order Birch-Murnaghan",
        "resolution": (
            "Retain BM3 because K0'=3.3(1), which is incompatible with a "
            "conventional BM2 fit where K0'=4 is fixed."
        ),
    }
    assert zircon["identifier"] == "zircon_hazen_1979_bm3_1"
    assert zircon["eos"] == {
        "type": "BM3",
        "parameters": {"V0": 260.79, "K0": 227.0, "K0_prime": 6.5},
        "model": "birch_murnaghan_3",
    }
    assert zircon["parameter_errors"]["V0"] == pytest.approx(0.04)
    assert zircon["fixed_parameters"] == ["K0_prime"]
    assert zircon["scientific_validation"]["status"] == "primary_source_validated"
    assert shen["scientific_validation"]["status"] == "primary_source_validated"
    assert shen["scientific_validation"]["primary_source_check"]["locations"] == [
        "Equation (4), page 8",
        "Table I and footnotes, page 4",
        "Table II, pages 9-12",
        "Section III.E, pages 8 and 12-14",
        "Experimental methods, page 2",
        "Supplemental Table S1, official Excel workbook",
    ]

    magnesite = get_material_document("magnesite")["eos_records"][0]
    molybdenum_carbide_document = get_material_document(
        "molybenum_carbide_mo2c"
    )
    molybdenum_carbide_records = {
        record["identifier"]: record
        for record in molybdenum_carbide_document["eos_records"]
    }
    molybdenum_carbide = molybdenum_carbide_records[
        "molybenum_carbide_mo2c_haines_2001_bm3_1"
    ]
    molybdenum_carbide_refit = molybdenum_carbide_records[
        "molybenum_carbide_mo2c_haines_2001_bm3_refit"
    ]
    platinum = get_material_document("platinum")["eos_records"][0]
    alumina = get_material_document("alumina")["eos_records"][0]
    cobalt = get_material_document("cobalt_hcp")["eos_records"][0]
    niobium = get_material_document("niobium")["eos_records"][0]
    tantalum = get_material_document("tantalum")["eos_records"][0]

    assert magnesite["eos"]["parameters"]["V0"] == pytest.approx(279.41)
    assert magnesite["parameter_errors"]["V0"] == pytest.approx(0.08)
    assert molybdenum_carbide["eos"]["parameters"]["V0"] == pytest.approx(148.9071)
    assert molybdenum_carbide["parameter_errors"]["V0"] == pytest.approx(0.049)
    assert molybdenum_carbide["experimental_pressure_range_gpa"] == [0.0, 46.0]
    assert molybdenum_carbide_refit["record_kind"] == "refit"
    assert molybdenum_carbide_refit["derived_from_record"] == (
        "molybenum_carbide_mo2c_haines_2001_bm3_1"
    )
    assert molybdenum_carbide_refit["fixed_parameters"] == ["V0"]
    assert molybdenum_carbide_refit["eos"]["parameters"] == pytest.approx(
        {
            "V0": 148.9071,
            "K0": 325.87419826054304,
            "K0_prime": 4.9092021387052664,
        }
    )
    assert molybdenum_carbide_refit["parameter_errors"]["V0"] is None
    assert molybdenum_carbide_refit["parameter_errors"]["K0"] == pytest.approx(
        9.50029103042738
    )
    assert molybdenum_carbide_refit["parameter_errors"][
        "K0_prime"
    ] == pytest.approx(0.6513930636044009)
    assert molybdenum_carbide_refit["fit_provenance"]["selection"] == {
        "predicate": "all Figure 2 markers",
        "included_rows": 16,
        "excluded_rows": 0,
    }
    assert molybdenum_carbide_refit["fit_provenance"]["objective"] == (
        "errors_in_variables"
    )
    assert molybdenum_carbide_document["datasets"][0]["used_by_eos_records"] == [
        "molybenum_carbide_mo2c_haines_2001_bm3_1",
        "molybenum_carbide_mo2c_haines_2001_bm3_refit",
    ]
    molybdenum_carbide_executable = Material.from_eosmat(
        molybdenum_carbide_document
    ).get_eos_record("molybenum_carbide_mo2c_haines_2001_bm3_refit")
    assert molybdenum_carbide_executable.eosmat_metadata["record_kind"] == "refit"
    assert platinum["identifier"] == "platinum_holmes_1989_vinet_1"
    assert platinum["eos"]["type"] == "Vinet"
    assert platinum["thermal"]["type"] == "LinearThermalPressure"
    assert platinum["experimental_pressure_range_gpa"] == [0.0, 550.0]
    assert alumina["fixed_parameters"] == []
    assert cobalt["eos"]["type"] == "Vinet"
    assert niobium["experimental_pressure_range_gpa"] == [0.0, 71.5]
    assert tantalum["eos"]["parameters"]["V0"] == pytest.approx(36.0835)
    assert tantalum["fixed_parameters"] == ["V0"]


def test_b4c_thermal_record_is_source_reproduced_and_diffraction_ready():
    document = get_material_document("b4c")
    records = {item["identifier"]: item for item in document["eos_records"]}
    record = records["b4c_somayazulu_2023_bm3_1"]
    berman = records["b4c_somayazulu_2023_berman_2"]
    refit = records["b4c_somayazulu_2023_berman_refit"]

    assert document["phase"] == "rhombohedral B4C"
    assert document["space_group"] == "R-3m"
    assert document["space_group_number"] == 166
    assert document["formula_units_per_cell"] == 9
    assert len(document["atom_sites"]) == 5
    assert all(site.get("wyckoff") for site in document["atom_sites"])

    # The asymmetric-site multiplicities and split occupancies must expand to
    # B36C9, i.e. nine B4C formula units in the conventional hexagonal cell.
    atom_counts = {"B": 0.0, "C": 0.0}
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        atom_counts[site["element"]] += multiplicity * site["occupancy"]
    assert atom_counts == pytest.approx({"B": 36.0, "C": 9.0})

    lattice = document["lattice"]
    cell_volume = math.sqrt(3.0) * lattice["a"] ** 2 * lattice["c"] / 2.0
    assert cell_volume == pytest.approx(328.38272251888424)
    assert cell_volume == pytest.approx(record["eos"]["parameters"]["V0"], abs=0.02)

    assert record["fixed_parameters"] == ["V0", "K0", "K0_prime"]
    assert record["experimental_temperature_range_k"] == [300.0, 2500.0]
    assert record["thermal"] == {
        "type": "MieGruneisenDebye",
        "model": "mie_gruneisen_debye",
        "debye_temperature_law": "integrated_gruneisen",
        "parameters": {
            "Tr": 298.0,
            "theta0": 1425.0,
            "gamma0": 0.8,
            "q": 2.1,
            "n": 5,
        },
        "parameter_errors": {
            "Tr": None,
            "theta0": None,
            "gamma0": None,
            "q": None,
            "n": None,
        },
        "fixed_parameters": ["Tr", "theta0", "gamma0", "n"],
    }

    material = Material.from_eosmat(document)
    executable = material.get_eos_record("b4c_somayazulu_2023_bm3_1")
    calculated_pressure = executable.pressure(289.1, 2023.0)
    assert calculated_pressure == pytest.approx(40.45288082692554)
    assert abs(calculated_pressure - 40.4) < 2.3

    assert berman["identifier"] == "b4c_somayazulu_2023_berman_2"
    assert berman["fixed_parameters"] == ["V0", "K0", "K0_prime"]
    assert berman["thermal"] == {
        "type": "AlphaKT",
        "model": "thermal_reference_state",
        "thermal_expansion_law": "linear_temperature",
        "reference_volume_law": "berman",
        "parameters": {
            "Tr": 298.0,
            "alpha0": 1.94e-5,
            "alpha1": 5.73e-10,
            "dK_dT": -0.008,
        },
        "parameter_errors": {
            "Tr": None,
            "alpha0": 1.6e-6,
            "alpha1": None,
            "dK_dT": 0.003,
        },
        "fixed_parameters": ["Tr", "alpha1"],
    }
    berman_executable = material.get_eos_record("b4c_somayazulu_2023_berman_2")
    assert berman_executable.pressure(289.1, 2023.0) == pytest.approx(43.38255199961772)

    assert refit["record_kind"] == "refit"
    assert refit["derived_from_record"] == "b4c_somayazulu_2023_berman_2"
    assert refit["fit_provenance"]["software"] == {
        "name": "EosFit7c",
        "version": "7.60",
        "version_date": "2021-05-13",
    }
    assert refit["fit_provenance"]["selection"]["included_rows"] == 41
    assert refit["thermal"]["parameters"] == {
        "Tr": 298.0,
        "alpha0": pytest.approx(1.81120e-5),
        "alpha1": pytest.approx(5.73e-10),
        "dK_dT": pytest.approx(-0.01311),
    }
    refit_executable = material.get_eos_record("b4c_somayazulu_2023_berman_refit")
    assert refit_executable.pressure(289.1, 2023.0) == pytest.approx(40.84280383878691)
    assert refit_executable.eosmat_metadata["record_kind"] == "refit"


def test_former_peak_only_materials_have_source_backed_crystal_models():
    expected_cell_counts = {
        "fe3s": {"Fe": 24.0, "S": 8.0},
        "feh2": {"Fe": 4.0, "H": 8.0},
        "iceviii": {"D": 16.0, "O": 8.0},
        "majorite": {"Mg": 32.0, "Si": 32.0, "O": 96.0},
        "mg7si2o8": {"Mg": 14.0, "Si": 4.0, "O": 28.0, "H": 12.0},
        "mgfe60o": {"Mg": 1.6, "Fe": 2.4, "O": 4.0},
        "perovskite_orthorhombic": {
            "Mg": 3.52,
            "Fe": 0.48,
            "Si": 4.0,
            "O": 12.0,
        },
        "phase_d": {"Mg": 1.11, "Si": 1.89, "O": 6.0, "H": 2.22},
        "sno2_cubic_27gpa": {"Sn": 4.0, "O": 8.0},
        "sno2_pa_3_at_48gpa": {"Sn": 4.0, "O": 8.0},
    }

    for material_id, expected in expected_cell_counts.items():
        document = get_material_document(material_id)
        assert document["space_group_number"]
        assert document["peaks"]
        assert document["atom_sites"]
        assert all(site.get("wyckoff") for site in document["atom_sites"])

        actual: dict[str, float] = {}
        for site in document["atom_sites"]:
            multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
            actual[site["element"]] = actual.get(site["element"], 0.0) + (
                multiplicity * site["occupancy"]
            )
        assert actual == pytest.approx(expected)

    fe3s = get_material_document("fe3s")
    assert fe3s["space_group"] == "I-4"
    assert fe3s["space_group_number"] == 82
    assert fe3s["formula_units_per_cell"] == 8
    assert "approximate Fe3P-type proxy" in fe3s["notes"]
    assert "not be presented as its original structure" in fe3s["notes"]


def test_b4c_primary_data_transcription_is_complete():
    document = get_material_document("b4c")
    assert len(document["datasets"]) == 1
    dataset = document["datasets"][0]
    assert dataset["identifier"] == "b4c_somayazulu_2023_table4_pvt"
    assert dataset["used_by_eos_records"] == [
        "b4c_somayazulu_2023_bm3_1",
        "b4c_somayazulu_2023_berman_2",
        "b4c_somayazulu_2023_berman_refit",
    ]
    names = [column["name"] for column in dataset["columns"]]
    data = np.array(
        [tuple(row) for row in dataset["rows"]],
        dtype=[(name, float) for name in names],
    )

    assert len(data) == 51
    assert np.count_nonzero(data["temperature_k"] > 300.0) == 41
    assert np.count_nonzero(data["temperature_k"] == 300.0) == 10
    representative = data[
        (data["volume_a3"] == 289.1) & (data["temperature_k"] == 2023.0)
    ]
    assert len(representative) == 1
    assert representative["pressure_gpa"][0] == pytest.approx(40.4)
    assert representative["pressure_sigma_gpa"][0] == pytest.approx(2.3)
    material = Material.from_eosmat(document)
    assert len(material.datasets) == 1
    assert material.to_eosmat()["datasets"] == document["datasets"]


def test_eosmat_dataset_validation_rejects_ragged_and_dangling_data():
    document = get_material_document("b4c")
    document["datasets"][0]["rows"][0].pop()
    with pytest.raises(ValueError, match="must match the column count"):
        validate_eosmat_document(document)

    document = get_material_document("b4c")
    document["datasets"][0]["used_by_eos_records"] = ["missing_record"]
    with pytest.raises(ValueError, match="references unknown EOS record"):
        validate_eosmat_document(document)


def test_eosmat_refit_validation_requires_resolvable_lineage_and_dataset():
    document = get_material_document("b4c")
    refit = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == "b4c_somayazulu_2023_berman_refit"
    )
    refit["derived_from_record"] = "missing_record"
    with pytest.raises(ValueError, match="derives from unknown EOS record"):
        validate_eosmat_document(document)

    document = get_material_document("b4c")
    refit = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == "b4c_somayazulu_2023_berman_refit"
    )
    refit["fit_provenance"]["dataset"] = "missing_dataset"
    with pytest.raises(ValueError, match="references unknown dataset"):
        validate_eosmat_document(document)

    document = get_material_document("b4c")
    refit = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == "b4c_somayazulu_2023_berman_refit"
    )
    del refit["fit_provenance"]
    with pytest.raises(ValueError, match="refit records require"):
        validate_eosmat_document(document)


def test_newly_validated_primary_records_retain_published_errors():
    ca_records = {
        record["identifier"]: record
        for record in get_material_document("ca_perovskite")["eos_records"]
    }
    shim = ca_records["ca_perovskite_shim_2000_bm3_1"]
    mao = ca_records["ca_perovskite_mao_1989_bm3_2"]
    cao_b1 = get_material_document("cao")["eos_records"][0]
    cao_b2 = get_material_document("cao_b2")["eos_records"][0]
    geo2 = get_material_document("geo2_rutile")["eos_records"][0]
    sno2 = get_material_document("sno2")["eos_records"][0]
    pbs = get_material_document("pbs_b1")["eos_records"][0]
    wadsleyite = get_material_document("wadsleyite")["eos_records"][0]
    jadeite = get_material_document("naalsi2o6")["eos_records"][0]
    perovskite = get_material_document("perovskite_orthorhombic")["eos_records"][0]
    kcl = next(
        record
        for record in get_material_document("kcl")["eos_records"]
        if record["identifier"] == "kcl_walker_2002_bm3_2"
    )

    # Shim et al. fixed V0 in the Table 2 fit but explicitly report its
    # uncertainty in the abstract, so it must not be discarded.
    assert shim["parameter_errors"] == {
        "V0": pytest.approx(0.05),
        "K0": pytest.approx(4.0),
        "K0_prime": pytest.approx(0.2),
    }
    assert mao["parameter_errors"] == {
        "V0": pytest.approx(0.08),
        "K0": pytest.approx(4.0),
        "K0_prime": None,
    }
    assert cao_b1["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(1.0),
        "K0_prime": pytest.approx(0.2),
    }
    assert cao_b2["parameter_errors"] == {
        "V0": pytest.approx(0.3321),
        "K0": pytest.approx(20.0),
        "K0_prime": pytest.approx(0.5),
    }
    assert geo2["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(5.0),
        "K0_prime": None,
    }
    assert sno2["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(2.0),
        "K0_prime": None,
    }
    assert pbs["parameter_errors"] == {
        "V0": pytest.approx(0.0421),
        "K0": pytest.approx(1.2),
        "K0_prime": pytest.approx(0.9),
    }
    assert wadsleyite["parameter_errors"] == {
        "V0": pytest.approx(0.02),
        "K0": pytest.approx(0.9),
        "K0_prime": pytest.approx(0.1),
    }
    assert wadsleyite["thermal"]["parameter_errors"] == {
        "Tr": None,
        "theta0": None,
        "gamma0": pytest.approx(0.02),
        "q": pytest.approx(0.1),
        "n": None,
    }
    assert jadeite["parameter_errors"] == {
        "V0": pytest.approx(0.08),
        "K0": pytest.approx(4.0),
        "K0_prime": None,
    }
    assert perovskite["parameter_errors"] == {
        "V0": pytest.approx(0.39),
        "K0": pytest.approx(6.0),
        "K0_prime": pytest.approx(0.4),
    }
    assert kcl["thermal"]["parameter_errors"] == {
        "Tr": None,
        "alpha_KT": pytest.approx(0.00009),
    }


@pytest.mark.parametrize(
    (
        "material_identifier",
        "record_identifier",
        "reference_volume",
        "expected_pressure",
        "bulk_modulus_error",
        "derivative_error",
    ),
    SHEN_SMITH_2026_RECORDS,
)
def test_shen_smith_2026_table_ii_vinet_regression_and_errors(
    material_identifier,
    record_identifier,
    reference_volume,
    expected_pressure,
    bulk_modulus_error,
    derivative_error,
):
    document = get_material_document(material_identifier)
    source_record = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == record_identifier
    )
    loaded_record = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).eos_records[0]

    # Equation (4) and Table II are independently reproduced at X=0.9,
    # where X=(V/V0)^(1/3). Some phases are outside their fitted interval at
    # this common compression, so validity checking is disabled for this
    # equation/parameter regression only.
    assert loaded_record.pressure(
        reference_volume * 0.9**3, 300.0, check_validity=False
    ) == pytest.approx(expected_pressure)
    assert source_record["eos"]["parameters"]["V0"] == pytest.approx(reference_volume)
    assert source_record["fixed_parameters"] == ["V0"]
    assert source_record["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(bulk_modulus_error),
        "K0_prime": pytest.approx(derivative_error),
    }
    assert "parameter_error_confidence" not in source_record

    # The article gives no confidence level or covariance. Peritheos therefore
    # treats the printed +/- values as standard errors and records the
    # independent-parameter assumption rather than fabricating covariance.
    uncertainty = loaded_record._uncertainty().parameter_uncertainty
    assert uncertainty.standard_errors["K0"] == pytest.approx(bulk_modulus_error)
    assert uncertainty.standard_errors["K0_prime"] == pytest.approx(derivative_error)
    prediction = loaded_record.pressure_with_uncertainty(
        reference_volume * 0.95**3,
        300.0,
        volume_sigma=0.01,
        check_validity=False,
    )
    assert prediction.standard_error > 0.0
    assert "parameter errors treated as mutually independent" in prediction.assumptions
    with pytest.raises(ValueError, match="isothermal EOS record"):
        loaded_record.pressure_with_uncertainty(
            reference_volume * 0.95**3,
            300.0,
            temperature_sigma=0.5,
            check_validity=False,
        )


@pytest.mark.parametrize(
    ("material_identifier", "record_identifier"),
    NINETY_FIVE_PERCENT_INTERVAL_RECORDS,
)
def test_published_95_percent_intervals_are_normalized_when_loaded(
    material_identifier, record_identifier
):
    document = get_material_document(material_identifier)
    source_record = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == record_identifier
    )
    loaded_record = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).eos_records[0]

    assert source_record["parameter_error_confidence"] == 0.95
    assert loaded_record.parameter_error_confidence == 0.95

    standard_errors = loaded_record._uncertainty().parameter_uncertainty.standard_errors
    for name, interval_half_width in source_record["parameter_errors"].items():
        if interval_half_width is None:
            continue
        loaded_name = name if name in standard_errors else f"rt_eos.{name}"
        assert standard_errors[loaded_name] == pytest.approx(
            interval_half_width / 1.95996398454
        )


def test_campbell_kcl_composite_has_explicit_primary_inputs_and_errors():
    document = get_material_document("kcl")
    identifier = "kcl_campbell_1991_bm2_1"
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]

    assert source["scientific_validation"]["status"] == "primary_source_validated"
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(0.8483 * 62.36),
        "K0": pytest.approx(28.7),
    }
    assert source["parameter_errors"] == {
        "V0": pytest.approx(0.0057 * 62.36),
        "K0": pytest.approx(0.6),
    }
    assert source["fixed_parameters"] == ["V0"]
    assert "Dewaele et al. (2012)" in source["parameter_provenance"]["V0"]
    assert record.pressure(
        record.reference_volume, 298.0, check_validity=False
    ) == pytest.approx(0.0, abs=1e-14)
    assert record._uncertainty().parameter_uncertainty.standard_errors == {
        "V0": pytest.approx(0.0057 * 62.36),
        "K0": pytest.approx(0.6),
    }


def test_munoz_1993_inn_uses_theoretical_murnaghan_cell():
    document = get_material_document("indium_nitride")
    source = document["eos_records"][0]
    record = Material.from_eosmat(document).eos_records[0]
    expected_volume = math.sqrt(3.0) / 2.0 * 3.483**2 * 5.7039

    assert source["identifier"] == "indium_nitride_munoz_1993_murnaghan_1"
    assert source["scientific_validation"]["status"] == "primary_source_validated"
    assert source["eos"] == {
        "type": "Murnaghan",
        "model": "murnaghan",
        "parameters": {
            "V0": pytest.approx(expected_volume),
            "K0": pytest.approx(166.0),
            "K0_prime": pytest.approx(3.8),
        },
    }
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": None,
        "K0_prime": None,
    }
    assert source["pressure_range_status"] == "theoretical"
    assert record.pressure(expected_volume) == pytest.approx(0.0, abs=1e-14)
    pressure = record.pressure(expected_volume * 0.9)
    assert record.volume(pressure) == pytest.approx(expected_volume * 0.9)


def test_formerly_unverified_migrated_reductions_are_corrected_or_removed():
    record_ids = {
        record["identifier"]
        for material_id in list_material_documents()
        for record in get_material_document(material_id)["eos_records"]
    }
    lithium = get_material_document("li_bcc")["eos_records"][0]
    majorite = get_material_document("majorite")
    phase_d = get_material_document("phase_d")["eos_records"]

    assert lithium["identifier"] == "li_hanfland_1999_vinet_1"
    assert lithium["eos"]["type"] == "Vinet"
    assert lithium["scientific_validation"]["status"] == "primary_source_validated"
    assert lithium["experimental_pressure_range_gpa"] == [0.0, 21.1]
    assert "feo_fei_1995_bm3_1" not in record_ids
    assert "tungsten_hixson_1992_bm3_1" not in record_ids
    assert "mgsio3" not in list_material_documents()
    assert majorite["formula_units_per_cell"] == 32
    assert majorite["space_group_number"] == 88
    assert majorite["eos_records"][0]["identifier"] == "majorite_yagi_1992_bm3_1"
    assert majorite["eos_records"][0]["eos"]["parameters"]["K0"] == 161.2
    assert [record["identifier"] for record in phase_d] == [
        "phase_d_ant_a_shieh_2000_bm2_1",
        "phase_d_ant_b_shieh_2000_bm2_1",
    ]
    assert [record["eos"]["parameters"]["V0"] for record in phase_d] == [
        88.12,
        87.191,
    ]
    assert all(
        record["scientific_validation"]["status"] == "primary_source_validated"
        for record in phase_d
    )


@pytest.mark.parametrize(
    ("material_id", "record_id", "volume", "published_pressure", "absolute_error"),
    [
        ("cscl", "cscl_campbell_1994_bm3_1", 3.5308**3, 28.7, 1.5),
        ("rbcl", "rbcl_b2_campbell_1994_bm3_1", 3.3301**3, 32.3, 2.0),
        (
            "fe3o4",
            "fe3o4_mao_1974_bm3_1",
            591.434826984 * 0.923,
            20.0,
            3.0,
        ),
        ("li_bcc", "li_hanfland_1999_vinet_1", 2.0 * 11.4371, 21.1, 0.1),
        ("majorite", "majorite_yagi_1992_bm3_1", 1435.7, 9.72, 0.4),
        (
            "mgfe60o",
            "mgfe60o_richet_1989_bm2_1",
            3.986**3,
            49.4,
            1.8,
        ),
        (
            "nis",
            "nis_campbell_1993_bm3_1",
            math.sqrt(3.0) / 2.0 * 3.2520**2 * 4.871,
            44.9,
            6.0,
        ),
        ("phase_d", "phase_d_ant_a_shieh_2000_bm2_1", 77.62, 24.6, 3.0),
        ("phase_d", "phase_d_ant_b_shieh_2000_bm2_1", 77.90, 22.0, 3.2),
        (
            "sno2_cubic_27gpa",
            "sno2_cubic_27gpa_ono_2000_bm3_1",
            119.85,
            25.64,
            0.6,
        ),
        ("sro", "sro_liu_1973_bm3_1", 137.388096 * 0.797, 34.05, 0.5),
        (
            "sro_b2",
            "sro_b2_sato_1981_bm2_1",
            28.0224 * 6.14 / 7.767,
            59.0,
            4.0,
        ),
    ],
)
def test_newly_validated_primary_table_regressions(
    material_id, record_id, volume, published_pressure, absolute_error
):
    document = get_material_document(material_id)
    record = Material.from_eosmat(document, record_identifiers=[record_id]).eos_records[
        0
    ]

    assert record.pressure(volume, check_validity=False) == pytest.approx(
        published_pressure, abs=absolute_error
    )
    calculated_pressure = record.pressure(volume, check_validity=False)
    assert record.volume(calculated_pressure, check_validity=False) == pytest.approx(
        volume
    )


def test_newly_validated_primary_errors_and_single_remaining_deferral():
    ono = get_material_document("sno2_cubic_27gpa")["eos_records"][0]
    magnetite = get_material_document("fe3o4")["eos_records"][0]
    mw60 = get_material_document("mgfe60o")["eos_records"][0]
    nis = get_material_document("nis")["eos_records"][0]
    sro_b2 = get_material_document("sro_b2")["eos_records"][0]
    deferred = [
        record["identifier"]
        for material_id in list_material_documents()
        for record in get_material_document(material_id)["eos_records"]
        if record["scientific_validation"]["status"] == "deferred"
    ]

    assert get_material_document("cscl")["phase"] == "B2 (CsCl-type), cubic Pm-3m"
    assert get_material_document("fe3o4")["phase"] == (
        "magnetite, cubic inverse-spinel Fd-3m"
    )
    assert get_material_document("nis")["phase"] == (
        "metastable NiAs-type, hexagonal P63/mmc"
    )
    assert get_material_document("sro")["phase"] == ("B1 (NaCl-type), cubic Fm-3m")
    assert get_material_document("sro_b2")["phase"] == ("B2 (CsCl-type), cubic Pm-3m")

    assert ono["parameter_errors"]["V0"] == pytest.approx(0.3)
    assert magnetite["parameter_errors"]["V0"] == pytest.approx(0.634133124)
    assert mw60["parameter_errors"]["V0"] == pytest.approx(0.0219872163)
    assert nis["parameter_errors"]["V0"] == pytest.approx(0.009596193736601715)
    assert sro_b2["parameter_errors"] == {
        "V0": pytest.approx(0.9127817589576549),
        "K0": pytest.approx(19.0),
    }
    assert deferred == []

    executable = Material.from_eosmat(
        get_material_document("sno2_cubic_27gpa")
    ).eos_records[0]
    prediction = executable.pressure_with_uncertainty(
        119.85,
        volume_sigma=0.02,
        check_validity=False,
    )
    assert prediction.standard_error > 0.0


def test_walker_2002_kcl_be1_and_reported_product_error_regression():
    document = get_material_document("kcl")
    source = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == "kcl_walker_2002_bm3_2"
    )
    record = Material.from_eosmat(
        document, record_identifiers=["kcl_walker_2002_bm3_2"]
    ).eos_records[0]

    # The preferred bold Table 3 row uses the staged fit; the italic row is the
    # separate simultaneous solution. V0 is therefore fixed for this record.
    assert source["fixed_parameters"] == ["V0"]
    assert "italic Table 3 row" in source["notes"]

    # Table 2 reports V=47.57(3) A^3 at 36.3(2) kbar and 23 degC.
    assert record.pressure(47.57, 296.15) == pytest.approx(3.63, abs=0.02)

    # Equation BE1 adds alpha0*K0*delta-T. The paper reports the identifiable
    # product as 0.0275(9) kbar/K = 0.00275(9) GPa/K.
    pressure_ambient = record.pressure(47.57, 296.15)
    pressure_600_c = record.pressure(47.57, 873.15)
    assert pressure_600_c - pressure_ambient == pytest.approx(0.00275 * 577.0)
    prediction = record.pressure_with_uncertainty(47.57, 873.15)
    assert prediction.standard_error == pytest.approx(0.00009 * 577.0)


@pytest.mark.parametrize(
    ("volume", "table_pressure_gpa"),
    [(48.145, 3.14), (29.827, 53.8), (22.065, 165.0)],
)
def test_dewaele_2012_kcl_b2_table_i_and_equation_2_regression(
    volume, table_pressure_gpa
):
    document = get_material_document("kcl")
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == "kcl_b2_dewaele_2012_vinet_3"
    )
    record = Material.from_eosmat(
        document, record_identifiers=[source["identifier"]]
    ).eos_records[0]

    assert source["default"] is True
    assert source["fixed_parameters"] == ["V0"]
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": None,
        "K0_prime": None,
    }
    # Table I observations have pressure uncertainties below 2%; the fit is
    # not expected to pass through every measured point exactly.
    assert record.pressure(volume, 300.0) == pytest.approx(
        table_pressure_gpa, rel=0.02, abs=0.1
    )

    # Equation (2) and Table V specify the exact additive thermal term.
    assert record.pressure(volume, 2000.0) - record.pressure(
        volume, 300.0
    ) == pytest.approx(0.00224 * 1700.0)


def test_dewaele_2012_kcl_b2_round_trip_validity_arrays_and_measurement_errors():
    record = Material.from_eosmat(
        get_material_document("kcl"),
        record_identifiers=["kcl_b2_dewaele_2012_vinet_3"],
    ).eos_records[0]

    volume = record.volume(100.0, 2000.0)
    assert record.pressure(volume, 2000.0) == pytest.approx(100.0)
    assert record.pressure([volume, volume], [2000.0, 2100.0]).shape == (2,)
    assert record.within_validity(volume, 2000.0)
    assert np.isfinite(record.pressure(volume, 7200.0))
    with pytest.raises(ValueError, match="outside the published calibration"):
        record.pressure(volume, 7200.0, check_validity=True)

    prediction = record.pressure_with_uncertainty(
        40.0,
        2000.0,
        volume_sigma=0.01,
        temperature_sigma=100.0,
        check_validity=False,
    )
    assert prediction.standard_error > 0.224
    assert "published parameter uncertainty not available" in prediction.assumptions


@pytest.mark.parametrize(
    ("lattice_ratio", "table_pressure_gpa"),
    [
        (0.995, 3.0),
        (0.980, 12.8),
        (0.970, 20.0),
        (0.960, 30.8),
    ],
)
def test_clendenen_1966_coo_table_iii_murnaghan_regression(
    lattice_ratio, table_pressure_gpa
):
    document = get_material_document("coo")
    source = document["eos_records"][0]
    record = Material.from_eosmat(document).eos_records[0]

    assert source["identifier"] == "coo_clendenen_1966_murnaghan_1"
    assert source["eos"] == {
        "type": "Murnaghan",
        "parameters": {
            "V0": pytest.approx(4.258**3),
            "K0": pytest.approx(190.5),
            "K0_prime": pytest.approx(3.9),
        },
        "model": "murnaghan",
    }
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": None,
        "K0_prime": None,
    }
    assert source["fixed_parameters"] == ["V0"]

    # Table III is a smoothed experimental table, not values generated exactly
    # by Equation (4). The paper says the residual is of experimental-error size.
    volume = source["eos"]["parameters"]["V0"] * lattice_ratio**3
    assert record.pressure(volume) == pytest.approx(table_pressure_gpa, abs=1.0)


def test_shim_2000_casio3_table_i_pressure_regression():
    document = get_material_document("ca_perovskite")
    record = Material.from_eosmat(
        document, record_identifiers=["ca_perovskite_shim_2000_bm3_1"]
    ).eos_records[0]

    # Table 1 reports V=35.909(55) A^3 at 89.6(3.0) GPa. The rounded values
    # happen to lie almost exactly on the published BM3 fit.
    assert record.pressure(35.909) == pytest.approx(89.6, abs=0.01)


def test_holmes_1989_platinum_equations_11_and_12_regression():
    document = get_material_document("platinum")
    identifier = "platinum_holmes_1989_vinet_1"
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]
    x = 0.9
    volume = 60.4000884 * x**3
    pressure_300 = (
        3.0 * 266.0 * (1.0 - x) / x**2 * math.exp(1.5 * (5.81 - 1.0) * (1.0 - x))
    )

    assert record.pressure(volume, 300.0) == pytest.approx(pressure_300)
    assert record.pressure(volume, 1000.0) == pytest.approx(
        pressure_300 + 0.0069426 * 700.0
    )


def test_luo_2023_mgo_complete_thermal_eos_and_primary_tables():
    document = get_material_document("mgo")
    identifier = "mgo_b1_luo_2023_vinet_thermal_5"
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]

    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(74.0741025123),
        "K0": pytest.approx(169.8),
        "K0_prime": pytest.approx(4.501),
    }
    assert source["thermal"]["model"] == "second_order_taylor_thermal_pressure"
    assert source["thermal"]["parameters"]["c0"] == pytest.approx(0.5096)

    # Tables II-III use the ambient initial volume from rho0=3.590 g/cm3,
    # whereas Equation (B2) uses the fitted zero-K V0K in eta.
    ambient_cell_volume = 74.5697677586
    table_volume = ambient_cell_volume * (1.0 - 0.42)
    assert record.pressure(table_volume, 8500.0) == pytest.approx(342.01, abs=1.5)
    pressure = record.pressure(table_volume, 8500.0)
    assert record.volume(pressure, 8500.0) == pytest.approx(table_volume)
    assert record.eos.temperature(pressure, table_volume) == pytest.approx(8500.0)

    shock = next(
        dataset
        for dataset in document["datasets"]
        if dataset["identifier"] == "mgo_luo_2023_table1_shock"
    )
    assert len(shock["rows"]) == 5
    for row in shock["rows"]:
        calculated_pressure = 3.590 * row[3] * row[5]
        assert calculated_pressure == pytest.approx(row[9], abs=0.6)

    grid = next(
        dataset
        for dataset in document["datasets"]
        if dataset["identifier"] == "mgo_luo_2023_tables2_3_pvt_grid"
    )
    resource = Path("peritheos/data").joinpath(grid["resource"]["path"])
    rows = list(csv.DictReader(resource.open(encoding="utf-8")))
    assert len(rows) == 576
    assert rows[0] == {
        "compression_ambient": "0.02",
        "temperature_k": "300",
        "pressure_gpa": "2.89",
        "source_table": "II",
    }
    assert rows[-1] == {
        "compression_ambient": "0.42",
        "temperature_k": "8500",
        "pressure_gpa": "342.01",
        "source_table": "III",
    }

    compression = np.array([float(row["compression_ambient"]) for row in rows])
    temperatures = np.array([float(row["temperature_k"]) for row in rows])
    published_pressures = np.array([float(row["pressure_gpa"]) for row in rows])
    volumes = ambient_cell_volume * (1.0 - compression)
    residuals = record.pressure(volumes, temperatures) - published_pressures
    assert np.sqrt(np.mean(residuals**2)) == pytest.approx(0.4840957067)
    assert np.max(np.abs(residuals)) == pytest.approx(1.4353026419)

    # Diagnostic only: Tables II-III are model output, not independent
    # observations, so these coefficients are not stored as another EOS.
    delta_eta = 1.0 - volumes / source["eos"]["parameters"]["V0"] - 0.02
    delta_temperature = temperatures - 300.0
    design = np.column_stack(
        [
            np.ones_like(delta_eta),
            delta_eta,
            delta_temperature,
            0.5 * delta_eta**2,
            0.5 * delta_temperature**2,
            0.5 * delta_eta * delta_temperature,
        ]
    )
    target_thermal_pressure = published_pressures - record.eos.rt_eos.pressure(volumes)
    coefficients, *_ = np.linalg.lstsq(design, target_thermal_pressure, rcond=None)
    assert coefficients == pytest.approx(
        [
            1.42397479,
            -17.5645093,
            0.00609260423,
            45.6965382,
            8.37791429e-8,
            0.00408491607,
        ]
    )
    diagnostic_residuals = design @ coefficients - target_thermal_pressure
    assert np.sqrt(np.mean(diagnostic_residuals**2)) == pytest.approx(0.4457892496)


@pytest.mark.parametrize(
    ("compression", "temperature", "table_pressure"),
    [
        (0.00, 3000.0, 19.28),
        (0.02, 500.0, 4.94),
        (0.10, 1000.0, 27.59),
        (0.20, 3000.0, 80.88),
        (0.34, 3000.0, 222.44),
    ],
)
def test_anderson_1989_gold_equation_29_table_v_regression(
    compression, temperature, table_pressure
):
    document = get_material_document("gold")
    identifier = "gold_anderson_1989_bm3_1"
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]
    volume = record.reference_volume * (1.0 - compression)

    # Top rows in Table V are Anderson et al.'s Equation (29), rounded to
    # 0.01 GPa. Bottom rows are the separate Heinz-Jeanloz comparison.
    assert record.pressure(volume, temperature) == pytest.approx(
        table_pressure, abs=0.005
    )


def test_anderson_1989_gold_partial_published_error_is_propagated():
    document = get_material_document("gold")
    identifier = "gold_anderson_1989_bm3_1"
    source_record = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]
    volume = record.reference_volume * 0.8
    prediction = record.pressure_with_uncertainty(volume, 3000.0)

    assert source_record["thermal"]["model"] == "log_volume_thermal_pressure"
    assert source_record["thermal"]["parameter_errors"] == {
        "Tr": None,
        "alpha_KT_ref": None,
        "dK_dT_V": pytest.approx(0.001),
    }
    assert prediction.standard_error == pytest.approx(
        math.log(1.0 / 0.8) * 2700.0 * 0.001
    )
    assert "parameter errors treated as mutually independent" in prediction.assumptions


def test_martinez_1996_aragonite_staged_bm2_thermal_metadata():
    document = get_material_document("aragonite")
    assert len(document["eos_records"]) == 1
    bm2_source = document["eos_records"][0]
    material = Material.from_eosmat(
        document, record_identifiers=["aragonite_martinez_1996_bm2_2"]
    )
    bm2 = material.eos_records[0]

    assert bm2_source["scientific_validation"]["status"] == ("primary_source_validated")
    assert bm2.reference_volume == pytest.approx(227.5)
    assert bm2.pressure(227.5) == pytest.approx(0.0, abs=1e-14)
    assert bm2.pressure(206.81, check_validity=False) == pytest.approx(7.14, abs=0.4)
    assert bm2.volume(bm2.pressure(213.08)) == pytest.approx(213.08)
    assert bm2._uncertainty().parameter_uncertainty.standard_errors == {
        "rt_eos.V0": pytest.approx(0.8),
        "rt_eos.K0": pytest.approx(3.48),
        "alpha0": pytest.approx(0.1e-5),
        "dK_dT": pytest.approx(0.002),
    }

    assert bm2_source["thermal"]["thermal_expansion_law"] == "constant"
    assert bm2_source["thermal"]["reference_volume_law"] == "linear_temperature"
    assert bm2_source["thermal"]["parameter_errors"] == {
        "Tr": None,
        "alpha0": pytest.approx(0.1e-5),
        "dK_dT": pytest.approx(0.002),
    }

    exported_thermal = material.to_eosmat()["eos_records"][0]["thermal"]
    assert exported_thermal["configuration"] == {
        "thermal_expansion_law": "constant",
        "reference_volume_law": "linear_temperature",
    }
    temperature = 873.0
    reference_volume = 227.5 * (1.0 + 6.5e-5 * (temperature - 298.0))
    assert bm2.pressure(
        reference_volume, temperature, check_validity=False
    ) == pytest.approx(0.0, abs=1e-13)
    prediction = bm2.pressure_with_uncertainty(
        216.88, temperature, check_validity=False
    )
    assert math.isfinite(prediction.standard_error)
    assert prediction.standard_error > 0.0

    # Table 6 printed isotherm fits independently reproduce the Table 7
    # temperature trends. The K0 slope depends mildly on use of its printed
    # fit errors, and the published -0.018(2) GPa/K lies between both results.
    table_temperatures = np.array([298, 373, 473, 573, 673, 773, 873, 973])
    table_volumes = np.array([227.5, 228.0, 229.6, 231.2, 232.9, 234.3, 235.6, 237.1])
    table_moduli = np.array([64.81, 68.09, 62.8, 60.7, 55.45, 54.8, 55.7, 54.8])
    modulus_errors = np.array([3.48, 4.6, 3.8, 4.7, 3.16, 4.09, 1.8, 2.51])
    delta_temperature = table_temperatures - 298.0
    volume_slope, volume_intercept = np.polyfit(delta_temperature, table_volumes, 1)
    unweighted_modulus_slope = np.polyfit(delta_temperature, table_moduli, 1)[0]
    weighted_modulus_slope = np.polyfit(
        delta_temperature, table_moduli, 1, w=1.0 / modulus_errors
    )[0]

    assert volume_slope / volume_intercept == pytest.approx(6.4835049e-5)
    assert unweighted_modulus_slope == pytest.approx(-0.0196859081)
    assert weighted_modulus_slope == pytest.approx(-0.0170213524)
    assert weighted_modulus_slope > -0.018 > unweighted_modulus_slope


def test_scott_2001_cementite_primary_bm3_parameters_and_errors():
    document = get_material_document("cementite")
    source = document["eos_records"][0]
    record = Material.from_eosmat(document).eos_records[0]

    assert source["scientific_validation"]["status"] == "primary_source_validated"
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(155.26),
        "K0": pytest.approx(175.4),
        "K0_prime": pytest.approx(5.1),
    }
    assert source["parameter_errors"] == {
        "V0": pytest.approx(0.14),
        "K0": pytest.approx(3.5),
        "K0_prime": pytest.approx(0.3),
    }
    assert source["parameter_error_confidence"] is None
    assert source["fixed_parameters"] == ["V0"]
    assert source["experimental_pressure_range_gpa"] == [0.0, 73.2]
    assert source["experimental_temperature_range_k"] == [300.0, 300.0]
    assert record.pressure(155.26, 300.0) == pytest.approx(0.0, abs=1e-14)
    pressure = record.pressure(130.0, 300.0)
    assert record.volume(pressure, 300.0) == pytest.approx(130.0)
    assert record._uncertainty().parameter_uncertainty.standard_errors == {
        "V0": pytest.approx(0.14),
        "K0": pytest.approx(3.5),
        "K0_prime": pytest.approx(0.3),
    }


def test_noguchi_1999_nio_primary_bm3_reference_volume_and_errors():
    document = get_material_document("nickel_oxide")
    source = document["eos_records"][0]
    record = Material.from_eosmat(document).eos_records[0]

    expected_volume = 0.75 * 4.177**3
    expected_error = 0.75 * 3.0 * 4.177**2 * 0.001

    assert source["scientific_validation"]["status"] == "primary_source_validated"
    assert source["eos"]["model"] == "birch_murnaghan_3"
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(expected_volume),
        "K0": pytest.approx(191.0),
        "K0_prime": pytest.approx(3.9),
    }
    assert source["parameter_errors"] == {
        "V0": pytest.approx(expected_error),
        "K0": None,
        "K0_prime": None,
    }
    assert source["fixed_parameters"] == ["V0"]
    assert source["experimental_pressure_range_gpa"] == [0.0, 147.6]
    assert source["temperature_ref"] == pytest.approx(300.0)
    assert record.pressure(record.reference_volume, 300.0) == pytest.approx(
        0.0, abs=1e-14
    )


def test_ross_1997_magnesite_table_i_pressure_regression():
    record = Material.from_eosmat(get_material_document("magnesite")).eos_records[0]

    # Table I reports V=272.33(4) A^3 at 3.09 GPa. This is a measured state,
    # so compare with the published fit within the experimental residual.
    assert record.pressure(272.33, check_validity=False) == pytest.approx(
        3.09, abs=0.01
    )


def test_frank_2004_ice_vii_table_i_pressure_regression():
    document = get_material_document("ice_vii")
    identifier = "ice_vii_frank_2004_bm3_2"
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]
    # Table I experiment 10: 8.315 cm^3/mol at 20.65(38) GPa. Convert the
    # molar volume to the two-formula-unit cubic cell used by .eosmat.
    volume = 8.315 * 1.0e24 * 2.0 / 6.02214076e23

    assert record.pressure(volume) == pytest.approx(20.65, abs=0.38)


def test_primary_audit_restores_model_inputs_omitted_during_migration():
    aluminum = get_material_document("aluminum")["eos_records"][1]
    silica = get_material_document("silica_cacl2")["eos_records"][0]
    ice = get_material_document("ice_vi")["eos_records"][0]

    assert aluminum["eos"]["parameters"]["n"] == 1.0
    assert aluminum["eos"]["parameters"]["Z"] == 13.0
    assert aluminum["thermal"]["parameters"]["n"] == 1.0
    assert silica["thermal"]["parameters"]["n"] == 3.0
    assert ice["thermal"]["parameters"]["Tr"] == 300.0
    assert {item["path"] for item in aluminum["audit_corrections"]} >= {
        "eos.parameters.n",
        "eos.parameters.Z",
        "thermal.parameters.n",
    }


@pytest.mark.parametrize(
    ("material", "volume", "reported_pressure", "tolerance"),
    [
        ("coesite", 522.61, 5.19, 0.05),
        ("zircon", 255.59, 4.81, 0.08),
    ],
)
def test_no_doi_primary_source_table_regressions(
    material, volume, reported_pressure, tolerance
):
    record = Material.from_eosmat(get_material_document(material)).eos_records[0]

    # These are measured P-V rows rather than values generated by the fit.
    assert record.pressure(volume, check_validity=False) == pytest.approx(
        reported_pressure, abs=tolerance
    )


def test_no_doi_primary_sources_use_stable_article_or_report_locators():
    for material in ("coesite", "lead_fcc", "zircon"):
        record = get_material_document(material)["eos_records"][0]
        validation = record["scientific_validation"]
        assert validation["status"] == "primary_source_validated"
        assert validation["primary_source_check"]["doi"] is None
        assert validation["primary_source_check"]["access_url"].startswith("https://")


@pytest.mark.parametrize(
    ("material", "volume", "temperature", "reported_pressure", "tolerance"),
    [
        ("ice_vi", 206.233, 340.7, 2.56, 0.02),
        ("ice_vii", 32.406, 300.6, 8.12, 0.06),
    ],
)
def test_bezacier_2014_table_i_pressure_regression(
    material, volume, temperature, reported_pressure, tolerance
):
    document = get_material_document(material)
    record_identifier = next(
        record["identifier"]
        for record in document["eos_records"]
        if "bezacier_2014" in record["identifier"]
    )
    record = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).eos_records[0]

    # Table I contains measured states, so comparison is to the reported fit
    # within its experimental residual rather than exact numeric identity.
    assert record.pressure(volume, temperature) == pytest.approx(
        reported_pressure, abs=tolerance
    )


def test_bezacier_2014_uncertainty_includes_state_and_published_fit_errors():
    record = Material.from_eosmat(get_material_document("ice_vi")).eos_records[0]

    prediction = record.pressure_with_uncertainty(
        206.233,
        340.7,
        volume_sigma=0.004,
        temperature_sigma=0.4,
    )

    assert prediction.value == pytest.approx(2.5511132321)
    assert prediction.standard_error > 0.0
    assert "parameter errors treated as mutually independent" in prediction.assumptions
    assert "state-variable errors treated as independent" in prediction.assumptions


def test_dioptas_010_shape_remains_accepted_without_peritheos_extensions():
    document = get_material_document("gold")
    document.pop("format")
    document.pop("identifier")
    document.pop("datasets", None)
    document["format_version"] = 2
    for record in document["eos_records"]:
        record.pop("identifier")
        record.pop("scientific_validation")

    validate_eosmat_document(document)


def test_peritheos_eos_only_material_does_not_require_crystal_structure():
    document = get_material_document("gold")
    for key in (
        "symmetry",
        "lattice",
        "formula_units_per_cell",
        "space_group",
        "space_group_number",
        "atom_sites",
        "source",
        "peaks",
    ):
        document.pop(key, None)

    validate_eosmat_document(document)


def test_material_document_returns_a_defensive_copy_and_reports_unknown_id():
    first = get_material_document("gold")
    first["name"] = "changed"
    assert get_material_document("gold")["name"] == "Gold"
    with pytest.raises(KeyError, match="Unknown material document"):
        get_material_document("not_a_material")


def test_eosmat_file_round_trip(tmp_path):
    document = get_material_document("mgo")
    path = tmp_path / "mgo.eosmat"

    save_eosmat(path, document)

    assert load_eosmat(path) == document


def test_normative_schema_is_bundled():
    schema = eosmat_schema()

    assert schema["$schema"].endswith("2020-12/schema")
    assert schema["properties"]["format_version"]["const"] == 3
    assert schema["properties"]["eos_records"]["type"] == "array"
    assert schema["properties"]["datasets"]["type"] == "array"
    assert schema["$defs"]["dataset"]["oneOf"] == [
        {"required": ["rows"]},
        {"required": ["resource"]},
    ]
    assert schema["additionalProperties"] is True
    assert (
        schema["$defs"]["eos_record"]["properties"]["source_lineage"]["type"] == "array"
    )
    thermal = schema["$defs"]["thermal"]["properties"]
    assert thermal["debye_temperature_law"]["default"] == "integrated_gruneisen"
    assert thermal["debye_temperature_law"]["enum"] == [
        "integrated_gruneisen",
        "variable_exponent",
    ]
    assert thermal["reference_volume_law"]["enum"] == [
        "integrated_expansivity",
        "linear_temperature",
        "berman",
    ]
    assert len(schema["$defs"]["equation"]["allOf"][0]["oneOf"]) == 11
    assert len(schema["$defs"]["thermal"]["allOf"][0]["oneOf"]) == 13


def test_normative_schema_validates_every_bundled_document():
    schema = eosmat_schema()
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)

    failures = []
    for identifier in list_material_documents():
        for error in validator.iter_errors(get_material_document(identifier)):
            failures.append(
                f"{identifier}:{'/'.join(map(str, error.path))}: {error.message}"
            )
    assert failures == []


def test_normative_schema_validates_every_executable_material_export():
    validator = Draft202012Validator(eosmat_schema())
    failures = []
    for material in list_materials():
        for error in validator.iter_errors(material.to_eosmat()):
            failures.append(
                f"{material.identifier}:{'/'.join(map(str, error.path))}: "
                f"{error.message}"
            )
    assert failures == []


def test_documented_complete_eosmat_example_is_structurally_valid():
    documentation = Path("docs/eosmat-schema.md").read_text(encoding="utf-8")
    example = re.search(
        r"## Complete EOS-only example.*?```json\n(.*?)\n```",
        documentation,
        flags=re.DOTALL,
    )
    assert example is not None

    document = json.loads(example.group(1))
    validate_eosmat_document(document)
    assert document["eos_records"][0]["thermal"]["debye_temperature_law"] == (
        "variable_exponent"
    )


def test_fei_2007_migrations_use_published_variable_exponent_debye_law():
    records = [
        record
        for identifier in list_material_documents()
        for record in get_material_document(identifier)["eos_records"]
        if isinstance(record["reference"], dict)
        and record["reference"].get("doi") == "10.1073/pnas.0609013104"
        and record.get("thermal") is not None
    ]

    assert {record["identifier"] for record in records} == {
        "gold_fei_2007_vinet_2",
        "neon_fcc_fei_2007_vinet_2",
    }
    assert {
        (
            record["thermal"]["type"],
            record["thermal"]["model"],
            record["thermal"]["debye_temperature_law"],
        )
        for record in records
    } == {("MieGruneisenDebye", "mie_gruneisen_debye", "variable_exponent")}
    assert all(
        record["migration_corrections"][0]["primary_reference"]["location"]
        == "Equation 3 and the definition immediately following it"
        for record in records
    )

    # Other MGD records retain the conventional integrated constant-q law.
    wadsleyite = get_material_document("wadsleyite")["eos_records"][0]
    assert wadsleyite["thermal"]["type"] == "MieGruneisenDebye"
    assert wadsleyite["thermal"]["model"] == "mie_gruneisen_debye"
    assert "debye_temperature_law" not in wadsleyite["thermal"]


@pytest.mark.parametrize(
    ("doi", "expected_law", "corrected"),
    [
        (
            "https://doi.org/10.1073/pnas.0609013104",
            "variable_exponent",
            True,
        ),
        (
            "10.0000/unrelated",
            None,
            False,
        ),
    ],
)
def test_importer_corrects_only_primary_identified_fei_thermal_records(
    tmp_path, doi, expected_law, corrected
):
    source = {
        "format_version": 2,
        "name": "Example",
        "formula": "X",
        "eos_records": [
            {
                "label": "Example record",
                "reference": {"authors": ["Author"], "year": 2007, "doi": doi},
                "eos": {"type": "Vinet", "parameters": {"V0": 1.0}},
                "parameter_errors": {},
                "fixed_parameters": [],
                "thermal": {
                    "type": "MieGruneisenDebye",
                    "parameters": {
                        "Tr": 300.0,
                        "theta0": 170.0,
                        "gamma0": 2.97,
                        "q": 0.6,
                        "n": 1,
                    },
                },
            }
        ],
    }
    path = tmp_path / "example.json"
    path.write_text(json.dumps(source), encoding="utf-8")

    record = migrate_document(path)["eos_records"][0]

    assert record["thermal"]["type"] == "MieGruneisenDebye"
    assert record["thermal"]["model"] == "mie_gruneisen_debye"
    assert record["thermal"].get("debye_temperature_law") == expected_law
    assert ("migration_corrections" in record) is corrected


def test_eosmat_debye_temperature_law_defaults_and_validation():
    document = get_material_document("gold")
    fei = next(
        record["thermal"]
        for record in document["eos_records"]
        if record["identifier"] == "gold_fei_2007_vinet_2"
    )
    fei.pop("debye_temperature_law")
    validate_eosmat_document(document)

    fei["debye_temperature_law"] = "unknown"
    with pytest.raises(ValueError, match="debye_temperature_law is invalid"):
        validate_eosmat_document(document)

    fei["debye_temperature_law"] = "variable_exponent"
    fei["type"] = "AlphaKT"
    fei["model"] = "thermal_reference_state"
    with pytest.raises(ValueError, match="requires MieGruneisenDebye"):
        validate_eosmat_document(document)


def test_eosmat_thermal_expansion_law_defaults_and_validation():
    document = get_material_document("aragonite")
    thermal = document["eos_records"][0]["thermal"]
    thermal.pop("reference_volume_law")
    thermal.pop("thermal_expansion_law")
    validate_eosmat_document(document)

    thermal["thermal_expansion_law"] = "unknown"
    with pytest.raises(ValueError, match="thermal_expansion_law is invalid"):
        validate_eosmat_document(document)

    thermal["thermal_expansion_law"] = "linear_temperature"
    with pytest.raises(ValueError, match="requires alpha1"):
        validate_eosmat_document(document)

    thermal["parameters"]["alpha1"] = 3.1e-9
    thermal["reference_volume_law"] = "linear_temperature"
    with pytest.raises(ValueError, match="requires constant thermal expansion"):
        validate_eosmat_document(document)

    thermal["reference_volume_law"] = "unknown"
    with pytest.raises(ValueError, match="reference_volume_law is invalid"):
        validate_eosmat_document(document)

    thermal["reference_volume_law"] = "berman"
    thermal["thermal_expansion_law"] = "constant"
    thermal["parameters"]["alpha1"] = 0.0
    with pytest.raises(ValueError, match="berman reference volume requires"):
        validate_eosmat_document(document)

    thermal["reference_volume_law"] = "integrated_expansivity"
    thermal["type"] = "MieGruneisenDebye"
    thermal["model"] = "mie_gruneisen_debye"
    with pytest.raises(ValueError, match="requires AlphaKT"):
        validate_eosmat_document(document)


def test_migration_manifest_does_not_claim_a_dioptas_data_license():
    root = resources.files("peritheos.data.materials")
    manifest = json.loads(root.joinpath("manifest.json").read_text(encoding="utf-8"))

    assert manifest["source"]["project"] == "Dioptas"
    assert manifest["source"]["version"] == "0.10.0"
    assert "license" not in manifest["source"]
    assert not root.joinpath("DIOPTAS_LICENSE.txt").is_file()
    assert manifest["materials"] == 116
    assert manifest["eos_records"] == 163
    assert manifest["scientific_validation"]["audit_date"] == "2026-09-04"


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda d: d.update(format="unknown"), "Unsupported material format"),
        (lambda d: d.update(format_version=2), "format 3"),
        (lambda d: d.update(name=None), "name must be a string"),
        (lambda d: d.update(identifier=3), "identifier must be a string"),
        (lambda d: d.update(lattice=[]), "lattice must be a JSON object"),
        (
            lambda d: d["lattice"].update(a=float("nan")),
            "lattice.a must be finite",
        ),
        (
            lambda d: d.update(formula_units_per_cell=0),
            "formula_units_per_cell must be greater than zero",
        ),
        (lambda d: d.update(peaks={}), "peaks must be a JSON array"),
        (
            lambda d: d["eos_records"][0].update(label=None),
            "label must be a string",
        ),
        (
            lambda d: d["eos_records"][0].update(reference=None),
            "reference must be a string or object",
        ),
        (
            lambda d: d["eos_records"][0]["eos"].update(type="unknown"),
            "unsupported",
        ),
        (
            lambda d: d["eos_records"][0]["eos"]["parameters"].pop("V0"),
            "requires V0",
        ),
        (
            lambda d: d["eos_records"][0].update(fixed_parameters={}),
            "fixed_parameters must be an array",
        ),
        (
            lambda d: d["eos_records"][0].update(
                experimental_pressure_range_gpa=[2, 1]
            ),
            "must be ordered",
        ),
    ],
)
def test_structural_validator_rejects_invalid_documents(mutation, message):
    document = copy.deepcopy(get_material_document("gold"))
    mutation(document)
    with pytest.raises(ValueError, match=message):
        validate_eosmat_document(document)


def test_validator_rejects_duplicate_record_id_and_multiple_defaults():
    document = get_material_document("gold")
    document["eos_records"][1]["identifier"] = document["eos_records"][0]["identifier"]
    with pytest.raises(ValueError, match="Duplicate EOS record identifier"):
        validate_eosmat_document(document)

    document = get_material_document("gold")
    document["eos_records"][0]["default"] = True
    document["eos_records"][1]["default"] = True
    with pytest.raises(ValueError, match="at most one default"):
        validate_eosmat_document(document)
