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
    spin_populations,
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


def test_solomatova_recovered_earlier_observations_are_complete_and_immutable():
    resources = {
        "mg065fe035o-chen-2012-table1-pv.csv": (
            36,
            "2c0a1be7f09c40bee120077216fd8b5a2f6f856a21a468d7387e5fd687e86eb1",
        ),
        "mg061fe039o-zhuravlev-2010-table5-pv.csv": (
            71,
            "7109f75899eb7025fa4233bbb955ad630412fe08ba0b80c08f1a5505ff4d5908",
        ),
        "mg061fe039o-fei-2007-table-s2-pv.csv": (
            77,
            "263f5c5ae3f780290827874235fe11442916eb0729f36cf161479bf1d6f81ccc",
        ),
        "mg090fe010o-marquardt-2009-table2-pv.csv": (
            29,
            "ac504ad2b8175dc1c0f199309cc1c143ec2744e4299305f699dd9a4410cb3b24",
        ),
        "mg083fe017o-lin-2005-figure2-300k-digitized.csv": (
            43,
            "39ec70f8a50dd5b3ae76484ee37fa8315ffab9329c80e576b8413721c04695ea",
        ),
    }
    for name, (count, checksum) in resources.items():
        path = ROOT / "peritheos/data/datasets" / name
        with path.open(newline="", encoding="utf-8") as stream:
            assert len(list(csv.DictReader(stream))) == count
        assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum

    zhuravlev_path = (
        ROOT / "peritheos/data/datasets/mg061fe039o-zhuravlev-2010-table5-pv.csv"
    )
    with zhuravlev_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert sum(row["used_in_solomatova_2016_fit"] == "1" for row in rows) == 54
    assert sum(row["used_in_zhuravlev_2010_hs_fit"] == "1" for row in rows) == 43
    assert sum(row["compression_path"] == "decompression" for row in rows) == 17

    fei_path = ROOT / "peritheos/data/datasets/mg061fe039o-fei-2007-table-s2-pv.csv"
    with fei_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert sum(row["compression_path"] == "compression" for row in rows) == 38
    assert sum(row["compression_path"] == "decompression" for row in rows) == 39


def test_zhuravlev_source_owned_high_spin_bm2_is_first_class_and_executable():
    document = json.loads(
        (ROOT / "peritheos/data/materials/mg061fe039o.eosmat").read_text(
            encoding="utf-8"
        )
    )
    identifier = "mg061fe039o_zhuravlev_2010_high_spin_bm2_5"
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )

    assert source["reference"]["doi"] == "10.1007/s00269-009-0347-6"
    assert source["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 77.4, "K0": 161},
    }
    assert source["parameter_errors"] == {"V0": 0.2, "K0": 3}
    assert source["experimental_pressure_range_gpa"] == [15.5, 70.9]
    assert "alternative_fits" in source["parameter_provenance"]

    executable = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).get_eos_record(identifier)
    assert executable.pressure(77.4) == pytest.approx(0.0)


def test_solomatova_three_level_populations_and_coupled_refits_are_reproducible():
    populations = spin_populations(np.array([-10.0, 0.0, 10.0]))
    assert populations[0, 2] > 0.999
    assert populations[2, 0] > 0.999
    assert populations.sum(axis=1) == pytest.approx(np.ones(3))

    diagnostics = reproduce()
    refits = diagnostics["coupled_refits"]
    assert refits["mg090_marq"]["ledger_status"] == "similar"
    assert refits["mg065_chen"]["ledger_status"] == "similar"
    assert refits["mg061_fei"]["ledger_status"] == "similar"
    assert refits["mg061_zhuravlev"]["ledger_status"] == "similar"
    assert refits["mg075_mao"]["ledger_status"] == "parity_not_achieved"
    assert refits["mg090_marq"]["observations"] == 29
    assert refits["mg065_chen"]["observations"] == 36
    assert refits["mg061_fei"]["observations"] == 38
    assert refits["mg061_zhuravlev"]["observations"] == 54
    assert refits["mg075_mao"]["observations"] == 42
    assert refits["mg090_marq"]["transition_pressure_50_percent_ls_gpa"] == (
        pytest.approx(53.821, abs=0.002)
    )
    assert refits["mg065_chen"]["transition_pressure_50_percent_ls_gpa"] == (
        pytest.approx(64.176, abs=0.002)
    )
    assert refits["mg061_fei"]["transition_pressure_50_percent_ls_gpa"] == (
        pytest.approx(56.561, abs=0.002)
    )
    assert refits["mg061_zhuravlev"]["transition_pressure_50_percent_ls_gpa"] == (
        pytest.approx(73.186, abs=0.002)
    )
    assert diagnostics["lin_normalized_shape_diagnostic"]["observations"] == 43
    assert diagnostics["lin_normalized_shape_diagnostic"]["status"] == (
        "normalization_anchored_diagnostic_only"
    )


def test_solomatova_material_provenance_distinguishes_data_from_constraints():
    documents = {document["identifier"]: document for document in _documents()}
    for material, dataset in {
        "mg090fe010o": "mg090fe010o_marquardt_2009_table2_pv",
        "mg065fe035o": "mg065fe035o_chen_2012_table1_pv",
        "mg061fe039o": "mg061fe039o_fei_2007_table_s2_pv",
    }.items():
        checks = [
            record["scientific_validation"]["primary_data_check"]
            for record in documents[material]["eos_records"]
            if "solomatova_2016" in record["identifier"]
            and dataset
            in record["scientific_validation"]["primary_data_check"].get(
                "dataset_identifiers", []
            )
        ]
        assert len(checks) == 2
        assert all(
            check["status"] == "complete_earlier_primary_table" for check in checks
        )

    zhuravlev_checks = [
        record["scientific_validation"]["primary_data_check"]
        for record in documents["mg061fe039o"]["eos_records"]
        if "solomatova_2016_zhuravlev" in record["identifier"]
    ]
    assert len(zhuravlev_checks) == 2
    assert all(
        check["status"] == "complete_earlier_primary_table"
        for check in zhuravlev_checks
    )

    lin_checks = [
        record["scientific_validation"]["primary_data_check"]
        for record in documents["mg083fe017o"]["eos_records"]
        if "solomatova_2016" in record["identifier"]
    ]
    assert all(
        check["status"] == "digitized_normalized_primary_figure" for check in lin_checks
    )
    mg40_findings = [
        record["scientific_validation"]["primary_data_check"]["finding"]
        for record in documents["mgfe60o"]["eos_records"]
        if "solomatova_2016" in record["identifier"]
    ]
    assert all(
        "constraints are not substituted for observations" in finding
        for finding in mg40_findings
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
