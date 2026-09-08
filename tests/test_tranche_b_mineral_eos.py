import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import (
    Material,
    get_material,
    list_eos_record_documents,
    validate_eosmat_document,
)
from peritheos.eosmat import load_eosmat
from scripts.reproduce_tranche_b_mineral_eos import (
    BM2,
    BM3,
    THERMAL_BM2,
    VINET,
    reproduce,
)

ROOT = Path(__file__).parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"
AUDITED_DOIS = {
    "10.1029/2001gl013523",
    "10.1038/srep15534",
    "10.1029/93jb02175",
    "10.1029/94gl00976",
    "10.1134/s0016702921080073",
    "10.1029/2001gl012910",
    "10.1029/94gl01592",
    "10.1038/s41598-020-66340-y",
    "10.1088/1674-0068/24/06/703-710",
    "10.1126/sciadv.1600427",
    "10.1029/2008gl036759",
    "10.1088/0953-8984/19/24/246103",
    "10.1038/s41467-018-07265-z",
    "10.1103/physrevlett.89.255507",
    "10.1103/physrevb.59.r14141",
}
ACCEPTED_DOIS = AUDITED_DOIS - {
    "10.1029/2001gl012910",
    "10.1029/94gl01592",
    "10.1103/physrevlett.89.255507",
    "10.1103/physrevb.59.r14141",
}
FILES = [
    "ca_perovskite.eosmat",
    "sio2_stv_andr.eosmat",
    "phase_h.eosmat",
    "wollastonite.eosmat",
    "pseudowollastonite.eosmat",
    "breyite.eosmat",
    "casio2o5_titanite.eosmat",
    "larnite.eosmat",
    "ca_perovskite_tetragonal.eosmat",
    "bridgmanite.eosmat",
    "fe088sio3_bridgmanite.eosmat",
    "klb1_mg_perovskite.eosmat",
    "klb1_ca_perovskite.eosmat",
    "klb1_ferropericlase.eosmat",
    "coesite.eosmat",
    "coesite_iv.eosmat",
    "coesite_v.eosmat",
]


def _records() -> list[tuple[dict, Material]]:
    found = []
    for filename in FILES:
        document = load_eosmat(MATERIALS / filename)
        validate_eosmat_document(document)
        if not document["eos_records"]:
            continue
        material = Material.from_eosmat(document)
        for record in document["eos_records"]:
            if record["reference"].get("doi", "").lower() in ACCEPTED_DOIS:
                found.append((record, material))
    return found


def test_twenty_four_primary_source_records_have_exact_parameters_or_reproduced_derivations():
    records = _records()
    assert len(records) == 24
    observed = {
        r["identifier"]: tuple(r["eos"]["parameters"].values()) for r, _ in records
    }
    assert observed == {**BM3, **BM2, **VINET}
    assert (
        sum(r["reference"]["doi"].lower() == "10.1029/2001gl013523" for r, _ in records)
        == 4
    )
    assert (
        sum(
            r["reference"]["doi"].lower() == "10.1134/s0016702921080073"
            for r, _ in records
        )
        == 7
    )
    assert (
        sum(r["reference"]["doi"].lower() == "10.1029/2008gl036759" for r, _ in records)
        == 4
    )


def test_all_accepted_curves_execute_and_anchor_zero_pressure():
    for record, material in _records():
        model = material.get_eos_record(record["identifier"])
        v0 = record["eos"]["parameters"]["V0"]
        if record["identifier"] in THERMAL_BM2:
            assert model.pressure(v0, 300.0) == pytest.approx(0.0, abs=1e-12)
            assert np.isfinite(model.pressure(0.9 * v0, 2000.0))
            assert model.pressure(0.9 * v0, 2000.0) > 0.0
        else:
            assert model.pressure(v0) == pytest.approx(0.0, abs=1e-12)
            assert np.isfinite(model.pressure(0.9 * v0))
            assert model.pressure(0.9 * v0) > 0.0


def test_phase_h_primary_table_is_complete_unchanged_and_reproduced():
    path = (
        ROOT / "peritheos/data/datasets/phase-h-tsuchiya-mookherjee-2015-table1-pv.csv"
    )
    with path.open(newline="", encoding="utf-8") as stream:
        assert len(list(csv.DictReader(stream))) == 12
    assert hashlib.sha256(path.read_bytes()).hexdigest() == (
        "4ad82028fa5483acc94f7d1218f1b13176c25f178ebe3750cdfb11309a0ccaf4"
    )
    result = reproduce()
    assert result["accepted_record_count"] == 24
    assert result["phase_h_table1"]["observations"] == 12
    assert result["phase_h_table1"]["published_bm3_pressure_rmse_gpa"] == pytest.approx(
        0.5462582, abs=1e-6
    )
    assert result["phase_h_table1"][
        "published_bm3_max_abs_pressure_residual_gpa"
    ] == pytest.approx(1.1329591, abs=1e-6)
    assert result["bykova_2018_table10"]["coesite_iv"]["observations"] == 3
    assert result["bykova_2018_table10"]["coesite_iv"][
        "pressure_rmse_gpa"
    ] == pytest.approx(0.080812, abs=1e-6)
    assert (
        result["bykova_2018_table10"]["coesite_v_held_candidate"]["observations"] == 1
    )
    assert (
        result["bykova_2018_table10"]["coesite_v_held_candidate"][
            "max_abs_pressure_residual_gpa"
        ]
        < 1e-5
    )


def test_bykova_combined_coesite_tables_and_v_held_candidate_are_preserved():
    source_tables = {
        "coesite-i-ii-cernok-2014-table1-pv.csv": (
            9,
            "733b94709bd150635b1fe5e7975a652055cc745e604dad723f2b06146b3ce575",
        ),
        "coesite-ii-iii-bykova-2018-table2-pv.csv": (
            15,
            "8aca860f04909fc827fc09ea6e24039d0364dcc21f3c774b07cb03cdcfe945e2",
        ),
    }
    for filename, (row_count, digest) in source_tables.items():
        path = ROOT / "peritheos" / "data" / "datasets" / filename
        with path.open(newline="", encoding="utf-8") as stream:
            assert len(list(csv.DictReader(stream))) == row_count
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest

    result = reproduce()
    combined = result["bykova_2018_combined_coesite_i_ii_iii"]
    assert combined["observations"] == 24
    assert combined["phase_counts"] == {
        "coesite-I": 7,
        "coesite-II": 6,
        "coesite-III": 11,
    }
    assert combined["pressure_range_gpa"] == [2.42, 36.9]
    assert combined["published_bm3_pressure_rmse_gpa"] == pytest.approx(1.2698145841)
    assert combined["partial_unweighted_refit"] == pytest.approx(
        {
            "V0": 542.2160125743,
            "K0": 126.3278193724,
            "K0_prime": 1.6951441000,
            "pressure_rmse_gpa": 0.6967701889,
            "max_abs_pressure_residual_gpa": 1.9967422491,
        }
    )

    sensitivity = result["bykova_2018_coesite_v_rounding_sensitivity"]
    assert sensitivity["anchor_pressure_interval_gpa"] == [56.5, 57.5]
    assert sensitivity["v0_interval_a3"] == pytest.approx(
        [426.7607834181, 428.0221644012]
    )
    assert sensitivity["maximum_abs_pressure_shift_26_64_gpa"] == pytest.approx(
        0.5244076465
    )

    coesite_v = load_eosmat(MATERIALS / "coesite_v.eosmat")
    assert coesite_v["eos_records"] == []
    with pytest.raises(KeyError, match="Unknown material 'coesite_v'"):
        get_material("coesite_v")
    assert "coesite_v_bykova_2018_am05_static_bm3_refit" not in (
        list_eos_record_documents()
    )
    candidate = coesite_v["source"]["held_eos_candidate"]
    assert candidate["status"] == "held_non_executable"
    assert candidate["missing_required_parameters"] == ["V0"]
    assert candidate["published_parameters"] == {
        "K0_gpa": 185.26,
        "K0_prime": 3.1,
    }
    assert candidate["diagnostic_reconstruction"]["not_for_quantitative_use"]
    resource = candidate["source_data"]["resource"]
    source_path = ROOT / "peritheos" / "data" / resource["path"]
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == resource["sha256"]


def test_hold_papers_do_not_create_production_records():
    found = []
    for path in MATERIALS.glob("*.eosmat"):
        document = json.loads(path.read_text(encoding="utf-8"))
        found.extend(
            record
            for record in document.get("eos_records", [])
            if record["reference"].get("doi", "").lower()
            in {
                "10.1029/2001gl012910",
                "10.1029/94gl01592",
                "10.1103/physrevlett.89.255507",
                "10.1103/physrevb.59.r14141",
            }
        )
    assert found == []


def test_audit_disposes_every_candidate_under_the_fifteen_dois_once():
    audit = (
        ROOT / "docs/literature-reproductions/tranche-b-mineral-eos-audit.md"
    ).read_text(encoding="utf-8")
    candidates = json.loads(
        (ROOT / "docs/data/litcurate-eos-candidates.json").read_text(encoding="utf-8")
    )["records"]
    same_doi = [
        row
        for row in candidates
        if row["publication"].get("doi", "").lower() in AUDITED_DOIS
    ]
    assert len(same_doi) == 93
    assert all(audit.count(row["identifier"]) == 1 for row in same_doi)
