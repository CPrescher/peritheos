import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]

AKBER_RECORDS = {
    "bridgmanite": ("bridgmanite_akber_knutson_2005_gga_bm3", 167.00, 235.0, 3.84),
    "mgsio3_post_perovskite": (
        "mgsio3_post_perovskite_akber_knutson_2005_gga_bm3",
        167.32,
        204.0,
        4.18,
    ),
    "mg09375al0125si09375o3_bridgmanite": (
        "mg09375al0125si09375o3_bridgmanite_akber_knutson_2005_gga_bm3",
        167.44,
        232.0,
        3.86,
    ),
    "mg09375al0125si09375o3_post_perovskite": (
        "mg09375al0125si09375o3_post_perovskite_akber_knutson_2005_gga_bm3",
        167.16,
        207.0,
        4.13,
    ),
    "al2o3_perovskite": (
        "al2o3_perovskite_akber_knutson_2005_gga_bm3",
        171.32,
        205.0,
        4.03,
    ),
    "al2o3_post_perovskite": (
        "al2o3_post_perovskite_akber_knutson_2005_gga_bm3",
        167.12,
        201.0,
        4.29,
    ),
    "mgal0125si0875o29375_bridgmanite": (
        "mgal0125si0875o29375_bridgmanite_akber_knutson_2005_gga_bm3",
        169.44,
        214.0,
        3.96,
    ),
    "mgal0125si0875o29375_post_perovskite": (
        "mgal0125si0875o29375_post_perovskite_akber_knutson_2005_gga_bm3",
        169.32,
        194.0,
        4.15,
    ),
    "mgal00625h00625si09375o3_bridgmanite": (
        "mgal00625h00625si09375o3_bridgmanite_akber_knutson_2005_gga_bm3",
        168.80,
        228.0,
        3.85,
    ),
    "mgal00625h00625si09375o3_post_perovskite": (
        "mgal00625h00625si09375o3_post_perovskite_akber_knutson_2005_gga_bm3",
        166.56,
        228.0,
        3.86,
    ),
}


def _load_script(name):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("material_identifier", AKBER_RECORDS)
def test_akber_knutson_records_preserve_table1_and_execute(material_identifier):
    identifier, v0, k0, k0_prime = AKBER_RECORDS[material_identifier]
    document = get_material_document(material_identifier)
    records = [
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    ]
    assert len(records) == 1
    record = records[0]
    assert record["reference"]["doi"] == "10.1029/2005GL023192"
    assert record["temperature_ref"] == 0
    assert record["pressure_range_status"] == "theoretical"
    assert record["eos"]["parameters"] == {"V0": v0, "K0": k0, "K0_prime": k0_prime}
    eos = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).get_eos_record(identifier)
    assert eos.pressure(v0) == pytest.approx(0.0, abs=1e-12)
    assert eos.pressure(0.65 * v0) > 200.0
    compressed = eos.volume(eos.pressure(0.65 * v0))
    assert compressed == pytest.approx(0.65 * v0, rel=2e-10)


def test_akber_reproduction_covers_all_ten_rows_and_source_envelope():
    module = _load_script("reproduce_akber_knutson_2005_al_perovskites.py")
    result = module.reproduce()
    assert set(result["records"]) == {value[0] for value in AKBER_RECORDS.values()}
    high = [row["pressure_at_0.65_v0_gpa"] for row in result["records"].values()]
    assert min(high) == pytest.approx(205.91219499, abs=1e-8)
    assert max(high) == pytest.approx(230.83306253, abs=1e-8)


def test_chantel_table3_model_and_complete_pure_table1_transcription():
    document = get_material_document("bridgmanite")
    identifier = "bridgmanite_chantel_2012_bm3_mgd"
    record = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )
    assert record["reference"]["doi"] == "10.1029/2012GL053075"
    assert record["eos"]["parameters"] == {
        "V0": pytest.approx(162.2014560815),
        "K0": 252.0,
        "K0_prime": 4.1,
    }
    assert record["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 901.0,
        "gamma0": 1.44,
        "q": 1.4,
        "n": 5,
    }

    resource = (
        ROOT
        / "peritheos/data/datasets/bridgmanite-chantel-2012-table1-density-velocity.csv"
    )
    assert (
        hashlib.sha256(resource.read_bytes()).hexdigest()
        == "ea823b81c0285ccf10a159a21802d5f17ad29394ea195afad8c1ab88c87cbfb5"
    )
    with resource.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 11
    assert rows[0]["density_g_cm3"] == "4.110"
    assert rows[-1]["temperature_k"] == "1200"
    assert sum(row["room_temperature_fit_included"] == "1" for row in rows) == 9

    eos = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).get_eos_record(identifier)
    assert eos.pressure(record["eos"]["parameters"]["V0"], 300.0) == pytest.approx(
        0.0, abs=1e-10
    )
    assert math.isfinite(eos.pressure(150.0, 1200.0))


def test_chantel_reproduction_and_withheld_source_rows():
    module = _load_script("reproduce_chantel_2012_bridgmanite.py")
    result = module.reproduce()
    assert result["observations"] == 11
    assert result["room_temperature_observations"] == 9
    assert result["accepted_table3_model"]["pressure_rmse_gpa"] == pytest.approx(
        0.6284631087, abs=2e-9
    )
    assert result["table2_acoustic_reconstruction"]["K_S"] == pytest.approx(
        246.06536677, abs=2e-7
    )
    assert result["table2_acoustic_reconstruction"]["K_S_prime"] == pytest.approx(
        4.57064356, abs=2e-7
    )

    accepted = []
    for path in (ROOT / "peritheos/data/materials").glob("*.eosmat"):
        document = json.loads(path.read_text(encoding="utf-8"))
        accepted.extend(
            record
            for record in document.get("eos_records", [])
            if record.get("reference", {}).get("doi", "").lower()
            == "10.1029/2012gl053075"
        )
    assert [record["identifier"] for record in accepted] == [
        "bridgmanite_chantel_2012_bm3_mgd"
    ]


@pytest.mark.parametrize(
    ("audit_file", "candidate_ids"),
    [
        (
            "caracas-cohen-2005-chemistry.md",
            [
                "litcurate_b6be4414ac2d99ec",
                "litcurate_2bc936756b968aaa",
                "litcurate_a3e008f1122a57ff",
                "litcurate_aaa863c524972351",
                "litcurate_b3fa9fccaff96ed0",
                "litcurate_8c3b3f3681b4aa1d",
                "litcurate_2f8d0ff887301bc7",
                "litcurate_30e8a07d1eea0ef0",
                "litcurate_b94bd4825cec6951",
                "litcurate_3bf5908c5d967ad9",
                "litcurate_9324d92b5b5ebba0",
                "litcurate_60d4e2334a2418a2",
                "litcurate_0f7a60cd0d32c1fd",
                "litcurate_77dc1b77253aa71f",
                "litcurate_1681bdcdc165180c",
                "litcurate_70cc2c2137aa02f1",
            ],
        ),
        (
            "akber-knutson-2005-al-perovskites.md",
            [
                "litcurate_8f807ef87788f85f",
                "litcurate_2bf137ddeb5434cd",
                "litcurate_5272b35caf831df6",
                "litcurate_bf1450f777d5bbe8",
                "litcurate_29261d2e7516b73b",
                "litcurate_02a65cdb4bb78819",
                "litcurate_028263f4d2c0e0fc",
                "litcurate_6e8b79e9dce41899",
                "litcurate_2dc493b3352b5541",
                "litcurate_0a5eeb36d988137b",
            ],
        ),
        (
            "chantel-2012-bridgmanite.md",
            [
                "litcurate_31bc34e8292357ee",
                "litcurate_079361f39d81f154",
                "litcurate_542dc7a1f73d9d73",
                "litcurate_259fceb6dd4123fc",
                "litcurate_646cbed7a410ad5d",
                "litcurate_67cb8f3a4139db77",
                "litcurate_c6e72c7743e81ed6",
                "litcurate_13f91e0f265496af",
                "litcurate_82fa467c469361c3",
                "litcurate_32212ca6614abcc8",
                "litcurate_59accf80318ce824",
                "litcurate_ed458727faf26675",
                "litcurate_25bcc1466760438f",
                "litcurate_94fc588faf4261ea",
                "litcurate_e28b7d2629a0bd84",
                "litcurate_b0e9cd7bc1f10123",
                "litcurate_4c8341508b0d1a6a",
                "litcurate_4de85b1c154b864a",
                "litcurate_781d39e2b88b3db0",
                "litcurate_adaf5cd6f8d92989",
                "litcurate_b2e503f20b0da628",
            ],
        ),
    ],
)
def test_every_candidate_row_has_exactly_one_audit_disposition(
    audit_file, candidate_ids
):
    text = (ROOT / "docs/literature-reproductions" / audit_file).read_text(
        encoding="utf-8"
    )
    assert all(text.count(identifier) == 1 for identifier in candidate_ids)


def test_caracas_candidates_remain_out_of_production_without_primary_table():
    for path in (ROOT / "peritheos/data/materials").glob("*.eosmat"):
        document = json.loads(path.read_text(encoding="utf-8"))
        assert all(
            record.get("reference", {}).get("doi", "").lower() != "10.1029/2005gl023164"
            for record in document.get("eos_records", [])
        )
