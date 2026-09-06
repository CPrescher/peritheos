import importlib.util
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]

RECORDS = {
    "bridgmanite": {
        "bridgmanite_li_zeng_2009_gga_bm3": ("BM3", 168.12, 234.8, 4.199),
        "bridgmanite_li_zeng_2009_gga_vinet": ("Vinet", 168.08, 235.7, 4.294),
        "bridgmanite_li_zeng_2009_gga_natural_strain3": (
            "NaturalStrain3",
            168.04,
            237.0,
            4.384,
        ),
    },
    "mgsio3_post_perovskite": {
        "mgsio3_post_perovskite_li_zeng_2009_gga_bm3": (
            "BM3",
            168.29,
            223.73,
            4.152,
        ),
        "mgsio3_post_perovskite_li_zeng_2009_gga_vinet": (
            "Vinet",
            168.05,
            224.2,
            4.406,
        ),
        "mgsio3_post_perovskite_li_zeng_2009_gga_natural_strain3": (
            "NaturalStrain3",
            167.75,
            226.7,
            4.684,
        ),
    },
}


def _load_reproduction():
    path = ROOT / "scripts/reproduce_li_zeng_2009_mgsio3.py"
    spec = importlib.util.spec_from_file_location("reproduce_li_zeng_2009_mgsio3", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("material_identifier", "identifier", "expected"),
    [
        (material, identifier, values)
        for material, records in RECORDS.items()
        for identifier, values in records.items()
    ],
)
def test_li_zeng_record_executes(material_identifier, identifier, expected):
    eos_type, v0, k0, k0_prime = expected
    document = get_material_document(material_identifier)
    record = next(
        row for row in document["eos_records"] if row["identifier"] == identifier
    )
    assert record["reference"]["doi"] == "10.1142/S0129183109014242"
    assert record["eos"]["type"] == eos_type
    assert record["eos"]["parameters"] == {
        "V0": v0,
        "K0": k0,
        "K0_prime": k0_prime,
    }
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": 4.0,
        "molar_mass_g_mol": 100.387,
    }
    eos = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).get_eos_record(identifier)
    assert eos.pressure(v0) == pytest.approx(0.0, abs=1e-10)
    assert eos.pressure(0.8 * v0) > 75.0


def test_li_zeng_reproduction_and_exact_doi_count():
    result = _load_reproduction().reproduce()
    expected_ids = {
        identifier for records in RECORDS.values() for identifier in records
    }
    assert set(result["records"]) == expected_ids
    pressures = [row["pressure_at_0.8_v0_gpa"] for row in result["records"].values()]
    assert min(pressures) == pytest.approx(79.5057322895, abs=1e-9)
    assert max(pressures) == pytest.approx(83.9029151590, abs=1e-9)

    doi_records = []
    for path in (ROOT / "peritheos/data/materials").glob("*.eosmat"):
        document = json.loads(path.read_text(encoding="utf-8"))
        doi_records.extend(
            record
            for record in document.get("eos_records", [])
            if record.get("reference", {}).get("doi", "").lower()
            == "10.1142/s0129183109014242"
        )
    assert {record["identifier"] for record in doi_records} == expected_ids


@pytest.mark.parametrize("material_identifier", RECORDS)
def test_li_zeng_materials_match_normative_schema(material_identifier):
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    document = get_material_document(material_identifier)
    assert list(Draft202012Validator(schema).iter_errors(document)) == []


def test_li_zeng_audit_enumerates_all_candidates():
    text = (ROOT / "docs/literature-reproductions/li-zeng-2009-mgsio3.md").read_text(
        encoding="utf-8"
    )
    candidate_ids = {
        "litcurate_7aff77eec160b0a9",
        "litcurate_bfcc5153b94a67c9",
        "litcurate_efaa92c0d651d646",
        "litcurate_bfbff5de9c3f67b5",
        "litcurate_45b95058dca79981",
        "litcurate_be47411b3936fcff",
        "litcurate_429a5c1961d33dda",
        "litcurate_0d351dc013942d75",
        "litcurate_a6db893a40f4491a",
        "litcurate_ea7e1400594f5b03",
    }
    assert all(identifier in text for identifier in candidate_ids)
