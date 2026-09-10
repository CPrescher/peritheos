import importlib.util
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_zhu_reproduction_and_schemas():
    path = ROOT / "scripts/reproduce_zhu_2025_pressure_standards.py"
    spec = importlib.util.spec_from_file_location("reproduce_zhu_2025", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert len(result["records"]) == 3
    expected = {
        "gold_ye_2017_vinet_300k": (5.897, 0.022, 28),
        "mgo_ye_2017_vinet_300k": (4.182, 0.019, 28),
        "platinum_ye_2017_vinet_300k": (5.230, 0.033, 39),
    }
    for identifier, (k0_prime, error, rows) in expected.items():
        refit = result["records"][identifier]
        assert refit["K0_prime"] == pytest.approx(k0_prime, abs=5e-4)
        assert refit["K0_prime_standard_error"] == pytest.approx(error, abs=5e-4)
        assert refit["observations"] == rows
        assert refit["rmse_pressure_difference_gpa"] < 0.6
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    validator = Draft202012Validator(schema)
    for material in ("platinum", "gold", "mgo"):
        document = get_material_document(material)
        assert list(validator.iter_errors(document)) == []
        ye_record = next(
            item
            for item in document["eos_records"]
            if item["identifier"] == f"{material}_ye_2017_vinet_300k"
        )
        assert ye_record["aliases"] == [f"{material}_zhu_2025_vinet_300k"]
        assert ye_record["reference"]["doi"] == "10.1002/2016JB013811"
        assert ye_record["scientific_validation"]["primary_data_check"]["status"] == (
            "bundled"
        )
        thermal_record = next(
            item
            for item in document["eos_records"]
            if item["identifier"] == f"{material}_zhu_2025_pvt"
        )
        assert thermal_record["equation_kind"] == "thermal"
        assert thermal_record["thermal"]["type"] == (
            "AsymptoticPowerLawMieGruneisenDebyeExcess"
        )
        assert thermal_record["fit_datasets"]
        assert result["thermal_records"][thermal_record["identifier"]]["parity"]

    assert sum(item[2] for item in module.THERMAL_DATASETS.values()) == 320
    assert result["thermal_reproduction_scope"]["status"] == ("thermal_refit_parity")
    expected_thermal_refits = {
        "gold_zhu_2025_pvt": (2.93599996, 2.62559912, 22),
        "platinum_zhu_2025_pvt": (2.75224510, 5.10989873, 85),
        "mgo_zhu_2025_pvt": (1.53118931, 1.43014746, 213),
    }
    for identifier, (gamma0, exponent, rows) in expected_thermal_refits.items():
        refit = result["thermal_records"][identifier]
        assert refit["gamma0"] == pytest.approx(gamma0, abs=5.0e-7)
        assert refit["b"] == pytest.approx(exponent, abs=5.0e-7)
        assert refit["observations"] == rows


def test_zhu_audit_enumerates_all_candidates():
    text = (
        ROOT / "docs/literature-reproductions/zhu-2025-pressure-standards.md"
    ).read_text(encoding="utf-8")
    assert "litcurate_b1a3d7a39cb0fae8" in text
    assert "litcurate_fdde71245d253100" in text
