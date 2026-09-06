import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
RECORD = "mgo_liquid_ghosh_karki_2016_3000k_bm3_1"


def test_ghosh_karki_liquid_mgo_record_executes():
    document = get_material_document("mgo_liquid")
    record = document["eos_records"][0]
    assert record["identifier"] == RECORD
    assert record["reference"]["doi"] == "10.1038/srep37269"
    assert record["temperature_ref"] == 3000.0
    assert record["eos"]["parameters"] == {
        "V0": 30.0,
        "K0": 16.6,
        "K0_prime": 6.14,
    }
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": 1.0,
        "molar_mass_g_mol": 40.304,
    }
    eos = Material.from_eosmat(document).get_eos_record(RECORD)
    assert eos.pressure(30.0) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(20.0)) == pytest.approx(20.0, rel=1e-10)


def test_ghosh_karki_reproduction_checkpoints():
    path = ROOT / "scripts/reproduce_ghosh_karki_2016_mgo_liquid.py"
    spec = importlib.util.spec_from_file_location("reproduce_mgo_liquid", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["temperature_k"] == 3000.0
    assert result["pressures_gpa"][0] == pytest.approx(0.0, abs=1e-12)
    assert result["pressures_gpa"][2] == pytest.approx(5.265502349081)
    assert result["pressures_gpa"][4] == pytest.approx(22.757226429494)
