import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "mgo_hama_suito_1996_qsm_static_vinet"


def test_hama_suito_record_preserves_coefficients_and_executes():
    document = get_material_document("mgo")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == IDENTIFIER
    )
    assert record["reference"]["doi"] == "10.1088/0953-8984/8/1/008"
    assert record["eos"] == {
        "type": "Vinet",
        "model": "vinet",
        "parameters": {
            "V0": pytest.approx(73.34965396218291),
            "K0": 157.0,
            "K0_prime": 4.37,
        },
    }
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": 4.0,
        "molar_mass_g_mol": 40.304,
    }
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    assert eos.pressure(record["eos"]["parameters"]["V0"]) == pytest.approx(
        0.0, abs=1e-12
    )
    volume = 0.75 * record["eos"]["parameters"]["V0"]
    assert eos.volume(eos.pressure(volume)) == pytest.approx(volume, rel=1e-10)


def test_hama_suito_reproduction_script():
    path = ROOT / "scripts/reproduce_hama_suito_1996_mgo.py"
    spec = importlib.util.spec_from_file_location("reproduce_hama_suito", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["v0_angstrom3_conventional_cell"] == pytest.approx(73.34965396218291)
    assert result["pressures_gpa"]["1.00"] == pytest.approx(0.0, abs=1e-12)
    assert result["pressures_gpa"]["0.90"] == pytest.approx(20.7606898793)
    assert result["pressures_gpa"]["0.35"] == pytest.approx(1245.7155013613)
