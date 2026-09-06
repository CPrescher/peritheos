import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "mg092fe008sio3_bridgmanite_gong_2004_adiabatic_bm3"


def test_gong_2004_record_preserves_shock_reduced_fit():
    document = get_material_document("mg092fe008sio3_bridgmanite")
    record = document["eos_records"][0]
    assert record["identifier"] == IDENTIFIER
    assert record["reference"]["doi"] == "10.1029/2003GL019132"
    assert record["eos"]["parameters"] == {
        "V0": pytest.approx(163.13738187176514),
        "K0": 260.1,
        "K0_prime": 4.18,
    }
    assert record["parameter_errors"] == {"V0": None, "K0": 0.9, "K0_prime": 0.04}
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": 4.0,
        "molar_mass_g_mol": 102.9102,
    }
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    assert eos.pressure(record["eos"]["parameters"]["V0"]) == pytest.approx(
        0.0, abs=1e-12
    )
    assert eos.pressure(0.9 * record["eos"]["parameters"]["V0"]) == pytest.approx(
        34.1718057495
    )


def test_gong_2004_primary_table_and_conversion_reproduction():
    path = ROOT / "scripts/reproduce_gong_2004_mg092fe008sio3.py"
    spec = importlib.util.spec_from_file_location("reproduce_gong_2004", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["rows"] == 13
    assert result["derived_v0_a3_conventional_z4"] == pytest.approx(163.13738187176514)
    assert result["us_up_rms_residual_km_s"] < 0.18
    assert result["pressure_at_0_9_v0_gpa"] == pytest.approx(34.1718057495)
    assert result["pressure_at_0_8_v0_gpa"] == pytest.approx(92.7360907675)
