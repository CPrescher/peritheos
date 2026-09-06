import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "cairo3_post_perovskite_martin_2007_bm2"


def test_martin_2007_record_preserves_fit_and_executes():
    document = get_material_document("cairo3_post_perovskite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == IDENTIFIER
    )
    assert record["identifier"] == IDENTIFIER
    assert record["reference"]["doi"] == "10.2138/am.2007.2473"
    assert record["eos"]["parameters"] == {"V0": 226.632, "K0": 180.2, "K0_prime": 4.0}
    assert record["fixed_parameters"] == ["K0_prime"]
    assert record["parameter_errors"] == {"V0": 0.045, "K0": 2.8, "K0_prime": None}
    assert record["volume_basis"]["molar_mass_g_mol"] == 280.292
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    assert eos.pressure(226.632) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(0.9 * 226.632)) == pytest.approx(
        0.9 * 226.632, rel=1e-10
    )


def test_martin_2007_primary_table_reproduction():
    path = ROOT / "scripts/reproduce_martin_2007_cairo3.py"
    spec = importlib.util.spec_from_file_location("reproduce_martin_2007", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["rows"] == 15
    assert result["rms_pressure_residual_gpa"] < 0.6
    assert result["max_abs_pressure_residual_gpa"] < 1.7
    assert result["pressure_at_0_9_v0_gpa"] == pytest.approx(23.4442816389)
