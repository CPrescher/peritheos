import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "mgsio3_post_perovskite_liu_2010_lda_static_bm3"


def test_liu_2010_record_preserves_coefficients_and_executes():
    document = get_material_document("mgsio3_post_perovskite")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == IDENTIFIER
    )
    assert record["reference"]["doi"] == "10.1142/S0217984910022391"
    assert record["eos"]["parameters"] == {"V0": 163.3, "K0": 219.3, "K0_prime": 4.4}
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": 4.0,
        "molar_mass_g_mol": 100.387,
    }
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    assert eos.pressure(163.3) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(0.8 * 163.3)) == pytest.approx(
        0.8 * 163.3, rel=1e-10
    )


def test_liu_2010_reproduction_script():
    path = ROOT / "scripts/reproduce_liu_2010_mgsio3_post_perovskite.py"
    spec = importlib.util.spec_from_file_location("reproduce_liu_2010", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()["pressures_gpa"]
    assert result["1.00"] == pytest.approx(0.0, abs=1e-12)
    assert result["0.90"] == pytest.approx(29.1540797098)
    assert result["0.70"] == pytest.approx(172.8919440328)
