import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "ca_perovskite_tetragonal_chen_2018_bm2"


def test_chen_2018_record_preserves_fit_and_executes():
    document = get_material_document("ca_perovskite_tetragonal")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == IDENTIFIER
    )
    assert record["reference"]["doi"] == "10.2138/am-2018-6087"
    assert record["eos"]["parameters"] == {"V0": 185.2, "K0": 223.0, "K0_prime": 4.0}
    assert record["parameter_errors"] == {"V0": 0.4, "K0": 6.0, "K0_prime": None}
    assert record["fixed_parameters"] == ["K0_prime"]
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    assert eos.pressure(185.2) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(0.85 * 185.2)) == pytest.approx(
        0.85 * 185.2, rel=1e-10
    )


def test_chen_2018_reproduction_script():
    path = ROOT / "scripts/reproduce_chen_2018_casio3.py"
    spec = importlib.util.spec_from_file_location("reproduce_chen_2018", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["v0_conventional_cell_a3"] == pytest.approx(185.2)
    assert result["pressures_gpa"]["1.00"] == pytest.approx(0.0, abs=1e-12)
    assert result["pressures_gpa"]["0.90"] == pytest.approx(29.0126237818)
    assert result["pressures_gpa"]["0.80"] == pytest.approx(77.8232952174)
