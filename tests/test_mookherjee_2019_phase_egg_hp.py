import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
RECORD = "phase_egg_mookherjee_2019_bm3_hp_1"


def test_phase_egg_hp_record_is_distinct_and_executable():
    document = get_material_document("phase_egg")
    record = next(x for x in document["eos_records"] if x["identifier"] == RECORD)
    assert record["reference"]["doi"] == "10.2138/am-2019-6694"
    assert record["eos"]["parameters"] == {
        "V0": 207.74,
        "K0": 222.8,
        "K0_prime": 4.44,
    }
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": 4.0,
        "molar_mass_g_mol": 120.071,
    }
    eos = Material.from_eosmat(document, record_identifiers=[RECORD]).get_eos_record(
        RECORD
    )
    assert eos.pressure(207.74) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(190.0)) == pytest.approx(190.0, rel=1e-10)


def test_phase_egg_hp_official_table_and_reproduction():
    document = get_material_document("phase_egg")
    dataset = next(
        x
        for x in document["datasets"]
        if x["identifier"] == "phase_egg_mookherjee_2019_supplement_hp_pv"
    )
    assert dataset["rows"] == [
        [220.0, -11.3],
        [210.0, -2.4],
        [200.0, 9.2],
        [190.0, 24.3],
        [180.0, 44.0],
    ]

    path = ROOT / "scripts/reproduce_mookherjee_2019_phase_egg_hp.py"
    spec = importlib.util.spec_from_file_location("reproduce_phase_egg_hp", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["observations"] == 5
    assert result["published_pressure_rmse_gpa"] == pytest.approx(0.039665757737)
    assert result["published_max_abs_residual_gpa"] == pytest.approx(0.049096930304)
    assert result["diagnostic_unweighted_pv_fit"]["V0"] == pytest.approx(207.71787917)
