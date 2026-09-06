import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]

RECORDS = {
    "mgo": {
        "mgo_b1_zhang_bukowinski_1991_mpib_bm3": (74.60, 180.0, 4.04, 4.0, 40.304),
    },
    "mgo_b2": {
        "mgo_b2_zhang_bukowinski_1991_mpib_bm3": (18.07, 175.0, 4.20, 1.0, 40.304),
    },
    "sio2_stv_andr": {
        "stishovite_zhang_bukowinski_1991_mpib_bm3": (
            45.54,
            378.0,
            6.4,
            2.0,
            60.083,
        ),
    },
}


@pytest.mark.parametrize(
    ("material_identifier", "record_identifier", "expected"),
    [
        (material, identifier, expected)
        for material, records in RECORDS.items()
        for identifier, expected in records.items()
    ],
)
def test_zhang_bukowinski_records_preserve_coefficients_and_execute(
    material_identifier, record_identifier, expected
):
    v0, k0, k0_prime, formula_units, molar_mass = expected
    document = get_material_document(material_identifier)
    record = next(
        row for row in document["eos_records"] if row["identifier"] == record_identifier
    )
    assert record["reference"]["doi"] == "10.1103/PhysRevB.44.2495"
    assert record["eos"]["parameters"] == {
        "V0": v0,
        "K0": k0,
        "K0_prime": k0_prime,
    }
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": formula_units,
        "molar_mass_g_mol": molar_mass,
    }
    eos = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).get_eos_record(record_identifier)
    assert eos.pressure(v0) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(0.9 * v0)) == pytest.approx(0.9 * v0, rel=1e-10)


def test_zhang_bukowinski_reproduction_script():
    path = ROOT / "scripts/reproduce_zhang_bukowinski_1991_oxides.py"
    spec = importlib.util.spec_from_file_location("reproduce_zhang_bukowinski", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert set(result["pressure_at_0.9_v0_gpa"]) == {
        "MgO_B1",
        "MgO_B2",
        "SiO2_stishovite",
    }
    assert result["pressure_at_0.9_v0_gpa"]["MgO_B1"] == pytest.approx(23.469382932694)
    assert result["pressure_at_0.9_v0_gpa"]["MgO_B2"] == pytest.approx(23.016261789559)
    assert result["pressure_at_0.9_v0_gpa"]["SiO2_stishovite"] == pytest.approx(
        55.619668433812
    )
