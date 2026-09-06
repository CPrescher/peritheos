import csv
import hashlib
import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]

RECORDS = {
    "srsio3_cubic_perovskite": {
        "srsio3_cubic_xiao_2013_experimental_bm2_1": (49.18, 211.0, 1.0),
        "srsio3_cubic_xiao_2013_gga_bm2_2": (49.97, 207.7, 1.0),
    },
    "srsio3_6h_perovskite": {
        "srsio3_6h_xiao_2013_gga_bm2_1": (311.23, 183.8, 6.0),
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
def test_xiao_records_preserve_source_coefficients_and_execute(
    material_identifier, record_identifier, expected
):
    v0, k0, formula_units = expected
    document = get_material_document(material_identifier)
    record = next(
        row for row in document["eos_records"] if row["identifier"] == record_identifier
    )
    assert record["reference"]["doi"] == "10.2138/am.2013.4470"
    assert record["eos"]["parameters"] == {
        "V0": v0,
        "K0": k0,
        "K0_prime": 4.0,
    }
    assert record["fixed_parameters"] == ["K0_prime"]
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": formula_units,
        "molar_mass_g_mol": 163.702,
    }
    eos = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).get_eos_record(record_identifier)
    assert eos.pressure(v0) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(0.9 * v0)) == pytest.approx(0.9 * v0, rel=1e-10)


def test_xiao_table1_and_reproduction():
    resource = ROOT / "peritheos/data/datasets/srsio3-xiao-2013-table1-pv.csv"
    assert (
        hashlib.sha256(resource.read_bytes()).hexdigest()
        == "2cbdb31c1e4f5bd4ce8da0c4ea070055faafe8e05c11a53c0de0242e7ebbdb8c"
    )
    with resource.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 17
    assert sum(row["fit_included"] == "1" for row in rows) == 16
    assert rows[0]["exclusion_reason"] == "partial_amorphization_below_9_gpa"

    path = ROOT / "scripts/reproduce_xiao_2013_srsio3.py"
    spec = importlib.util.spec_from_file_location("reproduce_xiao", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert result["observations"] == 17
    assert result["selected_observations"] == 16
    assert result["experimental_pressure_rmse_gpa"] == pytest.approx(0.33655633847)
    assert result["experimental_max_abs_residual_gpa"] == pytest.approx(0.72327807014)
