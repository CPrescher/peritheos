"""Primary-source checks for Sakai et al. (2025) nine-material R-S scales."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_sakai_2025_nine_materials import PARAMETERS, reproduce

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = tuple(PARAMETERS)
HASHES = {
    "copper": "43b9be203ad89da7c91bcaf6333c890142becd691191bfd87de1cd317fbc9ca5",
    "rhenium": "ddf159d5e99f776f368cc977463e6eb219e501bde08da6fae41f6126628738a1",
    "platinum": "f3f2a2737dc29add29775e9973e32e37d4372fc3fa0b1aff846df8646259a6cf",
    "tungsten": "96efd4a74899ae8c6e3d17214c2a26eb32aa71057d36fe87940b6b82fe19dd01",
    "gold": "c6ce297beef3b23e6be063b8809599323a7412bb5aa835bcaeb4293f16624a65",
    "molybdenum": "1f3bcfeda9c29d1e08328868d491f13d8d105bcc44b977e28fa749844566968c",
    "mgo": "2593fde0a35bd11bc4dbea03fb86e7b9e88687f1510ea8c174c9162d2c2f0652",
    "nacl_b2": "241365ad28e7bcce5d378965c64b91ff09d8a21d92ef24b50bf7ff98fcb86d19",
    "iron": "3e9705f05ddf2bd674b6eae14d7de0281dd42afcf89e986a022ab1174867bfe2",
}


@pytest.mark.parametrize("material", MATERIALS)
def test_record_is_executable_and_uses_the_source_model(material):
    document = get_material_document(material)
    identifier = f"{material}_sakai_2025_rydberg_stacey_1"
    record = next(
        item for item in document["eos_records"] if item["identifier"] == identifier
    )

    assert record["reference"]["doi"] == "10.1038/s43246-025-00792-5"
    assert record["eos"]["type"] == "RydbergStacey"
    assert record["eos"]["model"] == "rydberg_stacey"
    model = (
        Material.from_eosmat(document, record_identifiers=[identifier])
        .eos_records[0]
        .eos
    )
    assert model.pressure(model.V0) == pytest.approx(0.0)
    assert model.bulk_modulus(model.V0) == pytest.approx(model.K0)


@pytest.mark.parametrize("material", MATERIALS)
def test_official_table_s9_dataset_checksum(material):
    path = (
        ROOT
        / "peritheos"
        / "data"
        / "datasets"
        / f"{material}-sakai-2025-table-s9-rs-eos-grid.csv"
    )
    assert hashlib.sha256(path.read_bytes()).hexdigest() == HASHES[material]


def test_official_derived_grid_is_reproduced_from_rounded_table2_coefficients():
    metrics = reproduce()
    assert set(metrics) == set(MATERIALS)
    assert all(values["rows"] == 26 for values in metrics.values())
    # Table S9 was generated with the authors' internal fitted precision, while
    # Table 2 publishes rounded coefficients.  The sub-0.11 GPa discrepancy is
    # therefore the expected rounding-limited reproduction tolerance.
    assert max(values["max_abs_gpa"] for values in metrics.values()) < 0.11
