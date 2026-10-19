"""Archive fixtures must exercise the same scientific records as the repository."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "crates/peritheos/tests/data"


@pytest.mark.parametrize(
    "material",
    [
        "brucite",
        "diopside",
        "forsterite",
        "ferrous_ringwoodite",
        "nacl_b2",
        "palladium",
        "rbcl",
        "rbbr",
        "rbi",
    ],
)
def test_crate_catalog_fixtures_preserve_current_material_bytes(material):
    assert (FIXTURES / f"{material}.eosmat").read_bytes() == (
        ROOT / f"peritheos/data/materials/{material}.eosmat"
    ).read_bytes()


@pytest.mark.parametrize(
    "fixture,candidate",
    [
        ("brucite-staged.eosmat", "brucite/brucite.eosmat"),
        (
            "brucite-diagnostic-archive.eosmat",
            "brucite/brucite-diagnostic-archive.eosmat",
        ),
        (
            "ferrous-ringwoodite-staged.eosmat",
            "ferrous-ringwoodite/ferrous_ringwoodite.eosmat",
        ),
    ],
)
def test_crate_staged_fixtures_preserve_qualified_candidate_bytes(fixture, candidate):
    assert (FIXTURES / fixture).read_bytes() == (
        ROOT / "output/material-candidates" / candidate
    ).read_bytes()
