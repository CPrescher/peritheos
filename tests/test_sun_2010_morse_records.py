import csv
import hashlib
import json
from importlib import resources
from pathlib import Path

import pytest

from scripts.reproduce_sun_2010_morse_eos import verify_table

DOI = "10.1515/zna-2010-1-202"
TABLE = "sun-2010-table2-morse-eos-parameters.csv"
TABLE_SHA256 = "e71613962b7df52559b8633546a7513b902ce1b5d8d4a7666e0476464123ca68"


def _rows():
    resource = resources.files("peritheos").joinpath("data", "datasets", TABLE)
    with resource.open("rb") as handle:
        payload = handle.read()
    assert hashlib.sha256(payload).hexdigest() == TABLE_SHA256
    return list(csv.DictReader(payload.decode().splitlines()))


def test_sun_table_2_transcription_and_equation_invariants():
    rows = _rows()
    assert len(rows) == 50
    assert rows[0]["solid"] == "n-H2"
    assert rows[23]["solid"] == "Brass"
    assert rows[-1]["solid"] == "Nd"
    assert float(rows[22]["sms4_k0_gpa"]) == pytest.approx(147.3)
    assert float(rows[22]["sms4_k0_prime"]) == pytest.approx(5.8813)

    result = verify_table()
    assert result["rows"] == 50
    assert result["curves"] == 150
    assert result["largest_p0_gpa"] < 1.0e-12
    assert result["largest_k0_error_gpa"] < 1.0e-12
    assert result["largest_k0_prime_error"] < 1.0e-7


def test_sun_parameter_table_is_a_non_catalog_benchmark_fixture():
    records = []
    material_root = Path(__file__).parents[1] / "peritheos" / "data" / "materials"
    for path in material_root.glob("*.eosmat"):
        document = json.loads(path.read_text())
        records.extend(document.get("eos_records", []))

    assert not any(
        record.get("reference", {}).get("doi", "").lower() == DOI for record in records
    )
