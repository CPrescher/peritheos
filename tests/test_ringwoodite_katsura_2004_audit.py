import csv
import runpy
from pathlib import Path

import pytest

from peritheos import get_material_document

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "docs" / "data" / "ringwoodite-katsura-2004-table2-pvt.csv"
SCRIPT = ROOT / "scripts" / "reproduce_katsura_2004_ringwoodite.py"


def test_katsura_2004_table2_transcription_and_no_production_record():
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 127
    assert rows[0] == {
        "run": "M038",
        "temperature_k": "1500",
        "pressure_gpa": "21.84",
        "pressure_uncertainty_gpa": "0.09",
        "normalized_volume": "0.9300",
        "normalized_volume_uncertainty": "0.0003",
    }
    assert rows[33] == {
        "run": "M086",
        "temperature_k": "300",
        "pressure_gpa": "-0.01",
        "pressure_uncertainty_gpa": "0.03",
        "normalized_volume": "1.0001",
        "normalized_volume_uncertainty": "0.0002",
    }
    assert rows[-1]["run"] == "M092"
    assert rows[-1]["pressure_gpa"] == "19.74"
    assert rows[-1]["normalized_volume"] == "0.9443"

    document = get_material_document("ringwoodite")
    dois = {record["reference"]["doi"].lower() for record in document["eos_records"]}
    assert "10.1029/2004jb003094" not in dois
    assert {record["identifier"] for record in document["eos_records"]} == {
        "ringwoodite_meng_1994_bm3_1"
    }


def test_katsura_2004_reproduction_documents_hidden_normalization_blocker():
    reproduction = runpy.run_path(str(SCRIPT))["reproduce"]()

    assert reproduction["primary_data"] == {
        "rows": 127,
        "pressure_range_gpa": [-0.01, 23.18],
        "temperature_range_k": [300.0, 2000.0],
        "volume_ratio_range": [0.9106, 1.0001],
        "room_temperature_rows": 17,
    }
    chemical = reproduction["published_mgd_with_chemical_n_7"]
    assert chemical["pressure_rmse_gpa"] == pytest.approx(1.8899360537)
    assert chemical["heated_pressure_rmse_gpa"] == pytest.approx(2.0134196233)
    assert chemical["maximum_absolute_pressure_residual_gpa"] == pytest.approx(
        3.2636162190
    )

    alternate = reproduction["published_mgd_with_n_7_variable_exponent_theta"]
    assert alternate["pressure_rmse_gpa"] == pytest.approx(1.9048891592)

    bm3 = reproduction["unweighted_300_304_k_bm3_refit"]
    assert bm3["K0_prime"] == pytest.approx(4.8379027256)
    assert bm3["K0_prime_standard_error"] == pytest.approx(0.2297162097)

    thermal = reproduction["unweighted_mgd_refit_with_chemical_n_7"]
    assert thermal["theta0"] == pytest.approx(1035.4122650)
    assert thermal["gamma0"] == pytest.approx(1.3996316539)
    assert thermal["q"] == pytest.approx(2.8265592264, abs=1.0e-5)
    assert thermal["pressure_rmse_gpa"] == pytest.approx(0.2962327776)

    diagnostic = reproduction["n_5_hidden_normalization_diagnostic_only"]
    assert diagnostic["heated_pressure_rmse_gpa"] == pytest.approx(0.2421931692)
    assert reproduction["decision"].startswith("blocked:")
