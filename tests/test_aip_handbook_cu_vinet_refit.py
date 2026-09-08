import json
from pathlib import Path

import pytest

from scripts.reproduce_aip_handbook_cu_vinet_refit import (
    load_observations,
    refit,
)


def test_checked_cu_transcription_includes_corrected_printing_error():
    pressure, relative_volume = load_observations()
    assert len(pressure) == 25
    assert pressure[0] == pytest.approx(1.5)
    assert relative_volume[0] == pytest.approx(0.9900)
    assert pressure[-1] == pytest.approx(34.0)
    assert relative_volume[-1] == pytest.approx(0.849)


def test_cu_vinet_l8_refit_is_reproducible_and_improves_the_objective():
    result = refit()
    fitted = result["refit"]
    published = result["published"]

    assert result["observations"] == 25
    assert fitted["k0_gpa"] == pytest.approx(139.23431, abs=2.0e-5)
    assert fitted["k0_prime"] == pytest.approx(5.001127, abs=2.0e-6)
    assert fitted["rmse_gpa"] == pytest.approx(0.07015421, abs=1.0e-8)
    assert fitted["objective8_gpa8"] < published["objective8_gpa8"]


def test_cu_refit_is_a_separate_dataset_backed_catalog_record():
    root = Path(__file__).resolve().parents[1]
    path = root / "peritheos" / "data" / "materials" / "cu_handbook_1972_legacy.eosmat"
    document = json.loads(path.read_text(encoding="utf-8"))
    records = {record["identifier"]: record for record in document["eos_records"]}
    datasets = {dataset["identifier"]: dataset for dataset in document["datasets"]}

    fitted = records["cu_sun_2010_low_vn_refit"]
    assert len(records) == 1
    assert fitted["record_kind"] == "refit"
    assert "derived_from_record" not in fitted
    assert fitted["fit_provenance"]["dataset"] == "aip_handbook_table4d12_cu"
    assert fitted["fit_provenance"]["objective"] == (
        "sum((P_model_gpa - P_table_gpa)^8)"
    )
    assert fitted["eos"]["parameters"]["K0"] == pytest.approx(139.23431084)
    assert fitted["eos"]["parameters"]["K0_prime"] == pytest.approx(5.00112704)
    assert datasets["aip_handbook_table4d12_cu"]["used_by_eos_records"] == [
        fitted["identifier"]
    ]
