"""The installed static reconstruction cannot silently acquire thermal coverage."""

import hashlib
import json

import numpy as np
import pytest

from peritheos import (
    Material,
    get_eos_record,
    get_material,
    get_material_document,
    list_eos_records,
    load_eosmat,
    search_eos_records,
    validate_eosmat_document,
)
from peritheos.errors import MaterialError
from scripts.finger_1979_ringwoodite import (
    DATA,
    DATASET_ID,
    RECORD_ID,
    REPORT,
    ROOT,
    bm2,
    reproduce,
)
from scripts.reconstruction_plot_data import curve_comparison

THERMAL_ID = "ferrous_ringwoodite_armentrout_kavner_2011_bm2_mgd_candidate"
STAGE = ROOT / "docs/research/data/ferrous-ringwoodite-candidate.eosmat"


def test_static_only_catalog_preserves_source_and_separates_thermal_candidate():
    document = get_material_document("ferrous_ringwoodite")
    validate_eosmat_document(document)
    assert [r["identifier"] for r in document["eos_records"]] == [RECORD_ID]
    record = get_eos_record(RECORD_ID)
    material = get_material("ferrous_ringwoodite")
    assert material.default_record().identifier == RECORD_ID
    assert not record.is_thermal
    raw = document["eos_records"][0]
    assert raw["record_kind"] == "derived"
    assert raw["eos"]["parameters"] == {"V0": 558.7, "K0": 196.0}
    assert raw["parameter_errors"] == {"V0": 0.2, "K0": 6.0}
    assert raw["fixed_parameters"] == ["V0"]
    assert raw["validity"]["pressure_gpa"] == [0, 4]
    assert raw["validity"]["temperature_k"] == [296.15, 296.15]
    assert (
        raw["scientific_validation"]["numerical_reproduction"]["author_fit_reproduced"]
        is False
    )
    assert [s["identifier"] for s in document["datasets"]] == [DATASET_ID]
    staged = load_eosmat(STAGE)
    assert len(staged["eos_records"]) == 2
    assert any(r["identifier"] == THERMAL_ID for r in staged["eos_records"])
    assert THERMAL_ID not in {r.identifier for r in list_eos_records()}
    thermal = next(r for r in staged["eos_records"] if r["identifier"] == THERMAL_ID)
    assert thermal["scientific_validation"]["status"] == "not_reproduced"
    assert thermal["catalog_access"] == "explicit_selection"
    assert (
        Material.from_eosmat(material.to_eosmat()).to_eosmat() == material.to_eosmat()
    )


def test_independent_source_checks_uncertainties_and_native_inversion():
    result = reproduce()
    assert result == json.loads((ROOT / REPORT).read_text())
    assert result["source_sha256"] == hashlib.sha256(DATA.read_bytes()).hexdigest()
    assert result["source_checks"][0]["pressure_error_gpa"] is None
    assert result["source_checks"][0]["combined_coordinate_error_gpa"] is None
    assert all(
        r["within_combined_coordinate_error"] for r in result["source_checks"][1:]
    )
    assert result["diagnostic_K0_gpa"] == pytest.approx(196.2175873043)
    assert result["author_fit_reproduced"] is False
    record = get_eos_record(RECORD_ID)
    np.testing.assert_allclose(
        record.pressure(np.array([558.7, 553.5, 547.7]), 296.15),
        [0, 1.867379040453365, 4.0556076738036],
        atol=1e-10,
    )
    pressures = np.linspace(0, 4, 9)
    volumes = record.volume(pressures, 296.15, check_validity=True)
    np.testing.assert_allclose(bm2(volumes), pressures, atol=1e-10)
    for temperature in (296.14, 296.16, 300):
        with pytest.raises(MaterialError):
            record.volume(2, temperature, check_validity=True)
    for pressure in (-0.01, 4.01):
        with pytest.raises(MaterialError):
            record.volume(pressure, 296.15, check_validity=True)
    assert RECORD_ID in {
        r.identifier for r in search_eos_records(temperature_k=296.15, pressure_gpa=2)
    }
    assert RECORD_ID not in {
        r.identifier for r in search_eos_records(temperature_k=300)
    }
    assert RECORD_ID not in {
        r.identifier for r in search_eos_records(pressure_gpa=4.01)
    }


def test_stored_comparison_plots_use_the_named_two_compressed_states():
    document = load_eosmat(ROOT / "peritheos/data/materials/ferrous_ringwoodite.eosmat")
    comparison, reason = curve_comparison(document, document["eos_records"][0])
    assert comparison is not None, reason
    np.testing.assert_array_equal(comparison.series.pressure, [2.0, 4.0])
    np.testing.assert_array_equal(comparison.series.volume, [553.5, 547.7])
    assert comparison.source_rmse == pytest.approx(
        reproduce()["source_compressed_pressure_rms_gpa"]
    )
    assert comparison.refit_rmse == pytest.approx(
        reproduce()["diagnostic_compressed_pressure_rms_gpa"]
    )
