import numpy as np
import pytest

from peritheos import get_eos_record, get_material_document
from peritheos.eosmat import validate_eosmat_document
from peritheos.materials import Material
from scripts.reproduce_sakai_2018_rhenium import (
    DATASET,
    RECORD,
    fit_points,
    ledger_outcome,
    load_points,
    refit_figure10,
    reproduce,
)


def test_sakai_2018_published_pressure_and_rounding_boundary():
    result = reproduce()
    assert result["catalog_matches_independent_expression"]
    assert result["reported_value_within_rounding_bounds"]
    assert result["catalog_pressure_gpa"] == pytest.approx(464.0, abs=3.0)
    assert result["catalog_minus_reported_gpa"] != pytest.approx(0, abs=1.0)


def test_sakai_2018_source_metadata_and_interchange():
    document = get_material_document("rhenium")
    validate_eosmat_document(document)
    source = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert document["formula_units_per_cell"] == 2
    assert source["eos"]["parameters"] == {
        "V0": 29.47,
        "K0": 358.0,
        "K0_prime": 4.8,
    }
    assert source["fixed_parameters"] == ["V0"]
    assert source["parameter_errors"] == {"K0": 10.0, "K0_prime": 0.2}
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert (
        source["scientific_validation"]["primary_data_check"]["status"] == "plot_only"
    )
    method = source["pressure_calibration"]["methods"][0]
    assert method["reference"]["doi"] == "10.1103/PhysRevB.80.104114"
    assert "reference_eos_record" not in method
    restored = Material.from_eosmat(document).to_eosmat()
    record = next(r for r in restored["eos_records"] if r["identifier"] == RECORD)
    assert record["scientific_validation"] == source["scientific_validation"]
    assert record["pressure_calibration"] == source["pressure_calibration"]


def test_sakai_2018_reference_inversion_and_extrapolation():
    record = get_eos_record(RECORD)
    assert record.pressure(29.47) == pytest.approx(0, abs=1e-10)
    pressures = np.array([130.0, 200.0, 280.0, 310.0])
    volumes = record.volume(pressures)
    np.testing.assert_allclose(record.pressure(volumes), pressures, rtol=1e-10)
    assert record.validity.pressure_gpa == (130.0, 310.0)
    assert record.pressure(18.65) > record.validity.pressure_gpa[1]
    with pytest.raises(ValueError, match="validity"):
        record.pressure(18.65, check_validity=True)


def test_vector_transcription_preserves_series_and_overlapping_markers():
    rows = load_points()
    assert len(rows) == 58
    assert len({r["pdf_curve_index"] for r in rows}) == 58
    for series, count in (
        ("rp01_dewaele_pt", 26),
        ("rp01_yokoo_pt", 26),
        ("micro17_yokoo_pt", 6),
    ):
        assert sum(r["series"] == series for r in rows) == count
    first = next(r for r in rows if r["series"] == "rp01_yokoo_pt")
    assert float(first["pressure_gpa"]) == pytest.approx(138.52871095)
    assert float(first["volume_a3"]) == pytest.approx(23.31437662)
    document = get_material_document("rhenium")
    source = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert source["fit_datasets"] == [DATASET]


def test_independent_fit_recovers_published_coefficients_and_selection_sensitivity():
    fits = refit_figure10()
    fit = fits["rp01_yokoo_pt"]
    assert fit["observations"] == 26
    assert fit["parameters"]["K0"] == pytest.approx(358, abs=10)
    assert fit["parameters"]["K0_prime"] == pytest.approx(4.8, abs=0.2)
    assert fit["pressure_rmse_gpa"] == pytest.approx(3.1061686, abs=1e-5)
    # A remote starting point must reach the same optimum.
    rows = [r for r in load_points() if r["series"] == "rp01_yokoo_pt"]
    repeat = fit_points(
        [float(r["pressure_gpa"]) for r in rows],
        [float(r["volume_a3"]) for r in rows],
        initial=(600, 2),
    )
    assert repeat["parameters"] == pytest.approx(fit["parameters"], rel=1e-6)
    assert not fits["all_yokoo_pt"]["within_published_quoted_errors"]
    assert not fits["rp01_yokoo_pt_volume_objective"]["within_published_quoted_errors"]
    outcome = ledger_outcome({})
    assert outcome["status"] == "parity"
    assert outcome["parity_basis"] == "within_reported_parameter_errors"
    assert all(p["within_reported_error"] for p in outcome["parameters"])
    assert all(p["within_combined_2sigma"] is None for p in outcome["parameters"])
