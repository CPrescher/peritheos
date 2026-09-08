"""Source-data and reconstruction tests for the diamond thermal composites."""

from peritheos import get_material_document
from scripts.reproduce_diamond_thermal_composites import reproduce


def test_composite_cards_distinguish_derivation_and_h05_calibration():
    document = get_material_document("diamond")
    records = {record["identifier"]: record for record in document["eos_records"]}
    datasets = {dataset["identifier"]: dataset for dataset in document["datasets"]}

    for identifier in (
        "diamond_dewaele_2008_vinet_2",
        "diamond_correa_2008_dewaele_anchored",
        "diamond_benedict_2014_dewaele_anchored",
    ):
        method = records[identifier]["pressure_calibration"]["methods"][0]
        assert method["reference_calibration_record"] == "ruby_holzapfel_2005"
        assert method["reference"]["year"] == 2005

    for identifier in (
        "diamond_correa_2008_dewaele_anchored",
        "diamond_benedict_2014_dewaele_anchored",
    ):
        record = records[identifier]
        assert record["record_kind"] == "derived"
        assert record["derived_from_record"] == "diamond_dewaele_2008_vinet_2"
        assert "No coefficients are optimized" in record["derivation"]["method"]

    assert datasets["diamond_correa_2008_figure8_dft_md_vector_digitized"][
        "digitization"
    ]["method"].startswith("Marker centers extracted")
    assert datasets["diamond_benedict_2014_supplement_solid_dft_md"][
        "transcription"
    ]["arxiv_identifier"] == "1311.4577"


def test_correa_figure_8_validates_published_complete_pressure_model():
    result = reproduce()["correa_2008"]
    validation = result["source_model_pressure_validation"]

    assert validation["observations"] == 57
    assert validation["isochores"] == 3
    assert validation["absolute_pressure_rmse_gpa"] < 2.0
    assert validation["absolute_pressure_max_abs_residual_gpa"] < 7.0


def test_benedict_supplement_validates_pressure_and_caloric_model():
    result = reproduce()["benedict_2014"]
    pressure = result["source_model_pressure_validation"]
    energy = result["source_model_caloric_validation"]

    assert pressure["observations"] == 96
    assert pressure["isochores"] == 14
    assert pressure["absolute_pressure_rmse_gpa"] < 3.3
    assert pressure["absolute_pressure_max_abs_residual_gpa"] < 7.0
    assert energy["internal_energy_increment_comparisons"] == 82
    assert energy["internal_energy_increment_rmse_ev_per_atom"] < 0.062


def test_complete_composites_exactly_preserve_anchor_and_thermal_increments():
    result = reproduce()

    assert result["coefficient_optimization_performed"] is False
    for source in ("correa_2008", "benedict_2014"):
        identity = result[source]["complete_composite_identity"]
        assert identity["max_abs_composition_pressure_error_gpa"] < 2.0e-12
        assert identity["max_abs_reference_isotherm_error_gpa"] == 0.0
        assert identity["max_abs_thermal_energy_increment_error_j_per_mol"] < 3.0e-8
