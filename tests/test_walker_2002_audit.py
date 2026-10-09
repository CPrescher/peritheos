"""Walker primary ancestry, source disagreement, and unavailable replay limits."""

import json

import numpy as np
import pytest

from peritheos import (
    Material,
    get_eos_record,
    get_material_document,
    resolve_dataset_pressure,
)
from peritheos.errors import DatasetError
from scripts.audit_walker_2002 import (
    OUTPUT,
    b1_pressure_replay,
    bm3,
    check_report,
    load_rows,
    reproduce,
    values,
)


@pytest.mark.parametrize(
    "material,identifier,dataset,parameters,beta,error",
    [
        (
            "kcl",
            "kcl_walker_2002_bm3_2",
            "kcl_walker_2002_table2_pvt",
            {"V0": 53.53, "K0": 23.7, "K0_prime": 4.4},
            0.00275,
            0.00009,
        ),
        (
            "kcl_b1",
            "kcl_b1_walker_2002_bm3_linear_thermal",
            "kcl_b1_walker_2002_table1_pvt",
            {"V0": 37.50 * 4e24 / 6.02214076e23, "K0": 17.7, "K0_prime": 5},
            0.00195,
            0.00005,
        ),
    ],
)
def test_both_primary_attributions_preserve_coefficients_and_block_replay(
    material, identifier, dataset, parameters, beta, error
):
    document = get_material_document(material)
    record = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    calibration = record["pressure_calibration"]
    assert calibration["status"] == "partially_resolved"
    assert calibration["recalculation"]["status"] == "reference_eos_not_bundled"
    assert len(calibration["methods"]) == 1
    method = calibration["methods"][0]
    assert method["material"] == "NaCl (B1)"
    assert method["reference"]["doi"] == "10.1029/JB091iB05p04949"
    assert "page 806" in method["source_location"]
    assert "reference_eos_record" not in method
    assert record["eos"]["parameters"] == parameters
    assert record["thermal"]["parameters"] == {"Tr": 296.15, "alpha_KT": beta}
    assert record["parameter_errors"] == {"V0": None, "K0": None, "K0_prime": None}
    assert record["thermal"]["parameter_errors"]["alpha_KT"] == error
    round_trip = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).to_eosmat()["eos_records"][0]
    assert round_trip["pressure_calibration"] == calibration
    with pytest.raises(DatasetError, match="adjusted reference-isotherm BE2"):
        resolve_dataset_pressure(document, identifier, dataset)


def test_actual_temperature_and_volume_inputs_never_claim_complete_replay():
    result = reproduce()
    replay = result["nacl_pressure_replay"]
    assert replay["reference_temperature_k"] == 298.15
    assert replay["status"] == "conditional_run_normalized_replay"
    for number, filename in (
        (1, "kcl-walker-2002-table1-pvt.csv"),
        (2, "kcl-walker-2002-table2-pvt.csv"),
    ):
        raw, digest = load_rows(filename)
        inputs = replay[f"table{number}_rows"]
        assert len(inputs) == len(raw) == (30 if number == 1 else 39)
        assert result["source_resources"][f"table{number}_sha256"] == digest
        for source, row in zip(raw, inputs):
            temp = float(source["temperature_celsius"])
            assert row["temperature_k"] == pytest.approx(temp + 273.15)
            assert row["nacl_conventional_cell_volume_a3"] == pytest.approx(
                float(source["nacl_lattice_a_angstrom"]) ** 3
            )
            assert row["birch_abstract_thermal_increment_gpa"] == pytest.approx(
                0.00286 * (temp - 25)
            )
            assert row["outside_birch_abstract_temperature_range"] == (
                temp < 25 or temp > 500
            )
            if number == 1 and source["row_kind"] == "derived_reference":
                assert row["independently_replayed_pressure_gpa"] is None
            elif number == 2:
                assert isinstance(row["independently_replayed_pressure_gpa"], float)
                assert row["replay_kind"] == "conditional_run_normalized"
            else:
                assert isinstance(row["independently_replayed_pressure_gpa"], float)
                assert row["assumed_normalization_temperature_k"] == 296.15
                assert row["anchor_measured_temperature_k"] == 309.15
    assert replay["table2_rows"][7]["temperature_k"] == 297.15
    assert replay["table2_calibrant_reference_anchors"][0]["nacl_file"] == "r34439"
    assert replay["table2_calibrant_reference_anchors"][1]["nacl_file"] == "r35101"
    assert (
        max(
            abs(row["conditional_run_normalized_difference_gpa"])
            for row in replay["table2_rows"]
        )
        < 0.011
    )


def test_b1_two_reference_hypothesis_and_measured_temperature_sensitivity():
    replay = reproduce()["nacl_pressure_replay"]
    summary = replay["table1_replay"]
    assert summary["status"] == "conditional_reference_temperature_hypothesis"
    anchors = summary["calibrant_reference_anchors"]
    assert [a["calibrant_file"] for a in anchors] == ["r57689", "r57693"]
    assert [a["nacl_lattice_a_angstrom"] for a in anchors] == [5.6479, 5.6473]
    for kind, count, reference, rmse, maximum in (
        ("sample_observation", 22, "r57689", 0.0006997991, 0.0014577870),
        ("calibrant_spot_check", 5, "r57693", 0.0006382892, 0.0010323969),
    ):
        group = summary["nonzero_pressure_groups"][kind]
        assert group["observations"] == count
        conditional = group["assumed_23_celsius_normalization"]
        assert conditional["rmse_gpa"] == pytest.approx(rmse, abs=1e-10)
        assert conditional["max_absolute_residual_gpa"] == pytest.approx(
            maximum, abs=1e-10
        )
        sensitivity = group["measured_36_celsius_normalization"]
        assert 0.037 < sensitivity["rmse_gpa"] < 0.040
        assert 0.042 < sensitivity["max_absolute_residual_gpa"] < 0.047
        for index in group["source_row_indices"]:
            row = replay["table1_rows"][index]
            assert row["calibrant_reference_file"] == reference
            assert row["reported_pressure_gpa"] > 0
            assert row["replay_kind"] == "conditional_b1_reference_temperature"


def test_b1_imposed_zeros_keep_the_hypothesis_disagreement_visible():
    result = reproduce()
    replay = result["nacl_pressure_replay"]
    summary = replay["table1_replay"]
    assert summary["imposed_zero_source_row_indices"] == [1, 2]
    assert summary["unreplayed_source_row_indices"] == [0]
    for index in summary["imposed_zero_source_row_indices"]:
        row = replay["table1_rows"][index]
        assert row["reported_pressure_gpa"] == 0
        assert row["temperature_k"] == 309.15
        assert row["independently_replayed_pressure_gpa"] == pytest.approx(0.03718)
        assert row["conditional_reference_temperature_difference_gpa"] == pytest.approx(
            0.03718
        )
        assert row["measured_anchor_temperature_pressure_gpa"] == pytest.approx(
            0, abs=1e-10
        )
        assert row["reported_zero_is_imposed"]
        assert not row["included_in_nonzero_pressure_comparison"]
        assert row["replay_kind"] == "imposed_ambient_zero_diagnostic"
    # Calibrant spot checks still do not enter the KCl fit.
    assert result["b1_fit"]["observations"] == 23
    assert 1 in result["b1_fit"]["source_row_indices"]
    assert 2 not in result["b1_fit"]["source_row_indices"]


def test_b1_reference_selection_rejects_missing_or_misclassified_source():
    rows, _ = load_rows("kcl-walker-2002-table1-pvt.csv")
    with pytest.raises(ValueError, match="Expected one B1 ambient anchor r57693"):
        b1_pressure_replay([row for row in rows if row["spectrum"] != "r57693"])
    changed = [dict(row) for row in rows]
    changed[3]["row_kind"] = "unknown"
    with pytest.raises(ValueError, match="Unsupported B1 pressure row"):
        b1_pressure_replay(changed)


def test_b1_source_objective_and_figure_table_reference_disagreement():
    result = reproduce()["b1_fit"]
    assert result["observations"] == 23
    fitted = result["joint_fitted_parameters"]
    direct = result["direct_K0_alpha0_solver"]
    assert direct["success"]
    assert direct["K0"] == pytest.approx(fitted["K0"], abs=1e-7)
    assert direct["alpha_KT"] == pytest.approx(fitted["alpha_KT"], abs=1e-10)
    assert result["normal_equation_max_absolute_gradient"] < 1e-8
    assert result["joint_refit"]["rmse_gpa"] < result["published"]["rmse_gpa"]
    figure = result["figure1_reference_volume_diagnostic"]
    assert figure["conventional_cell_V0_a3"] == pytest.approx(249.080860076)
    assert round(figure["K0"], 1) == 17.7
    assert abs(figure["alpha_KT"] - 0.00195) < 0.00005
    assert figure["rmse_gpa"] == pytest.approx(result["joint_refit"]["rmse_gpa"])
    assert figure["rmse_gpa"] < result["table_reference_volume_diagnostic"]["rmse_gpa"]
    assert get_eos_record(
        "kcl_b1_walker_2002_bm3_linear_thermal"
    ).reference_volume == pytest.approx(249.080860076)
    literal = result["literal_printed_BE1_diagnostic"]
    assert abs(literal["K0"] - 17.7) > abs(fitted["K0"] - 17.7)


def test_b2_independent_staging_and_published_native_evaluation():
    result = reproduce()["b2_staged_fit"]["fits"]
    eight = result["eight_23_24_celsius_reference_rows"]
    seven = result["seven_exact_23_celsius_rows_sensitivity"]
    assert len(eight["cold_source_row_indices"]) == 8
    assert len(seven["cold_source_row_indices"]) == 7
    assert eight["K0"] == pytest.approx(23.77539804)
    assert eight["K0_prime"] == pytest.approx(4.41587367)
    assert eight["alpha_KT"] == pytest.approx(0.002766245675)
    rows, _ = load_rows("kcl-walker-2002-table2-pvt.csv")
    v = values(rows, "b2_kcl_cell_volume_a3")
    t = values(rows, "temperature_celsius") + 273.15
    native = get_eos_record("kcl_walker_2002_bm3_2").pressure(v, t)
    independent = bm3(v, 53.53, 23.7, 4.4) + 0.00275 * (t - 296.15)
    np.testing.assert_allclose(native, independent, rtol=0, atol=1e-12)


def test_generated_report_is_current():
    check_report(reproduce(), json.loads(OUTPUT.read_text(encoding="utf-8")))
