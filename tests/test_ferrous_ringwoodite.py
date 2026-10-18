import hashlib
import json
from collections import Counter

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, list_material_documents, load_eosmat, save_eosmat
from peritheos.eosmat import eosmat_schema
from peritheos.errors import MaterialError
from scripts.ferrous_ringwoodite_recover import (
    DATA,
    NORMALIZED,
    PREFIX,
    RAW,
    check_source_pdfs,
    normalized_csv,
    read_csv,
)
from scripts.ferrous_ringwoodite_reproduce import (
    CARD,
    CELL_TO_CM3_MOL,
    STAGE,
    STATIC_ID,
    THERMAL_ID,
    V0_CELL,
    V0_MOLAR,
    debye_energy,
    debye_energy_series,
    pressure,
    reproduce,
    theta,
)


@pytest.fixture(scope="module")
def result():
    return reproduce()


def candidate():
    return Material.from_eosmat(load_eosmat(CARD), require_primary_validation=False)


def test_staged_card_schema_qualification_and_lossless_roundtrip(tmp_path):
    document = load_eosmat(CARD)
    Draft202012Validator(eosmat_schema()).validate(document)
    assert document["identifier"] == "ferrous_ringwoodite"
    assert "ferrous_ringwoodite" not in list_material_documents()
    record = next(r for r in document["eos_records"] if r["identifier"] == THERMAL_ID)
    assert record["record_kind"] == "published"
    assert record["scientific_validation"]["status"] == "not_reproduced"
    assert record["catalog_access"] == "explicit_selection"
    assert record["default"] is False
    assert record["eos"]["parameters"]["K0"] == 202
    assert record["parameter_errors"] == {"V0": None, "K0": 4.0}
    assert record["thermal"]["parameter_errors"]["gamma0"] == 0.06
    assert record["thermal"]["parameters"]["n"] == 7
    assert record["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert (
        record["scientific_validation"]["numerical_reproduction"][
            "author_fit_reproduced"
        ]
        is False
    )
    validated = Material.from_eosmat(document)
    assert validated.default_record().identifier == STATIC_ID
    assert len(validated.eos_records) == 1
    with pytest.raises(MaterialError, match="not_reproduced"):
        Material.from_eosmat(document, record_identifiers=[record["identifier"]])
    path = tmp_path / "ferrous_ringwoodite.eosmat"
    save_eosmat(path, document)
    assert load_eosmat(path) == document
    material = candidate()
    restored = Material.from_eosmat(
        material.to_eosmat(), require_primary_validation=False
    )
    assert restored.to_eosmat() == material.to_eosmat()


def test_distinct_fe_structure_origin_cell_contents_and_reference_volume():
    document = load_eosmat(CARD)
    assert document["formula"] == "Fe2SiO4"
    assert document["formula_units_per_cell"] == 8
    assert document["source"]["origin_choice"] == 2
    assert document["space_group"] == "Fd-3m"
    assert document["space_group_number"] == 227
    contents = Counter()
    for site in document["atom_sites"]:
        contents[site["element"]] += int(site["wyckoff"][:-1]) * site["occupancy"]
    assert contents == {"Fe": 16, "Si": 8, "O": 32}
    assert document["atom_sites"][0]["x"] == 0.5
    assert document["atom_sites"][1]["x"] == 0.125
    assert document["atom_sites"][2]["x"] == 0.2409
    assert document["lattice"]["a"] ** 3 == pytest.approx(558.661848)
    thermal = next(r for r in document["eos_records"] if r["identifier"] == THERMAL_ID)
    assert thermal["eos"]["parameters"]["V0"] == pytest.approx(V0_CELL)
    assert V0_CELL * CELL_TO_CM3_MOL == pytest.approx(42.03)
    assert V0_CELL == pytest.approx(558.3396559465342)
    # Separate specimens and sources legitimately have different reference cells.
    assert document["lattice"]["a"] ** 3 != pytest.approx(V0_CELL)


def test_all_resources_and_primary_transcriptions_preserve_source_values():
    document = load_eosmat(CARD)
    for dataset in document["datasets"]:
        path = DATA.parent / dataset["resource"]["path"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
        rows = read_csv(path)
        assert len(rows) == dataset["row_count"]
        assert list(rows[0]) == [col["name"] for col in dataset["columns"]]
    raw = read_csv(RAW)
    assert len(raw) == 18
    assert raw[0]["measurement"] == "Fay008heat"
    assert raw[0]["ringwoodite_a_angstrom"] == "8.1573"
    assert raw[0]["ringwoodite_a_parentheses"] == "30"
    assert raw[0]["upstream_temperature_k"] == "300"
    assert raw[0]["downstream_temperature_k"] == "1292"
    assert raw[0]["source_density_column"] == "43.22"
    assert raw[0]["exclusion_footnote"] == "1"
    assert raw[-1]["ringwoodite_a_angstrom"] == "8.0785"
    assert raw[-1]["exclusion_footnote"] == "2"
    rows = read_csv(NORMALIZED)
    assert sum(int(row["included_by_table_s1"]) for row in rows) == 14
    assert rows[1]["ringwoodite_a_sd_angstrom"] == "0.0040"
    assert rows[1]["pressure_reported_uncertainty_gpa"] == "0.30"
    assert rows[1]["assumed_mean_sample_temperature_k"] == "1148.00"
    assert rows[1]["assumed_nacl_temperature_k"] == "724.00"
    assert NORMALIZED.read_text() == normalized_csv()
    liu = read_csv(DATA / f"{PREFIX}liu-2008-table1.csv")
    assert len(liu) == 29
    assert liu[0]["pressure_gpa"] == "3.73"
    assert liu[0]["density_g_cm3"] == "4.932"
    assert liu[-1]["temperature_k"] == "1073"
    assert liu[-1]["g_gpa"] == "69.4"
    suzuki = read_csv(DATA / f"{PREFIX}suzuki-1979-table2-calculated.csv")
    assert len(suzuki) == 10
    assert suzuki[-1]["volume_expansion_percent"] == "5.80"
    assert suzuki[-1]["volume_expansion_95ci_halfwidth_percent"] == "0.90"
    assert suzuki[-1]["volume_extrapolated"] == "1"


def test_optional_complete_source_pdf_comparison():
    sources = STAGE / "sources"
    if not (sources / "armentrout-2011-table-s1.pdf").is_file():
        pytest.skip(
            "Local publisher/author PDFs are audit evidence, not redistributable test fixtures"
        )
    result = check_source_pdfs(sources)
    assert result["table_s1_rows"] == 18
    assert result["liu_table1_rows"] == 29
    assert result["table_s1_all_fields_match"]
    assert result["liu_table1_all_cells_match"]
    if (sources / "stixrude-2005.pdf").is_file():
        assert result["stixrude_table_a1_all_values_errors_and_systematics_match"]
        assert result["stixrude_table_a1_parameters"] == 9


def test_independent_source_curve_checkpoints_and_law_ambiguity(result):
    assert result["source_equation_assessment"]["figure1_all_checkpoints_pass"]
    assert len(result["figure1_checkpoints"]) == 12
    for point in result["figure1_checkpoints"]:
        assert abs(point["difference"]) < point["tolerance"]
        # The source curve resolution cannot establish which theta law was used.
        assert (
            abs(
                point["variable_exponent_volume_ratio"]
                - point["calculated_volume_ratio"]
            )
            < point["tolerance"]
        )
    assert (
        result["source_equation_assessment"]["theta_law_uniquely_identified_by_figure"]
        is False
    )
    assert result["source_static_tradeoff"]["K0_at_Kprime_4_gpa"] == 201.5
    assert result["source_static_tradeoff"]["stored_preferred_K0_gpa"] == 202
    assert len(result["figure_s1_temperature_checkpoints"]) == 9
    assert all(
        r["rounded_mean_matches"] for r in result["figure_s1_temperature_checkpoints"]
    )
    hot = next(
        r
        for r in result["figure_s1_temperature_checkpoints"]
        if r["measurement"] == "Fay025heat"
    )
    assert hot["table_s1_arithmetic_mean_k"] == hot["figure_s1_temperature_k"] == 2515


def test_conditional_fits_do_not_claim_original_combined_fit(result):
    fit = result["conditional_gamma_refits"]
    assert fit["armentrout_2011_accepted_14"]["gamma0"] == pytest.approx(
        0.9402573001, abs=1e-7
    )
    assert fit["liu_2008_all_29"]["gamma0"] == pytest.approx(1.3040896821, abs=1e-7)
    assert fit["armentrout_14_pressure_error_weighted"]["gamma0"] == pytest.approx(
        0.9846885649, abs=1e-7
    )
    assert fit["armentrout_14_plus_liu_29"]["gamma0"] == pytest.approx(
        0.9825907351, abs=1e-7
    )
    supplemented = fit["armentrout_14_liu_29_suzuki_7_including_extrapolations"]
    assert supplemented["gamma0"] == pytest.approx(1.0344027677, abs=1e-7)
    assert supplemented["within_published_0_06_error_width"] is True
    assert (
        result["source_equation_assessment"]["exact_author_workflow_reproduced"]
        is False
    )
    assert (
        result["inconsistent_printed_density_column"][
            "maximum_difference_from_molar_volume_cm3_mol"
        ]
        > 4
    )
    assert result["primary_data"]["excluded_measurements"] == [
        "Fay008heat",
        "Fay014quench",
        "Fay019preheat",
        "Fay029quench",
    ]


def test_independent_quadrature_native_parity_and_public_cell_scaling(result):
    record = candidate().get_eos_record(THERMAL_ID)
    assert record.reference_volume == pytest.approx(V0_CELL)
    assert record.pressure(V0_CELL, 300) == pytest.approx(0, abs=1e-12)
    assert pressure(V0_MOLAR, 300) == pytest.approx(0, abs=1e-12)
    assert theta(V0_MOLAR, 1.08, 2, "integrated_gruneisen") == 685
    for point in result["independent_native_checkpoints"]:
        assert point["native_volume_a3"] == pytest.approx(point["volume_a3"], abs=3e-8)
        assert point["native_pressure_gpa"] == pytest.approx(
            point["pressure_gpa"], abs=3e-8
        )
    pressures = np.array([1.0, 10.0, 24.0])
    temperatures = np.array([300.0, 1500.0, 2000.0])
    volumes = record.volume(pressures, temperatures, check_validity=True)
    np.testing.assert_allclose(
        record.pressure(volumes, temperatures), pressures, atol=1e-8
    )
    for p, t in [(25, 2000), (10, 2515), (-1, 300)]:
        with pytest.raises(MaterialError, match="outside"):
            record.volume(p, t, check_validity=True)


def test_manifest_license_and_calibrant_scope():
    manifest = json.loads((DATA / "ferrous-ringwoodite-source.json").read_text())
    for filename, resource in manifest["datasets"].items():
        assert (
            hashlib.sha256((DATA / filename).read_bytes()).hexdigest()
            == resource["sha256"]
        )
    assert manifest["catalog_installed"] is False
    assert "contributor" in manifest["contributor_rights"].lower()
    record = next(
        r for r in load_eosmat(CARD)["eos_records"] if r["identifier"] == THERMAL_ID
    )
    calibration = record["pressure_calibration"]
    assert calibration["status"] == "partially_resolved"
    assert calibration["methods"][0]["reference"]["doi"] == "10.1063/1.371596"
    assert calibration["methods"][1]["kind"] == "ultrasonic"
    assert calibration["recalculation"]["status"] == "reference_eos_not_bundled"


def test_verified_static_reconstruction_preserves_published_bulk_modulus(result):
    document = load_eosmat(CARD)
    record = Material.from_eosmat(document).default_record()
    raw = next(r for r in document["eos_records"] if r["identifier"] == STATIC_ID)
    assert raw["record_kind"] == "derived"
    assert raw["scientific_validation"]["status"] == "primary_source_validated"
    assert raw["derivation"]["source_kind"] == "published_table"
    assert raw["eos"]["parameters"] == {"V0": 558.7, "K0": 196.0}
    assert raw["parameter_errors"] == {"V0": 0.2, "K0": 6.0}
    assert record.reference_temperature == 296.15
    assert record.pressure(558.7, 296.15) == 0
    static = result["static_1979_reconstruction"]
    assert static["conditional_unweighted_compressed_two_row_K0_gpa"] == pytest.approx(
        196.2175873043
    )
    assert static["author_fit_reproduced"] is False
    assert all(p["within_combined_coordinate_error"] for p in static["source_checks"])
    assert record.pressure(553.5, 296.15) == pytest.approx(1.867379040453365)
    assert record.pressure(547.7, 296.15) == pytest.approx(4.0556076738036)
    for p in [0.0, 2.0, 4.0]:
        volume = record.volume(p, 296.15, check_validity=True)
        assert record.pressure(volume, 296.15) == pytest.approx(p, abs=1e-8)
    with pytest.raises(MaterialError, match="isothermal"):
        record.pressure(553.5, 300)


def test_constant_q_theta_law_satisfies_thermodynamic_definition(result):
    audit = result["thermal_equation_followup"]
    for row in audit["identity_checks"]:
        assert row["integrated_gruneisen"] == pytest.approx(row["gamma"], abs=1e-8)
        expected_naive = row["gamma"] * (1 + 2 * np.log(row["volume_ratio"]))
        assert row["variable_exponent"] == pytest.approx(expected_naive, abs=1e-8)
    compressed = audit["identity_checks"][0]
    assert abs(compressed["variable_exponent"] - compressed["gamma"]) > 0.3


def test_author_bernoulli_series_and_published_residuals(result):
    # The author's series independently checks the normalized Debye integral,
    # despite the inverse prefactor printed in the adjacent prose definition.
    for temperature in (300, 700, 2000):
        assert debye_energy_series(temperature, 800, 20) == pytest.approx(
            debye_energy(temperature, 800), rel=1e-12
        )
    assert debye_energy(1e7, 685) / (3 * 7 * 8.31446261815324 * 1e7) == pytest.approx(
        1, abs=3e-5
    )
    with pytest.raises(ValueError, match="2\\*pi"):
        debye_energy_series(100, 800, 20)
    checks = result["thermal_equation_followup"]
    assert checks["series_checks"][-1]["maximum_pressure_difference_gpa"] < 1e-10
    assert all(row["passes"] for row in checks["figure9_residual_checks"])
    assert checks["figure9_residual_checks"][0]["measured_minus_model_gpa"] < -2.5


def test_upstream_debye_adoption_preserves_separate_uncertainty():
    rows = read_csv(DATA / f"{PREFIX}stixrude-2005-table-a1.csv")
    assert len(rows) == 9
    by_parameter = {r["parameter"]: r for r in rows}
    assert by_parameter["theta0"]["value"] == "685"
    assert by_parameter["theta0"]["reported_error"] == "100"
    assert by_parameter["theta0"]["from_systematics"] == "1"
    assert by_parameter["K0"]["value"] == "218"
    assert by_parameter["gamma0"]["value"] == "1.30"
    thermal = next(
        r for r in load_eosmat(CARD)["eos_records"] if r["identifier"] == THERMAL_ID
    )
    assert thermal["thermal"]["parameters"]["theta0"] == 685
    assert thermal["thermal"]["parameter_errors"]["theta0"] is None
    assert "theta0" in thermal["thermal"]["fixed_parameters"]
