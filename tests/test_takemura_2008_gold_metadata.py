"""Takemura Au readiness and conditional parity use observed rows and ruby physics."""

import json

import numpy as np
import pytest

from peritheos import (
    Material,
    get_material_document,
    get_pressure_calibration,
    recalculate_ruby_pressure,
)
from scripts.reproduce_takemura_2008_gold import (
    OUTPUT,
    check_saved,
    ledger_outcome,
    pressure,
    reproduce,
    volume,
)

RECORD = "gold_takemura_2008_vinet_6"
DATASET = "gold_takemura_2008_table3_compression"


def _source():
    document = get_material_document("gold")
    record = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    material = Material.from_eosmat(document, record_identifiers=[RECORD])
    return record, material.get_dataset(DATASET)


def test_ready_takemura_record_has_complete_observations_on_its_pressure_coordinate():
    record, data = _source()
    assert record["pressure_calibration"]["recalculation"]["status"] == "ready"
    assert record["scientific_validation"]["primary_data_check"][
        "dataset_identifiers"
    ] == [DATASET]
    assert data.checksum_verified
    assert data.resource.sha256 == (
        "04e2015f65f5203d89966ebf4b08eb3000d9a3fb7bfc7e3fd032e253e72b3858"
    )
    assert len(data) == 86
    assert np.unique(data["run"], return_counts=True)[1].tolist() == [
        13,
        21,
        6,
        10,
        9,
        9,
        18,
    ]
    for column in data.columns:
        assert column.role in {"flag", "value"}
        assert np.all(np.isfinite(data[column.name]))

    zha = data["ruby_pressure_zha_gpa"]
    dor = data["ruby_pressure_dorogokupets_gpa"]
    highest = np.argmax(dor)
    # Table III's paired coordinates describe the same observation, not two runs.
    assert (zha[highest], dor[highest], data["lattice_111_angstrom"][highest]) == (
        122.5,
        130.8,
        3.6699,
    )
    assert record["experimental_pressure_range_gpa"] == [dor.min(), dor.max()]
    assert record["eos"]["parameters"] == {"V0": 67.848, "K0": 167.0, "K0_prime": 5.9}
    assert record["fixed_parameters"] == ["V0", "K0"]
    assert record["parameter_errors"] == {"V0": None, "K0": None, "K0_prime": 0.1}


def test_all_takemura_rows_can_be_reduced_with_independent_published_ruby_equations():
    record, data = _source()
    identifier = record["pressure_calibration"]["methods"][0][
        "reference_calibration_record"
    ]
    assert identifier == "ruby_dorogokupets_oganov_2007"
    calibration = get_pressure_calibration(identifier)
    dor = data["ruby_pressure_dorogokupets_gpa"]
    # Analytic positive root of Dorogokupets-Oganov (2007), Eq. (16).
    shift = (np.sqrt(1 + 4 * 5.5 * dor / 1884) - 1) / (2 * 5.5)
    assert calibration.wavelength_ratio(dor) == pytest.approx(1 + shift, abs=1e-14)
    assert calibration.pressure_from_ratio(1 + shift) == pytest.approx(dor, abs=1e-10)

    # Independent Mao (1986) equation evaluates the recovered effective ratio.
    mao = recalculate_ruby_pressure(dor, identifier, "ruby_mao_1986")
    expected_mao = 1904 / 7.665 * ((1 + shift) ** 7.665 - 1)
    assert mao == pytest.approx(expected_mao, abs=1e-10)
    assert recalculate_ruby_pressure(mao, "ruby_mao_1986", identifier) == (
        pytest.approx(dor, abs=1e-10)
    )


def test_takemura_independent_equation_and_inverse_match_library_on_observed_states():
    record, data = _source()
    material = Material.from_eosmat(
        get_material_document("gold"), record_identifiers=[RECORD]
    )
    eos = material.eos_records[0]
    observed_volume = data["lattice_111_angstrom"] ** 3
    derivative = record["eos"]["parameters"]["K0_prime"]
    assert pressure(observed_volume, derivative) == pytest.approx(
        eos.pressure(observed_volume, check_validity=False), abs=1e-10
    )
    p = data["ruby_pressure_dorogokupets_gpa"]
    assert volume(p, derivative) == pytest.approx(eos.volume(p), abs=1e-8)


def test_takemura_parameter_parity_is_stable_to_pressure_coordinate_and_objective():
    report = reproduce()
    expected = {
        "dorogokupets_oganov": (5.9287095121, 5.9292840677),
        "zha": (5.5133984339, 5.5301250915),
    }
    assert report["status"] == "parity"
    assert report["parity_basis"] == "within_reported_parameter_errors"
    assert report["observations"] == 86
    assert report["author_fit_reproduced"] is False
    for name, values in expected.items():
        branch = report["branches"][name]
        for objective, expected_derivative in zip(
            ("equal_pressure_residuals", "equal_volume_residuals"), values
        ):
            fit = branch["fits"][objective]
            assert fit["K0_prime"] == pytest.approx(expected_derivative, abs=2e-7)
            assert abs(fit["K0_prime"] - branch["published_K0_prime"]) < 0.1
            assert fit["refit_error"] is None
            assert fit["reduced_chi_square"] is None
    record, _ = _source()
    outcome = ledger_outcome(record)
    assert outcome["parameters"][0]["within_reported_error"] is True
    assert outcome["parameters"][0]["within_combined_2sigma"] is None
    assert outcome["absolute_sigma"] is False
    assert outcome["rmse_gpa"] == pytest.approx(0.7778728515, abs=1e-9)


def test_takemura_dedicated_report_is_reproducible():
    check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), reproduce())


def test_takemura_common_ledger_uses_the_qualified_dedicated_result():
    from scripts.validate_primary_eos_refits import (
        DEFAULT_JSON,
        DEFAULT_MARKDOWN,
        _fit_record,
        render_markdown,
    )

    document = get_material_document("gold")
    record, _ = _source()
    dataset = next(d for d in document["datasets"] if d["identifier"] == DATASET)
    current = _fit_record(document, record, dataset)
    ledger = json.loads(DEFAULT_JSON.read_text(encoding="utf-8"))
    saved = next(r for r in ledger["records"] if r["record_identifier"] == RECORD)
    check_saved({key: saved[key] for key in current}, current)
    actual_row = next(
        line
        for line in DEFAULT_MARKDOWN.read_text(encoding="utf-8").splitlines()
        if f"[`{RECORD}`]" in line
    )
    assert actual_row in render_markdown(ledger).splitlines()
