"""Primary Fei Pt thermal equation and recovered joint-fit scope."""

import hashlib
import json

import numpy as np
import pytest

from peritheos import (
    Material,
    get_eos_record,
    get_material_document,
    search_eos_records,
)
from peritheos.materials import PT_FEI_2007
from scripts.reproduce_fei_2007_platinum import (
    DATA,
    OUTPUT,
    RECORD,
    ledger_outcome,
    observations,
    reproduce,
    source_pressure,
)


def test_canonical_thermal_scale_preserves_cold_projection_and_legacy_pressures():
    record = get_eos_record(RECORD)
    cold = get_eos_record("platinum_fei_2007_vinet_300k")
    assert get_eos_record("pt_fcc_fei_2007") is record
    assert record in search_eos_records(formula="Pt", thermal=True)
    volume = np.array([60.38, 55.0, 50.0])[:, None]
    temperature = np.array([300.0, 1000.0, 1873.0])[None, :]
    assert record.pressure(volume, 300.0) == pytest.approx(cold.pressure(volume))
    expected = source_pressure(volume, temperature)
    assert record.pressure(volume, temperature) == pytest.approx(expected, abs=1e-9)
    assert record.pressure(volume, temperature) == pytest.approx(
        PT_FEI_2007.pressure(volume, temperature), abs=1e-9
    )
    restored = Material.from_eosmat(
        get_material_document("platinum"), record_identifiers=[RECORD]
    ).eos_records[0]
    assert restored.pressure(volume, temperature) == pytest.approx(expected, abs=1e-9)
    assert restored.volume(expected, temperature) == pytest.approx(
        np.broadcast_to(volume, (3, 3)), rel=1e-8
    )
    assert restored.eos.temperature(
        expected, volume * restored.volume_scale
    ) == pytest.approx(np.broadcast_to(temperature, (3, 3)), rel=1e-8)


def test_au_rereduction_uses_paired_observations_and_keeps_original_pressure():
    volume, temperature, pressure, paired, reduced = observations()
    assert len(volume) == len(temperature) == len(pressure) == 78
    assert len(paired) == len(reduced) == 42
    assert pressure[36:] == pytest.approx(reduced)
    au_volume = np.array([float(row["gold_volume_a3"]) for row in paired])
    assert reduced == pytest.approx(
        get_eos_record("gold_fei_2007_vinet_2").pressure(au_volume, temperature[36:]),
        abs=1e-9,
    )
    original = np.array([float(row["pressure_gpa"]) for row in paired])
    assert np.max(abs(reduced - original)) > 0.5
    # Pt pressures do not define their own fitting target.
    assert np.max(abs(source_pressure(volume[36:], temperature[36:]) - reduced)) > 0.5


def test_published_uncertainties_and_temperature_law_remain_explicit():
    raw = next(
        r
        for r in get_material_document("platinum")["eos_records"]
        if r["identifier"] == RECORD
    )
    assert raw["thermal"]["debye_temperature_law"] == "variable_exponent"
    assert raw["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 230.0,
        "gamma0": 2.72,
        "q": 0.5,
        "n": 1.0,
    }
    assert raw["thermal"]["parameter_errors"]["gamma0"] == 0.03
    assert raw["thermal"]["parameter_errors"]["q"] == 0.5
    assert raw["parameter_error_confidence"] is None
    assert raw["parameter_covariance"] is None
    assert raw["experimental_temperature_range_k"] == [300.0, 1873.0]
    assert raw["pressure_calibration"]["methods"][1]["reference_eos_record"] == (
        "gold_fei_2007_vinet_2"
    )


def test_audit_reproduces_equations_without_claiming_joint_coefficient_parity():
    result = reproduce()
    saved = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert result["scope"] == saved["scope"]
    assert result["input_sha256"] == saved["input_sha256"]
    assert np.array(result["equation_grid"]["pressure_gpa"]) == pytest.approx(
        np.array(saved["equation_grid"]["pressure_gpa"]), abs=1e-9
    )
    for diagnostic in (
        "fixed_v0_joint_diagnostic",
        "free_v0_joint_sensitivity",
        "figure2_temperature_selection_sensitivity",
    ):
        assert result[diagnostic]["observations"] == saved[diagnostic]["observations"]
        assert result[diagnostic]["parameters"] == pytest.approx(
            saved[diagnostic]["parameters"], rel=1e-5, abs=1e-7
        )
        assert result[diagnostic]["rmse_gpa"] == pytest.approx(
            saved[diagnostic]["rmse_gpa"], rel=1e-7
        )
    for filename, expected in result["input_sha256"].items():
        assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == expected
    assert result["independent_equation_max_difference_gpa"] < 1e-9
    assert result["independent_au_reduction_max_difference_gpa"] < 1e-9
    fit = result["fixed_v0_joint_diagnostic"]
    assert fit["observations"] == 78
    assert fit["solver_success"]
    assert fit["rmse_gpa"] < fit["published_rmse_gpa"]
    assert fit["parameters"]["q"] < 0  # source prints +0.5(5)
    assert result["free_v0_joint_sensitivity"]["parameters"]["q"] < 0
    assert result["figure2_temperature_selection_sensitivity"]["observations"] == 71
    outcome = ledger_outcome({"identifier": RECORD})
    assert outcome["status"] == "parity_not_achieved"
    assert outcome["observations"] == 78
    assert all(p["within_combined_2sigma"] is None for p in outcome["parameters"])
