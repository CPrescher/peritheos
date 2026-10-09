"""Independent Ne equation/normalization and qualified recovered-input scope."""

import hashlib
import json

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.datasets import Dataset
from peritheos.eosmat import validate_eosmat_document
from scripts.reproduce_fei_2007_neon import (
    DATA,
    OUTPUT,
    PARAMETERS,
    RECORDS,
    SOURCE,
    observations,
    reproduce,
    source_pressure,
)


@pytest.mark.parametrize("branch", ["vinet", "bm3"])
def test_ne_thermal_branches_match_independent_quadrature_and_invert(branch):
    document = get_material_document("neon_fcc")
    validate_eosmat_document(document)
    material = Material.from_eosmat(document)
    record = material.get_eos_record(RECORDS[branch])
    assert material.default_record().identifier == RECORDS["vinet"]
    volumes = np.array([25.0, 31.0, 39.0, 54.0])[:, None]
    temperatures = np.array([293.0, 300.0, 1000.0, 2000.0])[None, :]
    calculated = source_pressure(volumes, temperatures, branch)
    assert record.pressure(volumes, temperatures) == pytest.approx(calculated, abs=1e-9)
    assert record.pressure(volumes, 300) == pytest.approx(
        source_pressure(volumes, 300, branch), abs=1e-9
    )
    assert record.volume(calculated, temperatures) == pytest.approx(
        np.broadcast_to(volumes, (4, 4)), rel=1e-8
    )
    # The four-atom cell and molar-energy basis must not silently differ by 4.
    thermal = record.pressure(39, 1000) - record.pressure(39, 300)
    assert thermal == pytest.approx(3.689432553826, abs=1e-8)
    raw = next(r for r in document["eos_records"] if r["identifier"] == RECORDS[branch])
    assert raw["eos"]["parameters"] == dict(
        zip(("V0", "K0", "K0_prime"), PARAMETERS[branch][:3])
    )
    assert raw["thermal"]["parameter_errors"]["q"] == 0.3
    assert raw["thermal"]["debye_temperature_law"] == "variable_exponent"
    assert raw["thermal"]["thermal_pressure_reference"] == "reference_temperature"
    assert raw["parameter_covariance"] is None
    assert raw["parameter_error_confidence"] is None


def test_source_rows_scales_and_uncertainty_types_remain_distinct():
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    for name, meta in source["files"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == meta["sha256"]
    v, t, p, sp, sv = observations()
    assert len(v) == len(t) == len(p) == 54
    assert len(sp) == len(sv) == 35  # Experimental errors exist for prior tables only.
    assert p[0] == 4.83  # Original Finger ruby pressure retained.
    assert p[14] == pytest.approx(10.047704716)  # Explicitly recalibrated Hemley.
    assert p[48] == pytest.approx(23.583191758)  # Digitized Fei Pt, not Ne inversion.
    assert t[48:].tolist() == [1000] * 6
    assert source["finger_1981"]["temperature_k"] == 293
    assert source["finger_1981"]["coexistence_observation"]["volume"] is None
    assert source["figure5_hot_digitization"]["bound_note"].startswith("Illustrative")
    metadata = get_material_document("neon_fcc")["datasets"][-2:]
    prior = Dataset.from_mapping(metadata[0]).as_pressure_volume()
    assert prior.pressure[0] == 4.83
    assert prior.pressure_uncertainty[0] == 0.05
    hot = Dataset.from_mapping(metadata[1]).as_pressure_volume()
    assert np.all(np.isnan(hot.pressure_uncertainty))
    assert np.all(np.isnan(hot.volume_uncertainty))
    material = Material.from_eosmat(get_material_document("neon_fcc"))
    assert (
        max(
            abs(
                material.get_eos_record(RECORDS["vinet"]).pressure(v[48:], 1000)
                - p[48:]
            )
        )
        > 1
    )


def test_fit_parity_and_curve_validation_preserve_exact_regression_limits():
    report = reproduce()
    saved = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert report["input_sha256"] == saved["input_sha256"]
    assert report["status"]["fit_reproduction"] == "parity"
    assert report["status"]["reproduction_outcome"] == "reproduced"
    assert report["status"]["exact_original_regression"] == "not_established"
    assert report["status"]["absolute_pressure_accuracy"] == "not_established"
    assert not report["scope"]["complete_author_selection"]
    assert not report["scope"]["paired_new_fei_calibrant_volumes_recovered"]
    assert report["hemley_ruby_recalculation_max_rounding_difference_gpa"] < 1e-8
    assert report["normalization"]["v0_molar_cm3_mol"] == pytest.approx(13.394294924873)
    for branch, result in report["results"].items():
        assert result["status"] == "parity"
        assert all(
            c["within_published_error_width"]
            and c["within_combined_2sigma"]
            and c["similar"]
            for c in result["parameter_comparisons"]
        )
        assert result["independent_equation_max_difference_gpa"] < 1e-9
        assert all(
            w < 1 for w in result["cold_individual_parameter_error_widths"].values()
        )
        assert result["hot_6_q_only"]["parameters"]["q"] == pytest.approx(0.6, abs=0.02)
        for key in (
            "cold_48_unweighted",
            "hot_6_q_only",
            "joint_54_unweighted",
            "source_tables_35_effective_error_weighted",
        ):
            assert result[key]["parameters"] == pytest.approx(
                saved["results"][branch][key]["parameters"], abs=1e-6
            )
            assert result[key]["solver_success"]
        assert result["cold_48_unweighted"]["observations"] == 48
        assert result["joint_54_unweighted"]["observations"] == 54
        assert (
            abs(
                result["source_tables_35_effective_error_weighted"]["parameters"][
                    "K0_prime"
                ]
                - result["source_tables_35_unweighted"]["parameters"]["K0_prime"]
            )
            > 0.5
        )
    # These 18 points were sampled on calculated source curves, not passed to
    # observations(), which has 54 separate measured/plot-digitized inputs.
    curve = report["curve_check"]
    assert len(curve["points"]) == 18
    assert all(s["samples"] == 6 for s in curve["summary"]["vinet"]["printed"].values())
    assert max(abs(row["difference_gpa"]) for row in curve["thermal_separation"]) < 0.16


@pytest.mark.parametrize("branch", ["vinet", "bm3"])
def test_refit_ledger_adapter_preserves_joint_thermal_scope(branch):
    from scripts.reproduce_fei_2007_neon import ledger_outcome

    record = next(
        r
        for r in get_material_document("neon_fcc")["eos_records"]
        if r["identifier"] == RECORDS[branch]
    )
    outcome = ledger_outcome(record)
    assert outcome["status"] == "parity"
    assert outcome["fit_kind"] == "diagnostic_joint_pvt"
    assert outcome["observations"] == 54
    assert set(outcome["free_parameters"]) == {"K0", "K0_prime", "q"}
    assert len(outcome["parameters"]) == 3
    assert outcome["solver_success"]
    assert "confidence" in outcome["parity_basis"]
