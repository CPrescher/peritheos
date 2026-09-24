"""Equal-relative-weight objective and identifiability disclosure checks."""

import json

import numpy as np
import pytest
from scipy.optimize import brentq

from scripts.refit_argon_xiao_2025 import (
    NAMES,
    OUTPUT,
    PUBLISHED,
    fit_case,
    predictions,
    recovered_data,
    residuals,
    sensitivity,
    spectrum,
)
from scripts.reproduce_argon_xiao_2025 import properties


def test_equal_weights_cover_each_recovered_input_without_property_balancing():
    data = recovered_data()
    assert len(data["volume"]) == 22
    assert len(data["bulk"]) == 5
    observed = np.concatenate([data["volume"], data["bulk"]])
    predicted = predictions(PUBLISHED, data, True)
    residual = residuals(np.zeros(9), data, True)
    assert residual == pytest.approx((predicted - observed) / observed, abs=1e-14)
    assert np.sqrt(np.mean(residual[:22] ** 2)) * 100 == pytest.approx(0.09517622713)
    # Separate source-equation implementation verifies predicted zero-P KT.
    independent = []
    for temperature in data["bulk_t"]:
        volume = brentq(lambda v: properties(v, temperature)["P"], 20, 25)
        independent.append(properties(volume, temperature)["K"] / 1000)
    assert predicted[22:] == pytest.approx(independent, rel=1e-8)
    assert (
        data["bulk_source"]["kind"]
        == "published_extrapolated_isothermal_bulk_modulus_estimates"
    )


def test_nine_parameter_sensitivity_does_not_assert_unique_recovery():
    data = recovered_data()
    result = spectrum(sensitivity(np.zeros(9), data, False))
    assert len(result["singular_values"]) == 9
    assert result["condition_number"] > 1e8
    assert result["rank_at_relative_singular_value_threshold"]["1e-06"] < 9
    fit = fit_case(data, False, [3, 4], starts=1, max_nfev=40)
    assert fit["free_parameters"] == ["theta0_k", "gamma0"]
    assert set(fit["fixed_parameters"]) == set(NAMES) - {"theta0_k", "gamma0"}
    assert fit["best_metrics"]["volume"]["relative_rms_percent"] < 0.04
    assert fit["solutions"][0]["optimizer_success"]


def test_report_preserves_all_multistarts_and_separates_original_fit():
    report = json.loads(OUTPUT.read_text())
    assert report["source_objective_reproduced"] is False
    assert report["fixed_v00_cm3_mol"] == 22.555
    assert report["parameter_uncertainties"] is None
    assert "Gas-phase" in report["excluded_parameter"]
    data = recovered_data()
    for name, case in report["cases"].items():
        include_bulk = "plus" in name
        assert case["input_entries"] == (27 if include_bulk else 22)
        for fit_name in ["full_nine_parameter_fit", "two_parameter_sensitivity_fit"]:
            fit = case[fit_name]
            assert len(fit["solutions"]) == fit["starts"] == 5
            selected = next(
                solution
                for solution in fit["solutions"]
                if solution["start_index"] == fit["best_start_index"]
            )
            assert fit["best_optimizer_success"] == selected["optimizer_success"]
            diagnostics = fit["best_thermodynamic_diagnostics"]
            assert diagnostics["passes_sampled_positive_cv_and_bulk"] == (
                diagnostics["minimum_cv_j_mol_k"] > 0
                and diagnostics["minimum_bulk_gpa"] > 0
            )
            if include_bulk and fit_name == "full_nine_parameter_fit":
                assert not diagnostics["passes_sampled_positive_cv_and_bulk"]
                assert diagnostics["negative_cv_temperatures_k"]

            assert fit["converged_starts"] == sum(
                solution["optimizer_success"] for solution in fit["solutions"]
            )

            assert len(fit["free_parameters"]) == (
                9 if fit_name.startswith("full") else 2
            )
            for solution in fit["solutions"]:
                parameters = np.array([solution["parameters"][name] for name in NAMES])
                residual = residuals(np.log(parameters / PUBLISHED), data, include_bulk)
                assert float(residual @ residual) == pytest.approx(
                    solution["metrics"]["all"]["sum_squared_relative_residuals"],
                    rel=1e-7,
                    abs=1e-12,
                )
            assert (
                fit["best_metrics"]["all"]["sum_squared_relative_residuals"]
                <= case["published_metrics"]["all"]["sum_squared_relative_residuals"]
            )
