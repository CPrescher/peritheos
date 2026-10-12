"""Complete-pair coverage, convention ambiguity and diagnostic Au recalibration."""

import csv
import hashlib
import io
import json

import numpy as np
import pytest

from scripts.audit_gold_speziale_reduction import (
    DATA,
    REPORT,
    RESIDUALS,
    observations,
    reproduce,
)
from scripts.check_numerical_archive import check_csv
from scripts.reproduce_fei_2007_gold import check_saved


@pytest.fixture(scope="module")
def result():
    return reproduce()


def test_complete_coverage_deduplicates_shared_heating_cycle(result):
    report, table = result
    rows = list(csv.DictReader(io.StringIO(table)))
    assert len(rows) == report["tabulated_pairs"] == 50
    assert report["unique_pairs"] == 48
    assert report["unique_hot_pairs"] == 38
    assert report["unique_cold_pairs"] == 10
    assert [r["duplicate_of"] for r in rows if r["duplicate_of"]] == [
        "Hirose2006:11",
        "Hirose2006:12",
    ]
    for duplicate in report["duplicates"]:
        point = next(
            r
            for r in rows
            if r["source"] == duplicate["source"]
            and int(r["source_row"]) == duplicate["source_row"]
        )
        source, number = duplicate["duplicate_of"].split(":")
        original = next(
            r for r in rows if r["source"] == source and r["source_row"] == number
        )
        for key in (
            "temperature_k",
            "mgo_volume_a3",
            "au_volume_a3",
            "reported_pressure_gpa",
        ):
            assert point[key] == original[key]


def test_cold_pressure_cannot_identify_gamma_or_theta(result):
    _, table = result
    for row in csv.DictReader(io.StringIO(table)):
        if float(row["temperature_k"]) != 300:
            continue
        values = [float(v) for k, v in row.items() if k.endswith("_300K_pressure_gpa")]
        assert max(values) - min(values) < 1e-12
    assert result[0]["native_equation_max_difference_gpa"] < 1e-9


def test_all_hot_pairs_support_gamma_but_do_not_identify_theta(result):
    conventions = result[0]["conventions"]
    correct = conventions["eq11_integrated_300K"]["unique_hot"]
    integrated = conventions["local_power_integrated_300K"]["unique_hot"]
    power = conventions["local_power_power_300K"]["unique_hot"]
    assert correct["inside_lattice_pressure_rounding_envelope"] == 0
    assert integrated["rmse_gpa"] < 0.018
    assert integrated["max_abs_residual_gpa"] < 0.054
    assert power["rmse_gpa"] < integrated["rmse_gpa"]
    assert power["inside_lattice_pressure_rounding_envelope"] == 38
    assert (
        conventions["local_power_integrated_298K"]["unique_hot"][
            "inside_lattice_pressure_rounding_envelope"
        ]
        == 38
    )


def test_recalibration_is_systematic_and_changes_gold_thermal_fits(result):
    refits = result[0]["diagnostic_au_refits"]
    fit2 = refits["Hirose2008_fit2"]["fits"]
    calibrated = fit2["published_speziale"]
    assert calibrated["pressure_target_shift"]["max_residual_gpa"] < 0
    assert (
        calibrated["refitted_curve_shift_from_reported_refit"]["max_residual_gpa"] < 0
    )
    assert (
        calibrated["refitted_curve_shift_from_reported_refit"]["max_abs_residual_gpa"]
        > 1.8
    )
    assert (
        fit2["local_power_candidate"]["refitted_curve_shift_from_reported_refit"][
            "max_abs_residual_gpa"
        ]
        < 0.04
    )
    assert set(refits) == {
        "Hirose2008_300K",
        "Hirose2008_fit1",
        "Hirose2008_fit2",
    }  # Fei EOS work is paused pending the author's reply.
    assert all(f["solver_success"] for r in refits.values() for f in r["fits"].values())
    assert result[0]["published_coefficients_changed"] is False
    assert result[0]["historical_reducer_identified"] is False


def test_source_hashes_archive_and_parenthetical_widths_are_preserved(result):
    report, residuals = result
    for name, digest in report["source_hashes"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    assert (
        next(
            r
            for r in observations()
            if r["source"] == "Hirose2006" and r["temperature_k"] == 300
        )["reported_pressure_width_gpa"]
        is None
    )
    saved = json.loads(REPORT.read_text())
    check_saved(saved, report)
    for name, refit in report["diagnostic_au_refits"].items():
        for target, fit in refit["fits"].items():
            for parameter, value in fit["parameters"].items():
                expected = saved["diagnostic_au_refits"][name]["fits"][target][
                    "parameters"
                ][parameter]
                assert np.isclose(expected, value, rtol=1e-4, atol=1e-16)
    numerical = {f"{name}_pressure_gpa" for name in report["conventions"]}
    check_csv(RESIDUALS.read_text(), residuals, numerical, atol=1e-8)
