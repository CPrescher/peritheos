"""Guard source coverage, phase distinctions and numerical recovery checks."""

import csv
import hashlib
import io
from collections import Counter
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from scripts.digitize_shim_2002_figure2 import (
    ADOPTED_V0_A3,
    PRESSURE_AXIS,
    RATIO_AXIS,
    linear,
)


def source_rows(material, identifier):
    document = get_material_document(material)
    dataset = next(d for d in document["datasets"] if d["identifier"] == identifier)
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_bytes()
    )
    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    rows = list(csv.DictReader(io.StringIO(payload.decode())))
    assert list(rows[0]) == [c["name"] for c in dataset["columns"]]
    return document, dataset, rows


def test_full_wang_table_preserves_runs_footnotes_and_printed_volume_errors():
    _, dataset, rows = source_rows(
        "ca_perovskite", "ca_perovskite_wang_1996_table1_pvt"
    )
    assert len(rows) == 66
    assert Counter(r["run"] for r in rows) == {"13": 34, "3": 15, "4": 12, "5": 5}
    assert [int(r["source_order"]) for r in rows] == list(range(1, 67))
    assert sum(float(r["temperature_k"]) > 306 for r in rows) == 52
    assert [rows[i]["pressure_gpa"] for i in [0, 34, 49, 61, 64, 65]] == [
        "12.04",
        "12.00",
        "11.61",
        "12.16",
        "12.01",
        "11.69",
    ]
    assert rows[32]["pressure_footnote"] == "dagger"
    assert rows[33]["pressure_footnote"] == "double_dagger"
    assert all(r["source_origin"] == "wang_weidner_1994" for r in rows[34:])
    for row in rows:
        central, digits = row["volume_printed_a3"][:-1].split("(")
        assert row["volume_a3"] == central
        assert float(row["volume_reported_uncertainty_a3"]) == int(digits) * 0.01
    error = next(c for c in dataset["columns"] if "reported_uncertainty" in c["name"])
    assert error["role"] == "uncertainty"
    assert dataset["used_by_eos_records"] == []


def test_wang_room_temperature_subset_keeps_12_fit_rows_and_corrects_two_glyphs():
    _, _, full = source_rows("ca_perovskite", "ca_perovskite_wang_1996_table1_pvt")
    _, dataset, subset = source_rows(
        "ca_perovskite", "ca_perovskite_wang_1996_table1_room_temperature"
    )
    assert len(subset) == 14
    assert [r["pressure_gpa"] for r in subset[-2:]] == ["1.13", "0.59"]
    assert sum(int(r["fit_included"]) for r in subset) == 12
    for original, row in zip(full[20:34], subset):
        for field in (
            "pressure_gpa",
            "temperature_k",
            "volume_a3",
            "differential_stress_gpa",
        ):
            assert original[field] == row[field]
        assert original["volume_reported_uncertainty_a3"] == row["volume_sigma_a3"]
    assert dataset["provenance"]["corrected_rows"][0]["previous_pressure_gpa"] == 1.137
    assert (
        next(c for c in dataset["columns"] if c["name"] == "volume_sigma_a3")["role"]
        == "uncertainty"
    )


def test_shim_vector_conversion_preserves_six_unique_states_and_cubic_point():
    document, dataset, rows = source_rows(
        "casio3_perovskite_tetragonal",
        "casio3_perovskite_tetragonal_shim_2002_figure2a_digitized",
    )
    assert len(rows) == 6
    assert Counter(r["source_refinement"] for r in rows) == {
        "cubic": 1,
        "tetragonal": 5,
    }
    assert len({(r["plot_x_pt"], r["plot_y_pt"]) for r in rows}) == 6
    assert rows[0]["source_refinement"] == "cubic"
    assert [float(r["pressure_gpa"]) for r in rows] == [
        19.7,
        24.2,
        25.2,
        26.6,
        36.1,
        45.8,
    ]
    for row in rows:
        assert float(row["plot_pressure_gpa"]) == pytest.approx(
            linear(float(row["plot_x_pt"]), PRESSURE_AXIS), abs=1e-6
        )
        ratio = linear(float(row["plot_y_pt"]), RATIO_AXIS)
        assert float(row["volume_ratio"]) == pytest.approx(ratio, abs=1e-8)
        assert float(row["volume_a3"]) == pytest.approx(ADOPTED_V0_A3 * ratio, abs=1e-6)
    roles = {c["name"]: c["role"] for c in dataset["columns"]}
    assert roles["pressure_plot_sigma_gpa"] == "standard_deviation"
    assert roles["pressure_digitization_bound_gpa"] == "bound"
    assert roles["volume_digitization_bound_a3"] == "bound"
    pressure_bound = next(
        c for c in dataset["columns"] if c["name"] == "pressure_digitization_bound_gpa"
    )
    assert pressure_bound["of"] == "plot_pressure_gpa"
    record = document["eos_records"][0]
    assert record["eos"]["parameters"] == {"V0": 45.58, "K0": 255.0, "K0_prime": 4.0}
    assert record["parameter_errors"] == {"V0": None, "K0": 5.0, "K0_prime": None}
    assert (
        record["scientific_validation"]["primary_data_check"]["status"]
        == "plot_digitized"
    )
    assert "fit_datasets" not in record  # Original numerical inputs remain unavailable.


def test_shim_plot_errors_support_published_modulus_but_weights_matter():
    document, _, rows = source_rows(
        "casio3_perovskite_tetragonal",
        "casio3_perovskite_tetragonal_shim_2002_figure2a_digitized",
    )
    p, v, sigma_p, sigma_v = [
        np.array([float(r[key]) for r in rows])
        for key in (
            "pressure_gpa",
            "volume_a3",
            "pressure_plot_sigma_gpa",
            "volume_plot_sigma_a3",
        )
    ]
    results = {}
    for weighting in ("unit", "pressure", "both"):

        def residual(k):
            model = BM3(V0=45.58, K0=k[0], K0_prime=4)
            sigma = 1.0 if weighting == "unit" else sigma_p
            if weighting == "both":
                step = 1e-4
                derivative = (model.pressure(v + step) - model.pressure(v - step)) / (
                    2 * step
                )
                sigma = np.sqrt(sigma_p**2 + (derivative * sigma_v) ** 2)
            return (model.pressure(v) - p) / sigma

        fit = least_squares(residual, [255.0])
        assert fit.success
        results[weighting] = fit.x[0]
    assert results["pressure"] == pytest.approx(254.7407, abs=1e-3)
    assert results["both"] == pytest.approx(255.0081, abs=1e-3)
    assert results["unit"] == pytest.approx(262.2289, abs=1e-3)
    check = document["eos_records"][0]["scientific_validation"][
        "figure_digitization_check"
    ]
    assert check["status"] == "qualified_numerical_check"
    assert abs(results["both"] - 255) < 5
