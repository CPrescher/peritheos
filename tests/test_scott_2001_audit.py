"""Scientific and preservation checks for the separate Scott diagnostic."""

import hashlib
import json

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.fitting import fit_rt_eos
from scripts.audit_scott_2001 import (
    CSV,
    EXTRACTION,
    REPORT,
    SOURCE,
    V0,
    linear_fit,
    pressure,
    replay,
    volume,
)
from scripts.validate_primary_eos_refits import _series


def original_arrays():
    rows = np.genfromtxt(CSV, delimiter=",", names=True)
    return rows, rows["volume_ratio"] * V0, rows["pressure_gpa"]


def test_independent_finite_strain_equation_and_volume_basis():
    # Eulerian finite strain form, independent of the direct power implementation.
    f = np.array([0.0, 0.025, 0.05, 0.08])
    v = V0 / (1 + 2 * f) ** 1.5
    expected = 3 * 175.4 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (5.1 - 4) * f)
    np.testing.assert_allclose(pressure(v, *SOURCE), expected, atol=2e-12)
    np.testing.assert_allclose(pressure(v / 4, *SOURCE, V0 / 4), expected, atol=2e-12)
    np.testing.assert_allclose(
        BM3(V0=V0, K0=175.4, K0_prime=5.1).pressure(v), expected, atol=2e-12
    )
    np.testing.assert_allclose(volume(expected, *SOURCE), v, atol=2e-10)


def test_original_data_and_published_coefficients_preserved():
    assert (
        hashlib.sha256(CSV.read_bytes()).hexdigest()
        == "bf198807ec30f0b978d73e6d31791e07677ba8c8d393ff706e5875d6039e5b19"
    )
    doc = get_material_document("cementite")
    record = next(
        r for r in doc["eos_records"] if r["identifier"] == "cementite_scott_2001_bm3_1"
    )
    for name, value in [("V0", 155.26), ("K0", 175.4), ("K0_prime", 5.1)]:
        assert record["eos"]["parameters"][name] == value
    dataset = next(
        d
        for d in doc["datasets"]
        if d["identifier"] == "cementite_scott_2001_figure3_digitized"
    )
    series = _series(doc, record, dataset)
    assert len(series.volume) == 21
    assert series.pressure_sigma is None and series.volume_sigma is None


def test_fixed_v0_global_linear_solution_matches_both_fitting_implementations():
    _, v, p = original_arrays()
    direct = replay(v, p)
    expected = [161.3414242272, 6.0441242841]
    np.testing.assert_allclose(linear_fit(v, p), expected, atol=2e-7)
    np.testing.assert_allclose(list(direct["parameters"].values()), expected, atol=2e-5)
    fits = [
        fit_rt_eos(
            BM3,
            v,
            p,
            initial={"K0": 175.4, "K0_prime": 5.1},
            fixed={"V0": V0},
            absolute_sigma=flag,
        )
        for flag in [True, False]
    ]
    np.testing.assert_allclose(
        fits[0].covariance, direct["covariance_residual_scaled"], rtol=2e-6
    )
    np.testing.assert_array_equal(fits[0].covariance, fits[1].covariance)
    assert not direct["central_pair_within_printed_widths"]


def test_visible_width_weighting_effect_is_separate_from_subset_effect():
    rows, v, p = original_arrays()
    mask = np.isfinite(rows["pressure_error_gpa"])
    assert mask.sum() == 13
    equal = linear_fit(v[mask], p[mask])
    weighted = linear_fit(v[mask], p[mask], rows["pressure_error_gpa"][mask])
    np.testing.assert_allclose(equal, [161.398747, 6.0930355], atol=2e-5)
    np.testing.assert_allclose(weighted, [171.982417, 5.27476475], atol=2e-5)
    assert np.all(np.abs(weighted - SOURCE) < [3.5, 0.3])
    assert abs(equal[0] - linear_fit(v, p)[0]) < 0.1


def test_separate_graph_extraction_replays_conditional_agreement():
    extraction = json.loads(EXTRACTION.read_text())
    report = json.loads(REPORT.read_text())
    assert (
        extraction["original_csv_sha256"]
        == hashlib.sha256(CSV.read_bytes()).hexdigest()
    )
    assert extraction["native_image_size"] == [2381, 3387]
    assert len(extraction["points"]) == 21
    v = np.array([r["volume_ratio"] for r in extraction["points"]]) * V0
    p = np.array([r["pressure_gpa"] for r in extraction["points"]])
    conditional = linear_fit(v[1:], p[1:], p[1:])
    np.testing.assert_allclose(conditional, [173.624746, 5.21386780], atol=3e-5)
    assert np.all(np.abs(conditional - SOURCE) < [3.5, 0.3])
    assert (
        report["drawn_curve_equation_check"]["source_comparison"]["published"][
            "max_abs_x_pixels"
        ]
        < 3
    )
    assert report["drawn_curve_equation_check"]["conditional_curve_fit"][
        "central_pair_within_printed_widths"
    ]
    # The exact ambient point has zero leverage on fixed-V0 equal-pressure fits.
    _, original_v, original_p = original_arrays()
    np.testing.assert_allclose(
        linear_fit(original_v[1:], original_p[1:]),
        linear_fit(original_v, original_p),
        atol=1e-10,
    )


@pytest.mark.parametrize("v0", [155.12, 155.4])
def test_ambient_error_range_alone_does_not_restore_pair(v0):
    _, v, p = original_arrays()
    pair = linear_fit(v, p, v0=v0)
    assert not np.all(np.abs(pair - SOURCE) <= [3.5, 0.3])
