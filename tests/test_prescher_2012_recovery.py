"""Scientific/source boundaries for the independent Prescher investigation."""

import json

import numpy as np
import pytest

from peritheos.eos.rt import BM3
from peritheos.fitting import fit_rt_eos
from scripts.recover_prescher_2012_fe3c import (
    TABLE,
    VECTORS,
    bm3,
    fit,
    polynomial_pressure_fit,
    selection,
)


def test_vector_markers_are_separate_from_all_19_published_rows():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    evidence = json.loads(VECTORS.read_text())
    assert len(table) == 19
    assert table["volume_a3"][0] == 152.98
    assert table["pressure_gpa"].min() == 1.8
    markers = evidence["markers"]
    assert len(markers) == 22
    assert sorted(
        r["matching_table_s1_row"]
        for r in markers
        if r["matching_table_s1_row"] is not None
    ) == list(range(1, 20))
    extra = [r for r in markers if r["matching_table_s1_row"] is None]
    assert [r["pressure_gpa"] for r in extra] == pytest.approx(
        [0.0019, 5.9312, 44.5251], abs=0.0001
    )
    for row in markers:
        assert row["pressure_gpa"] == pytest.approx(
            np.polyval(evidence["pressure_affine"], row["pdf_x"])
        )
        assert row["volume_a3"] == pytest.approx(
            np.polyval(evidence["volume_affine"], row["pdf_top"])
        )


@pytest.mark.parametrize("name", ["full_range", "below_22_gpa", "above_22_gpa"])
def test_direct_pressure_optimum_matches_linear_bm3_and_library(name):
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    mask = selection(name, table["pressure_gpa"])
    p, v = table["pressure_gpa"][mask], table["volume_a3"][mask]
    result = fit(p, v)
    params = result["parameters"]
    assert list(params.values()) == pytest.approx(
        polynomial_pressure_fit(p, v), rel=1e-6
    )
    model = BM3(**params)
    assert model.pressure(v) == pytest.approx(bm3(v, *params.values()), abs=1e-10)
    for mode, coordinate in (("equal_pressure", False), ("coordinate_widths", True)):
        direct = fit(
            p,
            v,
            mode,
            sp=table["pressure_error_gpa"][mask],
            sv=table["volume_error_a3"][mask],
        )
        native = fit_rt_eos(
            BM3,
            v,
            p,
            initial={"V0": 155, "K0": 160, "K0_prime": 6},
            pressure_sigma=table["pressure_error_gpa"][mask] if coordinate else None,
            volume_sigma=table["volume_error_a3"][mask] if coordinate else None,
            max_nfev=5000,
        )
        assert list(native.parameters.values()) == pytest.approx(
            list(direct["parameters"].values()), rel=2e-5
        )
        assert list(native.standard_errors.values()) == pytest.approx(
            direct["residual_scaled_errors"], rel=2e-3
        )


def test_ambient_marker_changes_the_full_range_optimum_without_a_fixed_v0():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    ambient = json.loads(VECTORS.read_text())["markers"][0]
    result = fit(
        np.r_[table["pressure_gpa"], ambient["pressure_gpa"]],
        np.r_[table["volume_a3"], ambient["volume_a3"]],
    )
    assert result["parameters"]["K0"] == pytest.approx(161.3573, abs=0.002)
    assert result["parameters"]["K0_prime"] == pytest.approx(5.90421, abs=0.0002)
    assert result["fixed_parameters"] == {}
    assert result["observations"] == 20
