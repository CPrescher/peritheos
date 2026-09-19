"""Primary Fischer tables, composition, independent curves and acceptance scope."""

import csv
import hashlib
import json
import re
from collections import Counter

import numpy as np
import pytest
from scipy.constants import Avogadro

from peritheos import Material, get_material_document, list_material_documents
from scripts.reproduce_fischer_2014_fesi import (
    ATOMS_PER_CELL,
    ROOT,
    STATIC,
    VINET_STATIC,
    load_data,
    pressure,
    reproduce,
)
from scripts.validate_primary_eos_refits import _fit_record, _series

MATERIALS = ("fesi_b20", "fe084si016_d03", "fe084si016_hcp", "fe073si027_d03")


def test_primary_tables_preserve_full_observations_markers_and_missing_errors():
    for table, expected, counts in (
        (1, 211, {"D03": 15, "hcp": 76, "hcp+B2": 96, "fcc+B2": 13, "fcc+hcp+B2": 11}),
        (2, 180, {"B20": 27, "B20+B2": 39, "B2": 114}),
    ):
        path = ROOT / f"peritheos/data/datasets/fischer-2014-table-s{table}-pvt.csv"
        rows = list(csv.DictReader(path.read_text().splitlines()))
        assert len(rows) == expected
        assert Counter(row["phase"] for row in rows) == counts
        assert rows[0]["pressure_uncertainty_gpa"] == "0.0"
        assert rows[0]["kbr_a_angstrom"] == ""
        assert all(row["kbr_a_angstrom"] for row in rows[1:])
        if table == 1:
            assert float(rows[0]["d03_a_angstrom"]) == 5.697387235035957
            flagged = next(r for r in rows if r["source_row"] == "209")
            assert flagged["hcp_a_uncertainty_angstrom"] == ""
            assert "hcp_a_uncertainty_angstrom=*" in flagged["source_flags"]
        else:
            assert float(rows[-1]["pressure_gpa"]) == 146.9148894617161
            flagged = next(r for r in rows if r["source_row"] == "62")
            assert flagged["b20_a_angstrom"] == flagged["b2_a_angstrom"] == ""
            assert "b2_a_angstrom=**" in flagged["source_flags"]
    rows = list(
        csv.DictReader(
            (ROOT / "peritheos/data/datasets/fischer-2012-table-s2-pvt.csv")
            .read_text()
            .splitlines()
        )
    )
    assert len(rows) == 185
    assert Counter(r["phase"] for r in rows) == {
        "D03": 98,
        "hcp+B2": 78,
        "B2_volume_only": 9,
    }
    assert float(rows[0]["d03_molar_atomic_volume_cm3_mol"]) == 6.799456088159677
    assert float(rows[-1]["kbr_temperature_k"]) == 1281.625


@pytest.mark.parametrize("material", MATERIALS)
def test_composition_cell_contents_and_source_structure_are_consistent(material):
    d = get_material_document(material)
    counts = Counter()
    for site in d["atom_sites"]:
        counts[site["element"]] += (
            int(re.match(r"\d+", site["wyckoff"])[0]) * site["occupancy"]
        )
    atoms = ATOMS_PER_CELL[material]
    si_fraction = (
        0.5
        if material == "fesi_b20"
        else 0.27
        if material == "fe073si027_d03"
        else 0.16
    )
    assert counts["Si"] == pytest.approx(atoms * si_fraction)
    assert counts["Fe"] == pytest.approx(atoms * (1 - si_fraction))
    assert d["space_group_number"] in (198, 225, 194)
    if material.endswith("d03"):
        assert "not a published refinement" in d["source"]["structure_location"]
    dataset = d["datasets"][0]
    assert (
        hashlib.sha256(
            (ROOT / "peritheos/data" / dataset["resource"]["path"]).read_bytes()
        ).hexdigest()
        == dataset["resource"]["sha256"]
    )


@pytest.mark.parametrize("material", MATERIALS)
def test_native_models_match_independent_atomic_normalization_and_inverse(material):
    d = get_material_document(material)
    m = Material.from_eosmat(d)
    rows = load_data(material)
    phase = material.rsplit("_", 1)[1]
    cell_v = np.array([float(r[phase + "_volume_a3"]) for r in rows])
    atomic_v = cell_v * Avogadro / 1e24 / ATOMS_PER_CELL[material]
    t = np.array([float(r["temperature_k"]) for r in rows])
    for stored, record in zip(d["eos_records"], m.eos_records):
        model = "bm3" if stored["eos"]["type"] == "BM3" else "vinet"
        expected = pressure(material, atomic_v, t, model)
        np.testing.assert_allclose(
            record.pressure(cell_v, t), expected, atol=3e-10, rtol=0
        )
        np.testing.assert_allclose(
            record.volume(expected, t), cell_v, atol=2e-8, rtol=0
        )
        p_low, p_high = stored["experimental_pressure_range_gpa"]
        t_low, t_high = stored["experimental_temperature_range_k"]
        p_mid, t_mid = (p_low + p_high) / 2, (t_low + t_high) / 2
        assert record.volume(p_mid, t_mid, check_validity=True) > 0
        with pytest.raises(ValueError, match="validity"):
            record.volume(p_high + 1, t_mid, check_validity=True)
        with pytest.raises(ValueError, match="validity|isothermal"):
            record.volume(p_mid, t_high + 1, check_validity=True)
        source_static = (STATIC if model == "bm3" else VINET_STATIC)[material]
        assert stored["eos"]["parameters"]["V0"] == pytest.approx(
            source_static[0] * ATOMS_PER_CELL[material] / 0.602214076
        )
        assert stored["eos"]["parameters"]["K0"] == source_static[1]
        assert stored["eos"]["parameters"]["K0_prime"] == source_static[2]
        assert stored["parameter_error_confidence"] is None
        assert stored["parameter_covariance"] is None
        assert stored["record_kind"] == "published"
        assert record.pressure(stored["eos"]["parameters"]["V0"], 300) == pytest.approx(
            0, abs=1e-10
        )
        if "thermal" in stored:
            assert stored["thermal"]["parameters"]["n"] == 1
            assert stored["thermal"]["debye_temperature_law"] == "integrated_gruneisen"


@pytest.mark.parametrize("material", MATERIALS)
def test_native_refit_selects_single_phase_and_retains_source_constraints(material):
    d = get_material_document(material)
    for record in d["eos_records"]:
        result = _fit_record(d, record, d["datasets"][0])
        assert result["status"] == "parity"
        assert (
            result["observations"]
            == {
                "fesi_b20": 27,
                "fe084si016_d03": 15,
                "fe084si016_hcp": 76,
                "fe073si027_d03": 66,
            }[material]
        )
        series = _series(d, record, d["datasets"][0])
        assert len(series.pressure) == result["observations"]
        if material == "fe073si027_d03":
            assert record["fixed_parameters"] == ["V0", "K0", "K0_prime"]
            assert np.all(series.temperature > 300)
            assert result["parameters"][0]["refit"] == pytest.approx(
                1.8612195, abs=2e-6
            )


def test_primary_observational_reproduction_and_b2_normalization_hold():
    result = reproduce()
    # Published observed-minus-calculated residual limits, rounded to 0.1 GPa.
    # Allow 0.15 GPa for source coefficient rounding and thermal q precision.
    for key, limits in (
        ("fesi_b20_bm3", [-0.9, 0.7]),
        ("fe084si016_d03_bm3", [-1.1, 0.5]),
        ("fe084si016_hcp_bm3", [-2.9, 2.7]),
    ):
        residual = result[key]["published_residual_range_gpa"]
        np.testing.assert_allclose(
            [-residual[1], -residual[0]], limits, atol=0.17, rtol=0
        )
    b2 = result["fesi_b2_bm3"]
    assert b2["observations"] == 114
    assert b2["published_rmse_gpa"] == pytest.approx(7.147574357, abs=2e-8)
    assert b2["parameters"]["gamma0"] > 2.5
    diagnostic = b2["unpublished_double_thermal_amplitude_diagnostic"]
    assert diagnostic["production_accepted"] is False
    assert diagnostic["residual_range_gpa"] == pytest.approx(
        [-5.999806268, 3.924682636], abs=2e-8
    )
    assert result["fe073si027_d03_vinet"]["static_observations"] == 32
    assert result["fe073si027_d03_vinet"]["thermal_observations"] == 66
    # No unsupported normalization or homogeneous phase-mixture is promoted.
    source_records = [
        r
        for name in list_material_documents()
        for r in get_material_document(name)["eos_records"]
        if r.get("reference", {}).get("doi") == "10.1002/2013JB010898"
    ]
    assert len(source_records) == 9
    published = [r for r in source_records if r["record_kind"] == "published"]
    assert len(published) == 7
    assert all(
        "b2" not in r["identifier"] or "b20" in r["identifier"] for r in published
    )
    assert {r["identifier"] for r in source_records if r["record_kind"] == "refit"} == {
        "fesi_b2_fischer_2014_bm3_refit",
        "fesi_b2_fischer_2014_vinet_refit",
    }
    audit = json.loads((ROOT / "docs/data/fischer-2014-source-audit.json").read_text())
    assert len(audit["supplement_workbooks"]) == 7
