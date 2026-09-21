import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_litasov_2007_superhydrous_phase_b import (
    CASES,
    DATA,
    DATASET,
    OUTPUT,
    PREFIX,
    bm3,
    ledger_outcome,
    observations,
    parameter_names,
    pressure,
    reproduce,
    rows,
)

ROOT = Path(__file__).resolve().parents[1]


def test_complete_primary_data_uncertainties_and_peak_assignments():
    doc = get_material_document("superhydrous_phase_b_lt")
    data = rows()
    assert [int(r["measurement_number"]) for r in data] == list(range(4, 70)) + [
        -1,
        -2,
        -3,
    ]
    assert sum(r["temperature_k"] == "300" for r in data) == 20
    assert sum(bool(r["au_volume_a3"]) for r in data) == 66
    assert float(data[0]["volume_a3"]) == 555.31
    assert float(data[22]["c_sigma_angstrom"]) == 0.012
    assert data[-1]["au_volume_a3"] == data[-1]["pressure_sigma_gpa"] == ""
    assert (
        doc["datasets"][0]["uncertainty"]["reported_as"]
        == "estimated_standard_deviation"
    )
    for dataset in doc["datasets"]:
        resource = ROOT / "peritheos/data" / dataset["resource"]["path"]
        assert (
            hashlib.sha256(resource.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
        assert len(resource.read_text().splitlines()) - 1 == dataset["row_count"]
    material = Material.from_eosmat(doc)
    assert len(material.get_dataset(DATASET)["temperature_k"]) == 69
    peaks = material.get_dataset(PREFIX + "table3")
    assert len(peaks["relative_intensity"]) == 30
    assert sum(peaks["phase_code"] == 1) == 10
    assert DATA.name == "superhydrous-phase-b-litasov-2007-table2.csv"


@pytest.mark.parametrize("suffix", CASES)
def test_published_models_reproduce_primary_observations_and_invert(suffix):
    data, (v, t, observed) = observations(suffix)
    record = get_eos_record(PREFIX + suffix)
    calculated = record.pressure(v, t)
    assert calculated == pytest.approx(pressure(suffix, v, t), abs=5e-10)
    assert record.volume(calculated, t) == pytest.approx(v, abs=1e-8)
    # Source pressures have 0.1–0.3 GPa estimated SDs; the complete source-table
    # comparison gives RMS <0.15 GPa, not just a reference-state identity.
    assert np.sqrt(np.mean((calculated - observed) ** 2)) < 0.15
    assert len(data) == (20 if suffix.startswith("bm3") else 69)
    pars = CASES[suffix][0]
    assert record.pressure(pars[0], 300) == pytest.approx(0, abs=1e-12)
    assert record.pressure(600.0, 300.0) == pytest.approx(
        bm3(600.0, *pars[:3]), abs=1e-10
    )
    stored = next(
        r
        for r in get_material_document("superhydrous_phase_b_lt")["eos_records"]
        if r["identifier"] == record.identifier
    )
    assert stored["record_kind"] == "published"
    assert stored["eos"]["parameters"] == dict(zip(["V0", "K0", "K0_prime"], pars[:3]))
    assert stored["experimental_pressure_range_gpa"] == [
        float(min(observed)),
        float(max(observed)),
    ]
    assert stored["fit_datasets"] == [DATASET]
    assert ledger_outcome(stored)["status"] == "similar"


def test_independent_refits_and_sign_controls_are_current(assert_audit_close):
    result = reproduce()
    saved = json.loads(OUTPUT.read_text())

    # Compare fitted coefficients rather than magnifying their numerical drift
    # by dividing through small published-minus-refitted differences. Check those
    # derived fields for internal consistency, and keep the physical checks below.
    def comparable(report):
        if isinstance(report, list):
            return [comparable(item) for item in report]
        if not isinstance(report, dict):
            return report
        derived = set()
        if {"parameter", "published", "refit"} <= report.keys():
            difference = report["refit"] - report["published"]
            assert report["difference"] == pytest.approx(difference, abs=1e-14)
            assert report["relative_difference"] == pytest.approx(
                abs(difference / report["published"]), abs=1e-14
            )
            derived = {"difference", "relative_difference"}
        return {k: comparable(v) for k, v in report.items() if k not in derived}

    assert_audit_close(comparable(result), comparable(saved), rel=1e-4, abs=1e-12)
    for suffix in CASES:
        fit = result[suffix]
        assert fit["rmse_gpa"] == pytest.approx(
            saved[suffix]["rmse_gpa"], rel=1e-8, abs=1e-10
        )
        _, (volume, temperature, _) = observations(suffix)

        def coefficients(report):
            fitted = {p["parameter"]: p["refit"] for p in report["parameters"]}
            return [
                fitted.get(name, value)
                for name, value in zip(parameter_names(suffix), CASES[suffix][0])
            ]

        np.testing.assert_allclose(
            pressure(suffix, volume, temperature, coefficients(fit)),
            pressure(suffix, volume, temperature, coefficients(saved[suffix])),
            rtol=0,
            atol=5e-6,
        )
        assert fit["solver_success"]
        assert all(p["similar"] for p in fit["parameters"])
        assert all(
            abs(p["difference"]) < 2 * p["published_error"] for p in fit["parameters"]
        )
    free = result["mgd_free_debye_a89"]

    def theta(fit):
        return next(p["refit"] for p in fit["parameters"] if p["parameter"] == "theta0")

    assert theta(free) == pytest.approx(553, abs=1)
    assert theta(free["literal_equation9_control"]) == pytest.approx(673, abs=1)
    assert result["elastic_debye_temperature_k"] == pytest.approx(860, abs=3)
    assert result["ambient_mean_volume_a3"] == pytest.approx(623.34, abs=1e-10)
    assert min(free["wrong_atom_count_controls_rmse_gpa"].values()) > 0.3
    for state in result["figure10_density_checks"]:
        assert state["scale_difference_g_cm3"] <= 0.03
        assert state["model_spread_g_cm3"] <= 0.02
    au = result["gold_anderson_diagnostic"]
    assert au["max_abs_residual_gpa"] == pytest.approx(0.0904154, abs=1e-6)
    assert au["rmse_gpa"] < 0.04


def test_corrected_structure_has_seventy_atoms_in_source_setting():
    doc = get_material_document("superhydrous_phase_b_lt")
    assert doc["formula_units_per_cell"] == 2
    assert doc["space_group_number"] == 34
    assert doc["lattice"]["a"] == 14.024
    counts = Counter()
    for site in doc["atom_sites"]:
        x, y, z = site["x"], site["y"], site["z"]
        positions = np.mod(
            [
                [x, y, z],
                [0.5 + x, 0.5 - y, 0.5 + z],
                [0.5 - x, 0.5 + y, 0.5 + z],
                [-x, -y, z],
            ],
            1,
        )
        multiplicity = len(np.unique(np.round(positions, 8), axis=0))
        assert multiplicity == int(site["wyckoff"][0])
        counts[site["element"]] += multiplicity * site["occupancy"]
    assert counts == {"Mg": 20, "Si": 6, "O": 36, "H": 8}
    sites = {s["label"]: s for s in doc["atom_sites"]}
    assert sites["O6"]["z"] == 0.007
    assert [sites["H2"][k] for k in ["x", "y", "z"]] == [0.45, 0.18, 0.27]
    archive = doc["structure_provenance"]["archive_resource"]
    assert (
        hashlib.sha256(
            (ROOT / "peritheos/data" / archive["path"]).read_bytes()
        ).hexdigest()
        == archive["sha256"]
    )
    assert "not rescaled" in doc["notes"]
    for record in doc["eos_records"]:
        calibration = record["pressure_calibration"]
        assert calibration["status"] == "partially_resolved"
        assert "reference_eos_record" not in calibration["methods"][0]
    assert (
        next(r for r in doc["eos_records"] if r["default"])["identifier"]
        == PREFIX + "htbm_a89"
    )


def test_incomplete_and_axial_alternatives_remain_outside_production():
    audit = json.loads(
        (
            ROOT / "docs/data/litasov-2007-superhydrous-phase-b-source-audit.json"
        ).read_text()
    )
    doc = get_material_document("superhydrous_phase_b_lt")
    assert audit["accepted_records"] == [r["identifier"] for r in doc["eos_records"]]
    assert len(audit["held_candidates"]) == 7
    assert all(r["reason"] for r in audit["held_candidates"])
    assert audit["table7_axial_fits"]["a"]["K0_gpa"] == [125.7, 0.3]
    assert all(
        r["reference"]["doi"] == "10.1016/j.pepi.2007.06.003"
        for r in doc["eos_records"]
    )
