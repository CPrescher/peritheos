"""Primary-data and unit checks without importing the native EOS extension."""

import hashlib
import json

import pytest

from scripts.reproduce_argon_barker_1987 import (
    DATA,
    REPORT,
    ROOT,
    audit,
    load_rows,
    molar_to_atomic_a3,
    pixel_to_source,
)


def test_molar_volume_and_axis_units():
    # One mole at 1 A^3/atom occupies 0.602214076 cm^3, independent SI check.
    assert molar_to_atomic_a3(0.602214076) == pytest.approx(1)
    assert pixel_to_source(150, 918) == (8, 0)
    assert pixel_to_source(958, 114) == (22, 450)
    assert pixel_to_source(554, 516) == (15, 225)


def test_liquid_table_and_pressure_diagnostic():
    rows = load_rows("argon-barker-1987-table1-liquid.csv")
    assert [(float(r["pv_over_nkt"]), float(r["u_over_nkt"])) for r in rows] == [
        (2.15, -7.21),
        (0.90, -7.64),
        (2.09, -7.21),
        (0.77, -7.64),
        (1.95, -7.21),
        (0.71, -7.62),
        (1.89, -7.20),
    ]
    result = audit()
    assert result["liquid_table"][0]["pressure_gpa"] == pytest.approx(0.06611, rel=1e-5)
    assert result["liquid_table"][2][
        "relative_pressure_difference_percent"
    ] == pytest.approx(-2.7906976744)
    assert result["status"] == "not_reproduced"
    assert not result["independent_model_reproduction"]
    assert not result["fit_performed"]
    assert result == json.loads(REPORT.read_text())


def test_figure_preserves_theory_and_experiment_separately():
    rows = load_rows("argon-barker-1987-figure5-selected.csv")
    assert sum(r["determination_method"] == "theoretical" for r in rows) == 21
    assert sum(r["determination_method"] == "experimental" for r in rows) == 11
    for r in rows:
        v, p = pixel_to_source(float(r["pixel_x"]), float(r["pixel_y"]))
        assert float(r["molar_volume_cm3_mol"]) == pytest.approx(v, abs=5e-7)
        assert float(r["pressure_gpa"]) == pytest.approx(p / 10, abs=5e-7)
        assert float(r["volume_atomic_a3"]) == pytest.approx(
            molar_to_atomic_a3(v), abs=5e-7
        )
        assert r["temperature_k"] == "298"
    report = audit()["figure5"]
    assert report["experimental_points_outside_sampled_curve"] == 2
    assert len(report["comparisons"]) == 9


def test_study_assets_and_reproduction_status():
    study = json.loads(
        (ROOT / "peritheos/data/studies/argon-barker-1987.json").read_text()
    )
    assert study["reproduction"]["status"] == "not_reproduced"
    assert study["reproduction"]["eos_record_identifiers"] == []
    assert study["source_audit"]["pdf"]["identity_verified"]
    assert study["source_audit"]["pdf"]["provenance_type"] == "user_supplied"
    for dataset in study["datasets"]:
        path = DATA / dataset["resource"]["path"].split("/")[-1]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
