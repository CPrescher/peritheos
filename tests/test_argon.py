"""Published argon curves, source tables and the historical hcp diagnostic."""

import csv
import hashlib
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, get_material_document
from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import Dewaele2006
from peritheos.errors import EosValidationError, MaterialError
from scripts.reproduce_argon import REPORT, ROOT, dewaele, reproduce, rows


class PythonDewaele(Dewaele2006):
    """Exercise the independent Python path."""


def test_fcc_identity_and_lossless_source_tables():
    doc = get_material_document("argon_fcc")
    schema = json.loads((ROOT / "peritheos/data/eosmat-v3.schema.json").read_text())
    assert not list(Draft202012Validator(schema).iter_errors(doc))
    assert (doc["formula_units_per_cell"], doc["space_group_number"]) == (4, 225)
    assert doc["lattice"]["a"] == 5.0868
    assert doc["atom_sites"][0]["wyckoff"] == "4a"
    for ds in doc["datasets"]:
        path = ROOT / "peritheos/data" / ds["resource"]["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == ds["resource"]["sha256"]
        with path.open() as stream:
            assert next(csv.reader(stream)) == [c["name"] for c in ds["columns"]]
    table = rows("argon-fcc-dewaele-2021-supplement.csv")
    assert len(table) == 288
    assert sum(r["selected_for_room_temperature_fit"] == "1" for r in table) == 95
    assert table[0]["argon_a_angstrom"] == "5.0865998"
    missing = next(r for r in table if r["sample"] == "ArNe12_017")
    assert missing["pressure_gpa"] == ""
    assert missing["gold_a_angstrom"] == "4.8579998"
    assert missing["selected_for_room_temperature_fit"] == "0"
    assert sum(r["sample"] == "Ar_cell4_027" for r in table) == 2
    assert min(float(r["temperature_k"]) for r in table) == 5.5
    assert (
        max(float(r["pressure_gpa"]) for r in table if r["pressure_gpa"]) == 113.73376
    )
    ono = rows("argon-fcc-ono-2020-table1.csv")
    assert len(ono) == 19
    assert [
        float(ono[-1][key])
        for key in (
            "pressure_gpa",
            "pressure_sigma_gpa",
            "volume_a3",
            "volume_sigma_a3",
        )
    ] == [136.7, 0.3, 47.26, 0.10]


def test_source_reproduction_and_diagnostics(assert_audit_close):
    result = reproduce()
    assert_audit_close(result, json.loads(REPORT.read_text()))
    for key in ("dewaele_2021", "ono_2020"):
        audit = result[key]
        assert audit["native_max_difference_gpa"] < 2e-11
        for point in audit["benchmarks"]:
            assert (
                abs(point["source_pressure_gpa"] - point["calculated_pressure_gpa"])
                < point["tolerance_gpa"]
            )
        assert all(f["solver_success"] for f in audit["refits"].values())
    assert result["dewaele_2021"]["low_temperature_native_max_difference_gpa"] < 2e-11
    assert result["dewaele_2021"]["published_rmse_gpa"] < 0.35
    assert result["ono_2020"]["published_rmse_gpa"] < 0.85
    assert result["ono_2020"]["refits"]["equal_pressure"][
        "parameters"
    ] == pytest.approx([184.5, 1.07, 8.02], rel=0.005)
    # Published Table 2 coefficients are retained despite alternative weights.
    doc = get_material_document("argon_fcc")
    assert doc["eos_records"][1]["eos"]["parameters"] == dict(
        V0=184.5, K0=1.07, K0_prime=8.02
    )


def test_absolute_zero_pressure_and_reference_increment_roundtrip():
    doc = get_material_document("argon_fcc")
    material = Material.from_eosmat(doc)
    record = material.get_eos_record("argon_fcc_dewaele_2021_vinet_mgd")
    model = record.eos
    fallback = PythonDewaele(
        model.rt_eos,
        **doc["eos_records"][0]["thermal"]["parameters"],
        thermal_pressure_reference="absolute_zero",
    )
    v = np.array([130.0, 90.0, 60.0])
    t = np.array([296.0, 150.0, 10.0])
    assert record.pressure(v, t) == pytest.approx(dewaele(v, t), abs=5e-11)
    assert fallback.pressure(v * record.volume_scale, t) == pytest.approx(
        record.pressure(v, t), abs=5e-11
    )
    assert record.volume(record.pressure(v, t), t) == pytest.approx(v, rel=1e-10)
    assert model.thermal_pressure_increment(
        v * record.volume_scale, 296
    ) == pytest.approx(0.0, abs=1e-13)
    assert model.thermal_pressure(model.rt_eos.V0, 296) > 0
    assert model.thermal_pressure_increment(
        v * record.volume_scale, t
    ) == pytest.approx(
        model._native.evaluate_array(
            "thermal_pressure_increment", v * record.volume_scale, t
        ),
        abs=5e-11,
    )
    restored = Material.from_eosmat(material.to_eosmat()).get_eos_record(
        record.identifier
    )
    assert restored.pressure(v, t) == pytest.approx(record.pressure(v, t), abs=5e-11)
    assert restored.eos.configuration_values() == {
        "thermal_pressure_reference": "absolute_zero"
    }
    with pytest.raises(MaterialError):
        record.pressure(40, 296, check_validity=True)


def test_absolute_baseline_rejects_invalid_choices_and_preserves_default():
    params = dict(
        Tr=296,
        theta0=93.3,
        gamma0=2.7,
        gamma_inf=0.5,
        beta=1,
        anharmonic_a=1e-5,
        anharmonic_m=1,
        electronic_e=1e-5,
        electronic_g=1,
        n=1,
    )
    reference = Vinet(2.288, 2.65, 7.423)
    for choice in ("reference_isentrope", "unknown", None):
        with pytest.raises(EosValidationError):
            Dewaele2006(reference, **params, thermal_pressure_reference=choice)
    default = Dewaele2006(reference, **params)
    absolute = Dewaele2006(
        reference, **params, thermal_pressure_reference="absolute_zero"
    )
    assert default.thermal_pressure(1.4, 296) == pytest.approx(0)
    assert absolute.thermal_pressure_increment(1.4, 500) == pytest.approx(
        default.thermal_pressure(1.4, 500)
    )
    for temperature in (0, -1, float("nan")):
        with pytest.raises(EosValidationError):
            absolute.pressure(1.4, temperature)


def test_wittlinger_diagnostic_compares_original_observation_residuals():
    audit = reproduce()["wittlinger_1997"]
    assert audit["published_rmse_gpa"] == pytest.approx(0.772262, abs=1e-6)
    fit = audit["errors_in_variables"]
    assert 0.44 < fit["original_coordinate_rmse_gpa"] < 0.46
    assert fit["adjusted_coordinate_rmse_gpa"] < 0.02
    assert audit["published_reduced_chi_square"] < 1
    doc = get_material_document("argon_hcp")
    record = doc["eos_records"][0]
    validation = record["scientific_validation"]
    assert validation["reproduction_status"] == "not_reproduced"
    assert validation["independent_numerical_check"]["status"] == "not_reproduced"
    assert audit["reproduction_status"] == "not_reproduced"
    from scripts.reproduce_argon import ledger_outcome

    outcome = ledger_outcome(record)
    assert outcome["status"] == "not_refittable"
    assert outcome["reproduction_status"] == "not_reproduced"
    assert "NOT REPRODUCED" in outcome["reason"]
    assert record["eos"]["parameters"] == {"V0": 78.0, "K0": 6.5}
    assert "8.5 GPa" in record["scientific_validation"]["usage_recommendation"]
    with pytest.raises(MaterialError):
        Material.from_eosmat(doc).eos_records[0].volume(50, check_validity=True)
