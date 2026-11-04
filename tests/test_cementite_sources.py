import copy
import csv
import hashlib
import json
from collections import Counter

import numpy as np
import pytest
from jsonschema import Draft202012Validator
from scipy.constants import Avogadro, R
from scipy.integrate import quad

from peritheos import (
    Material,
    convert_eosmat,
    eosmat_schema,
    get_material_document,
    list_eos_record_documents,
    load_eosmat,
)
from scripts.reproduce_cementite import ROOT, SCOTT_ID, qualify_scott, reproduce

DATA = ROOT / "peritheos/data/datasets"
RESEARCH = ROOT / "peritheos/data/research/cementite"


@pytest.fixture(scope="module")
def audit():
    return reproduce()


def test_scott_broad_uncertainty_overlap_cannot_establish_coefficient_parity(audit):
    outcome = audit["scott_2001"]
    assert outcome["status"] == "similar"
    assert outcome["source_parameter_widths_agreement"] is False
    assert not outcome["pressure_sigmas_used"]
    assert not outcome["volume_sigmas_used"]
    assert outcome["original_author_workflow_parity"] is False
    assert all(p["within_combined_2sigma"] for p in outcome["parameters"])
    original = copy.deepcopy(outcome)
    inflated = copy.deepcopy(outcome)
    inflated["status"] = "parity"
    for p in inflated["parameters"]:
        p["refit_error"] *= 100
    assert qualify_scott(inflated)["status"] == "similar"
    assert outcome == original
    assert -0.99 < outcome["K0_K0_prime_correlation"] < -0.98
    assert outcome["curve_comparison"]["max_abs_relative_volume_difference"] < 0.0031
    record = get_material_document("cementite", format_version=4)["eos_records"][0]
    assert record["parameters"]["V0"]["value"] == 155.26
    assert record["parameters"]["K0"]["value"] == 175.4
    assert record["parameters"]["K0_prime"]["value"] == 5.1
    assert record["evidence"]["reproduction"]["assessment"]["status"] == "similar"
    ledger = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    row = next(r for r in ledger["records"] if r["record_identifier"] == SCOTT_ID)
    assert row["status"] == "similar"
    assert row["parameters"] == outcome["parameters"]


def test_litasov_all_rows_errors_and_original_tokens_are_preserved():
    with (DATA / "cementite-litasov-2013-table1.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    raw = json.loads((DATA / "cementite-litasov-2013-table1-tokens.json").read_text())
    assert len(rows) == len(raw) == 81
    assert Counter(int(r["run"]) for r in rows) == {
        1: 29,
        2: 9,
        3: 8,
        4: 11,
        5: 8,
        6: 11,
        7: 5,
    }
    assert float(rows[0]["pressure_sigma_gpa"]) == 0.04
    assert float(rows[0]["volume_sigma_a3"]) == 0.01
    assert float(rows[-1]["volume_a3"]) == 148.59
    assert max(float(r["pressure_gpa"]) for r in rows) == 31
    assert sum(r["mgo_volume_a3"] == "" for r in rows) == 3
    for row, original in zip(rows, raw):
        assert int(row["source_row"]) == original["row"]
        assert float(row["volume_a3"]) == float(original["tokens"][-1].split("(")[0])
        if original["tokens"][1] is None:
            assert row["mgo_volume_a3"] == row["mgo_volume_sigma_a3"] == ""


def test_takahashi_all_rows_and_excel_error_precision_survive_loading():
    document = load_eosmat(
        ROOT / "peritheos/data/materials/cementite_nonmagnetic.eosmat"
    )
    material = Material.from_eosmat(document)
    dataset = material.get_dataset("cementite_takahashi_2019_table_s1")
    assert dataset.checksum_verified
    assert len(dataset) == 60
    assert dataset["volume_a3"][1] == 121.30
    assert dataset["volume_a3_error"][1] == 0.08
    assert dataset["temperature_k"].max() == 2600
    assert np.isnan(dataset["temperature_k_error"][2])
    assert np.isfinite(dataset["pressure_kbr_gpa"]).sum() == 50
    assert np.isfinite(dataset["pressure_mgo_gpa"]).sum() == 50
    assert np.isfinite(dataset["pressure_sio2_gpa"]).sum() == 10
    assert dataset["pressure_sio2_gpa"][-1] == 326.4
    assert dataset["sio2_volume_a3"][-1] == 57.30
    assert dataset["sio2_volume_a3_error"][-1] == 0.37
    raw = json.loads(
        (DATA / "cementite-takahashi-2019-table-s1-cells.json").read_text()
    )
    assert raw["rows"][1]["cells"]["P5"] == {"value": 121.3, "number_format": "0.00"}
    assert raw["rows"][1]["cells"]["Q5"]["value"] == "(8)"
    assert raw["rows"][-1]["cells"]["A63"]["value"] == "Fe3C040_43"


def independent_pressure(volume, temperature, values):
    """Direct source Equations 1-6, quadrature and SI cell-volume conversion."""
    v0, k0, kp, gamma0, q = values
    eta = (v0 / volume) ** (1 / 3)
    static = 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))
    gamma = gamma0 * (volume / v0) ** q
    theta = 314 * np.exp((gamma0 - gamma) / q)

    def energy(t):
        integral = quad(lambda x: x**3 / np.expm1(x), 0, theta / t)[0]
        return 9 * 4 * R * t * (t / theta) ** 3 * integral

    molar_m3 = volume * 1e-30 * Avogadro / 4
    return static + gamma * (energy(temperature) - energy(300)) / molar_m3 / 1e9


@pytest.mark.parametrize(
    "scale,values",
    [
        ("kbr", [148.8, 311.1, 3.4, 1.06, 1.92]),
        ("mgo", [147.3, 323, 3.43, 1.37, 0.98]),
    ],
)
def test_catalog_models_match_source_equations_and_roundtrip(scale, values):
    document = load_eosmat(
        ROOT / "peritheos/data/materials/cementite_nonmagnetic.eosmat"
    )
    Draft202012Validator(eosmat_schema(4)).validate(document)
    material = Material.from_eosmat(document)
    identifier = f"cementite_nonmagnetic_takahashi_2019_{scale}_bm3_mgd"
    record = material.get_eos_record(identifier)
    volume = np.array([120.0, 110.0, 94.0])
    temperature = np.array([300.0, 2000.0, 300.0])
    expected = [independent_pressure(v, t, values) for v, t in zip(volume, temperature)]
    pressure = record.pressure(volume, temperature)
    np.testing.assert_allclose(pressure, expected, rtol=2e-10)
    np.testing.assert_allclose(record.volume(pressure, temperature), volume, rtol=1e-10)
    # Exercise the Python thermal-pressure path; the RT backend remains native.
    record.eos._native = None
    np.testing.assert_allclose(
        record.pressure(volume, temperature), expected, rtol=2e-10
    )
    assert material.default_equilibrium_record().identifier.endswith("_kbr_bm3_mgd")
    assert identifier in list_eos_record_documents()
    restored = convert_eosmat(convert_eosmat(document, version=3), version=4)
    # Native v4 may omit its fixed implicit unit declaration when normalizing.
    # Compare the complete compatibility view, including effective units.
    assert convert_eosmat(restored, version=3) == convert_eosmat(document, version=3)


def test_takahashi_source_density_checks_do_not_promote_conditional_refits(audit):
    for result in audit["takahashi_2019"].values():
        assert result["solver_success"]
        assert result["coefficient_reconstruction"] == "similar"
        assert result["original_author_workflow_parity"] is False
        assert not all(
            p["within_printed_source_error_width"] for p in result["parameters"]
        )
        b = result["density_benchmark"]
        assert (
            abs(b["calculated_density_g_cm3"] - b["source_value_g_cm3"])
            < b["tolerance_g_cm3"]
        )


def test_source_checksums_and_prescher_missing_volumes_are_explicit():
    manifest = json.loads((RESEARCH / "sources.json").read_text())
    for filename, checksum in manifest["datasets"].items():
        assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == checksum
    prescher = manifest["sources"][1]
    assert prescher["year"] == 2012
    assert "unresolved" in prescher["reference_volume_status"]
    assert all("V0" not in r for r in prescher["published_bm3_fits"])
