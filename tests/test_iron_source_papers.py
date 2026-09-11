"""Primary-source import, equations, and independently published checkpoints."""

import hashlib
import json
from collections import Counter

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import Dewaele2006
from peritheos.errors import EosValidationError
from peritheos.materials import Material
from scripts.reproduce_iron_source_papers import (
    DATA,
    REPORT,
    belongs,
    pressure,
    reproduce,
    rows,
)


def source_records():
    for name in ("iron", "fe09ni01_hcp", "nacl_b2"):
        document = get_material_document(name)
        for r in document["eos_records"]:
            if belongs(r["identifier"]):
                yield document, r


def test_complete_primary_tables_and_unit_cell_geometry():
    y = rows("iron-yamazaki-2012-table-s1")
    s = rows("sakai-2014-table2-pvt")
    b = rows("iron-brown-2000-tables-i-ii")
    assert (len(y), len(s), len(b)) == (207, 104, 37)
    assert Counter(x["material"] for x in s) == {"iron": 27, "fe09ni01_hcp": 77}
    assert sum(float(x["temperature_k"]) > 300 for x in s) == 32
    assert y[0]["run"] == "M1067" and y[0]["pressure_gpa"] == "42.8"
    assert y[-1]["run"] == "M1230"
    for r in y + s:
        calculated = (
            np.sqrt(3) / 2 * float(r["a_angstrom"]) ** 2 * float(r["c_angstrom"])
        )
        # Printed a/c and V are separately rounded, including some 3-decimal a.
        assert calculated == pytest.approx(
            float(r["volume_a3_conventional_cell"]), abs=0.025
        )
    first = next(x for x in s if x["run"] == "FNC02_004a")
    assert first["pressure_p1_gpa"] == "24.6"
    assert first["pressure_p4_gpa"] == ""
    assert first["nacl_volume_a3"] == "120.9"
    hot = next(x for x in s if x["run"] == "F10N09_024")
    assert hot["pressure_p7_gpa"] == "123.5" and hot["pressure_p1_gpa"] == ""
    assert hot["mgo_volume_a3"] == "52.53"
    assert b[-1]["pressure_gpa"] == "442.1"
    assert b[-1]["density_g_cm3_error"] == "0.11"
    for document, _ in source_records():
        for ds in document["datasets"]:
            if ds["identifier"] in (
                "iron_yamazaki_2012_table_s1",
                "sakai_2014_table2_pvt",
                "iron_brown_2000_tables_i_ii",
            ):
                path = DATA / ds["resource"]["path"].split("/")[-1]
                assert (
                    hashlib.sha256(path.read_bytes()).hexdigest()
                    == ds["resource"]["sha256"]
                )
                assert list(rows(path.stem)[0]) == [c["name"] for c in ds["columns"]]


def test_all_source_equations_native_and_roundtrip():
    records = list(source_records())
    assert len(records) == 43
    for document, r in records:
        material = Material.from_eosmat(document)
        selected = material.get_eos_record(r["identifier"])
        v = (
            np.array([19.0, 18.0, 16.0])
            if r["equation_kind"] == "hugoniot"
            else np.array([18.0, 16.0, 14.0])
        )
        t = np.array([300.0, 1500.0, 2300.0]) if r.get("thermal") else 300.0
        expected = pressure(r, v, t)
        actual = selected.pressure(v, t) if r.get("thermal") else selected.pressure(v)
        assert actual == pytest.approx(expected, rel=2e-10)
        restored = Material.from_eosmat(material.to_eosmat()).get_eos_record(
            r["identifier"]
        )
        assert (
            restored.pressure(v, t) if r.get("thermal") else restored.pressure(v)
        ) == pytest.approx(actual)
        assert r["record_kind"] == "published" and "default_for" not in r
    alloy = get_material_document("fe09ni01_hcp")
    assert [(s["element"], s["occupancy"]) for s in alloy["atom_sites"]] == [
        ("Fe", 0.9),
        ("Ni", 0.1),
    ]


def test_zero_gamma_infinity_is_exact_power_law_and_rejects_negative():
    params = dict(
        Tr=300,
        theta0=1173,
        gamma0=3.2,
        gamma_inf=0.0,
        beta=0.8,
        anharmonic_a=3.7e-5,
        anharmonic_m=1.87,
        electronic_e=1.95e-4,
        electronic_g=1.339,
        n=1,
    )
    model = Dewaele2006(BM3(0.667, 202, 4.5), **params)
    vs = np.array([0.4, 0.5, 0.6])
    ratio = vs / 0.667
    assert model.gruneisen_parameter(vs) == pytest.approx(3.2 * ratio**0.8)
    assert model.characteristic_temperature(vs) == pytest.approx(
        1173 * np.exp(3.2 / 0.8 * (1 - ratio**0.8))
    )
    native = model.pressure(vs, 1500)
    del model._native
    assert model.pressure(vs, 1500) == pytest.approx(native, rel=1e-10)
    with pytest.raises(EosValidationError):
        Dewaele2006(BM3(0.667, 202, 4.5), **{**params, "gamma_inf": -0.01})


def test_source_checkpoints_and_explicit_regression_limits(assert_audit_close):
    report = reproduce()["records"]
    brown = report["iron_brown_2000_linear_hugoniot"]
    assert brown["status"] == "parity"
    assert round(brown["parameters"]["c0"], 3) == 3.935
    assert round(brown["parameters"]["s"], 3) == 1.578
    assert round(brown["quadratic_parameters"]["q"], 3) == -0.038
    assert report["iron_dewaele_2006_bm3"]["status"] == "similar"
    dub = report["iron_dubrovinsky_2000_bm3_thermal"]
    assert dub["status"] == "not_refittable"
    assert dub["checkpoint_bulk_modulus_211gpa_300k"] == pytest.approx(1110, rel=0.002)
    assert dub["checkpoint_average_alpha_202gpa_5200k"] == pytest.approx(
        9.13e-6, rel=0.01
    )
    for suffix, published in [
        ("type1", 14.02),
        ("type2", 13.87),
        ("type4_1", 13.88),
        ("type4_2", 13.87),
        ("type4_3", 13.94),
    ]:
        result = report[f"fe09ni01_hcp_sakai_2014_{suffix}_thermal"]
        # Source values printed to .01 A^3; rounding of coefficients also matters.
        assert result["checkpoint_329gpa_5000k_volume_a3"] == pytest.approx(
            published, abs=0.014
        )
        assert result["status"] == "parity_not_achieved"
    assert report["iron_yamazaki_2012_bm3_thermal"][
        "checkpoint_330gpa_6000k_density_g_cm3"
    ] == pytest.approx(13.12, abs=0.025)
    assert report["iron_yamazaki_2012_vinet_thermal"][
        "checkpoint_330gpa_6000k_density_g_cm3"
    ] == pytest.approx(13.18, abs=0.025)
    stored = json.loads(REPORT.read_text(encoding="utf-8"))
    replayed = reproduce()
    # Compare fitted parameters at solver precision. Their near-zero differences
    # from published values lose relative precision, so validate those against
    # each report's own coefficients before comparing the rest of the reports.
    for audit in (stored, replayed):
        for record in audit["records"].values():
            for comparison in record.get("parameter_comparisons", []):
                difference = comparison.pop("absolute_difference")
                assert difference == pytest.approx(
                    abs(comparison["fitted"] - comparison["published"]),
                    rel=1e-12,
                    abs=1e-12,
                )
    assert_audit_close(replayed, stored)
