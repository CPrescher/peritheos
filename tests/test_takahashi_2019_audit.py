"""Independent source-equation and objective checks for the scoped audit."""

import csv
import json

import numpy as np
import pytest
from scipy.constants import Avogadro, R
from scipy.integrate import quad

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from scripts.audit_takahashi_2019 import (
    DATA,
    ROOT,
    SOURCE,
    WIDTH,
    effective_variance_fit,
    pressure,
    solve,
    volume,
    workbook_audit,
)


def observations(scale):
    rows = list(
        csv.DictReader(
            (DATA / "cementite-takahashi-2019-table-s1.csv").read_text().splitlines()
        )
    )
    arrays = [
        np.array([float(r[f]) for r in rows])
        for f in ("volume_a3", "temperature_k", "volume_a3_error")
    ]
    v, t, sv = arrays
    p = np.array(
        [float(r[f"pressure_{scale}_gpa"] or r["pressure_sio2_gpa"]) for r in rows]
    )
    sp = np.array(
        [
            float(r[f"pressure_{scale}_gpa_error"] or r["pressure_sio2_gpa_error"])
            for r in rows
        ]
    )
    st = np.array([float(r["temperature_k_error"] or 0) for r in rows])
    return rows, v, t, p, sp, sv, st


def test_original_workbook_values_formats_and_error_widths():
    rows, *_ = observations("kbr")
    audit = workbook_audit(rows)
    for name in (
        "raw_cell_value_mismatches",
        "raw_cell_number_format_mismatches",
        "csv_central_value_mismatches",
        "csv_error_value_mismatches",
    ):
        assert audit[name] == []
    assert audit["formula_count"] == 0
    assert audit["cell_P5"]["number_format"] == "0.00"
    assert float(rows[1]["volume_a3_error"]) == 0.08


@pytest.mark.parametrize("scale", ["kbr", "mgo"])
def test_quadrature_si_units_and_both_peritheos_paths(scale):
    _, v, t, *_ = observations(scale)
    x = SOURCE[scale]
    factor = Avogadro * 1e-25 / 4
    model = MieGruneisenDebye(
        BM3(V0=x[0] * factor, K0=x[1], K0_prime=x[2]),
        Tr=300,
        theta0=314,
        gamma0=x[3],
        q=x[4],
        n=4,
    )
    direct = pressure(v, t, x)
    np.testing.assert_allclose(model.pressure(v * factor, t), direct, rtol=1e-11)
    model._native = None
    np.testing.assert_allclose(model.pressure(v * factor, t), direct, rtol=1e-11)
    np.testing.assert_allclose(volume(direct, t, x), v, atol=1e-9)
    vv, tt = 94, 6000
    gamma = x[3] * (vv / x[0]) ** x[4]
    theta = 314 * np.exp((x[3] - gamma) / x[4])

    def energy(temp):
        return (
            36
            * R
            * temp
            * (temp / theta) ** 3
            * quad(lambda a: a**3 / np.expm1(a), 0, theta / temp)[0]
        )

    thermal = gamma * (energy(tt) - energy(300)) / (vv * 1e-30 * Avogadro / 4) / 1e9
    assert abs(pressure(vv, tt, x) - pressure(vv, 300, x) - thermal) < 1e-10


@pytest.mark.parametrize("scale", ["kbr", "mgo"])
def test_coordinate_width_objective_changes_coefficient_parity(scale):
    _, v, t, p, sp, sv, st = observations(scale)
    source = SOURCE[scale]
    base = solve(lambda x: pressure(v, t, x) - p, source)
    weighted, iterations = effective_variance_fit(v, t, p, sp, sv, st, source)
    assert base.success and weighted.success
    assert iterations < 30
    assert not np.all(np.abs(base.x - source) <= WIDTH[scale])
    assert np.all(np.abs(weighted.x - source) <= WIDTH[scale])
    if scale == "kbr":
        # All central values and unscaled inverse-normal errors recover the
        # printed decimal precision, a stronger property than width parity.
        _, singular, vectors = np.linalg.svd(weighted.jac, full_matrices=False)
        scaled = vectors / singular[:, None]
        covariance = np.einsum("ij,ik->jk", scaled, scaled)
        errors = np.sqrt(np.diag(covariance))
        decimals = [1, 1, 2, 2, 2]
        np.testing.assert_allclose(
            [round(x, d) for x, d in zip(weighted.x, decimals)], source
        )
        np.testing.assert_allclose(
            [round(e, d) for e, d in zip(errors, decimals)], WIDTH[scale]
        )
    assert len(v) == 60 and max(t) == 2600
    np.testing.assert_array_equal(SOURCE[scale], source)


@pytest.mark.parametrize("scale", ["kbr", "mgo"])
def test_catalog_and_ledger_parity_keep_uncertainty_and_workflow_separate(scale):
    from scripts.reproduce_cementite import takahashi_ledger_outcome

    document = get_material_document("cementite_nonmagnetic", format_version=4)
    identifier = f"cementite_nonmagnetic_takahashi_2019_{scale}_bm3_mgd"
    record = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    ledger = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    row = next(r for r in ledger["records"] if r["record_identifier"] == identifier)
    replay = takahashi_ledger_outcome(record)
    assessment = record["evidence"]["reproduction"]["assessment"]
    assert replay["status"] == row["status"] == assessment["status"] == "parity"
    assert assessment["scope"] == "parameter_agreement"
    assert row["original_author_workflow_parity"] is False
    assert (
        record["evidence"]["reproduction"]["original_author_workflow_parity"] is False
    )
    assert row["absolute_sigma"] is True
    assert "Unscaled inverse normal" in row["covariance_convention"]
    assert "source confidence" in row["qualification"]
    expected = (
        "reproduced_at_printed_precision"
        if scale == "kbr"
        else "not_reproduced_at_printed_precision"
    )
    assert row["uncertainty_reproduction"]["status"] == expected
    assert (
        record["evidence"]["reproduction"]["uncertainty_reproduction"]
        == row["uncertainty_reproduction"]
    )
    assert row["previous_equal_pressure_reconstruction"]["status"] == "similar"
    assert assessment["parameters"] == row["parameters"]
    assert all(p["within_printed_source_error_width"] for p in row["parameters"])
    np.testing.assert_allclose(
        [p["refit"] for p in row["parameters"]],
        [p["refit"] for p in replay["parameters"]],
        rtol=1e-5,
    )
    np.testing.assert_allclose(
        np.diag(row["parameter_covariance"]),
        [p["refit_error"] ** 2 for p in row["parameters"]],
    )
    published = [record["parameters"][n]["value"] for n in ["V0", "K0", "K0_prime"]]
    published += [record["thermal"]["parameters"][n]["value"] for n in ["gamma0", "q"]]
    np.testing.assert_array_equal(published, SOURCE[scale])
