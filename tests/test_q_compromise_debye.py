"""Independent source-law, native fitting, uncertainty and interchange checks."""

import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    convert_eosmat,
    eosmat_schema,
    get_material,
    get_material_document,
)
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import QCompromiseDebye
from peritheos.errors import EosValidationError, UnsupportedOperationError
from peritheos.fitting import fit_thermal_eos
from peritheos.uncertainty import EOSUncertainty
from scripts.audit_garnet_angel_2022 import independent_pressure, reproduce, source

ROOT = Path(__file__).resolve().parents[1]


def model():
    p, _ = source("pyrope")
    return QCompromiseDebye(BM3(p[1] / 10, p[2], p[3]), 298, p[11], p[18], 20)


@pytest.mark.parametrize("name", ["pyrope", "grossular"])
def test_independent_quadrature_and_source_table_expansivity(name):
    report = reproduce()["records"][name]
    assert report["maximum_independent_difference_gpa"] < 2e-9
    # Table 1: reported alpha errors are 5e-8 and 3e-8 K^-1.
    assert report["alpha_298k_inverse_k"] == pytest.approx(
        report["table1_alpha_inverse_k"], abs=5e-8
    )
    p, _ = source(name)
    record = get_material(name).get_eos_record(
        f"{name}_angel_2022_bm3_q_compromise_debye"
    )
    cell = p[1] * 8e24 / 6.02214076e23 * 0.95
    assert record.pressure(cell, 1000) == pytest.approx(
        independent_pressure(p, p[1] * 0.95, 1000), abs=2e-9
    )
    assert record.volume(record.pressure(cell, 1000), 1000) == pytest.approx(
        cell, rel=1e-9
    )
    raw = get_material_document(name)
    for version in [3, 4]:
        document = convert_eosmat(raw, version=version)
        Draft202012Validator(eosmat_schema(version)).validate(document)
        restored = Material.from_eosmat(document)
        assert restored.get_eos_record(record.identifier).pressure(
            cell, 1000
        ) == pytest.approx(record.pressure(cell, 1000), abs=1e-10)
    assert record.parameter_covariance is not None
    assert np.linalg.eigvalsh(record.parameter_covariance)[0] > 0


def test_constant_thermal_pressure_derivative_shapes_and_inversion():
    m = model()
    v = m.rt_eos.V0 * np.array([0.85, 0.95, 1.0])[:, None]
    t = np.array([90.0, 298.0, 1000.0, 1700.0])[None, :]
    pressure = m.pressure(v, t)
    assert pressure.shape == (3, 4)
    assert np.isscalar(m.pressure(float(v[0, 0]), 1000))
    np.testing.assert_allclose(
        m.thermal_pressure(v, t),
        np.broadcast_to(m.thermal_pressure(m.rt_eos.V0, t), (3, 4)),
        atol=1e-12,
    )
    np.testing.assert_allclose(m.thermal_pressure(v, 298), 0, atol=0)
    np.testing.assert_allclose(
        m.bulk_modulus(v, t),
        np.broadcast_to(m.rt_eos.bulk_modulus(v), (3, 4)),
        rtol=2e-8,
    )
    np.testing.assert_allclose(
        m.volume(pressure, t), np.broadcast_to(v, (3, 4)), rtol=1e-9
    )
    np.testing.assert_allclose(
        m.temperature(pressure, v), np.broadcast_to(t, (3, 4)), atol=2e-5
    )
    fallback = copy.copy(m)
    del fallback._native
    np.testing.assert_allclose(fallback.pressure(v, t), pressure, atol=2e-9)
    for fn in [m.molar_heat_capacity_v, m.molar_heat_capacity_p]:
        with pytest.raises(UnsupportedOperationError):
            fn(m.rt_eos.V0, 1000)


@pytest.mark.parametrize(
    "name,value", [("Tr", 0), ("theta0", -1), ("n", 0), ("gamma0", float("nan"))]
)
def test_invalid_parameters(name, value):
    kwargs = dict(Tr=298, theta0=771.096, gamma0=1.1847, n=20)
    kwargs[name] = value
    with pytest.raises(EosValidationError):
        QCompromiseDebye(model().rt_eos, **kwargs)


@pytest.mark.parametrize(
    "volume,temperature", [(0, 300), (-1, 300), (1, 0), (1, float("inf"))]
)
def test_invalid_states(volume, temperature):
    with pytest.raises(EosValidationError):
        model().pressure(volume, temperature)


def test_native_fit_and_gamma_uncertainty_match_independent_sensitivity():
    m = model()
    p, _ = source("pyrope")
    v, t = np.meshgrid(p[1] * np.array([0.85, 0.95, 1]), [400.0, 700.0, 1000.0, 1500.0])
    pressure = np.array(
        [
            independent_pressure(p, float(a), float(b))
            for a, b in zip(v.ravel(), t.ravel())
        ]
    )
    fit = fit_thermal_eos(
        QCompromiseDebye,
        m.rt_eos,
        volume=v.ravel() / 10,
        temperature=t.ravel(),
        pressure=pressure,
        initial={"gamma0": 1.1},
        fixed={"Tr": 298, "theta0": p[11], "n": 20},
    )
    assert fit.parameters["gamma0"] == pytest.approx(p[18], abs=1e-8)
    sigma = 0.012
    prediction = EOSUncertainty(m, parameter_errors={"gamma0": sigma}).pressure(
        m.rt_eos.V0, 1000
    )
    sensitivity = independent_pressure(p, p[1], 1000) / p[18]
    assert prediction.standard_error == pytest.approx(
        abs(sensitivity) * sigma, rel=1e-7
    )


def test_source_and_structural_preservation():
    current = get_material_document("pyrope")
    native = convert_eosmat(current, version=4)
    assert (
        hashlib.sha256(
            json.dumps(native["eos_records"][0], sort_keys=True).encode()
        ).hexdigest()
        == "a1ef6a856836ca36f9940814beffde1b8a23c61379d0bcf27ccc95bdf152850c"
    )
    assert native["atom_sites"][-1]["x"] == 0.0329
    gross = get_material_document("grossular")
    assert gross["formula_units_per_cell"] == 8
    assert gross["lattice"]["a"] == 11.8418
    assert gross["atom_sites"][-1]["y"] == 0.5453
    assert gross["atom_sites"][-1]["z"] == 0.3487


@pytest.mark.parametrize("name", ["pyrope", "grossular"])
def test_source_covariance_propagates_in_public_cell_units(name):
    p, full = source(name)
    indices = [1, 2, 10, 17]
    covariance = full[np.ix_(indices, indices)]
    record = get_material(name).get_eos_record(
        f"{name}_angel_2022_bm3_q_compromise_debye"
    )
    np.testing.assert_array_equal(record.parameter_covariance, covariance)
    volume = p[1] * 0.95
    gradient = []
    for index in [2, 3, 11, 18]:
        step = abs(p[index]) * 1e-5
        plus, minus = dict(p), dict(p)
        plus[index] += step
        minus[index] -= step
        gradient.append(
            (
                independent_pressure(plus, volume, 1000)
                - independent_pressure(minus, volume, 1000)
            )
            / (2 * step)
        )
    expected = np.sqrt(np.asarray(gradient) @ covariance @ np.asarray(gradient))
    cell_volume = volume * 8e24 / 6.02214076e23
    prediction = record.pressure_with_uncertainty(cell_volume, 1000)
    assert prediction.standard_error == pytest.approx(expected, rel=1e-6)


def test_preservation_receipt_rejects_changes_to_new_or_old_scientific_values():
    from scripts.scientific_source_snapshot import (
        assert_snapshot_lineage,
        scientific_snapshot,
    )

    proof = json.loads(
        (ROOT / "docs/data/garnet-angel-2022/integration-preservation.json").read_text()
    )
    snapshot = scientific_snapshot()
    assert_snapshot_lineage(snapshot, proof["before_scientific_sha256"])
    for material in ["grossular", "pyrope"]:
        changed = copy.deepcopy(snapshot)
        document = next(d for d in changed["documents"] if d["identifier"] == material)
        document["eos_records"][0]["parameters"]["K0"]["value"] += 1
        with pytest.raises(AssertionError):
            assert_snapshot_lineage(changed, proof["before_scientific_sha256"])


def test_additive_source_binding_recovers_the_frozen_pyrope_document():
    from scripts.merge_source_bindings import accepts_merge_binding

    proof = json.loads(
        (ROOT / "docs/data/garnet-angel-2022/integration-preservation.json").read_text()
    )
    path = "peritheos/data/materials/pyrope.eosmat"
    binding = proof["source_bindings"]["release_document"][path]
    # Frozen release receipts bind the explicit source document, including
    # redundant unit declarations that a compact conversion may omit.
    current = json.loads((ROOT / path).read_text())
    for record in current["eos_records"]:
        record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)

    def fingerprint(document):
        return hashlib.sha256(
            json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()

    assert fingerprint(current) == binding["current_sha256"]
    assert current["eos_records"].pop()["identifier"] == (
        "pyrope_angel_2022_bm3_q_compromise_debye"
    )
    assert fingerprint(current) == binding["predecessor_sha256"]
    assert accepts_merge_binding(
        "release_document",
        path,
        binding["historical_sha256"],
        binding["current_sha256"],
    )
    assert not accepts_merge_binding(
        "release_document", path, binding["historical_sha256"], "0" * 64
    )
