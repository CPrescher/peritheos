"""Primary observations, model convention, structure and source-curve checks."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator
from scipy.constants import Avogadro

from peritheos import Material, get_material, get_material_document
from peritheos.eosmat import eosmat_schema, load_eosmat, save_eosmat
from peritheos.errors import EosValidationError, MaterialError
from scripts.recover_rubidium_halides_farla_2025 import (
    benchmark_rows,
    observation_rows,
    write_csv,
)
from scripts.reproduce_rubidium_halides_farla_2025 import (
    DATA,
    PREFIX,
    SOURCE,
    benchmark_diagnostics,
    debye_anchor_check,
    independent_fit,
    reference_pressure,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("compound", SOURCE)
def test_fixed_debye_temperatures_match_supplemental_elastic_inputs(compound):
    result = debye_anchor_check(compound)
    assert result["calculated_theta_k"] == pytest.approx(
        result["source_fixed_theta_k"], abs=0.04
    )
    assert result["calculated_average_speed_m_s"] == pytest.approx(
        result["source_average_speed_m_s"], abs=0.4
    )


@pytest.mark.parametrize("compound,count", [("rbcl", 125), ("rbbr", 126), ("rbi", 125)])
def test_complete_source_table_transcription_and_checksum(compound, count, tmp_path):
    raw = json.loads(
        (DATA / f"{PREFIX}-source-tables.json").read_text(encoding="utf-8")
    )
    assert raw["source_sha256"] == (
        "2c56d355b97e7fd9dfdc1062e32e0c01eaccef2927e51a92e38ee5de4dc430bf"
    )
    assert len(raw["tables"]) == 7
    rows = observation_rows(raw["tables"], compound)
    assert len(rows) == count
    assert sum(r["source_table"] == "S1" for r in rows) == count - 109
    assert {r["run"] for r in rows if r["source_table"] != "S1"} == {
        "BT792",
        "BT793",
        "BT795",
        "BT815",
    }
    for kind, normalized in (
        ("observations", rows),
        ("benchmarks", benchmark_rows(raw["tables"], compound)),
    ):
        restored = tmp_path / f"{kind}.csv"
        write_csv(restored, normalized)
        assert (
            restored.read_bytes()
            == (DATA / f"{PREFIX}-{compound}-{kind}.csv").read_bytes()
        )
    document = get_material_document(compound)
    material = get_material(compound)
    for metadata in document["datasets"]:
        if "farla_2025" not in metadata["identifier"]:
            continue
        path = ROOT / "peritheos/data" / metadata["resource"]["path"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == metadata["resource"]["sha256"]
        )
        dataset = material.get_dataset(metadata["identifier"])
        assert dataset.checksum_verified
        assert len(dataset) == (54 if "benchmarks" in metadata["identifier"] else count)


def test_source_irregularities_and_calibrants_are_retained():
    raw = json.loads(
        (DATA / f"{PREFIX}-source-tables.json").read_text(encoding="utf-8")
    )
    assert raw["tables"][0][3][4] == "1.19(0"
    br = observation_rows(raw["tables"], "rbbr")
    assert br[0]["pressure_gpa_source_token"] == "1.19(0"
    assert float(br[0]["pressure_uncertainty_gpa"]) == 0
    assert float(br[-18]["cscl_pressure_gpa"]) == 15.17  # BT795 end at 1800 K
    cl = observation_rows(raw["tables"], "rbcl")
    assert cl[0]["volume_a3"] == "58.6949"
    assert float(cl[0]["volume_uncertainty_a3"]) == pytest.approx(0.0050)
    first = next(r for r in cl if r["run"] == "BT792")
    assert first["cscl_volume_a3_source_token"] == "61.810(12)"
    assert first["pt_volume_a3_source_token"] == "59.847(12)"
    assert first["mo_volume_a3_source_token"] == "30.824(5)"
    assert first["temperature_uncertainty_k"] == "30"
    iodine = observation_rows(raw["tables"], "rbi")
    missing = [r for r in iodine if not r["volume_a3"]]
    assert len(missing) == 1
    assert missing[0]["temperature_k"] == "500"
    assert missing[0]["cscl_pressure_gpa"] == "7.04"
    assert missing[0]["mo_volume_a3"] == ""
    for compound in SOURCE:
        bench = benchmark_rows(raw["tables"], compound)
        assert len(bench) == 54
        assert sum(r["molten_source_flag"] for r in bench) == 1
        assert bench[-1]["pressure_gpa"] == "100"


@pytest.mark.parametrize("compound", SOURCE)
def test_schema_loading_reference_state_structure_and_metadata_round_trip(
    compound, tmp_path
):
    document = get_material_document(compound)
    Draft202012Validator(eosmat_schema()).validate(document)
    material = Material.from_eosmat(document)
    record = material.get_eos_record(f"{compound}_b2_farla_2025_bm3_mgd_cscl")
    assert record.determination_method == "hybrid"
    molar, theta, parameters = SOURCE[compound]
    assert record.reference_volume == pytest.approx(molar * 1e24 / Avogadro)
    assert record.pressure(record.reference_volume, 300) == pytest.approx(0, abs=1e-12)
    assert record.thermal_pressure_increment(record.reference_volume * 0.8, 300) == 0
    assert record.eos.theta0 == theta
    assert record.eos.n == 2
    volume = record.reference_volume * 0.8
    gamma = parameters[2] * 0.8 ** parameters[3]
    expected_theta = theta * np.exp((parameters[2] - gamma) / parameters[3])
    assert record.eos.characteristic_temperature(
        volume * record.volume_scale
    ) == pytest.approx(expected_theta)
    if compound != "rbcl":
        assert material.space_group_number == 221
        assert material.formula_units_per_cell == 1
        assert sum(s["multiplicity"] * s["occupancy"] for s in material.atom_sites) == 2
        assert material.lattice["a"] ** 3 == pytest.approx(record.reference_volume)
    else:
        assert material.default_record().identifier == "rbcl_b2_campbell_1994_bm3_1"
    restored = Material.from_eosmat(material.to_eosmat())
    assert restored.eosmat_metadata == material.eosmat_metadata
    restored_metadata = restored.get_eos_record(record.identifier).eosmat_metadata
    for key, value in record.eosmat_metadata.items():
        # Executable components gain normalized configuration/volume scales;
        # preserve all source provenance and scientific extension metadata.
        if key not in {"eos", "thermal", "volume"}:
            assert restored_metadata[key] == value
    assert restored.get_eos_record(record.identifier).pressure(
        50, 1000
    ) == pytest.approx(record.pressure(50, 1000), abs=1e-12)
    path = tmp_path / f"{compound}.eosmat"
    save_eosmat(path, document)
    assert load_eosmat(path) == document


@pytest.mark.parametrize("compound", SOURCE)
def test_all_source_checkpoints_native_equation_and_coefficient_precision(compound):
    audit = benchmark_diagnostics(compound)
    assert audit["rows"] == 54
    assert audit["max_native_reference_ratio_difference"] < 2e-12
    assert audit["outside_coefficient_rounding_envelope"] == 0
    # Preserve the scientific qualification: there is no universal last-digit parity.
    assert audit["outside_last_digit_rounding"] > 0
    in_range = [
        r
        for r in audit["states"]
        if r["pressure_gpa"] <= 20 and r["temperature_k"] <= 1500
    ]
    assert max(abs(r["difference"]) for r in in_range) < 0.00025


@pytest.mark.parametrize("compound", SOURCE)
def test_pressure_roundtrips_broadcasting_and_validity(compound):
    record = get_material(compound).get_eos_record(
        f"{compound}_b2_farla_2025_bm3_mgd_cscl"
    )
    p = np.array([[5.0], [10.0], [20.0]])
    t = np.array([300.0, 1000.0, 1800.0])
    v = record.volume(p, t, check_validity=True)
    assert v.shape == (3, 3)
    np.testing.assert_allclose(
        record.pressure(v, t, check_validity=True),
        np.broadcast_to(p, v.shape),
        atol=1e-9,
    )
    molar, theta, parameters = SOURCE[compound]
    np.testing.assert_allclose(
        reference_pressure(v, t, molar, theta, parameters),
        p + np.zeros_like(t),
        atol=1e-10,
    )
    with pytest.raises(MaterialError):
        record.volume(30, 1000, check_validity=True)
    with pytest.raises(MaterialError):
        record.volume(10, 2000, check_validity=True)
    with pytest.raises(EosValidationError):
        record.pressure(-1, 300)


@pytest.mark.parametrize("compound,count", [("rbcl", 122), ("rbbr", 123), ("rbi", 121)])
def test_independent_refit_discloses_selection_and_parity_limit(compound, count):
    fitted = independent_fit(compound)
    assert fitted["success"]
    assert fitted["rows"] == count
    assert not fitted["author_fit_parity"]
    assert not fitted["author_workflow_parity"]
    assert fitted["source_pressure_rmse_gpa"] < 0.17
    assert fitted["source_pressure_max_abs_gpa"] < 0.47
    assert fitted["fitted_pressure_rmse_gpa"] < fitted["source_pressure_rmse_gpa"]


def test_source_manifest_checks_every_bundled_product():
    manifest = json.loads((DATA / f"{PREFIX}-source.json").read_text(encoding="utf-8"))
    for product in manifest["products"]:
        path = ROOT / "peritheos/data" / product["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == product["sha256"]
    assert manifest["public_supplement"]["license"] == "CC BY 4.0"
