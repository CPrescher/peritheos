"""Source transcription, phase separation and independent Nisr EOS reproduction."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_nisr_2017_silica import (
    PUBLISHED,
    ledger_outcome,
    load_rows,
    pressure,
    reproduce,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = {
    "stishovite_nisr_2017_dry_bm3": "sio2_stv_andr",
    "hydrous_stishovite_nisr_2017_bm3": "hydrous_stishovite_nisr_2017",
    "hydrous_silica_cacl2_nisr_2017_bm2": "hydrous_silica_cacl2_nisr_2017",
}


def test_complete_supplement_preserves_errors_and_coexisting_phases():
    rows = load_rows()
    assert Counter(r["table"] for r in rows) == {"S1": 20, "S2": 7, "S3": 8, "S4": 8}
    assert Counter(r["phase"] for r in rows if r["table"] == "S3") == {
        "hydrous_stishovite": 4,
        "hydrous_cacl2": 4,
    }
    assert rows[0]["volume_a3"] == "46.612"
    assert rows[20]["volume_a3"] == "47.191"
    assert rows[-1]["volume_a3"] == "40.2441"
    assert rows[-1]["lattice_c_sigma_angstrom"] == "0.000497"
    assert [
        (r["table"], r["pressure_gpa"])
        for r in rows
        if float(r["volume_sigma_a3"]) == 0
    ] == [("S3", "38.334"), ("S4", "44.857")]
    assert not any("gold" in name for name in rows[0])
    path = ROOT / "peritheos/data/datasets/nisr-2017-silica-tables-s1-s4.csv"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = json.loads((path.parent / "nisr-2017-source-manifest.json").read_text())
    assert manifest["transcription"]["sha256"] == digest
    for mid in MATERIALS.values():
        dataset = next(
            d
            for d in get_material_document(mid)["datasets"]
            if d["identifier"] == "nisr_2017_silica_tables_s1_s4"
        )
        assert dataset["resource"]["sha256"] == digest
        assert dataset["row_count"] == 43
        ambient = next(
            d
            for d in get_material_document(mid)["datasets"]
            if d["identifier"] == "nisr_2017_corrected_table2_ambient"
        )
        assert ambient["rows"][0][3:7] == [2.6647, 0.0009, 46.569, 0.011]
        assert ambient["rows"][1][5:] == [47.198, 0.021, 3.2, 0.5]


@pytest.mark.parametrize("identifier", PUBLISHED)
def test_published_coefficients_and_native_evaluation(identifier):
    record = get_eos_record(identifier)
    document = get_material_document(MATERIALS[identifier])
    metadata = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    table, v0, k0, kp = PUBLISHED[identifier]
    assert metadata["eos"]["parameters"] == (
        {"V0": v0, "K0": k0} if table == "S4" else {"V0": v0, "K0": k0, "K0_prime": kp}
    )
    assert metadata["record_kind"] == "published"
    assert metadata["pressure_calibration"]["status"] == "partially_resolved"
    assert "reference_eos_record" not in metadata["pressure_calibration"]["methods"][0]
    assert record.pressure(v0) == pytest.approx(0, abs=1e-12)
    assert record.eos.bulk_modulus(v0) == pytest.approx(k0)
    grid = np.linspace(*metadata["experimental_pressure_range_gpa"], 7)
    v = record.volume(grid, 300, check_validity=True)
    assert record.pressure(v, 300) == pytest.approx(grid, abs=1e-8)
    assert pressure(v, v0, k0, kp) == pytest.approx(grid, abs=1e-8)
    restored = Material.from_eosmat(json.loads(json.dumps(document)))
    assert restored.get_eos_record(identifier).pressure(v) == pytest.approx(
        grid, abs=1e-8
    )
    with pytest.raises(ValueError):
        record.volume(grid[-1] + 1, 300, check_validity=True)
    with pytest.raises(ValueError):
        record.volume(grid[1], 500, check_validity=True)


@pytest.mark.parametrize("identifier", PUBLISHED)
def test_independent_source_refits_and_objective_sensitivity(identifier):
    result = reproduce()["fits"][identifier]
    expected = {
        "stishovite_nisr_2017_dry_bm3": {"K0": 312.1355083},
        "hydrous_stishovite_nisr_2017_bm3": {"K0": 256.9167267},
        "hydrous_silica_cacl2_nisr_2017_bm2": {"V0": 47.2677666, "K0": 284.2176101},
    }
    assert result["pressure"]["parameters"] == pytest.approx(
        expected[identifier], rel=2e-7
    )
    assert all(r["solver_success"] for r in result.values())
    assert result["volume_sigma"]["excluded_zero_sigma_pressures_gpa"] == (
        [44.857] if "cacl2" in identifier else []
    )
    doc = get_material_document(MATERIALS[identifier])
    metadata = next(r for r in doc["eos_records"] if r["identifier"] == identifier)
    outcome = ledger_outcome(metadata)
    assert outcome["status"] == "parity"
    assert all(p["similar"] for p in outcome["parameters"])
    # Check saved independent reproduction stays current with numeric tolerance.
    saved = json.loads(
        (ROOT / "docs/data/nisr-2017-silica-reproduction.json").read_text()
    )["fits"][identifier]
    for mode in result:
        assert result[mode]["parameters"] == pytest.approx(
            saved[mode]["parameters"], rel=2e-7
        )
    # Independent high-pressure observations; tolerances are two source-data
    # volume RMS residuals, explicitly not the tiny Rietveld coordinate sigmas.
    benchmark = result["pressure"]["high_pressure_benchmark"]
    assert (
        abs(benchmark["calculated_volume_a3"] - benchmark["source_volume_a3"])
        < 2 * result["pressure"]["published_volume_rmse_a3"]
    )


def test_hydrous_identity_transition_gap_and_diffraction_fallback():
    water = reproduce()["water_content_audit"]
    assert water["formula_implied_water_wt_percent"] == pytest.approx(2.81, abs=0.01)
    assert water["corrected_table2_volume_expansion_percent"] == pytest.approx(
        1.3507, abs=0.0001
    )
    for mid in list(MATERIALS.values())[1:]:
        doc = get_material_document(mid)
        assert doc["formula"] == "Si0.954O2H0.184"
        assert doc["formula_units_per_cell"] == 2
        assert doc["atom_sites"] == []
        assert "not measured" in doc["source"]["structure"]["fallback"]
        a, b, c = (doc["lattice"][key] for key in "abc")
        for h, k, ell, d, intensity in doc["peaks"]:
            assert d == pytest.approx(
                1 / np.sqrt((h / a) ** 2 + (k / b) ** 2 + (ell / c) ** 2), abs=1e-9
            )
            assert intensity == 1  # explicitly documented uniform display marker
        with pytest.raises(ValueError):
            get_eos_record(doc["eos_records"][0]["identifier"]).volume(
                35, check_validity=True
            )
    high = get_material_document("hydrous_silica_cacl2_nisr_2017")
    assert high["source"]["structure"]["reference_pressure_gpa"] == 47.57
    assert np.prod([high["lattice"][k] for k in "abc"]) == pytest.approx(
        41.5235, abs=0.0002
    )
