"""Original Sata Fe3C observations, source fit and separate pressure bases."""

import json

import numpy as np
import pytest
from jsonschema import validate

from peritheos import (
    Material,
    eosmat_schema,
    get_eos_record,
    get_material,
    get_material_document,
    load_eosmat,
)
from peritheos.errors import MaterialError
from scripts.audit_sata_2010_fe3c import (
    DATASET_ID,
    PEAKS,
    RECORD_ID,
    REPORT,
    ROOT,
    SELECTED_ID,
    SOURCE,
    TOKENS,
    VR,
    bm3_fit,
    check_manuscript,
    data_agreement,
    digest,
    modified_bulk,
    modified_pressure,
    observations,
    pressure,
    reproduce,
    selected_rows,
)
from scripts.reconstruction_plot_data import curve_comparison


def test_all_original_tokens_rows_errors_and_shared_li_selection():
    raw, rows = observations()
    assert len(rows) == len(raw) == 9
    assert [r["measurement_id"] for r in rows] == [
        29,
        45,
        56,
        83,
        111,
        130,
        151,
        172,
        193,
    ]
    assert [r["run"] for r in rows] == [1, 1, 2, 2, 2, 2, 2, 2, 2]
    assert raw[0]["mgo_a_raw"] == "3.9321(11)"
    assert raw[-1]["pressure_raw"] == "186.60(33)"
    assert rows[-1]["error_pressure_gpa"] == 0.33
    assert rows[0]["error_mgo_volume_a3"] == 0.05
    joint = selected_rows()
    assert len(joint) == 12
    assert [r["source_row"] for r in joint[9:]] == [12, 13, 14]
    assert [r["pressure_gpa"] for r in joint[9:]] == [27.2, 30.3, 30.5]
    assert [r["volume_a3"] for r in joint[9:]] == [138.12, 136.74, 136.68]
    assert all(r["run"] is None and r["measurement_id"] is None for r in joint[9:])
    manifest = json.loads(SOURCE.read_text())
    for path in [TOKENS, PEAKS]:
        assert digest(path) == manifest["files"][path.name]["sha256"]
    assert manifest["manuscript"]["version"].startswith("2009 manuscript")


def test_literal_mbm3_and_published_bm3_with_independent_joint_refit(
    assert_audit_close,
):
    report = reproduce()
    archived = json.loads(REPORT.read_text())
    assert_audit_close(report, archived, ignore_keys={"effective_variance"})
    current_effective = report["modified_bm3_replays"]["effective_variance"]
    archived_effective = archived["modified_bm3_replays"]["effective_variance"]
    # The finite-difference effective-variance optimum varies slightly across
    # supported SciPy versions: bound its coefficients and native GPa residuals
    # separately, retaining exact structure, decisions and archived values.
    assert_audit_close(
        current_effective,
        archived_effective,
        ignore_keys={"parameters_Pr_Kr_Kprime", "residual_gpa"},
    )
    np.testing.assert_allclose(
        current_effective["parameters_Pr_Kr_Kprime"],
        archived_effective["parameters_Pr_Kr_Kprime"],
        rtol=3e-5,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        current_effective["residual_gpa"],
        archived_effective["residual_gpa"],
        rtol=1e-5,
        atol=1e-5,
    )
    assert VR == pytest.approx(149.456)
    volume = np.array([r["volume_a3"] for r in selected_rows()])
    np.testing.assert_allclose(
        modified_pressure(volume, [0, 290, 3.76]), pressure(volume), atol=2e-12
    )
    assert modified_pressure(VR, [1.0, 290.0, 3.76]) == pytest.approx(1)
    assert modified_bulk(VR, [1.0, 290.0, 3.76]) == pytest.approx(290)
    fitted = report["modified_bm3_replays"]["equal_pressure"]
    np.testing.assert_allclose(
        fitted["parameters_Pr_Kr_Kprime"],
        [0.01170009, 290.407438, 3.76370512],
        atol=2e-5,
    )
    assert np.all(
        np.abs(np.array(fitted["parameters_Pr_Kr_Kprime"]) - [0, 290, 3.76])
        < [1.6, 13, 0.18]
    )
    assert report["source_replay"]["rmse_gpa"] == pytest.approx(1.10673314549)
    fixed = bm3_fit()
    np.testing.assert_allclose(
        fixed["parameters"], [290.49452790647, 3.76249517105], atol=1e-9
    )
    assert fixed["dof"] == 10
    assert np.all(np.linalg.eigvalsh(fixed["covariance"]) > 0)
    assert report["source_only_modified_bm3"]["parameters_Pr_Kr_Kprime"][1] > 296


def test_mgo_replay_retains_the_unresolved_source_reduction_difference():
    report = reproduce()["mgo_pressure_replay"]
    assert report["max_abs_gpa"] == pytest.approx(0.35263898908)
    assert report["parameters"] == {"V0": 74.71, "K0": 160.2, "K0_prime": 3.99}
    v = np.array([r["mgo_volume_a3"] for r in observations()[1]])
    native = get_eos_record("mgo_speziale_2001_bm3_2").pressure(v, 300)
    np.testing.assert_allclose(native, report["replayed_pressure_gpa"], atol=2e-10)
    assert observations()[1][-1]["pressure_gpa"] == 186.6
    assert report["replayed_pressure_gpa"][-1] != pytest.approx(186.6, abs=0.01)


def test_poor_data_agreement_preserves_conditional_parameter_parity():
    rows = selected_rows()
    agreement = data_agreement(rows)
    assert agreement["status"] == "poor"
    assert agreement["sata_9"]["outside_pressure_width"] == 7
    assert agreement["sata_9"]["outside_combined_width"] == 6
    assert agreement["all_12"]["outside_combined_width"] == 6
    worst = agreement["rows"][3]
    assert worst["pressure_residual_gpa"] == pytest.approx(-2.591676960)
    assert worst["combined_coordinate_width_gpa"] == pytest.approx(0.287225471)
    assert worst["absolute_residual_over_combined_width"] == pytest.approx(
        9.023144616, abs=1e-8
    )
    # Check the P/V propagation against an independently differenced slope.
    for row, calculated in zip(rows, agreement["rows"]):
        v, h = row["volume_a3"], 1e-4
        slope = (pressure(v + h) - pressure(v - h)) / (2 * h)
        width = np.hypot(row["error_pressure_gpa"], slope * row["error_volume_a3"])
        assert calculated["combined_coordinate_width_gpa"] == pytest.approx(
            width, rel=1e-8
        )


def test_installed_source_defaults_limits_structure_and_litasov_separation():
    material = get_material("cementite_nonmagnetic")
    assert (
        material.default_record().identifier
        == "cementite_nonmagnetic_takahashi_2019_kbr_bm3_mgd"
    )
    document = get_material_document("cementite_nonmagnetic")
    validate(document, eosmat_schema())
    r = next(r for r in document["eos_records"] if r["identifier"] == RECORD_ID)
    assert r["record_kind"] == "published" and not r["default"]
    assert r["eos"]["parameters"] == {"V0": VR, "K0": 290.0, "K0_prime": 3.76}
    assert r["parameter_errors"] == {"V0": None, "K0": 13.0, "K0_prime": 0.18}
    assert r["fixed_parameters"] == []
    assert r["published_reference_fit"]["Pr_error_gpa"] == 1.6
    assert r["pressure_calibration"]["recalculation"]["status"] == "not_possible"
    assert r["scientific_validation"]["status"] == "primary_source_validated"
    assert r["reproduction"]["assessment"]["scope"] == "parameter_agreement"
    assert r["reproduction"]["assessment"]["status"] == "parity"
    assert r["reproduction"]["assessment"]["fit"]["fixed_parameters"] == ["V0"]
    validation = r["scientific_validation"]
    assert validation["data_agreement_status"] == "poor"
    assert validation["usage_recommendation"] == "qualified_use"
    assert (
        "point errors are smaller than the overall experimental scatter"
        in validation["warning"]
    )
    assert r["reproduction"]["assessment"]["fit"]["data_agreement_status"] == "poor"
    assert not any("litasov" in r["identifier"] for r in document["eos_records"])
    source = get_eos_record(RECORD_ID)
    p = np.linspace(27.2, 186.6, 11)
    v = source.volume(p, 300, check_validity=True)
    np.testing.assert_allclose(source.pressure(v, 300), p, atol=1e-9)
    np.testing.assert_allclose(source.pressure(v, 300), pressure(v), atol=2e-10)
    for p_bad, t_bad in [(27.19, 300), (186.61, 300), (100, 298), (100, 1000)]:
        with pytest.raises(MaterialError):
            source.volume(p_bad, t_bad, check_validity=True)
    exported = material.to_eosmat()
    assert Material.from_eosmat(exported).to_eosmat() == exported
    selected = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    exported_source = selected.to_eosmat()
    assert (
        Material.from_eosmat(
            exported_source, record_identifiers=[RECORD_ID]
        ).to_eosmat()
        == exported_source
    )
    assert len(selected.get_dataset(DATASET_ID)) == 9
    assert len(selected.get_dataset(SELECTED_ID)) == 12
    assert np.all(np.isnan(selected.get_dataset(SELECTED_ID).values("run")[9:]))
    peaks = next(
        d
        for d in document["datasets"]
        if d["identifier"] == "cementite_sata_2010_table4_diffraction"
    )
    assert peaks["used_by_eos_records"] == [] and peaks["provenance_only"]
    assert (
        peaks["rows"][1][-1] == -0.00046
    )  # Retain source delta, not recomputed -0.00045.


def test_generated_curve_adapter_uses_all_twelve_original_scale_rows():
    document = load_eosmat(
        ROOT / "peritheos/data/materials/cementite_nonmagnetic.eosmat"
    )
    record = next(r for r in document["eos_records"] if r["identifier"] == RECORD_ID)
    comparison, reason = curve_comparison(document, record)
    assert comparison is not None, reason
    assert len(comparison.series.pressure) == 12
    assert comparison.source_rmse == pytest.approx(
        reproduce()["source_replay"]["rmse_gpa"]
    )
    assert comparison.refit_rmse == pytest.approx(bm3_fit()["rmse_gpa"])


def test_optional_manuscript_correspondence():
    path = ROOT / "output/sata-2010-fe3c-2026-11-01/sources/sata-2010-manuscript.pdf"
    if not path.is_file():
        pytest.skip("Local source manuscript is not package data")
    assert check_manuscript(path)["numerical_tokens"] == 72
