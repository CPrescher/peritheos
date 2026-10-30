"""Published Fe3C coefficients, adopted reference and scoped parity stay distinct."""

import copy
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    convert_eosmat,
    eosmat_schema,
    get_eos_record,
    get_material_document,
    load_eosmat,
    save_eosmat,
)
from scripts.prescher_2012_catalog import (
    IDS,
    RETIRED_IDS,
    ROOT,
    V0,
    ledger_outcome,
    selected_rows,
)
from scripts.reconstruction_plot_data import curve_comparison
from scripts.standardize_reconstruction_metadata import assessment_for_record


@pytest.mark.parametrize(
    "side,values,errors",
    [
        ("full", [161, 5.9], [2, 0.2]),
    ],
)
def test_published_parameters_and_explicit_adoption_survive_public_loading(
    side, values, errors
):
    document = get_material_document("cementite", format_version=4)
    record = next(r for r in document["eos_records"] if r["identifier"] == IDS[side])
    assert [record["parameters"][k]["value"] for k in ("K0", "K0_prime")] == values
    assert [record["parameters"][k]["error"] for k in ("K0", "K0_prime")] == errors
    assert record["parameters"]["V0"] == {"value": V0, "fixed": True}
    assert record["parameter_error_confidence"] is None
    adoption = record["evidence"]["parameter_provenance"]["V0"]
    assert adoption["value_kind"] == "adopted"
    assert adoption["source_value_published"] is False
    assert adoption["source_constraint_recovered"] is False
    assessment = record["evidence"]["reproduction"]["assessment"]
    assert assessment["status"] == "parity"
    assert assessment["scope"] == "parameter_agreement"
    assert "uncertainty is not reproduced" in assessment["basis"].lower()
    assert "confidence" in json.dumps(assessment).lower()
    # Exact serialization checks use the saved audit, not a fresh regression:
    # supported NumPy/SciPy builds can differ in insignificant optimizer bits.
    ledger = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    outcome = next(
        row
        for row in ledger["records"]
        if row["record_identifier"] == record["identifier"]
    )
    organized = assessment_for_record(
        record,
        {
            **outcome,
            "record_identifier": record["identifier"],
            "material_identifier": "cementite",
        },
    )
    assert assessment == organized
    replay = {p["parameter"]: p for p in assessment["parameters"]}
    for name in ("K0", "K0_prime"):
        assert replay[name]["within_published_printing_precision"] is True
        assert replay[name]["within_combined_2sigma"] is None
    assert outcome["author_workflow_reproduced"] is False
    assert replay["V0"]["compared_for_parity"] is False
    assert replay["V0"]["refit"] == V0
    assert replay["K0"]["refit_error"] == pytest.approx(2.5321927447246857)
    assert replay["K0_prime"]["refit_error"] == pytest.approx(0.24086704947871293)
    assert replay["K0"]["uncertainty_matches_printed_precision"] is False
    assert replay["K0"]["uncertainty_finding"] == "not_reproduced"
    assert replay["K0_prime"]["uncertainty_matches_printed_precision"] is True
    assert assessment["basis_code"] == (
        "conditional_K0_and_K0_prime_at_published_printed_precision"
    )
    assert [replay[k]["display_decimal_places"] for k in ("K0", "K0_prime")] == [0, 1]
    assert assessment["fit"]["fixed_values"] == {"V0": V0}
    assert assessment["fit"]["degrees_of_freedom"] == 17
    model = get_eos_record(IDS[side])
    pressures = np.linspace(*record["experimental_pressure_range_gpa"], 9)
    volumes = model.volume(pressures)
    assert model.pressure(volumes) == pytest.approx(pressures, abs=1e-8)
    assert model.pressure(V0) == pytest.approx(0, abs=1e-12)
    assert model.eos.bulk_modulus(V0) == pytest.approx(values[0])


def test_reconstruction_plots_use_the_named_inputs_and_both_actual_curves():
    document = get_material_document("cementite", format_version=4)
    for side in IDS:
        record = next(
            r for r in document["eos_records"] if r["identifier"] == IDS[side]
        )
        comparison, reason = curve_comparison(document, record)
        assert comparison is not None, reason
        p, v = selected_rows(side)
        assert comparison.series.pressure == pytest.approx(p)
        assert comparison.series.volume == pytest.approx(v)
        fit = record["evidence"]["reproduction"]["assessment"]["fit"]
        assert comparison.source_rmse == pytest.approx(fit["published_rmse_gpa"])
        assert comparison.refit_rmse == pytest.approx(fit["rmse_gpa"])


def test_complete_datasets_stay_public_and_split_fits_stay_in_the_audit():
    document = get_material_document("cementite", format_version=4)
    data = Material.from_eosmat(document)
    final = data.get_dataset("cementite_prescher_2012_table_s1")
    graphic = data.get_dataset("cementite_prescher_2012_figure1_digitized")
    assert [len(final), len(graphic)] == [19, 22]
    assert final.checksum_verified
    assert final["volume_a3"][15] == 133.25
    assert graphic["volume_a3"][0] == pytest.approx(155.08711984)
    ids = {r["identifier"] for r in document["eos_records"]}
    assert not (ids & RETIRED_IDS)
    assert {r for r in ids if r.startswith("cementite_prescher_2012_")} == {IDS["full"]}
    assert not any(
        d["identifier"].endswith("working_cells_with_table_pressures")
        for d in document["datasets"]
    )
    archived = json.loads(
        (
            ROOT / "peritheos/data/research/cementite/prescher-2012-retired-splits.json"
        ).read_text()
    )
    assert {r["identifier"] for r in archived["eos_records"]} == RETIRED_IDS
    assert archived["datasets"][0]["rows"][15][2] == 133.58
    assert all(
        r["reference"]["doi"] == "10.1103/PhysRevB.85.140402"
        for r in document["eos_records"]
        if r["identifier"].startswith("cementite_prescher_2012_")
    )
    # Scott's separate source record remains first; addition does not select a new default.
    assert document["eos_records"][0]["identifier"] == "cementite_scott_2001_bm3_1"


def test_native_and_compatibility_round_trips_preserve_coefficients_and_evidence(
    tmp_path,
):
    document = get_material_document("cementite", format_version=4)
    before = copy.deepcopy(document)
    for version in (3, 4):
        converted = convert_eosmat(document, version=version)
        Draft202012Validator(eosmat_schema(version=version)).validate(converted)
        path = tmp_path / f"cementite-v{version}.eosmat"
        save_eosmat(path, converted)
        reloaded = convert_eosmat(load_eosmat(path), version=4)
        assert reloaded["eos_records"] == before["eos_records"]
        assert reloaded["datasets"] == before["datasets"]
    assert document == before


def test_adopted_v0_is_the_saved_drawn_curve_recovery_not_an_ambient_measurement():
    recovery = json.loads(
        (ROOT / "docs/data/prescher-2012-fe3c-recovery.json").read_text()
    )
    assert V0 == recovery["source_drawn_curve_reconstruction"]["parameters"]["V0"]
    # Check the catalog full-range rounded pair on independent graphic vertices,
    # allowing the effect of half-last-digit coefficient precision plus drawing precision.
    graphic = json.loads(
        (
            ROOT
            / "peritheos/data/research/cementite/prescher-2012-figure1-vectors.json"
        ).read_text()
    )
    from scripts.recover_prescher_2012_fe3c import bm3

    rows = graphic["drawn_fit_curve"]
    assert rows
    # Source curve dictionaries preserve their original vector coordinates.
    volume = np.array([r["volume_a3"] for r in rows])
    pressure = np.array([r["pressure_gpa"] for r in rows])
    central = bm3(volume, V0, 161, 5.9)
    bounds = np.maximum.reduce(
        [
            abs(bm3(volume, V0, k, kp) - central)
            for k in (160.5, 161.5)
            for kp in (5.85, 5.95)
        ]
    )
    assert np.all(abs(central - pressure) < bounds + 0.04)


def test_published_comparison_precision_exposes_the_remaining_k0_error_mismatch():
    from scripts.report_numbers import value_error

    assert value_error(161, 2, 0) == ("161", "2")
    assert value_error(160.994011, 2.532193, 0) == ("161", "3")
    assert value_error(5.919478, 0.240867, 1) == ("5.9", "0.2")
    assert value_error(5.9, 0.2, 1) == ("5.9", "0.2")


def test_numbered_primary_fit_replays_the_preserved_ledger_at_recorded_precision():
    from scripts.validate_primary_eos_refits import check_ledger_values

    record = next(
        r
        for r in get_material_document("cementite", format_version=4)["eos_records"]
        if r["identifier"] == IDS["full"]
    )
    ledger = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    stored = next(
        row for row in ledger["records"] if row["record_identifier"] == IDS["full"]
    )
    replay = ledger_outcome(record)
    check_ledger_values({key: stored[key] for key in replay}, replay)
    expected = {row["parameter"]: row for row in stored["parameters"]}
    actual = {row["parameter"]: row for row in replay["parameters"]}
    for name in ("K0", "K0_prime"):
        for key in ("refit", "refit_error"):
            assert actual[name][key] == pytest.approx(
                expected[name][key], rel=1e-7, abs=1e-8
            )


@pytest.mark.parametrize(
    "changes",
    [
        {"assessment_basis": ""},
        {"assessment_qualifications": []},
        {"assessment_qualifications": ["Reviewed context", 2]},
    ],
)
def test_explicit_assessment_narrative_rejects_missing_or_nonprose_values(changes):
    record = next(
        r
        for r in get_material_document("cementite", format_version=4)["eos_records"]
        if r["identifier"] == IDS["full"]
    )
    row = {
        **ledger_outcome(record),
        "record_identifier": record["identifier"],
        "material_identifier": "cementite",
        **changes,
    }
    with pytest.raises(ValueError, match="basis|qualifications"):
        assessment_for_record(record, row)
