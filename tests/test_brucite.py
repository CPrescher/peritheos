import csv
import hashlib
import json
from collections import Counter

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    eosmat_schema,
    get_eos_record,
    get_material,
    load_eosmat,
    search_eos_records,
)
from peritheos.errors import MaterialError, MaterialLookupError
from scripts.brucite_reproduce import (
    COLUMNS,
    DATA,
    REPORT,
    STAGED,
    TABLES,
    bm3,
    conditional_continuous_pressure,
    normalized_csv,
    observations,
    parse_number,
    reference_volume,
    reproduce,
    source,
)
from scripts.brucite_replay_checks import check_report
from scripts.brucite_stage import audit_document, catalog_document, conditional_thermal_document, staged_document
from scripts.eosmat_io import read_material
from scripts.run_eosmat_audit import historical_reconstruction_view


@pytest.fixture(scope="module")
def candidate():
    document = load_eosmat(STAGED)
    return document, Material.from_eosmat(document, require_primary_validation=False)


def test_brucite_staged_only_schema_and_validation_gate(candidate):
    document, material = candidate
    Draft202012Validator(eosmat_schema()).validate(document)
    assert document == staged_document()
    assert get_material("brucite").default_record().identifier == (
        "brucite_fukui_2003_300k_bm3"
    )
    with pytest.raises(MaterialError, match="pending_primary_source_check"):
        Material.from_eosmat(document)
    assert len(material.eos_records) == 4
    assert material.default_record().identifier == "brucite_fukui_2003_300k_bm3"
    assert all(not r.is_thermal for r in material.eos_records)
    restored = Material.from_eosmat(
        material.to_eosmat(), require_primary_validation=False
    )
    assert restored.to_eosmat() == material.to_eosmat()
    for record in restored.eos_records:
        assert record.scientific_validation_status == "pending_primary_source_check"


def test_brucite_full_tables_errors_and_extrema(candidate):
    document, material = candidate
    expected_counts = [11, 37, 54]
    assert sum(expected_counts) == 102
    for (owner, table, filename), dataset, count in zip(
        TABLES, document["datasets"], expected_counts
    ):
        payload = (DATA / filename).read_bytes()
        assert payload == normalized_csv(owner, table).encode()
        assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
        rows = list(csv.DictReader(payload.decode().splitlines()))
        assert len(rows) == count == dataset["row_count"]
        assert list(rows[0]) == list(COLUMNS)
        assert all(r["temperature_error_k"] == "" for r in rows)
        loaded = material.get_dataset(dataset["identifier"])
        assert loaded.checksum_verified
        assert len(loaded) == count
    ambient = observations("fukui_2003", "table1")
    assert ambient[0]["temperature_k"] == ambient[1]["temperature_k"] == "299.3"
    assert ambient[2]["lattice_c_error_angstrom"] == "0.0025"
    assert ambient[8]["cell_volume_error_a3"] == "0.005"
    assert all(r["pressure_gpa"] == r["pressure_error_gpa"] == "" for r in ambient)
    rows = observations("fukui_2003", "table4")
    assert rows[0]["pressure_gpa"] == "-0.06"
    assert rows[0]["pressure_error_gpa"] == "0.13"
    assert rows[5]["lattice_c_error_angstrom"] == "0.0160"
    assert rows[19]["cell_volume_error_a3"] == "0.474"
    assert Counter(r["temperature_k"] for r in rows) == {
        "300": 13,
        "473": 8,
        "673": 8,
        "873": 8,
    }
    assert max(float(r["pressure_gpa"]) for r in rows) == 22.56
    xia = observations("xia_1998", "table6")
    assert xia[12]["temperature_k"] == "1040"
    assert xia[41]["cell_volume_a3"] == "37.39"  # OCR incorrectly reads 37.35.
    assert xia[-1]["pressure_gpa"] == "3.06"
    assert max(float(r["pressure_gpa"]) for r in xia) == 9.24
    assert all(r["pressure_error_gpa"] == "" for r in xia)
    assert parse_number("4.4424(160)") == ("4.4424", "0.0160")
    typed = material.get_dataset("brucite_fukui_2003_table4").as_pressure_volume()
    assert typed.pressure[0] == -0.06
    assert typed.pressure_uncertainty[0] == 0.13
    assert typed.volume_uncertainty[19] == 0.474
    assert typed.pressure_uncertainty_role == "standard_deviation"


def test_brucite_ambient_ordered_structure_is_separate(candidate):
    document, _ = candidate
    assert document["space_group_number"] == 164
    assert document["formula_units_per_cell"] == 1
    counts = Counter()
    for site in document["atom_sites"]:
        counts[site["element"]] += int(site["wyckoff"][:-1]) * site["occupancy"]
    assert counts == {"Mg": 1, "O": 2, "H": 2}
    assert document["atom_sites"][1]["x"] == 1 / 3
    assert document["atom_sites"][2]["z"] == 0.4108
    lattice = document["lattice"]
    volume = np.sqrt(3) / 2 * lattice["a"] ** 2 * lattice["c"]
    assert volume == pytest.approx(40.9437, abs=0.0001)
    assert volume != pytest.approx(40.917, abs=0.01)
    structure = document["source"]
    assert "Primary doctoral thesis" in structure["structure_source_kind"]
    cif = DATA.parent / structure["structure_dataset"]
    assert (
        hashlib.sha256(cif.read_bytes()).hexdigest()
        == structure["structure_dataset_sha256"]
    )
    assert "0.2189(2)" in cif.read_text()
    assert "0.4108(4)" in cif.read_text()
    assert structure["structure_reference"]["doi"] == "10.15083/0002012205"
    # Mg(OH)2 has 5 atoms per formula unit; no factor of 4/5/Avogadro belongs
    # in the BM3 cell-volume calculation or dimensionless expansion integral.
    assert source()["volume_basis"]["atoms_per_formula_unit"] == 5
    assert reference_volume(300) * 0.602214076 == pytest.approx(24.640793347692)


@pytest.mark.parametrize("temperature", [300, 473, 673, 873])
def test_brucite_isotherm_equation_derivative_and_inversion(candidate, temperature):
    document, material = candidate
    identifier = f"brucite_fukui_2003_{temperature}k_bm3"
    record = material.get_eos_record(identifier)
    data = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    params = data["eos"]["parameters"]
    v0 = params["V0"]
    volumes = v0 * np.array([1.0, 0.98, 0.90, 0.82])
    expected = bm3(volumes, v0, params["K0"], params["K0_prime"])
    actual = record.pressure(volumes, temperature)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(record.volume(actual, temperature), volumes, rtol=1e-10)
    assert record.pressure(v0, temperature) == pytest.approx(0, abs=1e-12)
    assert record.eos.bulk_modulus(v0) == pytest.approx(params["K0"])
    h = 1e-4
    derivative = (
        -volumes
        * (
            record.pressure(volumes + h, temperature)
            - record.pressure(volumes - h, temperature)
        )
        / (2 * h)
    )
    np.testing.assert_allclose(record.eos.bulk_modulus(volumes), derivative)
    with pytest.raises(MaterialError):
        record.volume(10, temperature + 1, check_validity=True)


def test_brucite_source_observations_are_independent_checkpoints(candidate):
    _, material = candidate
    static = material.get_eos_record("brucite_fukui_2003_300k_bm3")
    hot = material.get_eos_record("brucite_fukui_2003_473k_bm3")
    # Direct Table 4 observations, not generated EOS fixtures. At 22.56 GPa
    # the combined P/V observational error is about 0.25 GPa. At 5.95 GPa
    # the pressure error alone is 0.03 GPa. These test the printed curve,
    # outside the conservative ambient phase-review domain when necessary.
    assert static.pressure(32.058, 300) == pytest.approx(22.56, abs=0.30)
    assert hot.pressure(36.883, 473) == pytest.approx(5.95, abs=0.03)
    assert hot.pressure(33.626, 473) == pytest.approx(15.15, abs=0.18)
    static.volume(22.56, 300, check_validity=True)
    with pytest.raises(MaterialError):
        static.volume(22.57, 300, check_validity=True)


def test_brucite_conditional_agreement_does_not_claim_exact_reproduction():
    report = reproduce()
    stored = json.loads(REPORT.read_text(encoding="utf-8"))
    check_report(report, stored)
    assert (
        max(
            abs(x)
            for x in report["ambient_diagnostic_refits"]["unweighted"][
                "delta_over_reported_errors"
            ]
        )
        < 0.04
    )
    for isotherm in report["isotherms"]:
        assert isotherm["source_equation_validated"]
        assert not isotherm["author_workflow_reproduced"]
        weighted = isotherm["diagnostic_refits"]["pressure_sigma"]
        assert weighted["success"]
        assert max(abs(x) for x in weighted["delta_over_reported_errors"]) < 1
    # The high-T source curve has an unresolved systematic offset. Preserve
    # this evidence even though an assumed weighted refit agrees in coefficients.
    hot = report["isotherms"][-1]["published_curve"]
    assert 0.48 < hot["pressure_rmse_gpa"] < 0.49
    assert hot["max_absolute_effective_sigma_residual"] > 6.7
    assert report["metric_checks"]["xia_table6"][
        "rows_exceeding_printed_volume_error"
    ] == [1, 10, 14, 18, 21, 30, 50, 52]


def test_brucite_continuous_thermal_and_calibration_holds(candidate):
    document, _ = candidate
    for record in document["eos_records"]:
        calibration = record["pressure_calibration"]
        assert calibration["methods"][0]["reference_eos_record"] == (
            "gold_anderson_1989_bm3_1"
        )
        assert (
            calibration["recalculation"]["status"] == "missing_calibrant_observations"
        )
        assert "thermal" not in record
    src = source()["sources"]
    assert src["fukui_2003"]["parameters"]["dK0_prime_dT_per_k"]["value"] == 0.00571
    assert (
        src["fukui_2003"]["phase_scope"]["ambient_dehydration_temperature_k"] == 600.6
    )
    assert src["xia_1998"]["table7_volume"]["alpha1_per_k2"] == "3.6(13)e-8"
    assert src["fei_mao_1993"]["abstract_parameters"]["K0_300_gpa"] == 54.3
    # Prove why freezing Kprime is not the published continuous law. This
    # conditional surface is an audit helper, never a thermal card record.
    v = 34.0
    temperature = 873
    frozen = bm3(v, reference_volume(temperature), 41.8 - 0.0319 * 573, 6.6)
    assert abs(conditional_continuous_pressure(v, temperature) - frozen) > 3


def test_staged_conditional_surface_keeps_scientific_holds():
    doc = conditional_thermal_document()
    Draft202012Validator(eosmat_schema()).validate(doc)
    material = Material.from_eosmat(doc, require_primary_validation=False)
    record = material.get_eos_record("brucite_fukui_2003_conditional_thermal")
    assert record.scientific_validation_status == "pending_primary_source_check"
    assert record.pressure(34, 873) == pytest.approx(14.051454515344124)
    assert len(material.eos_records) == 5
    assert not record.is_default


def test_catalog_brucite_published_ranges_and_qualified_source_records():
    document = catalog_document()
    Draft202012Validator(eosmat_schema()).validate(document)
    material = Material.from_eosmat(document)
    assert len(material.eos_records) == 3
    assert len(document["eos_records"]) == 3
    for temperature, pmax in [(300, 22.56), (473, 16.52), (673, 16.09)]:
        record = material.get_eos_record(f"brucite_fukui_2003_{temperature}k_bm3")
        v = record.volume(pmax, temperature, check_validity=True)
        assert record.pressure(v, temperature, check_validity=True) == pytest.approx(
            pmax
        )
        with pytest.raises(MaterialError):
            record.volume(pmax + 0.01, temperature, check_validity=True)
    static = material.default_record()
    assert static.pressure(40.909, 300, check_validity=True) == pytest.approx(
        0, abs=0.06
    )
    archive = audit_document()
    assert len(archive["eos_records"]) == 12
    assert (
        archive["curation_disposition"]["thermal_eos_recommendation"]
        == "not_recommended"
    )
    for record in archive["eos_records"]:
        if record["identifier"] in {r.identifier for r in material.eos_records}:
            continue
        with pytest.raises(MaterialLookupError):
            get_eos_record(record["identifier"])
        with pytest.raises(MaterialError):
            Material.from_eosmat(archive, record_identifiers=[record["identifier"]])
        assert (
            "Not recommended" in record["scientific_validation"]["usage_recommendation"]
        )
    registered = read_material(DATA.parent / "materials/brucite.eosmat")
    assert all("assessment" in r["reproduction"] for r in registered["eos_records"])
    # The later ledger summary is additive to the archived scientific document.
    assert historical_reconstruction_view(registered) == document
    assert search_eos_records(formula="Mg(OH)2", pressure_gpa=22.56)
    assert not search_eos_records(formula="Mg(OH)2", thermal=True)
    assert not search_eos_records(formula="Mg(OH)2", temperature_k=873)
    assert load_eosmat(STAGED.with_name("brucite-diagnostic-archive.eosmat")) == archive


def test_brucite_figure_normalization_is_evidence_not_a_silent_correction():
    report = reproduce()
    audit = report["figure6c_audit"]
    assert audit["mean_implied_normalization_volume_a3"] == pytest.approx(
        42.70, abs=0.02
    )
    assert audit["printed_expansion_reference_volume_a3"] == pytest.approx(
        42.5159510616
    )
    assert all(
        abs(p["digitized_pressure_gpa"] - p["printed_BM3_pressure_gpa"]) < 0.1
        for p in audit["solid_curve"]
    )
    original = Material.from_eosmat(
        audit_document(), require_primary_validation=False
    ).get_eos_record("brucite_fukui_2003_873k_bm3")
    assert original.eosmat_metadata["eos"]["parameters"]["K0"] == 23.9
    assert original.eosmat_metadata["eos"]["parameters"]["V0"] == pytest.approx(
        reference_volume(873)
    )


def test_brucite_independent_refits_use_table_rows_and_retain_conditional_covariance():
    material = Material.from_eosmat(audit_document(), require_primary_validation=False)
    hot = material.get_eos_record("brucite_peritheos_fukui_2003_873k_bm3_refit")
    thermal = material.get_eos_record("brucite_peritheos_fukui_2003_thermal_refit")
    rows = observations("fukui_2003", "table4")
    for record in (hot, thermal):
        meta = record.eosmat_metadata
        assert meta["record_kind"] == "refit"
        assert not record.is_default
        covariance = np.array(meta["parameter_covariance"]["matrix"])
        np.testing.assert_allclose(covariance, covariance.T, atol=1e-12)
        assert np.linalg.eigvalsh(covariance).min() > 0
        assert "conditional" in meta["fit_provenance"]["covariance_scaling"]
        # Verify the joint covariance is usable by public uncertainty propagation,
        # including correlation between the reference moduli and thermal slopes.
        names = ["K0", "K0_prime"] + (
            ["dK_dT", "dKprime_dT"] if record.is_thermal else []
        )
        values = {
            **meta["eos"]["parameters"],
            **meta.get("thermal", {}).get("parameters", {}),
        }

        def pressure_at(parameters):
            dt = 573 if record.is_thermal else 0
            return bm3(
                35.5,
                reference_volume(873),
                parameters["K0"] + parameters.get("dK_dT", 0) * dt,
                parameters["K0_prime"] + parameters.get("dKprime_dT", 0) * dt,
            )

        gradient = []
        for name in names:
            h = 1e-7 if name.startswith("dK") else 1e-5
            gradient.append(
                (
                    pressure_at({**values, name: values[name] + h})
                    - pressure_at({**values, name: values[name] - h})
                )
                / (2 * h)
            )
        expected_sigma = np.sqrt(np.array(gradient) @ covariance @ np.array(gradient))
        assert record.pressure_with_uncertainty(
            35.5, 873
        ).standard_error == pytest.approx(expected_sigma, rel=1e-5)
        chosen = [
            r for r in rows if record.is_thermal or float(r["temperature_k"]) == 873
        ]
        residual = [
            record.pressure(float(r["cell_volume_a3"]), float(r["temperature_k"]))
            - float(r["pressure_gpa"])
            for r in chosen
        ]
        assert np.sqrt(np.mean(np.array(residual) ** 2)) == pytest.approx(
            meta["fit_provenance"]["statistics"]["pressure_rmse_gpa"], abs=1e-12
        )
        restored = Material.from_eosmat(
            material.to_eosmat(), require_primary_validation=False
        ).get_eos_record(record.identifier)
        assert (
            restored.eosmat_metadata["parameter_covariance"]
            == meta["parameter_covariance"]
        )
    assert (
        hot.eosmat_metadata["fit_provenance"]["statistics"]["pressure_rmse_gpa"] < 0.16
    )
    assert (
        thermal.eosmat_metadata["fit_provenance"]["statistics"]["pressure_rmse_gpa"]
        < 0.32
    )
    for p, t in [(5.95, 300), (10, 573), (15.48, 873)]:
        v = thermal.volume(p, t, check_validity=True)
        assert thermal.pressure(v, t) == pytest.approx(p)
    for p, t in [(5.94, 873), (15.49, 873), (10, 874), (0, 873)]:
        with pytest.raises(MaterialError):
            thermal.volume(p, t, check_validity=True)


def test_primary_structure_and_phase_sources_do_not_imply_universal_boundary():
    src = source()["sources"]
    structure = src["structure"]
    assert structure["refinement"]["Nobs"] == 173
    assert structure["refinement"]["O_H_angstrom"] == "0.915(2)"
    assert structure["refinement"]["motion_corrected_O_H_angstrom"] == 0.9528
    assert src["structure_zigan_1967_archive"]["sites"][2]["z"] == "0.4303(12)"
    assert src["duffy_1995_phase"]["primary_findings"]["new_modes_gpa"] == 4.4
    phase = src["wu_2026_phase"]
    assert phase["reported_boundary"]["slope_gpa_per_k"] == -0.008
    assert len(phase["source_inconsistencies"]) == 3
    assert source()["catalog_disposition"]["held_isotherms_k"] == []
    assert src["sun_2025_phase"]["doi"] == "10.1016/j.gsf.2024.101940"


def test_brucite_replay_tolerances_preserve_sources_and_small_coefficients():
    import copy

    saved = {
        "published_parameters": [41.8, 6.6],
        "fit": {"parameters": {"alpha1": 4.4e-8}},
        "source_rows": [1, 2],
        "success": True,
    }
    replay = copy.deepcopy(saved)
    replay["fit"]["parameters"]["alpha1"] *= 1 + 1e-6
    check_report(replay, saved)
    for field, value in (
        ("published_parameters", [41.8 + 1e-10, 6.6]),
        ("source_rows", [1, 3]),
        ("success", False),
    ):
        changed = {**saved, field: value}
        with pytest.raises(ValueError):
            check_report(changed, saved)
    replay["fit"]["parameters"]["alpha1"] = 0.0
    with pytest.raises(ValueError, match="calculated value"):
        check_report(replay, saved)


def test_brucite_refit_requires_curve_parity_within_coefficient_tolerance(monkeypatch):
    from scripts import brucite_reproduce as audit

    record = next(
        r for r in audit_document()["eos_records"]
        if r["identifier"] == "brucite_peritheos_fukui_2003_thermal_refit"
    )
    report = reproduce()
    # A small correlated-coefficient change may meet the scalar comparison but
    # must still be rejected when it changes the evaluated pressure curve.
    report["peritheos_refits"]["continuous_thermal"]["parameters"]["K0"] *= 1.00005
    monkeypatch.setattr(audit, "reproduce", lambda: report)
    with pytest.raises(ValueError, match="stale"):
        audit.ledger_outcome(record)
