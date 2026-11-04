"""Source selection and exact, reversible pressure-coordinate views."""

import copy
from dataclasses import replace

import numpy as np
import pytest
from scipy.optimize import brentq

from peritheos import (
    EosmatError,
    Material,
    PressureScaledEOS,
    ValidationError,
    get_eos_record,
    get_material,
    get_material_document,
    list_eos_records,
    list_reference_variants,
    recalculate_eos_pressure_scale,
    search_eos_records,
    search_materials,
    validate_eosmat_document,
)

SOURCE = "nacl_b2_sakai_2011_matsui_pt_bm3"
TARGET = "platinum_fei_2007_vinet_300k"
VARIANT = "nacl_b2_sakai_2011_fei_pt_bm3"


def test_discovery_lists_sources_and_preserves_explicit_reference_access():
    source_ids = {
        r.identifier
        for r in search_eos_records(author="Sakai", doi="10.1063/1.3573393")
    }
    assert source_ids == {SOURCE, SOURCE.replace("bm3", "vinet")}
    all_ids = {
        r.identifier
        for r in search_eos_records(
            doi="10.1063/1.3573393", include_reference_variants=True
        )
    }
    assert len(all_ids) == 8
    variants = list_reference_variants(SOURCE)
    assert {r.identifier for r in variants} == {
        *(f"nacl_b2_sakai_2011_{scale}_pt_bm3" for scale in ("fei", "do", "holmes")),
        "nacl_b2_sakai_2014_yokoo_pt_bm3",
    }
    assert get_eos_record(VARIANT) in variants
    assert search_materials(VARIANT) == ()
    assert (
        search_materials(VARIANT, include_reference_variants=True)[0].identifier
        == "nacl_b2"
    )
    material = get_material("nacl_b2")
    assert get_eos_record(VARIANT) not in material.source_records
    assert get_eos_record(VARIANT) in material.reference_variants
    assert get_eos_record(VARIANT) not in material.equilibrium_records
    assert get_eos_record(VARIANT).pressure_scale_variant_of == SOURCE


def test_pt_conversion_matches_independent_equations_not_an_alternative_fit():
    def vinet(volume, k, kp):
        x = (volume / 60.38) ** (1 / 3)
        return 3 * k * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))

    sample_volume = 23.0
    f = ((37.73 / sample_volume) ** (2 / 3) - 1) / 2
    source_p = 3 * 47.0 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (4.10 - 4) * f)
    pt_volume = brentq(lambda v: vinet(v, 273, 5.20) - source_p, 30, 60.38)
    expected = vinet(pt_volume, 277, 5.08)
    view = get_eos_record(SOURCE).on_pressure_scale(TARGET)
    assert view.pressure(sample_volume, check_validity=True) == pytest.approx(
        expected, abs=1e-9
    )
    assert view.volume(expected, check_validity=True) == pytest.approx(
        sample_volume, abs=1e-9
    )
    assert abs(expected - get_eos_record(VARIANT).pressure(sample_volume)) > 0.008
    result = view.recalculate(sample_volume)
    assert result.calibration_path == view.calibration_path
    assert result.intermediate_states[-1]["implied_standard_volume"] == pytest.approx(
        pt_volume
    )
    assert view.volume_unit == get_eos_record(SOURCE).volume_unit


@pytest.mark.parametrize(
    "source,target,pressure,temperature,strict",
    [
        (SOURCE, TARGET, [65.0, 75.0], 300.0, True),
        (SOURCE, SOURCE, [65.0, 75.0], 300.0, True),
        ("copper_dewaele_2004_vinet_1", "ruby_mao_1986", [20.0, 40.0], 298.0, True),
        (
            "forsterite_finkelstein_2014_bm3_1",
            "gold_sokolova_2013_holzapfel_4",
            [20.0, 40.0],
            298.0,
            False,
        ),
        (
            "kcl_b2_tateno_2019_vinet_4",
            "platinum_fei_2007_vinet_mgd",
            [30.0, 50.0],
            [1000.0, 1500.0],
            True,
        ),
        (
            "kcl_b2_chidester_2021_bm3_5",
            "platinum_fei_2007_vinet_mgd",
            [50.0, 70.0],
            [1200.0, 1500.0],
            True,
        ),
    ],
)
def test_general_routes_invert_the_same_sample_state(
    source, target, pressure, temperature, strict
):
    record = get_eos_record(source)
    volume = record.volume(pressure, temperature)
    view = PressureScaledEOS(source, target)
    converted = view.pressure(volume, temperature, check_validity=strict)
    expected = recalculate_eos_pressure_scale(
        source, volume, target, temperature, check_validity=strict
    )
    np.testing.assert_allclose(
        converted, expected.target_pressure_gpa, atol=1e-8, rtol=0
    )
    np.testing.assert_allclose(
        view.volume(converted, temperature, check_validity=strict),
        volume,
        atol=1e-8,
        rtol=0,
    )


def test_standard_temperature_is_distinct_from_source_sample_temperature():
    source = "kcl_b2_chidester_2021_bm3_5"
    target = "platinum_fei_2007_vinet_mgd"
    volume = get_eos_record(source).volume([50.0, 70.0], [1200.0, 1500.0])
    view = PressureScaledEOS(source, target, standard_temperature_k=1300.0)
    expected = recalculate_eos_pressure_scale(
        source, volume, target, [1200.0, 1500.0], standard_temperature_k=1300.0
    )
    np.testing.assert_allclose(
        view.pressure(volume, [1200.0, 1500.0]), expected.target_pressure_gpa
    )
    np.testing.assert_allclose(
        view.volume(expected.target_pressure_gpa, [1200.0, 1500.0]),
        volume,
        atol=1e-8,
        rtol=0,
    )


def test_unresolved_calibration_and_invalid_states_do_not_select_reference_fits():
    with pytest.raises(ValidationError, match="No executable"):
        PressureScaledEOS("nacl_b2_sakai_2011_holmes_pt_bm3", TARGET)
    with pytest.raises(ValueError, match="calibration/data envelope"):
        PressureScaledEOS(SOURCE, TARGET).pressure(16.0, check_validity=True)
    with pytest.raises(ValueError, match="isothermal"):
        PressureScaledEOS(SOURCE, TARGET).volume(70.0, 1500.0)
    for value in (-1, float("nan"), [300.0, 400.0]):
        with pytest.raises(ValidationError, match="positive scalar"):
            PressureScaledEOS(SOURCE, TARGET, standard_temperature_k=value)


def test_selected_reference_export_preserves_its_source_dependency():
    material = Material.from_eosmat(
        get_material_document("nacl_b2"), record_identifiers=[VARIANT]
    )
    assert {r.identifier for r in material.eos_records} == {SOURCE, VARIANT}
    assert material.eos_records[0].identifier == VARIANT
    validate_eosmat_document(material.to_eosmat())
    assert material.default_record().identifier == SOURCE


def test_view_binds_the_actual_loaded_coefficients_without_reloading_them():
    document = get_material_document("nacl_b2")
    raw = next(r for r in document["eos_records"] if r["identifier"] == SOURCE)
    raw["eos"]["parameters"]["K0"] = 48.0
    source = Material.from_eosmat(document, record_identifiers=[SOURCE]).eos_records[0]
    view = source.on_pressure_scale(TARGET)
    assert view.recalculate(23.0).source_pressure_gpa == pytest.approx(
        source.pressure(23.0)
    )
    assert (
        abs(
            view.recalculate(23.0).source_pressure_gpa
            - get_eos_record(SOURCE).pressure(23.0)
        )
        > 1
    )
    assert view.volume(view.pressure(23.0)) == pytest.approx(23.0, abs=1e-8)


def test_external_material_uses_its_declared_calibrant_route():
    source = replace(get_eos_record(SOURCE), identifier="external_sample_bm3")
    view = source.on_pressure_scale(TARGET)
    assert view.calibration_path[0] == source.identifier
    assert view.pressure(23.0) == pytest.approx(65.073242681114, abs=1e-8)
    assert view.volume(
        view.pressure([23.0, 22.5]), check_validity=True
    ) == pytest.approx([23.0, 22.5])
    missing = replace(source, pressure_calibration={})
    with pytest.raises(ValidationError, match="No executable"):
        missing.on_pressure_scale(TARGET)
    conflicting = replace(get_eos_record(SOURCE), pressure_calibration={})
    with pytest.raises(ValidationError, match="distinct source record identifier"):
        conflicting.on_pressure_scale(TARGET)


@pytest.mark.parametrize("bad_parent", [None, "", "unknown", VARIANT])
def test_reference_links_must_resolve_to_another_source(bad_parent):
    document = get_material_document("nacl_b2")
    variant = next(r for r in document["eos_records"] if r["identifier"] == VARIANT)
    variant["pressure_scale_variant_of"] = bad_parent
    with pytest.raises(EosmatError, match="pressure_scale_variant_of"):
        validate_eosmat_document(document)


def test_reference_variants_cannot_chain_change_model_or_become_defaults():
    base = get_material_document("nacl_b2")
    for mutation in (
        {"default": True},
        {"fixed_parameters": []},
        {"pressure_scale_variant_of": "nacl_b2_sakai_2011_do_pt_bm3"},
        {"pressure_scale_variant_of": SOURCE.replace("bm3", "vinet")},
    ):
        document = copy.deepcopy(base)
        variant = next(r for r in document["eos_records"] if r["identifier"] == VARIANT)
        variant.update(mutation)
        with pytest.raises(EosmatError):
            validate_eosmat_document(document)


def test_catalog_audit_preserves_distinct_experiments_and_all_coefficients():
    references = [
        r
        for r in list_eos_records(include_reference_variants=True)
        if r.is_reference_variant
    ]
    assert len(references) == 51
    assert len({r.pressure_scale_variant_of for r in references}) == 27
    for identifier in (
        "mgo_jacobsen_2008_bm3_kcl_mao1978",
        "mgo_jacobsen_2008_bm3_helium_mao1986",
        "akimotoite_wang_2004_t0133_bm3",
        "akimotoite_wang_2004_t0150_bm3",
        "gold_dewaele_2004_mao_ruby_vinet",
        "platinum_dewaele_2004_mao_ruby_vinet",
        "tantalum_dewaele_2004_mao_ruby_vinet",
    ):
        assert not get_eos_record(identifier).is_reference_variant
