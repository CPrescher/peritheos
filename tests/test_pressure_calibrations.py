import copy

import numpy as np
import pytest

import peritheos.eosmat as eosmat_module
import peritheos.pressure_calibrations as calibration_module
from peritheos import (
    DiamondRamanCalibration,
    find_common_pressure_calibration_routes,
    find_pressure_calibration_path,
    get_cross_calibration_edge,
    get_material_document,
    get_pressure_calibration,
    get_pressure_calibration_document,
    list_cross_calibration_edges,
    list_diamond_raman_calibrations,
    list_material_documents,
    list_pressure_calibrations,
    list_ruby_pressure_calibrations,
    list_ruby_xrd_bridges,
    list_xrd_pressure_standards,
    recalculate_diamond_raman_pressure,
    recalculate_eos_pressure_scale,
    recalculate_pressure_calibration_path,
    recalculate_ruby_pressure,
    recalculate_ruby_to_xrd_pressure,
    recalculate_ruby_with_measured_xrd_standard,
    recalculate_xrd_pressure_scale,
    validate_pressure_calibration_references,
    xrd_standard_pressure,
)
from peritheos.errors import EosmatError, MaterialLookupError, ValidationError
from peritheos.materials import Material

EXPECTED_RUBY_SCALES = {
    "ruby_mao_1978",
    "ruby_mao_1986",
    "ruby_dewaele_2004",
    "ruby_holzapfel_2005",
    "ruby_dorogokupets_oganov_2007",
    "ruby_shen_2020",
}

EXPECTED_DIAMOND_RAMAN_SCALES = {
    "diamond_raman_akahama_2006",
    "diamond_raman_eremets_2023",
}


def test_bundled_ruby_calibration_documents_are_executable():
    assert set(list_ruby_pressure_calibrations()) == EXPECTED_RUBY_SCALES
    assert set(list_pressure_calibrations()) == (
        EXPECTED_RUBY_SCALES | EXPECTED_DIAMOND_RAMAN_SCALES
    )
    for identifier in list_ruby_pressure_calibrations():
        document = get_pressure_calibration_document(identifier)
        calibration = get_pressure_calibration(identifier)
        assert document["kind"] == "ruby_fluorescence"
        assert calibration.identifier == identifier
        assert calibration.pressure_from_ratio(1.0) == pytest.approx(0.0)


def test_bundled_diamond_raman_calibrations_are_executable_and_reversible():
    assert set(list_diamond_raman_calibrations()) == EXPECTED_DIAMOND_RAMAN_SCALES
    for identifier in list_diamond_raman_calibrations():
        calibration = get_pressure_calibration(identifier)
        assert isinstance(calibration, DiamondRamanCalibration)
        pressure = np.array([0.0, 50.0, 150.0, 300.0])
        wavenumber = calibration.wavenumber_from_pressure(pressure)
        assert calibration.pressure_from_wavenumber(wavenumber) == pytest.approx(
            pressure
        )

    converted = recalculate_diamond_raman_pressure(
        [50.0, 150.0],
        "diamond_raman_akahama_2006",
        "diamond_raman_eremets_2023",
    )
    assert recalculate_diamond_raman_pressure(
        converted,
        "diamond_raman_eremets_2023",
        "diamond_raman_akahama_2006",
    ) == pytest.approx([50.0, 150.0])


@pytest.mark.parametrize("identifier", sorted(EXPECTED_RUBY_SCALES))
def test_ruby_scales_round_trip_pressure_ratio_and_wavelength(identifier):
    calibration = get_pressure_calibration(identifier)
    pressure = np.array([0.0, 1.0, 50.0, 150.0, 300.0])
    ratio = calibration.wavelength_ratio(pressure)
    assert calibration.pressure_from_ratio(ratio) == pytest.approx(pressure)

    wavelength = calibration.wavelength_from_pressure(pressure)
    assert calibration.pressure_from_wavelength(wavelength) == pytest.approx(pressure)


def test_pressure_recalculation_between_ruby_scales_is_reversible():
    source_pressure = np.array([0.0, 10.0, 50.0, 100.0])
    recalculated = recalculate_ruby_pressure(
        source_pressure,
        "ruby_mao_1986",
        "ruby_dorogokupets_oganov_2007",
    )
    assert recalculated[-1] > source_pressure[-1]
    assert recalculate_ruby_pressure(
        recalculated,
        "ruby_dorogokupets_oganov_2007",
        "ruby_mao_1986",
    ) == pytest.approx(source_pressure)


def test_ruby_only_apis_reject_diamond_calibration_records():
    with pytest.raises(ValidationError, match="must identify a ruby calibration"):
        recalculate_ruby_with_measured_xrd_standard(
            50.0,
            "diamond_raman_akahama_2006",
            "gold_fei_2007_vinet_2",
            60.0,
        )
    with pytest.raises(ValidationError, match="must identify a ruby calibration"):
        recalculate_ruby_to_xrd_pressure(
            50.0,
            "diamond_raman_akahama_2006",
            "gold_dewaele_2004_vinet_5",
            "gold_fei_2007_vinet_2",
        )


def test_ruby_shift_api_and_input_validation():
    calibration = get_pressure_calibration("ruby_mao_1986")
    shift = (
        calibration.wavelength_from_pressure(50.0) - calibration.reference_wavelength_nm
    )
    assert calibration.pressure_from_shift(shift) == pytest.approx(50.0)
    with pytest.raises(ValidationError):
        calibration.pressure_from_ratio(0.99)
    with pytest.raises(ValidationError):
        calibration.wavelength_ratio(-1.0)
    with pytest.raises(ValidationError):
        calibration.wavelength_ratio(np.inf)
    with pytest.raises(ValidationError):
        calibration.pressure_from_wavelength(700.0, reference_wavelength_nm=0.0)
    with pytest.raises(ValidationError):
        calibration.pressure_from_shift(1.0, reference_wavelength_nm=np.nan)
    with pytest.raises(ValidationError):
        calibration.wavelength_from_pressure(1.0, reference_wavelength_nm=-1.0)
    with pytest.raises(ValidationError):
        get_pressure_calibration("ruby_holzapfel_2005").wavelength_ratio(2_000.0)


def test_every_identified_ruby_use_links_to_an_executable_calibration():
    ruby_methods = []
    for material_identifier in list_material_documents():
        for record in get_material_document(material_identifier)["eos_records"]:
            for method in record["pressure_calibration"]["methods"]:
                if method["kind"] == "ruby_fluorescence":
                    ruby_methods.append(method)

    assert len(ruby_methods) == 33
    assert all(
        method["reference_calibration_record"] in EXPECTED_RUBY_SCALES
        for method in ruby_methods
    )
    validate_pressure_calibration_references()


def test_holzapfel_2005_reference_is_not_the_unrelated_aip_doi():
    document = get_pressure_calibration_document("ruby_holzapfel_2005")
    assert document["parameters"] == {
        "A_gpa": 1845.0,
        "B": 14.7,
        "C": 7.5,
        "reference_wavelength_nm": 694.24,
    }
    assert document["reference"]["doi"] == "10.1080/09511920500147501"


def test_xrd_standard_pressure_uses_bundled_reference_eos():
    assert "gold_dewaele_2004_vinet_5" in list_xrd_pressure_standards()
    gold = Material.from_eosmat(
        get_material_document("gold"),
        record_identifiers=["gold_dewaele_2004_vinet_5"],
    ).get_eos_record("gold_dewaele_2004_vinet_5")
    volumes = gold.volume(np.array([0.0, 50.0, 100.0]))
    assert xrd_standard_pressure("gold_dewaele_2004_vinet_5", volumes) == pytest.approx(
        [0.0, 50.0, 100.0]
    )


def test_ruby_pressures_can_be_replaced_by_paired_gold_observations():
    gold = Material.from_eosmat(
        get_material_document("gold"),
        record_identifiers=["gold_dewaele_2004_vinet_5"],
    ).get_eos_record("gold_dewaele_2004_vinet_5")
    source_pressure = np.array([10.0, 50.0, 90.0])
    gold_pressure = np.array([10.5, 52.0, 95.0])
    result = recalculate_ruby_with_measured_xrd_standard(
        source_pressure,
        "ruby_mao_1986",
        "gold_dewaele_2004_vinet_5",
        gold.volume(gold_pressure),
    )
    assert result.target_pressure_gpa == pytest.approx(gold_pressure)
    assert result.pressure_difference_gpa == pytest.approx(
        gold_pressure - source_pressure
    )
    assert result.source_wavelength_ratio == pytest.approx(
        get_pressure_calibration("ruby_mao_1986").wavelength_ratio(source_pressure)
    )


def test_xrd_scales_are_crossed_through_a_virtual_same_standard_volume():
    gold = Material.from_eosmat(
        get_material_document("gold"),
        record_identifiers=[
            "gold_dewaele_2004_vinet_5",
            "gold_fei_2007_vinet_2",
        ],
    )
    source = gold.get_eos_record("gold_dewaele_2004_vinet_5")
    target = gold.get_eos_record("gold_fei_2007_vinet_2")
    source_pressure = np.array([10.0, 50.0, 90.0])
    result = recalculate_xrd_pressure_scale(
        source_pressure,
        source.identifier,
        target.identifier,
    )
    virtual_volume = source.volume(source_pressure, 300.0)
    assert result.implied_standard_volume == pytest.approx(virtual_volume)
    assert result.target_pressure_gpa == pytest.approx(
        target.pressure(virtual_volume, 300.0)
    )
    assert result.pressure_difference_gpa == pytest.approx(
        np.asarray(result.target_pressure_gpa) - source_pressure
    )


def test_ruby_scale_reaches_fei_through_a_ruby_linked_gold_eos():
    source_pressure = np.array([25.0, 75.0, 120.0])
    result = recalculate_ruby_to_xrd_pressure(
        source_pressure,
        "ruby_mao_1986",
        "gold_dewaele_2004_vinet_5",
        "gold_fei_2007_vinet_2",
    )
    bridge_pressure = recalculate_ruby_pressure(
        source_pressure,
        "ruby_mao_1986",
        "ruby_dewaele_2004",
    )
    assert result.bridge_pressure_gpa == pytest.approx(bridge_pressure)
    assert np.all(np.asarray(result.implied_standard_volume) > 0.0)
    assert result.bridge_calibration_record == "ruby_dewaele_2004"


def test_sample_eos_normalization_selects_xrd_or_ruby_provenance_edge():
    b4c = Material.from_eosmat(
        get_material_document("b4c"),
        record_identifiers=["b4c_somayazulu_2023_bm3_1"],
    ).get_eos_record("b4c_somayazulu_2023_bm3_1")
    b4c_volumes = b4c.volume(np.array([20.0, 50.0]))
    xrd_result = recalculate_eos_pressure_scale(
        b4c.identifier,
        b4c_volumes,
        "mgo_sokolova_2013_holzapfel_4",
        standard_temperature_k=300.0,
    )
    assert "mgo_b1_tange_2009_vinet" in xrd_result.calibration_path

    tungsten = Material.from_eosmat(
        get_material_document("tungsten"),
        record_identifiers=["tungsten_dewaele_2004_vinet_2"],
    ).get_eos_record("tungsten_dewaele_2004_vinet_2")
    tungsten_volumes = tungsten.volume(np.array([20.0, 50.0]))
    ruby_result = recalculate_eos_pressure_scale(
        tungsten.identifier,
        tungsten_volumes,
        "gold_fei_2007_vinet_2",
        bridge_eos_record="gold_dewaele_2004_vinet_5",
    )
    assert ruby_result.source_calibration_record == "ruby_dewaele_2004"


def test_fei_gold_has_two_ruby_linked_gold_bridge_options():
    assert list_ruby_xrd_bridges("gold_fei_2007_vinet_2") == (
        "gold_dewaele_2004_vinet_5",
        "gold_takemura_2008_vinet_6",
    )


def test_calibration_lookup_rejects_unknown_duplicate_and_unsupported_kind(monkeypatch):
    with pytest.raises(MaterialLookupError):
        get_pressure_calibration_document("unknown")

    mao = get_pressure_calibration_document("ruby_mao_1986")
    monkeypatch.setattr(
        calibration_module,
        "pressure_calibration_library",
        lambda: {"calibrations": [mao, mao]},
    )
    with pytest.raises(ValidationError, match="Duplicate"):
        get_pressure_calibration_document("ruby_mao_1986")

    non_ruby = dict(mao, kind="other_optical_gauge")
    monkeypatch.setattr(
        calibration_module,
        "get_pressure_calibration_document",
        lambda _identifier: non_ruby,
    )
    with pytest.raises(ValidationError, match="unsupported kind"):
        get_pressure_calibration("not_ruby")


def test_xrd_recalculation_rejects_missing_standard_and_shape_mismatch():
    with pytest.raises(MaterialLookupError):
        xrd_standard_pressure("missing_eos", 10.0)
    with pytest.raises(ValidationError, match="broadcast-compatible"):
        recalculate_ruby_with_measured_xrd_standard(
            [10.0, 20.0],
            "ruby_mao_1986",
            "gold_dewaele_2004_vinet_5",
            [60.0, 58.0, 56.0],
        )

    with pytest.raises(ValidationError, match="same XRD standard"):
        recalculate_xrd_pressure_scale(
            50.0,
            "gold_dewaele_2004_vinet_5",
            "mgo_b1_tange_2009_vinet",
        )
    with pytest.raises(ValidationError, match="exactly one executable ruby"):
        recalculate_ruby_to_xrd_pressure(
            50.0,
            "ruby_dewaele_2004",
            "gold_fei_2007_vinet_2",
            "gold_dewaele_2004_vinet_5",
        )
    result = recalculate_eos_pressure_scale(
        "tungsten_dewaele_2004_vinet_2",
        30.0,
        "gold_fei_2007_vinet_2",
    )
    assert result.calibration_path[0] == "tungsten_dewaele_2004_vinet_2"
    assert result.calibration_path[-1] == "gold_fei_2007_vinet_2"


def test_recursive_graph_connects_kcl_ruby_diamond_platinum_and_gold():
    expected_edges = {
        "kcl_tateno_2019_to_platinum_sokolova_2013",
        "kcl_chidester_2021_to_platinum_dorogokupets_2007",
        "diamond_akahama_2006_to_platinum_holmes_1989",
        "diamond_eremets_2023_to_gold_fratanduono_2021",
        "sokolova_2013_gold_platinum_family",
    }
    assert set(list_cross_calibration_edges()) == expected_edges
    assert get_cross_calibration_edge(
        "kcl_chidester_2021_to_platinum_dorogokupets_2007"
    )["executable"]
    path = find_pressure_calibration_path(
        "kcl_b2_tateno_2019_vinet_4",
        "platinum_sokolova_2013_holzapfel_3",
    )
    assert [edge["identifier"] for edge in path] == [
        "kcl_tateno_2019_to_platinum_sokolova_2013"
    ]

    cases = [
        "kcl_b2_tateno_2019_vinet_4",
        "kcl_b2_chidester_2021_bm3_5",
        "ruby_shen_2020",
        "diamond_raman_akahama_2006",
        "diamond_raman_eremets_2023",
    ]
    for source in cases:
        path = find_pressure_calibration_path(source, "gold_fei_2007_vinet_2")
        assert path
        result = recalculate_pressure_calibration_path(
            50.0,
            source,
            "gold_fei_2007_vinet_2",
            300.0,
        )
        assert result.target_pressure_gpa > 0.0
        assert result.calibration_path[0] == source
        assert result.calibration_path[-1] == "gold_fei_2007_vinet_2"


def test_common_route_discovery_finds_and_ranks_shared_xrd_targets():
    sources = (
        "forsterite_finkelstein_2014_bm3_1",
        "ca_perovskite_shim_2000_bm3_1",
        "tungsten_dewaele_2004_vinet_2",
    )
    targets = (
        "gold_fei_2007_vinet_2",
        "gold_sokolova_2013_holzapfel_4",
    )
    common = find_common_pressure_calibration_routes(
        sources,
        target_nodes=targets,
    )

    assert {result.target_node for result in common} == set(targets)
    assert common == tuple(
        sorted(
            common,
            key=lambda result: (
                result.maximum_path_length,
                result.total_edge_count,
                result.target_node,
            ),
        )
    )
    for result in common:
        assert tuple(result.paths) == sources
        for source, path in result.paths.items():
            if path:
                assert path[0]["source_node"] == source
                assert path[-1]["target_node"] == result.target_node
        assert result.maximum_path_length == max(map(len, result.paths.values()))
        assert result.total_edge_count == sum(map(len, result.paths.values()))


@pytest.mark.parametrize(
    ("sources", "targets", "match"),
    [
        ((), None, "at least one"),
        ("ruby_mao_1986", None, "not a string"),
        (("unknown",), None, "Unknown pressure-calibration source"),
        (("ruby_mao_1986",), (), "at least one"),
        (("ruby_mao_1986",), ("unknown",), "Unknown pressure-calibration target"),
    ],
)
def test_common_route_discovery_validates_node_sets(sources, targets, match):
    with pytest.raises(ValidationError, match=match):
        find_common_pressure_calibration_routes(sources, target_nodes=targets)


def test_chidester_edge_restores_average_pt_temperature_and_is_reversible():
    edge_identifier = "kcl_chidester_2021_to_platinum_dorogokupets_2007"
    pt_identifier = "platinum_dorogokupets_oganov_2007_vinet_4"
    source_temperature = 1500.0
    expected_surface_temperature = (4.0 * source_temperature - 295.0) / 3.0
    expected_pt_temperature = 0.97 * expected_surface_temperature

    forward = recalculate_pressure_calibration_path(
        100.0,
        "kcl_b2_chidester_2021_bm3_5",
        pt_identifier,
        source_temperature,
    )
    assert forward.edge_identifiers == (edge_identifier,)
    assert forward.target_pressure_gpa == pytest.approx(100.0)
    assert forward.intermediate_states[0]["source_temperature_k"] == pytest.approx(
        source_temperature
    )
    assert forward.intermediate_states[0]["target_temperature_k"] == pytest.approx(
        expected_pt_temperature
    )

    reverse = recalculate_pressure_calibration_path(
        forward.target_pressure_gpa,
        pt_identifier,
        "kcl_b2_chidester_2021_bm3_5",
        expected_pt_temperature,
    )
    assert reverse.target_pressure_gpa == pytest.approx(100.0)
    assert reverse.intermediate_states[0]["target_temperature_k"] == pytest.approx(
        source_temperature
    )

    gold = recalculate_pressure_calibration_path(
        100.0,
        "kcl_b2_chidester_2021_bm3_5",
        "gold_fei_2007_vinet_2",
        source_temperature,
    )
    assert gold.edge_identifiers[0] == edge_identifier
    assert "platinum_sokolova_2013_holzapfel_3" in gold.calibration_path
    assert gold.target_pressure_gpa == pytest.approx(103.66816140389605)

    with pytest.raises(ValidationError, match="requires temperature_k"):
        recalculate_pressure_calibration_path(
            100.0,
            "kcl_b2_chidester_2021_bm3_5",
            pt_identifier,
        )


def test_eosmat_rejects_invalid_calibration_links_and_missing_targets(monkeypatch):
    document = get_material_document("gold")
    method = document["eos_records"][0]["pressure_calibration"]["methods"][0]
    method["reference_calibration_record"] = ""
    with pytest.raises(EosmatError, match="must be a non-empty string"):
        eosmat_module.validate_eosmat_document(document)

    method["reference_calibration_record"] = "ruby_mao_1986"
    method["kind"] = "other"
    with pytest.raises(EosmatError, match="requires a ruby_fluorescence method"):
        eosmat_module.validate_eosmat_document(document)

    monkeypatch.setattr(eosmat_module, "list_pressure_calibrations", lambda: ())
    with pytest.raises(EosmatError, match="references missing pressure calibration"):
        eosmat_module.validate_pressure_calibration_references()


def test_eosmat_rejects_other_malformed_pressure_calibration_fields():
    base = get_material_document("gold")

    documents = []
    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["methods"][0]["kind"] = "invalid"
    documents.append(document)

    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["methods"][0][
        "source_location"
    ] = ""
    documents.append(document)

    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["methods"][0]["reference"] = 3
    documents.append(document)

    document = copy.deepcopy(base)
    method = document["eos_records"][0]["pressure_calibration"]["methods"][0]
    method["reference_eos_record"] = ""
    documents.append(document)

    document = copy.deepcopy(base)
    method = document["eos_records"][0]["pressure_calibration"]["methods"][0]
    method["kind"] = "equation_of_state"
    method.pop("reference_calibration_record")
    method.pop("reference")
    documents.append(document)

    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["recalculation"]["status"] = (
        "invalid"
    )
    documents.append(document)

    for document in documents:
        with pytest.raises(EosmatError):
            eosmat_module.validate_eosmat_document(document)
    (list_cross_calibration_edges,)
    (list_diamond_raman_calibrations,)
    (recalculate_diamond_raman_pressure,)
    (recalculate_pressure_calibration_path,)
