import json
from dataclasses import replace

import numpy as np
import pytest

import peritheos.materials as materials_module
from peritheos import get_material_document
from peritheos.materials import (
    AG_DEWAELE_2008,
    AG_SOKOLOVA_2013,
    AL_SOKOLOVA_2013,
    AU_DORFMAN_2012,
    AU_FEI_2007,
    AU_SOKOLOVA_2013,
    CBN_DATCHI_2007,
    CU_SOKOLOVA_2013,
    DEFERRED_EOS_RECORDS,
    DIAMOND_BENEDICT_2014,
    DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED,
    DIAMOND_CORREA_2008,
    DIAMOND_CORREA_2008_DEWAELE_ANCHORED,
    DIAMOND_DEWAELE_2008,
    DIAMOND_SOKOLOVA_2013,
    FEI_2007_EOS_RECORDS,
    KBR_B1_DEWAELE_2012,
    KBR_B2_DEWAELE_2012,
    KCL_B1_DEWAELE_2012,
    KCL_B2_DEWAELE_2012,
    LIF_B1_DEWAELE_2019,
    MGO_SOKOLOVA_2013,
    MGO_TANGE_2009,
    MO_DORFMAN_2012,
    MO_SOKOLOVA_2013,
    NACL_B1_DEWAELE_2019,
    NACL_B2_DEWAELE_2019,
    NACL_B2_DORFMAN_2012,
    NACL_B2_FEI_2007,
    NB_SOKOLOVA_2013,
    NE_DORFMAN_2012,
    NE_FEI_2007,
    NI_DEWAELE_2008,
    PT_DORFMAN_2012,
    PT_FEI_2007,
    PT_SOKOLOVA_2013,
    RE_HCP_ANZELLINI_2014,
    SOKOLOVA_2013_EOS_RECORDS,
    TA_SOKOLOVA_2013,
    W_SOKOLOVA_2013,
    EOSRecord,
    HugoniotRecord,
    Material,
    get_eos_record,
    get_material,
    list_eos_records,
    list_materials,
)

DORFMAN_EOS_RECORDS = (
    AU_DORFMAN_2012,
    PT_DORFMAN_2012,
    MO_DORFMAN_2012,
    NACL_B2_DORFMAN_2012,
    NE_DORFMAN_2012,
)


def test_catalog_listing_lookup_and_material_filter():
    records = list_eos_records()
    materials = list_materials()

    assert len(records) == 226
    assert len(materials) == 141
    assert all(isinstance(item, EOSRecord) for item in records)
    assert all(isinstance(item, Material) for item in materials)
    assert [item.identifier for item in records] == sorted(
        item.identifier for item in records
    )
    assert [item.identifier for item in materials] == sorted(
        item.identifier for item in materials
    )
    canonical_mgo = get_eos_record("mgo_b1_tange_2009_vinet")
    assert canonical_mgo.identifier == MGO_TANGE_2009.identifier
    assert canonical_mgo is not MGO_TANGE_2009
    assert {item.identifier for item in list_eos_records(formula="au")} == {
        "gold_anderson_1989_bm3_1",
        "gold_dewaele_2004_vinet_5",
        "gold_fei_2007_vinet_2",
        "gold_fratanduono_2021_vinet_7",
        "gold_shen_2026_vinet_3",
        "gold_sokolova_2013_holzapfel_4",
        "gold_takemura_2008_vinet_6",
    }
    assert get_material("mgo_b1").eos_records == (
        MGO_TANGE_2009,
        MGO_SOKOLOVA_2013,
    )
    assert len(get_material("diamond").eos_records) == 7
    assert get_material("au_fcc").get_eos_record("au_fcc_fei_2007") is AU_FEI_2007
    assert list_materials(formula="Au") == (get_material("gold"),)
    assert list_eos_records(formula="missing") == ()
    assert list_materials(formula="missing") == ()
    with pytest.raises(KeyError, match="Unknown EOS record"):
        get_eos_record("missing")
    with pytest.raises(KeyError, match="Unknown material"):
        get_material("missing")
    with pytest.raises(KeyError, match="available"):
        get_material("au_fcc").get_eos_record("missing")


def test_catalog_metadata_is_explicit_and_json_safe():
    document = get_material("mgo_b1").to_dict()
    record = document["eos_records"][0]

    assert document["format"] == "peritheos.material"
    assert document["format_version"] == 3
    assert document["identifier"] == "mgo_b1"
    assert all("roles" not in item for item in document["eos_records"])
    assert document["units"] == {
        "pressure": "GPa",
        "temperature": "K",
        "volume": "angstrom^3/conventional_unit_cell",
    }
    assert record["reference"]["doi"] == "10.1029/2008JB005813"
    assert record["eos"]["model"] == "vinet"
    assert record["eos"]["parameters"]["K0"] == 160.63
    assert record["thermal"]["model"] == ("asymptotic_power_law_mie_gruneisen_debye")
    assert "Table 4" in record["parameter_provenance"]["thermal_correction"]["gamma0"]
    assert record["scientific_validation"]["status"] == "primary_source_validated"
    json.dumps(document)
    for material in list_materials():
        json.dumps(material.to_dict(), allow_nan=False)


def test_material_document_records_actual_model_volume_basis():
    molar_record = get_material("au_fcc").to_dict()["eos_records"][1]
    cell_record = get_material("kcl_b2").to_dict()["eos_records"][0]

    assert molar_record["volume"]["model_unit"] == "J bar^-1 mol^-1"
    assert cell_record["volume"]["model_unit"] == ("angstrom^3/conventional_unit_cell")


def test_material_document_groups_scales_without_merging_equation_components():
    document = get_material("au_fcc").to_dict()

    assert document["formula"] == "Au"
    assert len(document["eos_records"]) == 3
    assert document["eos_records"][0]["eos"]["model"] == "holzapfel"
    assert (
        document["eos_records"][0]["thermal"]["model"]
        == "multi_oscillator_gruneisen_thermal_pressure"
    )
    assert document["eos_records"][1]["thermal"]["model"] == "mie_gruneisen_debye"
    assert document["eos_records"][1]["thermal"]["configuration"] == {
        "debye_temperature_law": "variable_exponent"
    }
    assert "thermal" not in document["eos_records"][2]
    assert document["eos_records"][1]["parameter_errors"]["K0_prime"] == 0.02
    assert document["eos_records"][1]["thermal"]["parameter_errors"]["gamma0"] == 0.03


def test_material_rejects_empty_duplicate_or_incompatible_records():
    gold = get_material("au_fcc")
    with pytest.raises(ValueError, match="identifier"):
        replace(gold, identifier="")
    with pytest.raises(ValueError, match="at least one"):
        replace(gold, eos_records=())
    with pytest.raises(ValueError, match="unique"):
        replace(gold, eos_records=(AU_FEI_2007, AU_FEI_2007))
    with pytest.raises(ValueError, match="match its material"):
        replace(gold, eos_records=(PT_FEI_2007,))


@pytest.mark.parametrize("material", list_materials())
def test_material_document_json_round_trip_reconstructs_catalog_material(material):
    payload = json.loads(json.dumps(material.to_dict(), allow_nan=False))
    loaded = Material.from_dict(payload)

    assert loaded.identifier == material.identifier
    assert loaded.formula == material.formula
    for loaded_record, record in zip(loaded.eos_records, material.eos_records):
        assert loaded_record.identifier == record.identifier
        assert type(loaded_record.eos) is type(record.eos)
        assert loaded_record.eos.parameter_values() == pytest.approx(
            record.eos.parameter_values()
        )
        loaded_configuration = (
            {}
            if isinstance(loaded_record, HugoniotRecord)
            else loaded_record.eos.configuration_values()
        )
        configuration = (
            {}
            if isinstance(record, HugoniotRecord)
            else record.eos.configuration_values()
        )
        assert loaded_configuration == configuration
        pressure = record.pressure(
            0.9 * record.reference_volume,
            record.reference_temperature,
            check_validity=False,
        )
        assert loaded_record.pressure(
            0.9 * loaded_record.reference_volume,
            loaded_record.reference_temperature,
            check_validity=False,
        ) == pytest.approx(pressure)


def test_eosmat_crystallography_is_preserved_while_selected_eos_is_executable():
    source = get_material_document("gold")
    source["source"] = {"kind": "cif", "text": "data_gold\n_cell_length_a 4.0862"}
    source["peaks"] = [[1, 1, 1, 2.3592, 100.0]]
    source["dioptas_extension"] = {"reflection_profile": "pseudo-Voigt"}
    identifier = "gold_fei_2007_vinet_2"

    material = Material.from_eosmat(
        source,
        record_identifiers=(identifier,),
    )
    document = material.to_eosmat()

    assert material.symmetry == source["symmetry"]
    assert dict(material.lattice) == source["lattice"]
    assert material.formula_units_per_cell == source["formula_units_per_cell"]
    assert material.space_group == source["space_group"]
    assert material.space_group_number == source["space_group_number"]
    assert [dict(site) for site in material.atom_sites] == source["atom_sites"]
    assert document["symmetry"] == source["symmetry"]
    assert document["lattice"] == source["lattice"]
    assert document["atom_sites"] == source["atom_sites"]
    assert document["space_group"] == source["space_group"]
    assert document["space_group_number"] == source["space_group_number"]
    assert document["source"] == source["source"]
    assert document["peaks"] == source["peaks"]
    assert document["dioptas_extension"] == source["dioptas_extension"]
    assert document.get("aliases", []) == source.get("aliases", [])
    record = material.get_eos_record(identifier)
    pressure = record.pressure(
        0.9 * record.reference_volume, 1000.0, check_validity=False
    )
    assert np.isfinite(pressure)


def test_eosmat_datasets_follow_selected_eos_records():
    source = get_material_document("gold")

    with_digitized_data = Material.from_eosmat(
        source, record_identifiers=("gold_fei_2007_vinet_2",)
    ).to_eosmat()
    assert [dataset["identifier"] for dataset in with_digitized_data["datasets"]] == [
        "gold_fei_2007_figure1_digitized"
    ]

    with_data = Material.from_eosmat(
        source, record_identifiers=("gold_dewaele_2004_vinet_5",)
    ).to_eosmat()
    assert [dataset["identifier"] for dataset in with_data["datasets"]] == [
        "gold_dewaele_2004_table1_compression"
    ]
    assert with_data["datasets"][0]["used_by_eos_records"] == [
        "gold_dewaele_2004_vinet_5"
    ]


def test_eosmat_record_provenance_and_extensions_survive_executable_round_trip():
    source = get_material_document("gold")
    identifier = "gold_fei_2007_vinet_2"
    source_record = next(
        record for record in source["eos_records"] if record["identifier"] == identifier
    )
    source_record["record_extension"] = {"reviewer": "laboratory"}
    source_record["eos"]["equation_extension"] = {"convention": "published"}
    source_record["thermal"]["thermal_extension"] = {"energy_zero": "Tr"}
    source_record.setdefault("validity", {})["validity_extension"] = (
        "joint P-T envelope"
    )
    source_record["parameter_provenance"] = {
        "reference_isotherm": {"K0": "Table 1"},
        "thermal_correction": {"gamma0": "Table 1"},
        "additional": {"volume_basis": "conventional cell"},
    }

    output = Material.from_eosmat(source, record_identifiers=(identifier,)).to_eosmat()[
        "eos_records"
    ][0]

    for key in (
        "default",
        "reference",
        "fixed_parameters",
        "parameter_errors",
        "parameter_provenance",
        "scientific_validation",
        "record_extension",
    ):
        assert output[key] == source_record[key]
    assert output["eos"]["equation_extension"] == {"convention": "published"}
    assert output["thermal"]["thermal_extension"] == {"energy_zero": "Tr"}
    assert output["validity"]["validity_extension"] == "joint P-T envelope"
    assert output["thermal"]["type"] == source_record["thermal"]["type"]


def test_legacy_snapshot_without_debye_law_uses_integrated_default():
    legacy_material = Material(
        identifier="diamond",
        formula=DIAMOND_DEWAELE_2008.material,
        phase=DIAMOND_DEWAELE_2008.phase,
        cell_contents=DIAMOND_DEWAELE_2008.cell_contents,
        eos_records=(DIAMOND_DEWAELE_2008,),
    )
    payload = legacy_material.to_snapshot_dict()
    record = next(
        item
        for item in payload["eos_records"]
        if item["identifier"] == DIAMOND_DEWAELE_2008.identifier
    )
    thermal = record["equation"]["thermal_correction"]
    assert thermal.pop("configuration") == {
        "debye_temperature_law": "integrated_gruneisen"
    }

    loaded = (
        Material.from_dict(payload).get_eos_record(DIAMOND_DEWAELE_2008.identifier).eos
    )

    assert loaded.debye_temperature_law == "integrated_gruneisen"


def test_material_document_loader_uses_model_registry_not_implementation_string():
    payload = json.loads(json.dumps(get_material("au_fcc").to_snapshot_dict()))
    reference = payload["eos_records"][1]["equation"]["reference_isotherm"]
    reference["implementation"] = "untrusted.module.ArbitraryCode"
    assert Material.from_dict(payload).identifier == "au_fcc"

    reference["model"] = "unknown_model"
    with pytest.raises(ValueError, match="Invalid EOS record"):
        Material.from_dict(payload)


def test_material_document_round_trips_component_qualified_covariance():
    covariance_record = replace(
        AU_FEI_2007,
        parameter_errors=None,
        parameter_covariance=((0.04, 0.001), (0.001, 0.09)),
        covariance_parameters=("rt_eos.K0_prime", "gamma0"),
    )
    material = replace(get_material("au_fcc"), eos_records=(covariance_record,))
    payload = json.loads(json.dumps(material.to_dict()))
    covariance = payload["eos_records"][0]["parameter_covariance"]
    assert covariance["parameter_order"] == ["rt_eos.K0_prime", "gamma0"]

    loaded = Material.from_dict(payload).eos_records[0]
    assert loaded.parameter_covariance == covariance_record.parameter_covariance
    assert loaded.covariance_parameters == covariance_record.covariance_parameters


def test_material_document_loader_rejects_unvalidated_composition():
    payload = get_material("au_fcc").to_dict()
    validation = payload["eos_records"][1]["scientific_validation"]
    validation["status"] = "pending_primary_source_check"
    with pytest.raises(ValueError, match="pending_primary_source_check"):
        Material.from_dict(payload)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("format", "other", "Unsupported material format"),
        ("format_version", 1, "format 3"),
    ],
)
def test_material_document_loader_rejects_wrong_format(field, value, message):
    payload = get_material("au_fcc").to_dict()
    payload[field] = value
    with pytest.raises(ValueError, match=message):
        Material.from_dict(payload)


@pytest.mark.parametrize("record", DORFMAN_EOS_RECORDS)
def test_dorfman_reference_state_identity_and_round_trip(record):
    assert record.pressure(
        record.reference_volume, check_validity=False
    ) == pytest.approx(0.0, abs=1.0e-13)

    pressure = sum(record.validity.pressure_gpa) / 2.0
    volume = record.volume(pressure)
    assert record.pressure(volume) == pytest.approx(pressure, rel=1.0e-10)


def test_tange_reference_state_identity_and_array_round_trip():
    assert MGO_TANGE_2009.pressure(74.698, 300.0) == pytest.approx(0.0, abs=1e-13)
    volumes = 74.698 * np.array([0.99, 0.90, 0.70])
    temperatures = np.array([500.0, 1500.0, 3000.0])

    pressures = MGO_TANGE_2009.pressure(volumes, temperatures)
    recovered = MGO_TANGE_2009.volume(pressures, temperatures)

    assert isinstance(pressures, np.ndarray)
    assert pressures.shape == (3,)
    assert np.allclose(recovered, volumes, rtol=1.0e-10)


@pytest.mark.parametrize("record", SOKOLOVA_2013_EOS_RECORDS)
def test_sokolova_reference_state_identity_and_round_trip(record):
    assert record.pressure(record.reference_volume, 298.15) == pytest.approx(
        0.0, abs=1.0e-12
    )

    volume = record.volume(100.0, 1500.0)
    assert record.pressure(volume, 1500.0) == pytest.approx(100.0, rel=1.0e-9)


@pytest.mark.parametrize("record", FEI_2007_EOS_RECORDS)
def test_fei_reference_state_identity_and_round_trip(record):
    assert record.pressure(
        record.reference_volume, 300.0, check_validity=False
    ) == pytest.approx(0.0, abs=1.0e-12)

    lower, upper = record.validity.pressure_gpa
    pressure = lower + 0.5 * (upper - lower)
    temperature = sum(record.validity.temperature_k) / 2.0
    volume = record.volume(pressure, temperature)
    assert record.pressure(volume, temperature) == pytest.approx(pressure, rel=1.0e-10)


@pytest.mark.parametrize("standard", [AU_SOKOLOVA_2013, AU_FEI_2007])
def test_added_thermal_families_broadcast_arrays(standard):
    pressures = np.array([20.0, 50.0, 80.0])
    temperatures = np.array([500.0, 1000.0, 1800.0])
    volumes = standard.volume(pressures, temperatures)
    recovered = standard.pressure(volumes, temperatures)

    assert volumes.shape == (3,)
    assert np.allclose(recovered, pressures, rtol=1.0e-10)


def test_added_thermal_family_calibration_range_guards_are_opt_in():
    assert np.isfinite(MGO_SOKOLOVA_2013.volume(401.0, 2000.0))
    assert np.isfinite(NE_FEI_2007.volume(50.0, 1100.0))
    with pytest.raises(ValueError, match="outside the published calibration"):
        MGO_SOKOLOVA_2013.volume(401.0, 2000.0, check_validity=True)
    with pytest.raises(ValueError, match="outside the published calibration"):
        NE_FEI_2007.volume(50.0, 1100.0, check_validity=True)


@pytest.mark.parametrize(
    ("standard", "V0", "K0", "K0_prime", "theta1", "theta2", "delta", "t"),
    [
        (MGO_SOKOLOVA_2013, 1.1248, 160.3, 4.10, 748.0, 401.0, -0.235, 0.301),
        (DIAMOND_SOKOLOVA_2013, 0.3414, 441.5, 3.90, 1561.0, 684.0, -0.506, 1.085),
        (AL_SOKOLOVA_2013, 0.998, 72.8, 4.51, 381.0, 202.0, -0.242, -0.958),
        (CU_SOKOLOVA_2013, 0.7112, 133.5, 5.32, 296.0, 169.0, -0.07, 1.401),
        (AG_SOKOLOVA_2013, 1.025, 100.0, 6.15, 199.0, 115.0, 0.178, 2.210),
        (AU_SOKOLOVA_2013, 1.0215, 167.0, 5.90, 179.5, 83.0, 0.134, 0.087),
        (PT_SOKOLOVA_2013, 0.9091, 275.0, 5.35, 177.0, 143.0, 0.167, -0.343),
        (NB_SOKOLOVA_2013, 1.0828, 170.5, 3.65, 302.0, 134.0, -0.326, -0.763),
        (TA_SOKOLOVA_2013, 1.0861, 191.0, 3.83, 254.0, 101.0, -0.101, -0.148),
        (MO_SOKOLOVA_2013, 0.9369, 260.0, 4.20, 353.0, 222.0, -0.802, -0.791),
        (W_SOKOLOVA_2013, 0.9552, 308.0, 4.12, 309.0, 172.0, -0.686, -0.591),
    ],
)
def test_sokolova_published_parameters(
    standard, V0, K0, K0_prime, theta1, theta2, delta, t
):
    eos = standard.eos
    assert eos.rt_eos.V0 == V0
    assert eos.rt_eos.K0 == K0
    assert eos.rt_eos.K0_prime == K0_prime
    assert eos.QE1o == theta1
    assert eos.QE2o == theta2
    assert eos.delta == delta
    assert eos.t == t
    assert "Table 1" in standard.parameter_provenance["rt_eos.V0"]
    assert "Table 4" in standard.parameter_provenance["rt_eos.K0"]


@pytest.mark.parametrize("standard", SOKOLOVA_2013_EOS_RECORDS)
def test_sokolova_records_use_scientific_fit_year_and_cite_2013(standard):
    assert standard.identifier.endswith("_sokolova_2013")
    assert standard.reference.year == 2013
    assert standard.reference.doi == "10.1016/j.rgg.2013.01.005"
    assert any("Table 1" in note and "Table 4" in note for note in standard.notes)
    assert any("Sokolova et al. (2016)" in note for note in standard.notes)


def test_sokolova_2016_material_aliases_are_removed():
    with pytest.raises(KeyError, match="Did you mean"):
        get_eos_record("diamond_sokolova_2016")
    assert not hasattr(materials_module, "DIAMOND_SOKOLOVA_2016")


def test_sokolova_mgo_records_2016_anharmonic_correction():
    assert "2016" in MGO_SOKOLOVA_2013.parameter_provenance["a_0"]
    assert "corrects" in MGO_SOKOLOVA_2013.parameter_provenance["a_0"]


@pytest.mark.parametrize(
    ("ratio", "expected"),
    [
        (0.95, 9.14084),
        (0.90, 20.9517),
        (0.85, 36.2876),
        (0.80, 56.3146),
        (0.75, 82.667),
        (0.70, 117.677),
        (0.65, 164.746),
        (0.60, 228.95),
    ],
)
def test_sokolova_mgo_figure2_spreadsheet_regression(ratio, expected):
    pressure = MGO_SOKOLOVA_2013.pressure(
        MGO_SOKOLOVA_2013.reference_volume * ratio,
        300.0,
        check_validity=False,
    )
    assert pressure == pytest.approx(expected, abs=0.00051)


def test_sokolova_mgo_effective_atomic_number_and_inactive_terms():
    assert MGO_SOKOLOVA_2013.eos.rt_eos.Z == pytest.approx(10.34)
    assert MGO_SOKOLOVA_2013.eos.e_0 == 0.0
    assert "inactive" in MGO_SOKOLOVA_2013.parameter_provenance["e_0"]


def test_fei_table1_parameters_and_published_errors():
    assert AU_FEI_2007.eos.rt_eos.K0 == 167.0
    assert AU_FEI_2007.eos.rt_eos.K0_prime == 6.0
    assert PT_FEI_2007.eos.gamma0 == 2.72
    assert NACL_B2_FEI_2007.eos.theta0 == 290.0
    assert NE_FEI_2007.eos.rt_eos.K0 == 1.16
    assert NE_FEI_2007.parameter_errors["rt_eos.K0"] == 0.14
    neon_records = {
        record["identifier"]: record
        for record in get_material_document("neon_fcc")["eos_records"]
    }
    assert neon_records["neon_fcc_fei_2007_bm3_1"]["fixed_parameters"] == ["V0"]
    assert neon_records["neon_fcc_fei_2007_vinet_2"]["fixed_parameters"] == ["V0"]
    for identifier in (
        "neon_fcc_fei_2007_bm3_1",
        "neon_fcc_fei_2007_vinet_2",
    ):
        record = neon_records[identifier]
        assert record["fit_datasets"] == [
            "neon_fei_2007_figure5_digitized",
            "neon_hemley_1989_table1_fei_recalculated",
        ]
        assert [item["source"] for item in record["fit_data_sources"]] == [
            "Fei et al. (2007), this study",
            "Hemley et al. (1989), Fei reference 45",
            "Finger et al. (1981), Fei reference 47",
        ]
        assert record["fit_data_sources"][-1]["availability"] == "not_yet_bundled"


@pytest.mark.parametrize(
    "standard",
    [
        LIF_B1_DEWAELE_2019,
        NACL_B1_DEWAELE_2019,
        NACL_B2_DEWAELE_2019,
        KCL_B1_DEWAELE_2012,
        KCL_B2_DEWAELE_2012,
        KBR_B1_DEWAELE_2012,
        KBR_B2_DEWAELE_2012,
        CBN_DATCHI_2007,
        DIAMOND_DEWAELE_2008,
        NI_DEWAELE_2008,
        AG_DEWAELE_2008,
    ],
)
def test_added_standards_reference_identity_and_round_trip(standard):
    assert standard.pressure(
        standard.reference_volume,
        standard.reference_temperature,
        check_validity=False,
    ) == pytest.approx(0.0, abs=1.0e-12)

    lower, upper = standard.validity.pressure_gpa
    pressure = lower + 0.4 * (upper - lower)
    temperature = sum(standard.validity.temperature_k) / 2.0
    volume = standard.volume(pressure, temperature)
    assert standard.pressure(volume, temperature) == pytest.approx(
        pressure, rel=1.0e-10
    )


@pytest.mark.parametrize(
    ("pressure", "volume_per_formula_unit"),
    [(0.722, 16.181), (53.8, 11.3973), (109.3, 9.74253)],
)
def test_lif_dewaele_2019_table4_regression(pressure, volume_per_formula_unit):
    calculated = LIF_B1_DEWAELE_2019.pressure(
        4.0 * volume_per_formula_unit, check_validity=False
    )
    assert calculated == pytest.approx(pressure, abs=0.4)


@pytest.mark.parametrize(
    ("pressure", "lattice_parameter", "temperature"),
    [(4.53, 3.5534, 298.0), (34.7, 3.4893, 750.0), (54.8, 3.4515, 900.0)],
)
def test_diamond_dewaele_2008_table1_regression(
    pressure, lattice_parameter, temperature
):
    calculated = DIAMOND_DEWAELE_2008.pressure(
        lattice_parameter**3, temperature, check_validity=False
    )
    assert calculated == pytest.approx(pressure, abs=0.5)


@pytest.mark.parametrize(
    ("volume_per_atom", "temperature", "expected_pressure"),
    [
        (5.7034, 300.0, 5.097381463860107),
        (5.4, 1000.0, 34.72410455615707),
        (4.654270411587497, 3000.0, 150.0),
    ],
)
def test_diamond_benedict_2014_library_regression(
    volume_per_atom, temperature, expected_pressure
):
    volume_per_cell = 8.0 * volume_per_atom
    calculated = DIAMOND_BENEDICT_2014.pressure(volume_per_cell, temperature)

    assert calculated == pytest.approx(expected_pressure, rel=5.0e-11)
    assert DIAMOND_BENEDICT_2014.volume(calculated, temperature) == pytest.approx(
        volume_per_cell, rel=1.0e-11
    )


@pytest.mark.parametrize(
    ("volume_per_atom", "temperature", "expected_pressure"),
    [
        (5.785, 300.0, 4.710783105788316),
        (5.571, 300.0, 20.018830109151665),
        (4.43, 5000.0, 202.62811519774186),
        (3.21, 9000.0, 729.9463078775337),
    ],
)
def test_diamond_correa_2008_library_regression(
    volume_per_atom, temperature, expected_pressure
):
    volume_per_cell = 8.0 * volume_per_atom
    calculated = DIAMOND_CORREA_2008.pressure(volume_per_cell, temperature)

    assert calculated == pytest.approx(expected_pressure, rel=5.0e-11)
    assert DIAMOND_CORREA_2008.volume(calculated, temperature) == pytest.approx(
        volume_per_cell, rel=1.0e-11
    )


@pytest.mark.parametrize(
    ("anchored", "absolute"),
    [
        (DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED, DIAMOND_BENEDICT_2014),
        (DIAMOND_CORREA_2008_DEWAELE_ANCHORED, DIAMOND_CORREA_2008),
    ],
)
@pytest.mark.parametrize(("volume", "temperature"), [(40.0, 298.0), (35.0, 3000.0)])
def test_diamond_simulation_models_can_be_anchored_to_dewaele(
    anchored, absolute, volume, temperature
):
    reference_pressure = DIAMOND_DEWAELE_2008.pressure(volume, 298.0)
    expected = reference_pressure + (
        absolute.pressure(volume, temperature) - absolute.pressure(volume, 298.0)
    )

    assert anchored.pressure(volume, 298.0) == pytest.approx(reference_pressure)
    assert anchored.pressure(volume, temperature) == pytest.approx(expected)
    assert anchored.volume(expected, temperature) == pytest.approx(volume)


def test_diamond_benedict_cold_reference_is_not_zero_total_pressure():
    record = DIAMOND_BENEDICT_2014

    assert record.reference_volume == pytest.approx(8.0 * 5.7034)
    assert record.pressure(record.reference_volume, 300.0) == pytest.approx(
        5.097381463860107
    )
    assert "motionless-ion 0 K cold curve" in " ".join(record.notes)


def test_diamond_benedict_record_dac_volume_pair_uses_cell_volumes():
    record = DIAMOND_BENEDICT_2014
    expected_temperature = 2000.0
    f_dac = 0.25
    heated_volume = 8.0 * 5.0
    reference_pressure = record.pressure(heated_volume, 300.0)
    thermal_increment = (
        record.pressure(heated_volume, expected_temperature) - reference_pressure
    )
    ambient_pressure = (
        record.pressure(heated_volume, expected_temperature) - f_dac * thermal_increment
    )
    ambient_volume = record.volume(ambient_pressure, 300.0)

    recovered = record.temperature_from_volumes(
        ambient_volume,
        heated_volume,
        f_dac=f_dac,
    )

    assert recovered == pytest.approx(expected_temperature, rel=1.0e-11)


def test_diamond_benedict_record_dac_forward_state_uses_cell_volume():
    record = DIAMOND_BENEDICT_2014
    cold_pressure = 60.0
    temperature = 2000.0
    f_dac = 0.25

    heated_volume = record.volume_with_dac_confinement(
        cold_pressure,
        temperature,
        f_dac=f_dac,
        check_validity=True,
    )
    thermal_increment = record.thermal_pressure_increment(
        heated_volume, temperature, check_validity=True
    )
    confinement_pressure = record.dac_thermal_pressure(
        heated_volume,
        temperature,
        f_dac=f_dac,
        check_validity=True,
    )

    assert heated_volume > 8.0 * 4.0
    assert confinement_pressure == pytest.approx(f_dac * thermal_increment, rel=1.0e-12)
    assert record.pressure(heated_volume, temperature) == pytest.approx(
        cold_pressure + confinement_pressure,
        rel=1.0e-11,
    )


def test_isothermal_record_rejects_dac_forward_methods():
    with pytest.raises(ValueError, match="isothermal"):
        CBN_DATCHI_2007.volume_with_dac_confinement(60.0, 2000.0, f_dac=0.25)
    with pytest.raises(ValueError, match="isothermal"):
        CBN_DATCHI_2007.thermal_pressure_increment(40.0, 2000.0)


@pytest.mark.parametrize(
    ("standard", "pressure", "volume_per_atom", "tolerance"),
    [
        (AG_DEWAELE_2008, 1.8, 16.7256, 0.2),
        (AG_DEWAELE_2008, 104.0, 11.803, 0.3),
        (NI_DEWAELE_2008, 1.5, 10.8504, 0.1),
        (NI_DEWAELE_2008, 65.0, 8.8219, 1.0),
    ],
)
def test_dewaele_2008_metals_table2_regression(
    standard, pressure, volume_per_atom, tolerance
):
    calculated = standard.pressure(4.0 * volume_per_atom, check_validity=False)
    assert calculated == pytest.approx(pressure, abs=tolerance)


def test_kcl_kbr_thermal_pressure_is_exactly_the_published_linear_term():
    kcl_volume = 40.0
    kbr_volume = 50.0
    assert KCL_B2_DEWAELE_2012.pressure(
        kcl_volume, 2000.0, check_validity=False
    ) - KCL_B2_DEWAELE_2012.pressure(
        kcl_volume, 300.0, check_validity=False
    ) == pytest.approx(0.00224 * 1700.0)
    assert KBR_B2_DEWAELE_2012.pressure(
        kbr_volume, 3000.0, check_validity=False
    ) - KBR_B2_DEWAELE_2012.pressure(
        kbr_volume, 300.0, check_validity=False
    ) == pytest.approx(0.00222 * 2700.0)


@pytest.mark.parametrize(
    ("ratio", "temperature", "expected"),
    [
        (1.00, 500.0, 1.06),
        (1.00, 3000.0, 16.76),
        (0.99, 300.0, 1.65),
        (0.95, 1500.0, 16.49),
        (0.90, 3000.0, 38.58),
        (0.80, 2000.0, 68.91),
        (0.65, 300.0, 168.81),
    ],
)
def test_tange_table5_printed_pressure_regression(ratio, temperature, expected):
    """Tange et al. (2009), Table 5, Fit3-Vinet columns, in GPa."""
    pressure = MGO_TANGE_2009.pressure(
        74.698 * ratio, temperature, check_validity=False
    )

    assert pressure == pytest.approx(expected, abs=0.0051)


def test_validity_checks_distinguish_extrapolation_from_evaluation():
    assert AU_DORFMAN_2012.within_validity(50.0, 300.0)
    assert AU_DORFMAN_2012.within_calibration_range(50.0, 300.0)
    assert not AU_DORFMAN_2012.within_validity(AU_DORFMAN_2012.reference_volume, 300.0)
    validity = AU_DORFMAN_2012.within_validity(np.array([50.0, 67.85]), 300.0)
    assert np.array_equal(validity, [True, False])

    with pytest.raises(ValueError, match="outside the published calibration"):
        AU_DORFMAN_2012.pressure(AU_DORFMAN_2012.reference_volume, check_validity=True)
    assert AU_DORFMAN_2012.pressure(AU_DORFMAN_2012.reference_volume) == pytest.approx(
        0.0
    )


@pytest.mark.parametrize("volume", [0.0, -1.0, np.nan, np.inf])
def test_eos_record_rejects_invalid_volume(volume):
    with pytest.raises(ValueError):
        MGO_TANGE_2009.pressure(volume, 1000.0, check_validity=False)


@pytest.mark.parametrize("temperature", [299.0, 301.0, np.nan, -1.0])
def test_isothermal_eos_record_rejects_unsupported_temperature(temperature):
    with pytest.raises(ValueError):
        AU_DORFMAN_2012.pressure(50.0, temperature)


def test_static_zero_kelvin_isothermal_record_accepts_only_zero_temperature():
    record = Material.from_eosmat(
        get_material_document("phase_egg"),
        record_identifiers=["phase_egg_mookherjee_2019_bm3_lp_1"],
    ).eos_records[0]

    assert record.reference_temperature == 0.0
    assert record.pressure(196.0) == pytest.approx(14.72595528543105)
    assert record.pressure(196.0, 0.0) == pytest.approx(14.72595528543105)
    with pytest.raises(ValueError, match="isothermal 0 K"):
        record.pressure(196.0, 1.0)
    with pytest.raises(ValueError, match="greater than zero"):
        record.pressure(196.0, -1.0)


def test_pressure_uncertainty_propagates_parameters_volume_and_temperature():
    volume = 74.698 * 0.9
    baseline = MGO_TANGE_2009.pressure_with_uncertainty(volume, 3000.0)
    measured = MGO_TANGE_2009.pressure_with_uncertainty(
        volume,
        3000.0,
        volume_sigma=0.01,
        temperature_sigma=50.0,
        full_covariance=True,
    )

    assert measured.value == pytest.approx(MGO_TANGE_2009.pressure(volume, 3000.0))
    assert measured.standard_error > baseline.standard_error
    assert measured.covariance.shape == (1, 1)
    assert any("state-variable" in item for item in measured.assumptions)
    assert any("independent" in item for item in measured.assumptions)


def test_isothermal_uncertainty_propagates_volume_but_rejects_temperature_sigma():
    prediction = AU_DORFMAN_2012.pressure_with_uncertainty(50.0, volume_sigma=0.01)
    assert prediction.value == pytest.approx(AU_DORFMAN_2012.pressure(50.0))
    assert prediction.standard_error > 0.0

    with pytest.raises(ValueError, match="isothermal EOS record"):
        AU_DORFMAN_2012.pressure_with_uncertainty(50.0, temperature_sigma=5.0)


def test_measurement_uncertainty_works_when_publication_has_no_parameter_errors():
    volume = KCL_B1_DEWAELE_2012.volume(1.0)
    prediction = KCL_B1_DEWAELE_2012.pressure_with_uncertainty(
        volume, volume_sigma=0.01
    )

    assert prediction.standard_error > 0.0
    assert "published parameter uncertainty not available" in prediction.assumptions


def test_sokolova_measurement_uncertainty_is_state_only():
    volume = MGO_SOKOLOVA_2013.volume(100.0, 1500.0)
    prediction = MGO_SOKOLOVA_2013.pressure_with_uncertainty(
        volume,
        1500.0,
        volume_sigma=0.01,
        temperature_sigma=20.0,
    )
    assert prediction.standard_error > 0.0
    assert "published parameter uncertainty not available" in prediction.assumptions


def test_volume_uncertainty_is_returned_in_public_unit_cell_volume():
    pressure = 50.0
    prediction = MGO_TANGE_2009.volume_with_uncertainty(
        pressure,
        1500.0,
        pressure_sigma=0.2,
        temperature_sigma=20.0,
        full_covariance=True,
    )

    assert isinstance(prediction.value, float)
    assert isinstance(prediction.standard_error, float)
    assert prediction.value == pytest.approx(MGO_TANGE_2009.volume(pressure, 1500.0))
    assert prediction.standard_error > 0.0
    assert prediction.covariance.shape == (1, 1)


def test_isothermal_volume_uncertainty_temperature_guard():
    with pytest.raises(ValueError, match="isothermal EOS record"):
        AU_DORFMAN_2012.volume_with_uncertainty(50.0, temperature_sigma=5.0)


@pytest.mark.parametrize(
    ("reported_pressure", "a", "c"),
    [
        (0.65, 2.7619, 4.4590),
        (36.9, 2.6847, 4.3322),
        (100.0, 2.5947, 4.1847),
        (144.0, 2.5537, 4.1147),
    ],
)
def test_anzellini_rhenium_table3_lattice_data(reported_pressure, a, c):
    volume = np.sqrt(3.0) * a**2 * c / 2.0
    calculated = RE_HCP_ANZELLINI_2014.pressure(volume, check_validity=False)

    # Table III gives rounded lattice parameters and reports pressure
    # uncertainty rising to about 2 GPa at 150 GPa.
    assert calculated == pytest.approx(reported_pressure, abs=3.0)


def test_anzellini_95_percent_parameter_intervals_are_not_treated_as_sigma():
    uncertainty = RE_HCP_ANZELLINI_2014._uncertainty().parameter_uncertainty

    assert RE_HCP_ANZELLINI_2014.parameter_error_confidence == 0.95
    assert uncertainty.standard_errors["K0"] == pytest.approx(8.0 / 1.95996398454)
    assert uncertainty.standard_errors["K0_prime"] == pytest.approx(
        0.17 / 1.95996398454
    )


def test_deferred_catalog_records_evidence_gaps():
    assert any(item.material == "Re" for item in DEFERRED_EOS_RECORDS)
    assert all(item.references and item.reason for item in DEFERRED_EOS_RECORDS)
