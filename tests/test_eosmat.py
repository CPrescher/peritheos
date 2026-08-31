import copy
import json
import math
import re
from importlib import resources
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    EOSMAT_FORMAT,
    EOSMAT_FORMAT_VERSION,
    Material,
    eosmat_schema,
    get_material_document,
    list_material_documents,
    list_materials,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
)
from scripts.import_dioptas_eos_database import migrate_document

NINETY_FIVE_PERCENT_INTERVAL_RECORDS = (
    ("alumina", "alumina_dewaele_2013_vinet_1"),
    ("aluminum", "aluminum_dewaele_2004_vinet_1"),
    ("copper", "copper_dewaele_2004_vinet_1"),
    ("diamond", "diamond_dewaele_2008_vinet_2"),
    ("lif_b1", "lif_b1_dewaele_2019_vinet_1"),
    ("nickel", "nickel_dewaele_2008_vinet_1"),
    ("rhenium", "rhenium_anzellini_2014_vinet_1"),
    ("silicon_v", "silicon_v_anzellini_2019_vinet_1"),
    ("silicon_vii", "silicon_vii_anzellini_2019_vinet_1"),
    ("silicon_x", "silicon_x_anzellini_2019_vinet_1"),
    ("silver", "silver_dewaele_2008_vinet_1"),
    ("titanium_alpha", "titanium_alpha_dewaele_2015_vinet_1"),
    ("titanium_omega", "titanium_omega_dewaele_2015_vinet_1"),
    ("tungsten", "tungsten_dewaele_2004_vinet_2"),
)


def test_complete_migrated_dioptas_library_is_bundled_and_valid():
    identifiers = list_material_documents()
    documents = [get_material_document(identifier) for identifier in identifiers]

    assert len(identifiers) == 116
    assert len(set(identifiers)) == 116
    assert sum(len(document["eos_records"]) for document in documents) == 147
    assert all(document["eos_records"] for document in documents)
    assert all(document["format"] == EOSMAT_FORMAT for document in documents)
    assert all(
        document["format_version"] == EOSMAT_FORMAT_VERSION for document in documents
    )
    assert all(document["identifier"] in identifiers for document in documents)
    assert all(
        document["units"]
        == {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        }
        for document in documents
    )
    json.dumps(documents, allow_nan=False)


def test_migrated_records_have_completed_primary_source_audit():
    records = [
        record
        for identifier in list_material_documents()
        for record in get_material_document(identifier)["eos_records"]
    ]

    assert len({record["identifier"] for record in records}) == 147
    statuses = [record["scientific_validation"]["status"] for record in records]
    assert set(statuses) == {"primary_source_validated", "deferred"}
    assert statuses.count("primary_source_validated") == 116
    assert statuses.count("deferred") == 31
    assert all(
        record["scientific_validation"]["audit_date"] == "2026-08-31"
        for record in records
    )
    assert all(
        record["scientific_validation"]["primary_source_check"] for record in records
    )
    assert all(
        len({item["path"] for item in record.get("audit_corrections", [])})
        == len(record.get("audit_corrections", []))
        for record in records
    )
    assert {
        record["scientific_validation"]["migration_source"]["version"]
        for record in records
    } == {"0.10.0"}
    assert {
        record["scientific_validation"]["migration_source"]["commit"]
        for record in records
    } == {"5a8bfd81d10bfab3499039603380aae34576d60a"}
    assert all(record["eos"]["model"] for record in records)
    assert all(
        not record.get("thermal") or record["thermal"]["model"] for record in records
    )


def test_primary_source_audit_report_covers_every_migrated_record():
    report = json.loads(
        resources.files("peritheos.data")
        .joinpath("primary-source-audit.json")
        .read_text(encoding="utf-8")
    )
    bundled_ids = {
        record["identifier"]
        for identifier in list_material_documents()
        for record in get_material_document(identifier)["eos_records"]
    }

    assert report["summary"] == {
        "records": 147,
        "deferred": 31,
        "primary_source_validated": 116,
    }
    assert {entry["record"] for entry in report["records"]} == bundled_ids
    assert len(report["records"]) == len(bundled_ids)


def test_every_primary_validated_migrated_record_is_executable():
    failures = []
    checked = 0
    for identifier in list_material_documents():
        document = get_material_document(identifier)
        for record in document["eos_records"]:
            if record["scientific_validation"]["status"] != "primary_source_validated":
                continue
            checked += 1
            try:
                Material.from_eosmat(
                    document, record_identifiers=[record["identifier"]]
                )
            except (TypeError, ValueError) as error:
                failures.append(f"{record['identifier']}: {error}")

    assert checked == 116
    assert failures == []


def test_primary_audit_records_corrections_and_known_source_limitations():
    graphite = get_material_document("graphite")["eos_records"][0]
    b4c = get_material_document("b4c")["eos_records"][0]
    zircon = get_material_document("zircon")["eos_records"][0]
    shen = next(
        record
        for record in get_material_document("gold")["eos_records"]
        if record["identifier"] == "gold_shen_2026_vinet_3"
    )

    assert graphite["parameter_errors"]["V0"] == pytest.approx(0.02)
    assert graphite["audit_corrections"][0]["primary_reference"]["location"] == (
        "page 12599, Murnaghan-fit paragraph"
    )
    assert b4c["scientific_validation"]["reported_inconsistencies"][0] == {
        "field": "eos.order",
        "abstract": "third-order Birch-Murnaghan",
        "figure_1_caption": "second-order Birch-Murnaghan",
        "resolution": (
            "Retain BM3 because K0'=3.3(1), which is incompatible with a "
            "conventional BM2 fit where K0'=4 is fixed."
        ),
    }
    assert zircon["identifier"] == "zircon_hazen_1979_bm3_1"
    assert zircon["eos"] == {
        "type": "BM3",
        "parameters": {"V0": 260.79, "K0": 227.0, "K0_prime": 6.5},
        "model": "birch_murnaghan_3",
    }
    assert zircon["parameter_errors"]["V0"] == pytest.approx(0.04)
    assert zircon["fixed_parameters"] == ["K0_prime"]
    assert zircon["scientific_validation"]["status"] == "primary_source_validated"
    assert shen["scientific_validation"]["status"] == "deferred"
    assert "24 April 2027" in shen["scientific_validation"]["note"]

    magnesite = get_material_document("magnesite")["eos_records"][0]
    molybdenum_carbide = get_material_document("molybenum_carbide_mo2c")["eos_records"][
        0
    ]
    platinum = get_material_document("platinum")["eos_records"][0]
    alumina = get_material_document("alumina")["eos_records"][0]
    cobalt = get_material_document("cobalt_hcp")["eos_records"][0]
    niobium = get_material_document("niobium")["eos_records"][0]
    tantalum = get_material_document("tantalum")["eos_records"][0]

    assert magnesite["eos"]["parameters"]["V0"] == pytest.approx(279.41)
    assert magnesite["parameter_errors"]["V0"] == pytest.approx(0.08)
    assert molybdenum_carbide["eos"]["parameters"]["V0"] == pytest.approx(148.9071)
    assert molybdenum_carbide["parameter_errors"]["V0"] == pytest.approx(0.049)
    assert molybdenum_carbide["experimental_pressure_range_gpa"] == [0.0, 46.0]
    assert platinum["identifier"] == "platinum_holmes_1989_vinet_1"
    assert platinum["eos"]["type"] == "Vinet"
    assert platinum["thermal"]["type"] == "LinearThermalPressure"
    assert platinum["experimental_pressure_range_gpa"] == [0.0, 550.0]
    assert alumina["fixed_parameters"] == []
    assert cobalt["eos"]["type"] == "Vinet"
    assert niobium["experimental_pressure_range_gpa"] == [0.0, 71.5]
    assert tantalum["eos"]["parameters"]["V0"] == pytest.approx(36.0835)
    assert tantalum["fixed_parameters"] == ["V0"]


def test_newly_validated_primary_records_retain_published_errors():
    shim, mao = get_material_document("ca_perovskite")["eos_records"]
    cao_b1 = get_material_document("cao")["eos_records"][0]
    cao_b2 = get_material_document("cao_b2")["eos_records"][0]
    geo2 = get_material_document("geo2_rutile")["eos_records"][0]
    sno2 = get_material_document("sno2")["eos_records"][0]
    pbs = get_material_document("pbs_b1")["eos_records"][0]
    wadsleyite = get_material_document("wadsleyite")["eos_records"][0]
    jadeite = get_material_document("naalsi2o6")["eos_records"][0]
    perovskite = get_material_document("perovskite_orthorhombic")["eos_records"][0]
    kcl = next(
        record
        for record in get_material_document("kcl")["eos_records"]
        if record["identifier"] == "kcl_walker_2002_bm3_2"
    )

    # Shim et al. fixed V0 in the Table 2 fit but explicitly report its
    # uncertainty in the abstract, so it must not be discarded.
    assert shim["parameter_errors"] == {
        "V0": pytest.approx(0.05),
        "K0": pytest.approx(4.0),
        "K0_prime": pytest.approx(0.2),
    }
    assert mao["parameter_errors"] == {
        "V0": pytest.approx(0.08),
        "K0": pytest.approx(4.0),
        "K0_prime": None,
    }
    assert cao_b1["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(1.0),
        "K0_prime": pytest.approx(0.2),
    }
    assert cao_b2["parameter_errors"] == {
        "V0": pytest.approx(0.3321),
        "K0": pytest.approx(20.0),
        "K0_prime": pytest.approx(0.5),
    }
    assert geo2["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(5.0),
        "K0_prime": None,
    }
    assert sno2["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(2.0),
        "K0_prime": None,
    }
    assert pbs["parameter_errors"] == {
        "V0": pytest.approx(0.0421),
        "K0": pytest.approx(1.2),
        "K0_prime": pytest.approx(0.9),
    }
    assert wadsleyite["parameter_errors"] == {
        "V0": pytest.approx(0.02),
        "K0": pytest.approx(0.9),
        "K0_prime": pytest.approx(0.1),
    }
    assert wadsleyite["thermal"]["parameter_errors"] == {
        "Tr": None,
        "theta0": None,
        "gamma0": pytest.approx(0.02),
        "q": pytest.approx(0.1),
        "n": None,
    }
    assert jadeite["parameter_errors"] == {
        "V0": pytest.approx(0.08),
        "K0": pytest.approx(4.0),
        "K0_prime": None,
    }
    assert perovskite["parameter_errors"] == {
        "V0": pytest.approx(0.39),
        "K0": pytest.approx(6.0),
        "K0_prime": pytest.approx(0.4),
    }
    assert kcl["thermal"]["parameter_errors"] == {
        "Tr": None,
        "alpha_KT": pytest.approx(0.00009),
    }


@pytest.mark.parametrize(
    ("material_identifier", "record_identifier"),
    NINETY_FIVE_PERCENT_INTERVAL_RECORDS,
)
def test_published_95_percent_intervals_are_normalized_when_loaded(
    material_identifier, record_identifier
):
    document = get_material_document(material_identifier)
    source_record = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == record_identifier
    )
    loaded_record = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).eos_records[0]

    assert source_record["parameter_error_confidence"] == 0.95
    assert loaded_record.parameter_error_confidence == 0.95

    standard_errors = loaded_record._uncertainty().parameter_uncertainty.standard_errors
    for name, interval_half_width in source_record["parameter_errors"].items():
        if interval_half_width is None:
            continue
        loaded_name = name if name in standard_errors else f"rt_eos.{name}"
        assert standard_errors[loaded_name] == pytest.approx(
            interval_half_width / 1.95996398454
        )


def test_walker_2002_kcl_be1_and_reported_product_error_regression():
    document = get_material_document("kcl")
    record = Material.from_eosmat(
        document, record_identifiers=["kcl_walker_2002_bm3_2"]
    ).eos_records[0]

    # Table 2 reports V=47.57(3) A^3 at 36.3(2) kbar and 23 degC.
    assert record.pressure(47.57, 296.15) == pytest.approx(3.63, abs=0.02)

    # Equation BE1 adds alpha0*K0*delta-T. The paper reports the identifiable
    # product as 0.0275(9) kbar/K = 0.00275(9) GPa/K.
    pressure_ambient = record.pressure(47.57, 296.15)
    pressure_600_c = record.pressure(47.57, 873.15)
    assert pressure_600_c - pressure_ambient == pytest.approx(0.00275 * 577.0)
    prediction = record.pressure_with_uncertainty(47.57, 873.15)
    assert prediction.standard_error == pytest.approx(0.00009 * 577.0)


def test_shim_2000_casio3_table_i_pressure_regression():
    document = get_material_document("ca_perovskite")
    record = Material.from_eosmat(
        document, record_identifiers=["ca_perovskite_shim_2000_bm3_1"]
    ).eos_records[0]

    # Table 1 reports V=35.909(55) A^3 at 89.6(3.0) GPa. The rounded values
    # happen to lie almost exactly on the published BM3 fit.
    assert record.pressure(35.909) == pytest.approx(89.6, abs=0.01)


def test_holmes_1989_platinum_equations_11_and_12_regression():
    document = get_material_document("platinum")
    identifier = "platinum_holmes_1989_vinet_1"
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]
    x = 0.9
    volume = 60.4000884 * x**3
    pressure_300 = (
        3.0 * 266.0 * (1.0 - x) / x**2 * math.exp(1.5 * (5.81 - 1.0) * (1.0 - x))
    )

    assert record.pressure(volume, 300.0) == pytest.approx(pressure_300)
    assert record.pressure(volume, 1000.0) == pytest.approx(
        pressure_300 + 0.0069426 * 700.0
    )


def test_ross_1997_magnesite_table_i_pressure_regression():
    record = Material.from_eosmat(get_material_document("magnesite")).eos_records[0]

    # Table I reports V=272.33(4) A^3 at 3.09 GPa. This is a measured state,
    # so compare with the published fit within the experimental residual.
    assert record.pressure(272.33, check_validity=False) == pytest.approx(
        3.09, abs=0.01
    )


def test_frank_2004_ice_vii_table_i_pressure_regression():
    document = get_material_document("ice_vii")
    identifier = "ice_vii_frank_2004_bm3_2"
    record = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]
    # Table I experiment 10: 8.315 cm^3/mol at 20.65(38) GPa. Convert the
    # molar volume to the two-formula-unit cubic cell used by .eosmat.
    volume = 8.315 * 1.0e24 * 2.0 / 6.02214076e23

    assert record.pressure(volume) == pytest.approx(20.65, abs=0.38)


def test_primary_audit_restores_model_inputs_omitted_during_migration():
    aluminum = get_material_document("aluminum")["eos_records"][1]
    silica = get_material_document("silica_cacl2")["eos_records"][0]
    ice = get_material_document("ice_vi")["eos_records"][0]

    assert aluminum["eos"]["parameters"]["n"] == 1.0
    assert aluminum["eos"]["parameters"]["Z"] == 13.0
    assert aluminum["thermal"]["parameters"]["n"] == 1.0
    assert silica["thermal"]["parameters"]["n"] == 3.0
    assert ice["thermal"]["parameters"]["Tr"] == 300.0
    assert {item["path"] for item in aluminum["audit_corrections"]} >= {
        "eos.parameters.n",
        "eos.parameters.Z",
        "thermal.parameters.n",
    }


@pytest.mark.parametrize(
    ("material", "volume", "reported_pressure", "tolerance"),
    [
        ("coesite", 522.61, 5.19, 0.05),
        ("zircon", 255.59, 4.81, 0.08),
    ],
)
def test_no_doi_primary_source_table_regressions(
    material, volume, reported_pressure, tolerance
):
    record = Material.from_eosmat(get_material_document(material)).eos_records[0]

    # These are measured P-V rows rather than values generated by the fit.
    assert record.pressure(volume, check_validity=False) == pytest.approx(
        reported_pressure, abs=tolerance
    )


def test_no_doi_primary_sources_use_stable_article_or_report_locators():
    for material in ("coesite", "lead_fcc", "zircon"):
        record = get_material_document(material)["eos_records"][0]
        validation = record["scientific_validation"]
        assert validation["status"] == "primary_source_validated"
        assert validation["primary_source_check"]["doi"] is None
        assert validation["primary_source_check"]["access_url"].startswith("https://")


@pytest.mark.parametrize(
    ("material", "volume", "temperature", "reported_pressure", "tolerance"),
    [
        ("ice_vi", 206.233, 340.7, 2.56, 0.02),
        ("ice_vii", 32.406, 300.6, 8.12, 0.06),
    ],
)
def test_bezacier_2014_table_i_pressure_regression(
    material, volume, temperature, reported_pressure, tolerance
):
    document = get_material_document(material)
    record_identifier = next(
        record["identifier"]
        for record in document["eos_records"]
        if "bezacier_2014" in record["identifier"]
    )
    record = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).eos_records[0]

    # Table I contains measured states, so comparison is to the reported fit
    # within its experimental residual rather than exact numeric identity.
    assert record.pressure(volume, temperature) == pytest.approx(
        reported_pressure, abs=tolerance
    )


def test_bezacier_2014_uncertainty_includes_state_and_published_fit_errors():
    record = Material.from_eosmat(get_material_document("ice_vi")).eos_records[0]

    prediction = record.pressure_with_uncertainty(
        206.233,
        340.7,
        volume_sigma=0.004,
        temperature_sigma=0.4,
    )

    assert prediction.value == pytest.approx(2.5511132321)
    assert prediction.standard_error > 0.0
    assert "parameter errors treated as mutually independent" in prediction.assumptions
    assert "state-variable errors treated as independent" in prediction.assumptions


def test_dioptas_010_shape_remains_accepted_without_peritheos_extensions():
    document = get_material_document("gold")
    document.pop("format")
    document.pop("identifier")
    document["format_version"] = 2
    for record in document["eos_records"]:
        record.pop("identifier")
        record.pop("scientific_validation")

    validate_eosmat_document(document)


def test_peritheos_eos_only_material_does_not_require_crystal_structure():
    document = get_material_document("gold")
    for key in (
        "symmetry",
        "lattice",
        "formula_units_per_cell",
        "space_group",
        "space_group_number",
        "atom_sites",
        "source",
        "peaks",
    ):
        document.pop(key, None)

    validate_eosmat_document(document)


def test_material_document_returns_a_defensive_copy_and_reports_unknown_id():
    first = get_material_document("gold")
    first["name"] = "changed"
    assert get_material_document("gold")["name"] == "Gold"
    with pytest.raises(KeyError, match="Unknown material document"):
        get_material_document("not_a_material")


def test_eosmat_file_round_trip(tmp_path):
    document = get_material_document("mgo")
    path = tmp_path / "mgo.eosmat"

    save_eosmat(path, document)

    assert load_eosmat(path) == document


def test_normative_schema_is_bundled():
    schema = eosmat_schema()

    assert schema["$schema"].endswith("2020-12/schema")
    assert schema["properties"]["format_version"]["const"] == 3
    assert schema["properties"]["eos_records"]["type"] == "array"
    assert schema["additionalProperties"] is True
    thermal = schema["$defs"]["thermal"]["properties"]
    assert thermal["debye_temperature_law"]["default"] == "integrated_gruneisen"
    assert thermal["debye_temperature_law"]["enum"] == [
        "integrated_gruneisen",
        "variable_exponent",
    ]
    assert len(schema["$defs"]["equation"]["allOf"][0]["oneOf"]) == 10
    assert len(schema["$defs"]["thermal"]["allOf"][0]["oneOf"]) == 8


def test_normative_schema_validates_every_bundled_document():
    schema = eosmat_schema()
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)

    failures = []
    for identifier in list_material_documents():
        for error in validator.iter_errors(get_material_document(identifier)):
            failures.append(
                f"{identifier}:{'/'.join(map(str, error.path))}: {error.message}"
            )
    assert failures == []


def test_normative_schema_validates_every_executable_material_export():
    validator = Draft202012Validator(eosmat_schema())
    failures = []
    for material in list_materials():
        for error in validator.iter_errors(material.to_eosmat()):
            failures.append(
                f"{material.identifier}:{'/'.join(map(str, error.path))}: "
                f"{error.message}"
            )
    assert failures == []


def test_documented_complete_eosmat_example_is_structurally_valid():
    documentation = Path("docs/eosmat-schema.md").read_text(encoding="utf-8")
    example = re.search(
        r"## Complete EOS-only example.*?```json\n(.*?)\n```",
        documentation,
        flags=re.DOTALL,
    )
    assert example is not None

    document = json.loads(example.group(1))
    validate_eosmat_document(document)
    assert document["eos_records"][0]["thermal"]["debye_temperature_law"] == (
        "variable_exponent"
    )


def test_fei_2007_migrations_use_published_variable_exponent_debye_law():
    records = [
        record
        for identifier in list_material_documents()
        for record in get_material_document(identifier)["eos_records"]
        if isinstance(record["reference"], dict)
        and record["reference"].get("doi") == "10.1073/pnas.0609013104"
        and record.get("thermal") is not None
    ]

    assert {record["identifier"] for record in records} == {
        "gold_fei_2007_vinet_2",
        "neon_fcc_fei_2007_vinet_2",
    }
    assert {
        (
            record["thermal"]["type"],
            record["thermal"]["model"],
            record["thermal"]["debye_temperature_law"],
        )
        for record in records
    } == {("MieGruneisenDebye", "mie_gruneisen_debye", "variable_exponent")}
    assert all(
        record["migration_corrections"][0]["primary_reference"]["location"]
        == "Equation 3 and the definition immediately following it"
        for record in records
    )

    # Other MGD records retain the conventional integrated constant-q law.
    wadsleyite = get_material_document("wadsleyite")["eos_records"][0]
    assert wadsleyite["thermal"]["type"] == "MieGruneisenDebye"
    assert wadsleyite["thermal"]["model"] == "mie_gruneisen_debye"
    assert "debye_temperature_law" not in wadsleyite["thermal"]


@pytest.mark.parametrize(
    ("doi", "expected_law", "corrected"),
    [
        (
            "https://doi.org/10.1073/pnas.0609013104",
            "variable_exponent",
            True,
        ),
        (
            "10.0000/unrelated",
            None,
            False,
        ),
    ],
)
def test_importer_corrects_only_primary_identified_fei_thermal_records(
    tmp_path, doi, expected_law, corrected
):
    source = {
        "format_version": 2,
        "name": "Example",
        "formula": "X",
        "eos_records": [
            {
                "label": "Example record",
                "reference": {"authors": ["Author"], "year": 2007, "doi": doi},
                "eos": {"type": "Vinet", "parameters": {"V0": 1.0}},
                "parameter_errors": {},
                "fixed_parameters": [],
                "thermal": {
                    "type": "MieGruneisenDebye",
                    "parameters": {
                        "Tr": 300.0,
                        "theta0": 170.0,
                        "gamma0": 2.97,
                        "q": 0.6,
                        "n": 1,
                    },
                },
            }
        ],
    }
    path = tmp_path / "example.json"
    path.write_text(json.dumps(source), encoding="utf-8")

    record = migrate_document(path)["eos_records"][0]

    assert record["thermal"]["type"] == "MieGruneisenDebye"
    assert record["thermal"]["model"] == "mie_gruneisen_debye"
    assert record["thermal"].get("debye_temperature_law") == expected_law
    assert ("migration_corrections" in record) is corrected


def test_eosmat_debye_temperature_law_defaults_and_validation():
    document = get_material_document("gold")
    fei = document["eos_records"][1]["thermal"]
    fei.pop("debye_temperature_law")
    validate_eosmat_document(document)

    fei["debye_temperature_law"] = "unknown"
    with pytest.raises(ValueError, match="debye_temperature_law is invalid"):
        validate_eosmat_document(document)

    fei["debye_temperature_law"] = "variable_exponent"
    fei["type"] = "AlphaKT"
    fei["model"] = "thermal_reference_state"
    with pytest.raises(ValueError, match="requires MieGruneisenDebye"):
        validate_eosmat_document(document)


def test_migration_manifest_and_dioptas_license_are_bundled():
    root = resources.files("peritheos.data.materials")
    manifest = json.loads(root.joinpath("manifest.json").read_text(encoding="utf-8"))
    license_text = root.joinpath("DIOPTAS_LICENSE.txt").read_text(encoding="utf-8")

    assert manifest["source"]["project"] == "Dioptas"
    assert manifest["source"]["version"] == "0.10.0"
    assert manifest["materials"] == 116
    assert manifest["eos_records"] == 147
    assert "Copyright (c) 2021-2026 Clemens Prescher" in license_text


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda d: d.update(format="unknown"), "Unsupported material format"),
        (lambda d: d.update(format_version=2), "format 3"),
        (lambda d: d.update(name=None), "name must be a string"),
        (lambda d: d.update(identifier=3), "identifier must be a string"),
        (lambda d: d.update(lattice=[]), "lattice must be a JSON object"),
        (
            lambda d: d["lattice"].update(a=float("nan")),
            "lattice.a must be finite",
        ),
        (
            lambda d: d.update(formula_units_per_cell=0),
            "formula_units_per_cell must be greater than zero",
        ),
        (lambda d: d.update(peaks={}), "peaks must be a JSON array"),
        (
            lambda d: d["eos_records"][0].update(label=None),
            "label must be a string",
        ),
        (
            lambda d: d["eos_records"][0].update(reference=None),
            "reference must be a string or object",
        ),
        (
            lambda d: d["eos_records"][0]["eos"].update(type="unknown"),
            "unsupported",
        ),
        (
            lambda d: d["eos_records"][0]["eos"]["parameters"].pop("V0"),
            "requires V0",
        ),
        (
            lambda d: d["eos_records"][0].update(fixed_parameters={}),
            "fixed_parameters must be an array",
        ),
        (
            lambda d: d["eos_records"][0].update(
                experimental_pressure_range_gpa=[2, 1]
            ),
            "must be ordered",
        ),
    ],
)
def test_structural_validator_rejects_invalid_documents(mutation, message):
    document = copy.deepcopy(get_material_document("gold"))
    mutation(document)
    with pytest.raises(ValueError, match=message):
        validate_eosmat_document(document)


def test_validator_rejects_duplicate_record_id_and_multiple_defaults():
    document = get_material_document("gold")
    document["eos_records"][1]["identifier"] = document["eos_records"][0]["identifier"]
    with pytest.raises(ValueError, match="Duplicate EOS record identifier"):
        validate_eosmat_document(document)

    document = get_material_document("gold")
    document["eos_records"][0]["default"] = True
    document["eos_records"][1]["default"] = True
    with pytest.raises(ValueError, match="at most one default"):
        validate_eosmat_document(document)
