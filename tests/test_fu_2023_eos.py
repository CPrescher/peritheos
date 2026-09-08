import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
CA_RECORD = "ca_perovskite_fu_2023_bm3_mgd_refit"
CA_REFIT_RECORD = (
    "ca_perovskite_fu_2023_candidate_data_unweighted_bm3_mgd_refit"
)


def test_fu_2023_reproduction_and_schemas():
    path = ROOT / "scripts/reproduce_fu_2023_eos.py"
    spec = importlib.util.spec_from_file_location("reproduce_fu_2023", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()["records"]
    assert len(result) == 4
    ca = result[module.CA_RECORD]
    assert ca["analytical_checkpoints"]["fit_observations"] is False
    assert ca["coefficient_audit"]["status"] == "parity_not_achieved"
    assert ca["coefficient_audit"]["candidate_fit_observations"] == 174
    assert ca["coefficient_audit"]["fit_outputs"] == 242
    refit = result[module.CA_REFIT_RECORD]
    assert refit["record_kind"] == "refit"
    assert refit["objective"] == "unweighted_absolute_gpa"
    assert refit["fit_observations"] == 174
    assert refit["fit_outputs"] == 242
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    validator = Draft202012Validator(schema)
    for material in ("mg088fe010al014si090o3_bridgmanite", "ca_perovskite"):
        assert list(validator.iter_errors(get_material_document(material))) == []


def test_fu_2023_casio3_audit_keeps_observations_and_checkpoints_distinct():
    audit = json.loads(
        (ROOT / "docs/data/fu-2023-casio3-refit-audit.json").read_text(encoding="utf-8")
    )
    assert audit["source_rows"] == {
        "greaux_2019_figure3a_tetragonal_excluded": 13,
        "greaux_2019_figure3b_cubic": 34,
        "sun_2016_excluded_above_2200_k": 4,
        "sun_2016_figure_s3_temperature_subset": 140,
        "sun_2016_table1_total": 144,
    }
    assert audit["fit_observation_count"] == 140 + 34
    assert audit["fit_output_count"] == 140 + 3 * 34
    assert audit["weighting_disclosure"] == "not published by Fu et al."
    assert audit["registered_refit"] == {
        "objective": "unweighted_absolute_gpa",
        "record_identifier": CA_REFIT_RECORD,
        "status": "opt-in Peritheos refit; not a reproduction of Fu's undisclosed regression protocol",
    }
    assert {item["objective"] for item in audit["sensitivity_fits"]} == {
        "unweighted_absolute_gpa",
        "reported_sigma",
        "propagated_pv_sigma_no_temperature",
        "equal_observable_groups",
    }
    assert all(item["success"] for item in audit["sensitivity_fits"])
    assert all(
        source["redistributed"] is False
        for source in audit["sources"].values()
        if source.get("role", "").startswith("candidate fit observations")
    )

    document = get_material_document("ca_perovskite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == CA_RECORD
    )
    assert "fit_datasets" not in record
    check = record["scientific_validation"]["primary_data_check"]
    assert check["status"] == "externally_recovered_not_redistributed"
    assert check["observation_count"] == 174

    refit = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == CA_REFIT_RECORD
    )
    diagnostic = next(
        item
        for item in audit["sensitivity_fits"]
        if item["objective"] == "unweighted_absolute_gpa"
    )
    assert refit["record_kind"] == "refit"
    assert refit["derived_from_record"] == CA_RECORD
    assert refit["fit_provenance"]["weights"] == []
    assert refit["eos"]["parameters"] == pytest.approx(
        {
            "V0": diagnostic["parameters"]["V0"],
            "K0": diagnostic["parameters"]["K0"],
            "K0_prime": 4.0,
        }
    )
    assert refit["thermal"]["parameters"]["gamma0"] == pytest.approx(
        diagnostic["parameters"]["gamma0"]
    )
    assert refit["thermal"]["parameters"]["q"] == pytest.approx(
        diagnostic["parameters"]["q"]
    )


def test_fu_2023_casio3_reconstructed_pressure_equation_matches_catalog_model():
    path = ROOT / "scripts/audit_fu_2023_casio3_refit.py"
    spec = importlib.util.spec_from_file_location("audit_fu_2023_casio3", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    volume = np.array([38.0])
    temperature = np.array([1200.0])
    reconstructed = module._predictions(module.PUBLISHED, volume, temperature)[0][0]
    model = Material.from_eosmat(
        get_material_document("ca_perovskite"), record_identifiers=[CA_RECORD]
    ).get_eos_record(CA_RECORD)
    assert reconstructed == pytest.approx(model.pressure(38.0, 1200.0), rel=2e-12)


def test_fu_2023_audit_enumerates_all_candidates():
    text = (ROOT / "docs/literature-reproductions/fu-2023-eos.md").read_text(
        encoding="utf-8"
    )
    identifiers = {
        "litcurate_5661c5f835e32b74",
        "litcurate_efb817fde2eec8ae",
        "litcurate_1347395c4b56d29a",
        "litcurate_a70856250f14a2ee",
        "litcurate_7a71575938f3f9fb",
        "litcurate_a2591f9f56f627a0",
        "litcurate_438c7d5970e1cd0f",
        "litcurate_a6ab9ce1494d1f92",
        "litcurate_b601d5b21bd6f007",
        "litcurate_a3d1728dde300a01",
        "litcurate_6f0c0a415873bd49",
        "litcurate_dc5880a87d2a1522",
        "litcurate_f150ba8e2eb3c5db",
        "litcurate_7e1c43165d898ecb",
        "litcurate_57faf342bab1540f",
        "litcurate_c570565a82d28af4",
        "litcurate_62c01575e6a31000",
        "litcurate_8f536afdfadf777c",
        "litcurate_b46a77f6e6173d23",
    }
    assert len(identifiers) == 19
    assert all(identifier in text for identifier in identifiers)
