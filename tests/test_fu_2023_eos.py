import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator
from scipy.optimize import least_squares

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
CA_RECORD = "ca_perovskite_fu_2023_bm3_mgd_refit"
CA_REFIT_RECORD = "ca_perovskite_fu_2023_candidate_data_unweighted_bm3_mgd_refit"


def test_bundled_fu_observations_preserve_sources_and_conversions(fu_audit_script):
    document = get_material_document("ca_perovskite")
    dataset = next(
        d
        for d in document["datasets"]
        if d["identifier"] == "ca_perovskite_fu_2023_candidate_composite_external"
    )
    resource = dataset["resource"]
    path = ROOT / "peritheos/data" / resource["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == resource["sha256"]
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        assert reader.fieldnames == [c["name"] for c in dataset["columns"]]
    assert len(rows) == 174
    sun = [r for r in rows if r["source_id"] == "sun_2016_table1"]
    greaux = [r for r in rows if r["source_id"] == "greaux_2019_figure3b"]
    assert len(sun) == 140 and len(greaux) == 34
    assert max(float(r["temperature_k"]) for r in sun) == 2200
    assert min(float(r["pressure_gpa"]) for r in rows) == 12.05
    assert max(float(r["pressure_gpa"]) for r in rows) == 151.8
    assert all(r["bulk_modulus_gpa"] == r["shear_modulus_gpa"] == "" for r in sun)
    raw_sun, raw_greaux = fu_audit_script.load_bundled_observations()
    assert len(raw_sun) == 144 and len(raw_greaux) == 34
    assert raw_greaux[0]["pressure_nacl_gpa"] == 14.48
    assert raw_greaux[-1]["temperature_k"] == 1700
    for row, source in zip(greaux, raw_greaux):
        rho = source["density_g_cm3"]
        assert float(row["pressure_gpa"]) == source["pressure_nacl_gpa"]
        assert float(row["volume_a3"]) == pytest.approx(116.162 / (0.602214076 * rho))
        assert float(row["volume_standard_deviation_a3"]) == pytest.approx(
            float(row["volume_a3"]) * source["density_sigma_g_cm3"] / rho
        )
        assert float(row["shear_modulus_gpa"]) == pytest.approx(
            rho * source["vs_km_s"] ** 2
        )
        assert float(row["bulk_modulus_gpa"]) == pytest.approx(
            rho * (source["vp_km_s"] ** 2 - 4 / 3 * source["vs_km_s"] ** 2)
        )
    assert dataset["used_by_eos_records"] == [CA_RECORD, CA_REFIT_RECORD]
    published = next(r for r in document["eos_records"] if r["identifier"] == CA_RECORD)
    assert published["supporting_datasets"] == [dataset["identifier"]]
    assert "fit_datasets" not in published  # Fu's actual selection is undisclosed.


def test_bundled_fu_data_reproduce_unweighted_refit(fu_audit_script):
    sun, greaux = fu_audit_script.load_bundled_observations()
    result = fu_audit_script.fit_variant(
        fu_audit_script._arrays(sun, greaux), "unweighted_absolute_gpa"
    )
    audit = json.loads((ROOT / "docs/data/fu-2023-casio3-refit-audit.json").read_text())
    expected = next(
        f
        for f in audit["sensitivity_fits"]
        if f["objective"] == "unweighted_absolute_gpa"
    )
    assert result["success"]
    assert result["residual_sum_squares"] == pytest.approx(
        expected["residual_sum_squares"], rel=1e-8
    )
    assert result["parameters"] == pytest.approx(expected["parameters"], rel=2e-5)


@pytest.fixture
def fu_audit_script():
    path = ROOT / "scripts/audit_fu_2023_casio3_refit.py"
    spec = importlib.util.spec_from_file_location("audit_fu_2023_casio3", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_refit_standard_errors_use_residual_degrees_of_freedom(fu_audit_script):
    # Seven independent linear coefficients with known sensitivities, plus
    # three residuals orthogonal to the design: RSS=14, dof=10-7=3.
    slopes = np.arange(1.0, 8.0)
    fit = least_squares(
        lambda x: np.r_[slopes * (x - 1.0), [1.0, 2.0, 3.0]], np.zeros(7)
    )
    uncertainty = fu_audit_script._fit_uncertainty(fit)
    expected_errors = np.sqrt(14.0 / 3.0) / slopes
    assert uncertainty["degrees_of_freedom"] == 3
    assert uncertainty["standard_errors"] == pytest.approx(
        dict(zip(fu_audit_script.PARAMETERS, expected_errors))
    )
    np.testing.assert_allclose(
        uncertainty["parameter_covariance"]["matrix"], np.diag(expected_errors**2)
    )


@pytest.mark.parametrize("problem", ["rank", "dof", "bounds", "convergence"])
def test_refit_rejects_unsupported_standard_errors(fu_audit_script, problem):
    fit = least_squares(lambda x: np.r_[x, [1.0]], np.ones(7))
    if problem == "rank":
        fit.jac[:, -1] = 0.0
    elif problem == "dof":
        fit.fun = fit.fun[:7]
        fit.jac = fit.jac[:7]
    elif problem == "bounds":
        fit.active_mask[0] = 1
    else:
        fit.success = False
    with pytest.raises(ValueError, match="parameter errors require"):
        fu_audit_script._fit_uncertainty(fit)


def test_fu_refit_covariance_is_available_through_public_api(fu_audit_script):
    audit = json.loads(
        (ROOT / "docs/data/fu-2023-casio3-refit-audit.json").read_text(encoding="utf-8")
    )
    fit = next(
        item
        for item in audit["sensitivity_fits"]
        if item["objective"] == "unweighted_absolute_gpa"
    )
    assert set(fit["standard_errors"]) == set(fu_audit_script.PARAMETERS)
    assert all(value > 0 for value in fit["standard_errors"].values())
    assert fit["degrees_of_freedom"] == 242 - 7
    joint = np.asarray(fit["parameter_covariance"]["matrix"])
    assert np.linalg.eigvalsh(joint).min() > 0
    document = get_material_document("ca_perovskite")
    stored = next(
        r for r in document["eos_records"] if r["identifier"] == CA_REFIT_RECORD
    )
    expected = json.loads(json.dumps(stored))
    # Regeneration must preserve the committed record, including parameter names
    # used by the thermal API and the marginal block of the full joint fit.
    fu_audit_script.update_refit_record(document, audit)
    assert stored == expected
    model = Material.from_eosmat(
        document, record_identifiers=[CA_REFIT_RECORD]
    ).get_eos_record(CA_REFIT_RECORD)
    assert model.covariance_parameters == ("rt_eos.V0", "rt_eos.K0", "gamma0", "q")
    names = ("V0", "K0", "gamma0", "q")
    indices = [fit["parameter_covariance"]["parameter_order"].index(n) for n in names]
    np.testing.assert_allclose(
        model.parameter_covariance, joint[np.ix_(indices, indices)]
    )
    assert dict(model.parameter_errors) == pytest.approx(
        {
            api_name: fit["standard_errors"][name]
            for api_name, name in zip(model.covariance_parameters, names)
        }
    )
    for component in (stored, stored["thermal"]):
        for name in component["fixed_parameters"]:
            assert component["parameter_errors"][name] is None
    for temperature in (300.0, 1200.0):
        prediction = model.pressure_with_uncertainty(38.0, temperature)
        assert np.isfinite(prediction.standard_error)
        assert prediction.standard_error > 0
    assert (
        stored["fit_provenance"]["weighting_description"]
        == fu_audit_script.EQUAL_WEIGHTING
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
        source["redistributed"] is True
        for source in audit["sources"].values()
        if source.get("role", "").startswith("candidate fit observations")
    )

    document = get_material_document("ca_perovskite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == CA_RECORD
    )
    assert "fit_datasets" not in record
    check = record["scientific_validation"]["primary_data_check"]
    assert check["status"] == "bundled"
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
