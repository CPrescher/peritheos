"""EOS-record, serialization, fitting, and uncertainty tests for Hugoniots."""

import copy
from dataclasses import replace

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import eosmat_schema, get_material_document, validate_eosmat_document
from peritheos.errors import EosmatError, MaterialError, UnsupportedOperationError
from peritheos.fitting import fit_linear_us_up
from peritheos.hugoniot import LinearUsUpHugoniot
from peritheos.materials import (
    EOSRecord,
    HugoniotBranchDomain,
    HugoniotInitialState,
    HugoniotRecord,
    HugoniotVolumeBasis,
    LiteratureReference,
    Material,
    ValidityRange,
)


@pytest.fixture
def hugoniot_record():
    return HugoniotRecord(
        identifier="x_alpha_shock_1",
        name="Synthetic alpha-phase principal Hugoniot",
        material="X",
        phase="alpha",
        cell_contents="1 formula unit per conventional unit cell",
        eos=LinearUsUpHugoniot(V0=10.0, rho0=8.0, c0=4.0, s=1.5),
        reference_temperature=300.0,
        reference=LiteratureReference(
            authors="Example et al.",
            year=2026,
            title="Synthetic test data",
            doi="",
            locations=("test fixture",),
        ),
        validity=ValidityRange(
            pressure_gpa=(0.0, 200.0),
            temperature_k=(300.0, 300.0),
            volume_ratio=(0.65, 1.0),
        ),
        parameter_provenance={
            "V0": "initial state",
            "rho0": "initial state",
            "c0": "linear Us-up fit",
            "s": "linear Us-up fit",
            "P0": "initial state",
        },
        parameter_errors={"c0": 0.04, "s": 0.02},
        loading_path="principal",
        branch_kind="untransformed",
        initial_state=HugoniotInitialState(
            phase="alpha",
            temperature_k=300.0,
            pressure_gpa=0.0,
            density_g_cm3=8.0,
            material_identifier="x_alpha",
        ),
        volume_basis=HugoniotVolumeBasis(
            formula_units=1.0, molar_mass_g_mol=48.17712608
        ),
        branch_domain=HugoniotBranchDomain(
            particle_velocity_km_s=(0.0, 3.0),
            kind="phase_stability",
            boundary_status="reported_exactly",
        ),
    )


def test_hugoniot_is_an_eos_record_with_path_specific_methods(hugoniot_record):
    assert hugoniot_record.equation_kind == "hugoniot"
    assert hugoniot_record.is_hugoniot
    assert not hugoniot_record.is_isothermal
    volume = 8.0
    assert hugoniot_record.pressure(volume) > 0.0
    assert hugoniot_record.shock_velocity(volume) > 0.0
    assert hugoniot_record.particle_velocity(volume) > 0.0
    assert hugoniot_record.density(volume) == pytest.approx(10.0)
    assert hugoniot_record.specific_internal_energy_change(volume) > 0.0
    assert hugoniot_record.tangent_modulus(volume) > 0.0
    state = hugoniot_record.state_from_particle_velocity(1.0)
    assert state.volume == pytest.approx(
        hugoniot_record.eos.volume_from_particle_velocity(1.0)
    )
    assert state.pressure == pytest.approx(
        hugoniot_record.eos.pressure_from_particle_velocity(1.0)
    )
    assert hugoniot_record.pressure(volume, temperature=300.0) == pytest.approx(
        hugoniot_record.pressure(volume)
    )
    with pytest.raises(MaterialError, match="initial-state metadata"):
        hugoniot_record.pressure(volume, temperature=301.0)


def test_equilibrium_record_rejects_hugoniot_only_operation():
    from peritheos.eos.rt import BM3

    record = EOSRecord(
        identifier="equilibrium",
        name="Equilibrium",
        material="X",
        phase="alpha",
        cell_contents="cell",
        eos=BM3(10.0, 100.0, 4.0),
        reference_temperature=300.0,
        reference=LiteratureReference("A", 2026, "T", "", ("test",)),
        validity=ValidityRange((0.0, 10.0), (300.0, 300.0)),
        parameter_provenance={},
    )
    assert not hasattr(record, "shock_velocity")


def test_hugoniot_eosmat_round_trip_and_schema(hugoniot_record):
    material = Material(
        identifier="x_alpha",
        name="Synthetic X",
        formula="X",
        phase="alpha",
        cell_contents="1 formula unit per conventional unit cell",
        eos_records=(hugoniot_record,),
        formula_units_per_cell=1.0,
    )
    document = material.to_eosmat()
    record_document = document["eos_records"][0]
    assert record_document["equation_kind"] == "hugoniot"
    assert record_document["loading_path"] == "principal"
    assert record_document["branch_kind"] == "untransformed"
    assert record_document["initial_state"] == {
        "phase": "alpha",
        "temperature_k": 300.0,
        "pressure_gpa": 0.0,
        "density_g_cm3": 8.0,
        "material_identifier": "x_alpha",
    }
    assert record_document["volume_basis"] == {
        "formula_units": 1.0,
        "molar_mass_g_mol": 48.17712608,
        "kind": "formula_units",
    }
    assert record_document["branch_domain"] == {
        "particle_velocity_km_s": [0.0, 3.0],
        "kind": "phase_stability",
        "boundary_status": "reported_exactly",
        "notes": [],
    }
    assert "temperature_ref" not in record_document
    assert "reference_isotherm" not in record_document["parameter_provenance"]
    Draft202012Validator(eosmat_schema()).validate(document)
    restored = Material.from_eosmat(document)
    assert restored.eos_records[0].is_hugoniot
    assert isinstance(restored.eos_records[0], HugoniotRecord)
    assert restored.hugoniot_records == restored.eos_records
    assert restored.equilibrium_records == ()
    assert restored.eos_records[0].pressure(8.0) == pytest.approx(
        hugoniot_record.pressure(8.0)
    )

    invalid = copy.deepcopy(document)
    invalid["eos_records"][0]["thermal"] = copy.deepcopy(
        document["eos_records"][0]["eos"]
    )
    assert list(Draft202012Validator(eosmat_schema()).iter_errors(invalid))
    with pytest.raises(EosmatError, match="cannot have a thermal component"):
        validate_eosmat_document(invalid)

    invalid = copy.deepcopy(document)
    invalid["eos_records"][0]["initial_state"]["density_g_cm3"] = 7.9
    with pytest.raises(EosmatError, match="must match eos.parameters.rho0"):
        validate_eosmat_document(invalid)

    invalid = copy.deepcopy(document)
    invalid["eos_records"][0]["volume_basis"]["formula_units"] = 2.0
    with pytest.raises(EosmatError, match="must match formula_units_per_cell"):
        validate_eosmat_document(invalid)

    invalid = copy.deepcopy(document)
    invalid["eos_records"][0]["volume_basis"]["molar_mass_g_mol"] = 40.0
    with pytest.raises(EosmatError, match="same mass basis"):
        validate_eosmat_document(invalid)

    invalid = copy.deepcopy(document)
    invalid["eos_records"][0]["loading_path"] = "reshock"
    assert list(Draft202012Validator(eosmat_schema()).iter_errors(invalid))
    with pytest.raises(EosmatError, match="loading_path is invalid"):
        validate_eosmat_document(invalid)

    invalid = copy.deepcopy(document)
    del invalid["formula_units_per_cell"]
    assert list(Draft202012Validator(eosmat_schema()).iter_errors(invalid))
    with pytest.raises(EosmatError, match="formula_units_per_cell is required"):
        validate_eosmat_document(invalid)

    with pytest.raises(UnsupportedOperationError, match="Snapshot-v2"):
        material.to_snapshot_dict()


def test_hugoniot_eosmat_validation_rejects_inconsistent_metadata(hugoniot_record):
    material = Material(
        identifier="x_alpha",
        name="Synthetic X",
        formula="X",
        phase="alpha",
        cell_contents="1 formula unit per conventional unit cell",
        eos_records=(hugoniot_record,),
        formula_units_per_cell=1.0,
    )
    source = material.to_eosmat()

    def reject(path, value, message, *, delete=False):
        document = copy.deepcopy(source)
        target = document
        for key in path[:-1]:
            target = target[key]
        if delete:
            del target[path[-1]]
        else:
            target[path[-1]] = value
        with pytest.raises(EosmatError, match=message):
            validate_eosmat_document(document)

    record = ("eos_records", 0)
    parameters = (*record, "eos", "parameters")
    initial_state = (*record, "initial_state")
    volume_basis = (*record, "volume_basis")
    branch_domain = (*record, "branch_domain")
    invalid_cases = [
        ((*record, "equation_kind"), "unknown", "equation_kind is invalid", False),
        (
            (*record, "equation_kind"),
            "isothermal",
            "equation_kind does not match",
            False,
        ),
        ((*record, "default_for"), "unknown", "default_for is invalid", False),
        (
            (*record, "default_for"),
            "equilibrium",
            "default_for does not match",
            False,
        ),
        ((*parameters, "c0"), None, "parameters requires c0", True),
        ((*parameters, "c0"), 0.0, "must be greater than zero", False),
        ((*record, "branch_kind"), "unknown", "branch_kind is invalid", False),
        ((*initial_state, "phase"), "", "phase must be a non-empty string", False),
        (
            (*initial_state, "material_identifier"),
            "",
            "material_identifier must be a non-empty string",
            False,
        ),
        (
            (*initial_state, "eos_record_identifier"),
            "",
            "eos_record_identifier must be a non-empty string",
            False,
        ),
        ((*initial_state, "phase"), "beta", "untransformed branches", False),
        ((*record, "branch_kind"), "transformed", "distinct precursor", False),
        ((*initial_state, "pressure_gpa"), 1.0, "pressure_gpa must match", False),
        ((*record, "temperature_ref"), 301.0, "temperature_ref must match", False),
        ((*volume_basis, "kind"), "atoms", "kind must be 'formula_units'", False),
        ((*volume_basis, "formula_units"), 2.0, "formula_units must match", False),
        ((*record, "loading_path"), "precompressed", "require positive P0", False),
        (
            (*branch_domain, "particle_velocity_km_s"),
            [-1.0, 3.0],
            "particle velocity must be non-negative",
            False,
        ),
        ((*branch_domain, "kind"), "unknown", "branch_domain.kind is invalid", False),
        (
            (*branch_domain, "boundary_status"),
            "unknown",
            "boundary_status is invalid",
            False,
        ),
        ((*branch_domain, "notes"), [1], "notes must contain strings", False),
    ]
    for path, value, message, delete in invalid_cases:
        reject(path, value, message, delete=delete)


def test_transformed_record_keeps_explicit_precursor_and_mass_basis(hugoniot_record):
    transformed = HugoniotRecord(
        **{
            name: getattr(hugoniot_record, name)
            for name in hugoniot_record.__dataclass_fields__
            if name
            not in {
                "phase",
                "cell_contents",
                "eos",
                "loading_path",
                "branch_kind",
                "initial_state",
                "volume_basis",
                "equation_kind",
            }
        },
        phase="beta",
        cell_contents="2 formula units per conventional unit cell",
        eos=hugoniot_record.eos.with_parameters(V0=20.0),
        loading_path="principal",
        branch_kind="transformed",
        initial_state=HugoniotInitialState(
            phase="alpha",
            temperature_k=300.0,
            pressure_gpa=0.0,
            density_g_cm3=8.0,
            material_identifier="x_alpha",
        ),
        volume_basis=HugoniotVolumeBasis(
            formula_units=2.0, molar_mass_g_mol=48.17712608
        ),
    )
    material = Material(
        identifier="x_beta",
        formula="X",
        phase="beta",
        cell_contents="2 formula units per conventional unit cell",
        eos_records=(transformed,),
        formula_units_per_cell=2.0,
    )
    document = material.to_eosmat()["eos_records"][0]
    assert document["loading_path"] == "principal"
    assert document["branch_kind"] == "transformed"
    assert document["initial_state"]["phase"] == "alpha"
    assert document["volume_basis"]["formula_units"] == 2.0


def test_hugoniot_record_rejects_implicit_or_inconsistent_path_metadata(
    hugoniot_record,
):
    with pytest.raises(MaterialError, match="require HugoniotRecord"):
        EOSRecord(
            **{
                name: getattr(hugoniot_record, name)
                for name in EOSRecord.__dataclass_fields__
                if name != "equation_kind"
            }
        )
    with pytest.raises(MaterialError, match="loading_path must be"):
        HugoniotRecord(
            **{
                name: getattr(hugoniot_record, name)
                for name in hugoniot_record.__dataclass_fields__
                if name not in {"path_kind", "loading_path", "equation_kind"}
            },
            loading_path="reshock",
        )
    with pytest.raises(MaterialError, match="density_g_cm3 must match"):
        HugoniotRecord(
            **{
                name: getattr(hugoniot_record, name)
                for name in hugoniot_record.__dataclass_fields__
                if name not in {"initial_state", "equation_kind"}
            },
            initial_state=HugoniotInitialState(
                phase="alpha",
                temperature_k=300.0,
                pressure_gpa=0.0,
                density_g_cm3=7.9,
                material_identifier="x_alpha",
            ),
        )


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"phase": ""}, "phase must not be empty"),
        ({"material_identifier": ""}, "material_identifier must not be empty"),
        ({"eos_record_identifier": ""}, "eos_record_identifier must not be empty"),
        ({"temperature_k": 0.0}, "temperature_k must be positive and finite"),
        ({"pressure_gpa": np.inf}, "pressure_gpa must be finite"),
        ({"density_g_cm3": 0.0}, "density_g_cm3 must be positive and finite"),
    ],
)
def test_hugoniot_initial_state_rejects_invalid_metadata(changes, message):
    values = {
        "phase": "alpha",
        "temperature_k": 300.0,
        "pressure_gpa": 0.0,
        "density_g_cm3": 8.0,
        "material_identifier": "x_alpha",
    }
    values.update(changes)
    with pytest.raises(MaterialError, match=message):
        HugoniotInitialState(**values)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"kind": "atoms"}, "kind must be 'formula_units'"),
        ({"formula_units": 0.0}, "formula_units must be positive and finite"),
        ({"molar_mass_g_mol": np.inf}, "molar_mass_g_mol must be positive and finite"),
    ],
)
def test_hugoniot_volume_basis_rejects_invalid_metadata(changes, message):
    values = {"formula_units": 1.0, "molar_mass_g_mol": 48.17712608}
    values.update(changes)
    with pytest.raises(MaterialError, match=message):
        HugoniotVolumeBasis(**values)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        (
            {"particle_velocity_km_s": (2.0, 1.0)},
            "particle velocities must be finite, non-negative, and ordered",
        ),
        ({"kind": "unknown"}, "branch-domain kind must be one of"),
        ({"boundary_status": "unknown"}, "boundary_status must be one of"),
    ],
)
def test_hugoniot_branch_domain_rejects_invalid_metadata(changes, message):
    values = {
        "particle_velocity_km_s": (0.0, 3.0),
        "kind": "phase_stability",
        "boundary_status": "reported_exactly",
    }
    values.update(changes)
    with pytest.raises(MaterialError, match=message):
        HugoniotBranchDomain(**values)


def test_hugoniot_record_rejects_all_inconsistent_metadata(hugoniot_record):
    from peritheos.eos.rt import BM3

    legacy = replace(
        hugoniot_record,
        path_kind="principal_hugoniot",
        loading_path="",
        branch_kind="",
    )
    assert legacy.loading_path == "principal"
    assert legacy.branch_kind == "untransformed"

    invalid_cases = [
        (
            {"eos": BM3(10.0, 100.0, 4.0), "equation_kind": "isothermal"},
            "eos must be a Hugoniot model",
        ),
        ({"path_kind": "reshock"}, "path_kind must be one of"),
        (
            {"path_kind": "principal_hugoniot", "loading_path": "precompressed"},
            "path_kind conflicts with loading_path",
        ),
        (
            {"path_kind": "principal_hugoniot", "branch_kind": "transformed"},
            "path_kind conflicts with branch_kind",
        ),
        ({"branch_kind": "unknown"}, "branch_kind must be one of"),
        ({"initial_state": None}, "initial_state must be HugoniotInitialState"),
        ({"volume_basis": None}, "volume_basis must be HugoniotVolumeBasis"),
        ({"branch_domain": None}, "branch_domain must be HugoniotBranchDomain"),
        (
            {"initial_state": replace(hugoniot_record.initial_state, phase="beta")},
            "untransformed Hugoniot branch",
        ),
        ({"branch_kind": "transformed"}, "transformed Hugoniot branch"),
        ({"loading_path": "precompressed"}, "requires positive P0"),
        (
            {"initial_state": replace(hugoniot_record.initial_state, pressure_gpa=1.0)},
            "pressure_gpa must match",
        ),
        (
            {
                "initial_state": replace(
                    hugoniot_record.initial_state, temperature_k=301.0
                )
            },
            "temperature_k must match",
        ),
    ]
    for changes, message in invalid_cases:
        with pytest.raises(MaterialError, match=message):
            replace(hugoniot_record, **changes)


def test_hugoniot_record_optional_domain_checks_and_vector_results(hugoniot_record):
    volumes = np.array([8.0, 9.0])
    pressures = hugoniot_record.pressure(volumes, check_domain=False)
    assert np.all(hugoniot_record.volume(pressures, check_domain=False) > 0.0)
    assert np.all(hugoniot_record.shock_velocity(volumes, check_domain=False) > 0.0)
    assert np.all(hugoniot_record.particle_velocity(volumes, check_domain=False) > 0.0)
    assert np.all(hugoniot_record.density(volumes, check_domain=False) > 0.0)
    assert np.all(
        hugoniot_record.specific_internal_energy_change(volumes, check_domain=False)
        >= 0.0
    )
    assert np.all(hugoniot_record.tangent_modulus(volumes, check_domain=False) > 0.0)
    assert hugoniot_record.state_from_particle_velocity(
        np.array([0.5, 1.0]), check_domain=False
    ).volume.shape == (2,)
    assert hugoniot_record.within_calibration_range(volumes).all()

    pressure = float(hugoniot_record.pressure(8.0))
    prediction = hugoniot_record.volume_with_uncertainty(pressure)
    assert prediction.value == pytest.approx(8.0)
    with pytest.raises(MaterialError, match="Temperature uncertainty"):
        hugoniot_record.volume_with_uncertainty(pressure, temperature_sigma=1.0)


def test_material_rejects_inconsistent_hugoniot_metadata(hugoniot_record):
    def material(record=hugoniot_record, **changes):
        values = {
            "identifier": "x_alpha",
            "formula": "X",
            "phase": "alpha",
            "cell_contents": "1 formula unit per conventional unit cell",
            "eos_records": (record,),
            "formula_units_per_cell": 1.0,
        }
        values.update(changes)
        return Material(**values)

    invalid_cases = [
        ({"formula_units_per_cell": 0.0}, "positive and finite"),
        ({"formula_units_per_cell": None}, "require formula_units_per_cell"),
        ({"formula_units_per_cell": 2.0}, "volume basis must match"),
        (
            {
                "record": replace(
                    hugoniot_record,
                    initial_state=replace(
                        hugoniot_record.initial_state,
                        material_identifier="other_alpha",
                    ),
                )
            },
            "must reference its containing material",
        ),
        ({"space_group_number": 231}, "space_group_number must be between"),
    ]
    for changes, message in invalid_cases:
        record = changes.pop("record", hugoniot_record)
        with pytest.raises(MaterialError, match=message):
            material(record, **changes)

    first = replace(hugoniot_record, identifier="first", is_default=True)
    second = replace(hugoniot_record, identifier="second", is_default=True)
    with pytest.raises(MaterialError, match="at most one default hugoniot record"):
        material(eos_records=(first, second))


def test_record_checks_branch_domain_by_default(hugoniot_record):
    with pytest.raises(MaterialError, match="outside the declared Hugoniot branch"):
        hugoniot_record.state_from_particle_velocity(3.1)
    extrapolated = hugoniot_record.state_from_particle_velocity(3.1, check_domain=False)
    assert extrapolated.particle_velocity == pytest.approx(3.1)
    outside_volume = hugoniot_record.eos.volume_from_particle_velocity(3.1)
    with pytest.raises(MaterialError, match="outside the declared Hugoniot branch"):
        hugoniot_record.pressure(outside_volume)
    assert hugoniot_record.pressure(outside_volume, check_domain=False) > 0.0


def test_mass_basis_rejects_unscaled_transformed_precursor_volume(hugoniot_record):
    with pytest.raises(MaterialError, match="same mass basis"):
        replace(
            hugoniot_record,
            phase="beta",
            cell_contents="2 formula units per conventional unit cell",
            branch_kind="transformed",
            initial_state=replace(
                hugoniot_record.initial_state,
                material_identifier="x_alpha_precursor",
            ),
            volume_basis=replace(hugoniot_record.volume_basis, formula_units=2.0),
        )


def test_precompressed_transformed_branch_is_representable(hugoniot_record):
    model = hugoniot_record.eos.with_parameters(V0=20.0, P0=5.0)
    record = replace(
        hugoniot_record,
        phase="beta",
        cell_contents="2 formula units per conventional unit cell",
        eos=model,
        loading_path="precompressed",
        branch_kind="transformed",
        initial_state=HugoniotInitialState(
            phase="alpha",
            material_identifier="x_alpha_precompressed",
            temperature_k=300.0,
            pressure_gpa=5.0,
            density_g_cm3=8.0,
        ),
        volume_basis=replace(hugoniot_record.volume_basis, formula_units=2.0),
    )
    assert record.loading_path == "precompressed"
    assert record.branch_kind == "transformed"


def test_defaults_are_scoped_by_equation_kind(hugoniot_record):
    from peritheos.eos.rt import BM3

    equilibrium = EOSRecord(
        identifier="equilibrium",
        name="Equilibrium",
        material="X",
        phase="alpha",
        cell_contents="1 formula unit per conventional unit cell",
        eos=BM3(10.0, 100.0, 4.0),
        reference_temperature=300.0,
        reference=LiteratureReference("A", 2026, "T", "", ("test",)),
        validity=ValidityRange((0.0, 100.0), (300.0, 300.0)),
        parameter_provenance={},
        is_default=True,
    )
    material = Material(
        identifier="x_alpha",
        formula="X",
        phase="alpha",
        cell_contents="1 formula unit per conventional unit cell",
        eos_records=(equilibrium, replace(hugoniot_record, is_default=True)),
        formula_units_per_cell=1.0,
    )
    assert material.default_record() is equilibrium
    assert material.default_equilibrium_record() is equilibrium
    assert material.default_hugoniot_record() is material.hugoniot_records[0]
    serialized = material.to_eosmat()
    assert [record["default_for"] for record in serialized["eos_records"]] == [
        "equilibrium",
        "hugoniot",
    ]
    validate_eosmat_document(serialized)


def test_sesame_derived_hugoniot_requires_structured_derivation(hugoniot_record):
    material = Material(
        identifier="x_alpha",
        formula="X",
        phase="alpha",
        cell_contents="1 formula unit per conventional unit cell",
        eos_records=(hugoniot_record,),
        formula_units_per_cell=1.0,
    )
    document = material.to_eosmat()
    record = document["eos_records"][0]
    record["record_kind"] = "derived"
    with pytest.raises(EosmatError, match="derived records require derivation"):
        validate_eosmat_document(document)
    record["derivation"] = {
        "source_kind": "sesame_table",
        "source_identifier": "SESAME 1234, release 2026-01",
        "method": "sampled principal Hugoniot and fit Us = c0 + s up by OLS",
        "sampling_domain": {"particle_velocity_km_s": [0.0, 3.0]},
        "access_and_licensing": "Only fitted coefficients are redistributed.",
    }
    validate_eosmat_document(document)


def test_linear_us_up_ordinary_least_squares_recovers_parameters():
    particle = np.linspace(0.2, 3.0, 20)
    shock = 4.2 + 1.45 * particle
    result = fit_linear_us_up(
        particle,
        shock,
        V0=10.0,
        rho0=8.0,
    )
    assert result.success
    assert result.parameters["c0"] == pytest.approx(4.2, rel=1.0e-8)
    assert result.parameters["s"] == pytest.approx(1.45, rel=1.0e-8)
    assert np.allclose(result.adjusted_particle_velocity, particle)
    assert result.model.pressure(8.0) > 0.0
    serialized = result.to_dict()
    assert serialized["model"]["class"] == "LinearUsUpHugoniot"
    assert serialized["parameters"]["c0"] == pytest.approx(4.2)
    assert '"adjusted_particle_velocity"' in result.to_json()


def test_bundled_single_phase_hugoniot_records():
    mgo_document = get_material_document("mgo")
    mgo = Material.from_eosmat(mgo_document)
    mgo_hugoniot = mgo.get_eos_record("mgo_b1_duffy_ahrens_1995_hugoniot_5")
    assert mgo.phase == "B1 (periclase)"
    assert mgo_hugoniot.loading_path == "principal"
    assert mgo_hugoniot.branch_kind == "untransformed"
    assert mgo_hugoniot.eos.c0 == pytest.approx(6.87)
    assert mgo_hugoniot.eos.s == pytest.approx(1.24)
    assert mgo_hugoniot.state_from_particle_velocity(0.5228).pressure == pytest.approx(
        14.0, abs=0.01
    )
    assert mgo_hugoniot.state_from_particle_velocity(3.3768).pressure == pytest.approx(
        133.0, abs=0.01
    )

    nio_document = get_material_document("nickel_oxide")
    nio = Material.from_eosmat(nio_document)
    nio_hugoniot = nio.get_eos_record("nickel_oxide_noguchi_1999_linear_hugoniot_2")
    assert nio.phase == "rhombohedral B1"
    assert nio_hugoniot.loading_path == "principal"
    assert nio_hugoniot.branch_kind == "untransformed"
    assert nio_hugoniot.eos.c0 == pytest.approx(5.3549770536)
    assert nio_hugoniot.eos.s == pytest.approx(1.2137281492)
    assert nio_hugoniot.covariance_parameters == ("c0", "s")
    assert np.allclose(
        np.asarray(nio_hugoniot.parameter_covariance),
        [
            [0.0028212816024, -0.0014193469559],
            [-0.0014193469559, 0.00087384760293],
        ],
    )


def test_linear_us_up_supports_weighting_and_particle_velocity_errors():
    particle = np.linspace(0.0, 3.0, 16)
    shock = 4.0 + 1.5 * particle
    shock[-1] += 0.1
    result = fit_linear_us_up(
        particle,
        shock,
        V0=10.0,
        rho0=8.0,
        shock_velocity_sigma=np.full_like(particle, 0.03),
        particle_velocity_sigma=np.full_like(particle, 0.01),
        absolute_sigma=True,
    )
    assert result.success
    assert result.covariance.shape == (2, 2)
    assert np.all(np.isfinite(result.adjusted_particle_velocity))


def test_hugoniot_parameter_uncertainty_propagates(hugoniot_record):
    prediction = hugoniot_record.pressure_with_uncertainty(8.0)
    assert prediction.value == pytest.approx(hugoniot_record.pressure(8.0))
    assert prediction.standard_error > 0.0

    fit = fit_linear_us_up(
        np.linspace(0.2, 3.0, 20),
        4.2 + 1.45 * np.linspace(0.2, 3.0, 20),
        V0=10.0,
        rho0=8.0,
    )
    propagated = fit.eos_uncertainty().evaluate("shock_velocity", 8.0)
    assert propagated.value == pytest.approx(fit.model.shock_velocity(8.0))
