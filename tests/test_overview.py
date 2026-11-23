"""Display-ready overviews, their text rendering and the show/list commands."""

import json
from importlib import resources
from pathlib import Path

import jsonschema
import numpy as np
import pytest

import peritheos._eosmat_audit as audit_module
from peritheos import (
    get_eos_record,
    get_material,
    get_material_document,
    list_materials,
    material_overviews,
    write_material_overviews,
)
from peritheos.__main__ import main
from peritheos._overview import (
    CHECK_EXPLANATIONS,
    COMPARED_PHRASES,
    error_convention_text,
    format_fitted_value_error,
    format_range,
    format_value_error,
    phase_text,
)
from peritheos.errors import MaterialError, MaterialLookupError
from peritheos.materials import _material_from_eosmat

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"


@pytest.fixture(scope="module")
def schema():
    text = (
        resources.files("peritheos.data")
        .joinpath("overview-v1.schema.json")
        .read_text(encoding="utf-8")
    )
    document = json.loads(text)
    jsonschema.Draft202012Validator.check_schema(document)
    return jsonschema.Draft202012Validator(document)


@pytest.fixture(scope="module")
def overviews():
    return list(material_overviews())


@pytest.fixture
def no_sidecar_reads(monkeypatch):
    calls = []

    def refuse(material):
        calls.append(material)
        raise AssertionError(f"audit sidecar of {material!r} was read")

    monkeypatch.setattr(audit_module, "load_bundled_audit", refuse)
    return calls


def stored_records():
    """``{record identifier: stored v5 record}`` from the bundled files."""
    records = {}
    for path in MATERIALS.glob("*.eosmat"):
        for record in json.loads(path.read_text(encoding="utf-8"))["eos_records"]:
            records[record["identifier"]] = record
    return records


# --- data -------------------------------------------------------------------


def test_every_overview_is_json_and_schema_valid(overviews, schema):
    assert [item["identifier"] for item in overviews] == [
        material.identifier for material in list_materials()
    ]
    for material in overviews:
        assert json.loads(json.dumps(material, allow_nan=False)) == material
        schema.validate(material)
        for record in material["records"]:
            schema.validate(record)
            assert record["overview_version"] == 1
            assert record["pressure_scale"] is None


def test_written_index_validates(tmp_path, schema):
    path = write_material_overviews(tmp_path / "overview.json")
    index = json.loads(path.read_text(encoding="utf-8"))
    schema.validate(index)
    assert index["overview_version"] == 1
    assert len(index["materials"]) == len(list_materials())


def test_check_values_equal_the_stored_block(overviews):
    stored = stored_records()
    for material in overviews:
        for record in material["records"]:
            check = stored[record["identifier"]]["check"]
            for key in (
                "reproduction",
                "compared",
                "why",
                "warning",
                "source_check",
                "source_checked_on",
                "source_note",
            ):
                assert record["check"][key] == check.get(key), (record, key)
            assert record["audit"] == check["audit"]
            display = record["check"]["display"]
            assert display["compared"] == COMPARED_PHRASES[check["compared"]]
            assert (
                display["explanation"]
                == CHECK_EXPLANATIONS[(check["reproduction"], check["compared"])]
            )


def test_parameters_equal_the_stored_values(overviews):
    stored = stored_records()
    for material in overviews:
        for record in material["records"]:
            source = stored[record["identifier"]]
            components = {
                "model": source["parameters"],
                "thermal_model": source.get("thermal", {}).get("parameters", {}),
            }
            listed = {(row["component"], row["name"]) for row in record["parameters"]}
            # Every coefficient with a value is listed. Value-less annotations
            # (an error-only entry) carry no coefficient and must be empty.
            for component, parameters in components.items():
                for name, value in parameters.items():
                    if isinstance(value, dict) and "value" not in value:
                        assert value.get("error") is None and not value.get("fixed")
                        continue
                    assert (component, name) in listed, (record["identifier"], name)
            assert len(listed) == len(record["parameters"])
            for row in record["parameters"]:
                value = components[row["component"]][row["name"]]
                if not isinstance(value, dict):
                    value = {"value": value}
                assert row["value"] == value.get("value")
                assert row["error"] == value.get("error")
                assert row["fixed"] is bool(value.get("fixed", False))


def test_default_flags_agree_with_default_record(overviews):
    by_id = {item["identifier"]: item for item in overviews}
    for material in list_materials():
        overview = by_id[material.identifier]
        default = material.default_record()
        assert overview["default_record"] == (
            None if default is None else default.identifier
        )
        if default is not None:
            assert overview["records"][0]["identifier"] == default.identifier
        selected = (
            material.default_equilibrium_record(),
            material.default_hugoniot_record(),
        )
        for record in overview["records"]:
            executable = material.get_eos_record(record["identifier"])
            chosen = any(executable is item for item in selected)
            assert record["is_default"] is (chosen and executable.is_default)
            assert record["selected_when_no_default"] is (
                chosen and not executable.is_default
            )
            # A record's own overview resolves the same selection.
            own = executable.overview()
            assert own["is_default"] is record["is_default"]
            assert own["selected_when_no_default"] is record["selected_when_no_default"]


def test_hold_and_designated_defaults():
    palladium = get_material("palladium").overview()
    assert palladium["default_record"] is None
    assert palladium["default_records"] == {"equilibrium": None, "hugoniot": None}
    assert palladium["default_hold"].startswith("Frost experimental fit")
    assert not any(record["is_default"] for record in palladium["records"])
    mgo = get_material("mgo").overview()
    assert mgo["default_records"]["equilibrium"] == {
        "identifier": "mgo_speziale_2001_bm3_2",
        "label": "Speziale et al. (2001), BM3 + variable-q Debye",
    }
    assert mgo["default_records"]["hugoniot"]["identifier"] == (
        "mgo_b1_duffy_ahrens_1995_hugoniot_5"
    )
    assert mgo["default_hold"] is None


def test_overview_never_reads_the_sidecar_by_default(no_sidecar_reads):
    material = _material_from_eosmat(
        get_material_document("mgo", format_version=None, include_audit=False),
        require_primary_validation=False,
        bundled_audit=True,
    )
    material.overview()
    material.summary()
    for record in material.eos_records:
        record.overview()
    assert no_sidecar_reads == []


def test_include_audit_adds_the_pressure_scale():
    record = get_eos_record("mgo_speziale_2001_bm3_2")
    stored = json.loads(
        (ROOT / "peritheos/data/audits/mgo.audit.json").read_text(encoding="utf-8")
    )["records"][record.identifier]["pressure_calibration"]
    scale = record.overview(include_audit=True)["pressure_scale"]
    assert scale["status"] == stored["status"]
    assert [method["kind"] for method in scale["methods"]] == [
        method["kind"] for method in stored["methods"]
    ]
    assert scale["methods"][0]["reference_calibration_record"] == "ruby_mao_1986"
    assert scale["methods"][0]["citation"].startswith("Mao et al. (1986)")


def test_record_overview_fields():
    overview = get_eos_record("mgo_speziale_2001_bm3_2").overview()
    assert overview["citation"] == (
        "Speziale et al. (2001). J. Geophys. Res. 106, 515-528. "
        "doi:10.1029/2000JB900318"
    )
    assert overview["model"] == {
        "identifier": "BM3",
        "name": "Birch-Murnaghan, third order",
    }
    assert overview["thermal_model"]["identifier"] == "VariableQMieGruneisenDebye"
    assert overview["category"] == "equilibrium"
    assert overview["report_url"].endswith(
        "/reconstruction-reports/mgo/#record-mgo_speziale_2001_bm3_2"
    )
    units = {row["name"]: row["unit"] for row in overview["parameters"]}
    assert units["V0"] == "angstrom^3/conventional_unit_cell"
    assert units["K0"] == "GPa" and units["K0_prime"] == "1"
    assert units["theta0"] == "K" and units["q1"] is None
    ranges = overview["ranges"]
    assert ranges["validity"]["pressure_gpa"] == [0.0, 203.0]
    assert ranges["experimental"]["pressure_gpa"] == [0.0, 52.0]
    assert ranges["effective"]["source"] == "validity"
    assert ranges["effective"]["pressure_gpa"] == [0.0, 203.0]
    assert ranges["effective"]["temperature_k"] == [300.0, 3663.0]
    # A user file has no report link.
    compatibility = get_material("mgo_b1").eos_records[0].overview()
    assert compatibility["report_url"] is None
    assert compatibility["parameters"][0]["fixed"] is None


@pytest.mark.parametrize(
    ("identifier", "temperature"),
    [
        ("mgo_speziale_2001_bm3_2", 300.0),
        ("mgo_shen_2026_vinet_3", None),
        ("palladium_fedotenko_2020_bm3_1", None),
    ],
)
def test_effective_pressure_range_is_what_evaluation_enforces(identifier, temperature):
    record = get_eos_record(identifier)
    lower, upper = record.overview()["ranges"]["effective"]["pressure_gpa"]
    for bound in (lower, upper):
        record.volume(bound, temperature, check_validity=True)
    for outside, side in ((upper + 0.5, "above the upper"), (lower - 0.5, "below")):
        with pytest.raises(MaterialError, match=f"Pressure .* is {side}"):
            record.volume(outside, temperature, check_validity=True)


def test_effective_temperature_range_is_what_evaluation_enforces():
    record = get_eos_record("mgo_speziale_2001_bm3_2")
    lower, upper = record.overview()["ranges"]["effective"]["temperature_k"]
    volume = record.volume(20.0, 1000.0)
    record.pressure(volume, upper, check_validity=True)
    with pytest.raises(MaterialError, match="Temperature 3700 K is above the upper"):
        record.pressure(volume, 3700.0, check_validity=True)
    with pytest.raises(MaterialError, match=f"below the lower bound {lower:g} K"):
        record.pressure(volume, 250.0, check_validity=True)


# --- number display ---------------------------------------------------------


def test_display_rounds_only_long_computed_numbers():
    assert format_value_error(1.09175271088e-06, None) == ("1.09175e-06", None)
    assert format_value_error(74.71, 0.01) == ("74.71", "0.01")
    # Stored short values print exactly as stored, even beside a wider error.
    assert format_value_error(1.65, 0.4) == ("1.65", "0.4")
    # A computed value is rounded to its error's decimal place.
    assert format_value_error(8.62870429904, 0.05) == ("8.63", "0.05")
    assert format_value_error(8.62870429904, 0.0123456789) == ("8.629", "0.012")
    assert format_range([8.62870429904, 140.545132104]) == "8.629-140.5"
    assert format_range([None, 40.0]) == "<=40"
    assert format_range([-8.834, 220.5]) == "-8.834 to 220.5"


@pytest.mark.parametrize(
    ("value", "error", "texts"),
    [
        # More than three significant digits in the error (a numerical
        # output): the authoring-guide rule (two digits for a leading 1 or 2,
        # else one), the value to match.
        (1.8112e-05, 2.1554e-06, ("1.81e-05", "2.2e-06")),
        (25.447259, 0.091348, ("25.45", "0.09")),
        (24.9234, 0.1125, ("24.92", "0.11")),
        (1.00321, 0.09991, ("1.00", "0.10")),
        (5.3549770536, 0.0531157378, ("5.35", "0.05")),
        # An error larger than the value; a coarser stored value keeps its text.
        (4e-06, 2.5512e-05, ("4e-06", "2.6e-05")),
        (3e-05, 1.2345e-06, ("3e-05", "1.2e-06")),
        (0.0, 0.12345, ("0", "0.12")),
        # Three or fewer significant digits: value and error as stored, so a
        # published pair reads as printed.
        (-0.01311, 0.00419, ("-0.01311", "0.00419")),
        (1162.0, 204.0, ("1162", "204")),
        (1.07, 1.33, ("1.07", "1.33")),
        (600.0, 474.0, ("600", "474")),
        (0.003, 0.0123, ("0.003", "0.0123")),
        (45.58, 0.04, ("45.58", "0.04")),
        (1.65, 0.4, ("1.65", "0.4")),
        (232.0, 7.0, ("232", "7")),
        # A long computed value beside a short error: the stored-number rule.
        (8.62870429904, 0.05, ("8.63", "0.05")),
        # Without an error the stored-number rule applies.
        (1.09175271088e-06, None, ("1.09175e-06", None)),
        (None, None, ("n/a", None)),
    ],
)
def test_only_numerical_errors_are_rounded_for_display(value, error, texts):
    assert format_fitted_value_error(value, error) == texts


def test_error_rounding_changes_only_the_text():
    record = get_eos_record("b4c_somayazulu_2023_berman_refit")
    stored = {row["name"]: row for row in record.overview()["parameters"]}
    assert stored["alpha0"]["value"] == 1.8112e-05
    assert stored["alpha0"]["error"] == 2.1554e-06
    text = record.summary()
    assert "alpha0 = 1.81e-05 +/- 2.2e-06" in text
    assert "dK/dT  = -0.01311 +/- 0.00419" in text
    # Fixed coefficients and coefficients without an error are unchanged.
    assert "alpha1 = 5.73e-10 fixed" in text
    shen = get_eos_record("mgo_shen_2026_vinet_3").summary()
    assert "V0  = 74.636 fixed" in shen


def test_published_errors_are_not_rounded_further():
    """A published error may be a paper's printed digits, however many."""
    record = get_eos_record("fes_vi_morard_2026_bm3_mgd")
    assert record.overview()["record_kind"] in (None, "published")
    stored = {row["name"]: row for row in record.overview()["parameters"]}
    assert (stored["K0"]["value"], stored["K0"]["error"]) == (115.5, 27.91)
    text = record.summary()
    assert "  Parameters, errors as published: " in text
    assert "115.5 +/- 27.91" in text


def test_refit_errors_are_labelled_as_the_refit_s_own():
    refit = get_eos_record(
        "ca_perovskite_wang_1996_unweighted_bm3_linear_thermal_refit"
    )
    assert refit.overview()["record_kind"] == "refit"
    text = refit.summary()
    assert "  Parameters, errors from the Peritheos refit: " in text
    assert "errors as published" not in text
    published = get_eos_record("mgo_speziale_2001_bm3_2").summary()
    assert "  Parameters, errors as published: " in published
    derived = get_eos_record("nickel_oxide_noguchi_1999_linear_hugoniot_2").summary()
    assert "  Parameters, errors as stored for this derived record: " in derived
    diagnostic = get_eos_record("argon_fcc_chen_2010_bm3_reported_constants")
    assert diagnostic.overview()["record_kind"] == "diagnostic"
    assert (
        "  Parameters, errors as stored for this diagnostic record: "
        in diagnostic.summary()
    )
    unstated = get_eos_record("akimotoite_siersch_2021_bm3_1")
    assert unstated.overview()["record_kind"] is None
    assert "  Parameters, errors as published: " in unstated.summary()
    lines = get_material("ca_perovskite").summary().splitlines()

    def row(identifier):
        return next(line for line in lines if f" {identifier} " in line)

    assert row(refit.identifier).endswith("  Peritheos refit: Wang et al. 1996")
    assert row("ca_perovskite_wang_1996_preferred_bm3").endswith("  Wang et al. 1996")


def test_covariance_line_uses_display_names():
    text = get_eos_record(
        "ca_perovskite_wang_1996_unweighted_bm3_linear_thermal_refit"
    ).summary()
    assert "    Covariance reported for: V0, K0, alpha_KT." in text
    assert "rt_eos." not in text


def test_a_coefficient_without_a_value_reads_not_stated():
    text = get_eos_record("diamond_correa_2008_double_debye_log_moment_5").summary()
    assert "n/a" not in text
    assert any(
        line.split()[:3] == ["Tr", "=", "not"] and "not stated" in line
        for line in text.splitlines()
    )


def test_phase_from_the_crystal_system_prints_in_lower_case():
    assert phase_text("CUBIC") == "cubic"
    assert phase_text("B1 (periclase)") == "B1 (periclase)"
    material = get_material("ca_perovskite")
    assert material.phase == material.symmetry == "CUBIC"  # stored data unchanged
    assert "  Formula: CaSiO3    Phase: cubic" in material.summary()
    record = material.eos_records[0]
    assert "; CaSiO3, cubic\n" in record.summary()
    assert record.overview()["material"]["phase"] == "CUBIC"


def test_one_formula_unit_is_singular():
    material = get_material("ca_perovskite")
    assert material.formula_units_per_cell == 1
    assert "Cell contents: 1 formula unit per conventional unit cell" in (
        material.summary()
    )


# --- text snapshots -----------------------------------------------------------


MGO_DEFAULT = """\
mgo_speziale_2001_bm3_2
  Speziale et al. (2001), BM3 + variable-q Debye
  Material: Periclase B1 (mgo); MgO, B1 (periclase)
  Default: the designated equilibrium record
  Model: Birch-Murnaghan, third order (BM3) + VariableQMieGruneisenDebye; \
reference temperature 300 K
  Check: reproduced (refit vs published fit). A refit of the source data gave \
the published coefficients.
  Parameters, errors as published: confidence level not stated; used as 1 sigma
    Birch-Murnaghan, third order (BM3):
      V0     = 74.71 +/- 0.01                     angstrom^3/conventional_unit_cell
      K0     = 160.2 fixed                        GPa
      K0'    = 3.99 +/- 0.01                      dimensionless
    VariableQMieGruneisenDebye:
      Tr     = 300 fixed                          K
      theta0 = 773 fixed                          K
      gamma0 = 1.524 fixed (published +/- 0.025)  dimensionless
      q0     = 1.65 fixed (published +/- 0.4)
      q1     = 11.8 +/- 0.2
      n      = 2 fixed                            dimensionless
    No documented unit for: q0, q1 (see the model's equation page).
  Ranges (P in GPa, T in K):
    Validity range (enforced by check_validity=True): P 0-203, T 300-3663
    Experimental data range: P 0-52
  Reference: Speziale et al. (2001). J. Geophys. Res. 106, 515-528. \
doi:10.1029/2000JB900318
    Quasi-hydrostatic compression of magnesium oxide to 52 GPa: Implications \
for the pressure-volume-temperature equation of state
  Source check: verified (2026-10-09). Reference isotherm from the constrained \
BM3 fit of Section 4.2; theta0 and gamma0 from Section 5.1, q0 and q1 from \
Section 5.4 (Equations 10-11).
  Report: https://peritheos.readthedocs.io/en/latest/reconstruction-reports/mgo/\
#record-mgo_speziale_2001_bm3_2"""

FROST_CLOSE = """\
palladium_frost_2023_vinet_1
  Frost et al. (2023), Vinet
  Material: Palladium (palladium); Pd, cubic
  Model: Vinet; reference temperature 300 K
  Check: close (refit vs published fit). K0': refit 5.6037 vs published 5.473 \
(+2.4%), outside the combined error.
  Warning: Not recommended for quantitative use. Default selection on hold \
pending review. Never chosen automatically; select it explicitly. The authors' \
fitting objective and covariance convention are unstated, so the published \
fit is not fully reproduced.
  Parameters, errors as published: confidence level not stated; used as 1 sigma
    Vinet:
      V0  = 58.678 +/- 0.073  angstrom^3/conventional_unit_cell
      K0  = 189.3 +/- 3       GPa
      K0' = 5.473 +/- 0.063   dimensionless
  Ranges (P in GPa, T in K):
    Experimental data range (enforced by check_validity=True): P 0-182"""

SHEN_WARNING = """\
mgo_shen_2026_vinet_3
  Shen and Smith (2026), Cu-referenced Vinet
  Material: Periclase B1 (mgo); MgO, B1 (periclase)
  Model: Vinet; reference temperature 300 K
  Check: disagrees (refit vs published fit). K0': refit 3.954 vs published 4.08 \
(-3.1%), outside the combined error.
  Warning: Not recommended for quantitative use. Never chosen automatically; \
select it explicitly. The Table II fit is not reproduced from the published \
data; the published curve leaves residuals up to about +2.8 GPa.
  Parameters, errors as published: confidence level not stated; used as 1 sigma
    Vinet:
      V0  = 74.636 fixed   angstrom^3/conventional_unit_cell
      K0  = 161.9 +/- 0.7  GPa
      K0' = 4.08 +/- 0.03  dimensionless"""

PALLADIUM = """\
Palladium (palladium)
  Formula: Pd    Phase: cubic
  Cell contents: 4 formula units per conventional unit cell
  Volume unit: angstrom^3/conventional_unit_cell
  Default: none; selection on hold. Frost experimental fit reproduction remains \
unresolved. No replacement equilibrium default has been audited and selected. \
The separate Baty 0 K DFT record must not be selected automatically as a \
replacement room-temperature EOS.
  7 EOS records (* default, ! warning; P in GPa and T in K are the range \
enforced by check_validity=True, - none stated):
     identifier                      model  P (GPa)  T (K)  check       reference
     palladium_baty_2024_bm3_dft_2   BM3    -        -      untested    Baty et al. 2024
     palladium_guigue_2020_vinet_1   Vinet  5-65     -      untested    Guigue et al. 2020
     palladium_fedotenko_2020_bm3_1  BM3    0-31.2   -      reproduced  Fedotenko et al. 2020
  !  palladium_fedotenko_2020_bm2_2  BM2    0-31.2   -      reproduced  Fedotenko et al. 2020
  !  palladium_frost_2023_vinet_1    Vinet  0-182    -      close       Frost et al. 2023
  !  palladium_frost_2023_bm3_2      BM3    0-182    -      close       Frost et al. 2023
  !  palladium_baty_2024_bm3_1       BM3    0-80     -      disagrees   Baty et al. 2024
  Warnings:
    ! Not recommended for quantitative use. Never chosen automatically; select \
it explicitly.
      palladium_fedotenko_2020_bm2_2, palladium_baty_2024_bm3_1
    ! Not recommended for quantitative use. Default selection on hold pending \
review. Never chosen automatically; select it explicitly.
      palladium_frost_2023_vinet_1, palladium_frost_2023_bm3_2"""

NO_DESIGNATED_DEFAULT = """\
Al2O3 post-perovskite model crystal (al2o3_post_perovskite)
  Formula: Al2O3    Phase: hypothetical CaIrO3-type Cmcm post-perovskite model crystal
  Cell contents: Four Al2O3 formula units on the conventional Cmcm EOS basis; \
the source end-member calculation used an 80-atom supercell.
  Volume unit: angstrom^3/conventional_unit_cell
  Aliases: CCM 100 mol.% Al2O3 post-perovskite
  Family: post_perovskite
  Default: none designated; the first record is used \
(al2o3_post_perovskite_akber_knutson_2005_gga_bm3)
  1 EOS record (* default, ! warning; P in GPa and T in K are the range \
enforced by check_validity=True, - none stated):
    identifier                                        model  P (GPa)          \
T (K)  check     reference
    al2o3_post_perovskite_akber_knutson_2005_gga_bm3  BM3    -8.834 to 220.5  \
0      untested  Akber-Knutson et al. 2005"""


def test_record_summary_snapshots():
    assert get_eos_record("mgo_speziale_2001_bm3_2").summary() == MGO_DEFAULT
    assert (
        get_eos_record("palladium_frost_2023_vinet_1")
        .summary()
        .startswith(FROST_CLOSE + "\n")
    )
    assert (
        get_eos_record("mgo_shen_2026_vinet_3")
        .summary()
        .startswith(SHEN_WARNING + "\n")
    )


def test_material_summary_snapshots():
    assert get_material("palladium").summary() == PALLADIUM
    assert get_material("al2o3_post_perovskite").summary() == NO_DESIGNATED_DEFAULT


def test_summary_uses_no_attribute_names():
    for identifier in ("mgo_speziale_2001_bm3_2", "gold_dewaele_2004_vinet_5"):
        text = get_eos_record(identifier).summary()
        for name in ("parameter_covariance", "experimental_pressure_range_gpa"):
            assert name not in text
    text = get_eos_record("gold_dewaele_2004_vinet_5").summary()
    if get_eos_record("gold_dewaele_2004_vinet_5").parameter_covariance is not None:
        assert "Covariance reported for:" in text


# --- lookup messages ----------------------------------------------------------


def test_deferred_record_lookup_says_it_is_deferred():
    for lookup in (
        lambda: get_eos_record("feo_campbell_2009_bm3_mgd"),
        lambda: get_material("feo").get_eos_record("feo_campbell_2009_bm3_mgd"),
    ):
        with pytest.raises(MaterialLookupError, match="is deferred") as error:
            lookup()
        assert "get_material_document('feo')" in str(error.value)
        assert "Unknown" not in str(error.value)


def test_branch_domain_error_names_the_bound():
    record = get_eos_record("mgo_b1_duffy_ahrens_1995_hugoniot_5")
    upper = record.overview()["ranges"]["effective"]["particle_velocity_km_s"][1]
    with pytest.raises(MaterialError, match="above the upper bound") as error:
        record.state_from_particle_velocity(upper + 0.5)
    assert f"{upper:.6g} km/s" in str(error.value)


# --- command line ---------------------------------------------------------------


def test_show_prints_the_summary(capsys):
    assert main(["show", "periclase"]) == 0
    assert capsys.readouterr().out == get_material("mgo").summary() + "\n"
    assert main(["show", "MgO", "mgo_speziale_2001_bm3_2"]) == 0
    assert capsys.readouterr().out == MGO_DEFAULT + "\n"
    assert main(["show", "mgo_speziale_2001_bm3_2"]) == 0
    assert capsys.readouterr().out == MGO_DEFAULT + "\n"


def test_show_json_prints_the_overview(capsys):
    assert main(["show", "mgo", "mgo_speziale_2001_bm3_2", "--json"]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed == get_eos_record("mgo_speziale_2001_bm3_2").overview()
    assert main(["show", "palladium", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["default_hold"]


def test_show_reports_lookup_errors(capsys):
    assert main(["show", "periclas"]) == 2
    error = capsys.readouterr().err
    assert "Unknown material 'periclas'" in error and "mgo (Periclase B1)" in error
    assert main(["show", "mgo", "speziale"]) == 2
    assert "Did you mean: mgo_speziale_2001_bm3_2?" in capsys.readouterr().err


def test_list_prints_one_line_per_material(capsys):
    assert main(["list"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split()[:3] == ["identifier", "formula", "name"]
    assert len(lines) == len(list_materials()) + 1
    assert main(["list", "periclase"]) == 0
    rows = capsys.readouterr().out.splitlines()[1:]
    mgo = next(row for row in rows if row.startswith("mgo "))
    assert "37 records" in mgo and "Speziale et al. (2001)" in mgo
    assert np.all([row.split()[0] for row in rows])


@pytest.mark.parametrize(
    ("confidence", "text"),
    [
        (None, "confidence level not stated; used as 1 sigma"),
        # 0.68 is not 1 sigma: propagation divides by ndtri(0.84) = 0.9945.
        (0.68, "68% half-widths (0.994 sigma)"),
        (0.6827, "1 sigma"),
        (0.682689492137, "1 sigma"),
        (0.6826894921370859, "1 sigma"),
        (0.95, "95% half-widths (1.96 sigma)"),
        (0.9545, "2 sigma"),
        (0.9544997361036416, "2 sigma"),
        (0.9, "90% half-widths (1.64 sigma)"),
        (0.99, "99% half-widths (2.58 sigma)"),
        (0.6826, "68.26% half-widths (0.9998 sigma)"),
        # Never 100% or 0%.
        (0.99999, "99.999% half-widths (4.42 sigma)"),
        (0.9999999999999999, "99.99999999999999% half-widths (8.29 sigma)"),
        (1e-7, "1e-05% half-widths (1.25e-07 sigma)"),
    ],
)
def test_summary_states_how_printed_errors_relate_to_sigma(confidence, text):
    assert error_convention_text(confidence) == text


@pytest.mark.parametrize("confidence", [0.68, 0.6827, 0.95, 0.9545, 0.9])
def test_the_printed_sigma_factor_is_the_propagation_factor(confidence):
    from scipy.special import ndtri

    text = error_convention_text(confidence)
    factor = 1.0 if text == "1 sigma" else 2.0 if text == "2 sigma" else None
    if factor is None:
        factor = float(text.split("(")[1].split()[0])
    assert factor == pytest.approx(ndtri(0.5 * (1 + confidence)), rel=3e-3)
