"""EOSMAT v5 ``check`` fields, the ``default_records`` map and record layout."""

import copy
import json
import re
from pathlib import Path

import pytest

from peritheos import Material, convert_eosmat, get_material, get_material_document
from peritheos._eosmat_audit import (
    SOURCE_CHECK_VALUES,
    WHY_MAX_LENGTH,
    _is_fixed,
    _parameter_name,
    check_block,
    compared_scope,
    reproduction_why,
)
from peritheos._eosmat_format import format_audit, format_eosmat
from peritheos._overview import format_number
from peritheos.eosmat import list_material_documents, validate_eosmat_document
from peritheos.errors import EosmatError
from tests.catalog_state import added_records

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
AUDITS = ROOT / "peritheos/data/audits"

#: Close/disagrees records whose stored comparison table names no failed
#: criterion (and, for close, no coefficient without a combined-error test):
#: no ``why``.
RECORDS_WITHOUT_WHY = {
    "bridgmanite_dorogokupets_2015_298k_rydberg_stacey",
    "cementite_scott_2001_bm3_1",
    "fe_fcc_campbell_2009_bm3_mgd",
    "mg075fe025o_solomatova_2016_hs_bm3_reference_1",
    "mg075fe025o_solomatova_2016_ls_bm3_reference_2",
    "mgsio3_post_perovskite_dorogokupets_2015_298k_rydberg_stacey",
}


def record(status=None, scope="parameter_agreement", rows=None, **extra):
    evidence = {
        "scientific_validation": {
            "status": "primary_source_validated",
            "audit_date": "2026-01-02",
            "primary_source_check": {"finding": "Table 2 values match."},
        }
    }
    if status is not None:
        assessment = {"status": status, "scope": scope, "basis": "Refit of Table 1."}
        if rows is not None:
            assessment["parameters"] = rows
        evidence["reproduction"] = {"assessment": assessment}
    return {
        "identifier": "example_bm3",
        "model": "BM3",
        "temperature_ref": 300,
        "parameters": {"V0": 74.7, "K0": 160, "K0_prime": 4},
        "reference": "Fictional example; no experimental validation.",
        "evidence": evidence,
        **extra,
    }


def document(*records, **extra):
    return {
        "format_version": 5,
        "identifier": "example",
        "formula": "MgO",
        "eos_records": list(records),
        **extra,
    }


def row(parameter, published, refit, **flags):
    return {
        "parameter": parameter,
        "published": published,
        "published_error": None,
        "refit": refit,
        "refit_error": None,
        "difference": refit - published,
        "relative_difference": abs(refit - published) / published,
        **flags,
    }


@pytest.mark.parametrize(
    ("status", "scope", "compared"),
    [
        ("parity", "parameter_agreement", "published_fit"),
        ("similar", "parameter_agreement", "published_fit"),
        ("not_reproduced", "parameter_agreement", "published_fit"),
        ("parity", "independent_refit", "peritheos_refit"),
        ("parity", "source_equation", "source_equation"),
        ("source_reconstruction", "source_equation", "source_equation"),
        ("not_refittable", "parameter_agreement", "none"),
        ("not_refittable", "source_equation", "none"),
        ("not_assessed", "parameter_agreement", "none"),
        (None, None, "none"),
    ],
)
def test_compared_is_a_fixed_mapping_of_status_and_scope(status, scope, compared):
    source = record(status, scope)
    assert compared_scope(source) == compared
    loaded = Material.from_eosmat(document(source)).eos_records[0]
    assert loaded.compared == compared
    assert loaded.source_check == "verified"
    assert loaded.source_note == "Table 2 values match."


def test_bundled_compared_counts_and_source_check_values():
    compared = {}
    source = {}
    added = added_records()
    for material in list_material_documents():
        for stored in json.loads((MATERIALS / f"{material}.eosmat").read_bytes())[
            "eos_records"
        ]:
            check = stored["check"]
            assert not {"finding", "source"} & set(check)
            if stored["identifier"] in added:
                continue  # counted below: curated before the revision ledger
            compared[check["compared"]] = compared.get(check["compared"], 0) + 1
            source[check["source_check"]] = source.get(check["source_check"], 0) + 1
    assert compared == {
        "published_fit": 474,
        "peritheos_refit": 16,
        "source_equation": 19,
        "none": 178,
    }
    assert source == {"verified": 666, "unresolved": 18, "deferred": 3}


def test_source_check_is_a_fixed_renaming_kept_out_of_the_sidecar():
    assert SOURCE_CHECK_VALUES == {
        "primary_source_validated": "verified",
        "pending_primary_source_check": "pending",
        "not_reproduced": "unresolved",
        "deferred": "deferred",
    }
    for path in AUDITS.glob("*.audit.json"):
        for entry in json.loads(path.read_bytes())["records"].values():
            status = entry["evidence"]["scientific_validation"]["status"]
            assert status in SOURCE_CHECK_VALUES


def with_errors(source, published, refit):
    return {**source, "published_error": published, "refit_error": refit}


def thermal(model, parameters):
    return {"thermal": {"model": model, "parameters": parameters}}


#: One case per ``why`` template; crates/peritheos/tests/eosmat_audit.rs
#: checks the same cases against the Rust reader.
WHY_CASES = {
    # Every failing coefficient with both values, unit and percentage.
    "values": (
        record(
            "similar",
            rows=[
                row("V0", 162.4, 162.41, similar=True, within_combined_2sigma=True),
                row("K0", 256, 263.7307, similar=True, within_combined_2sigma=False),
            ],
        ),
        "Refit vs published: K0 263.7 vs 256 GPa (+3.0%) differs by more than "
        "twice the combined error.",
    ),
    # The second spelling of the table; three failed criteria in one row.
    "three_criteria": (
        record(
            "not_reproduced",
            rows=[
                {
                    "parameter": "rt_eos.K0_prime",
                    "published": 4.91,
                    "published_sigma": 0.12,
                    "refit": 4.183456,
                    "refit_sigma": None,
                    "similar": False,
                    "within_combined_2sigma": False,
                    "within_reported_error": False,
                }
            ],
        ),
        "Refit vs published: K0' 4.183 vs 4.91 (-14.8%) differs by more than "
        "twice the combined error, the agreement limit and the published error.",
    ),
    # Volume in Å³; published values print as stored, the refit with at most
    # one more decimal.
    "volume": (
        record(
            "similar",
            rows=[row("V0", 74.7, 75.0, similar=True, within_combined_2sigma=False)],
        ),
        "Refit vs published: V0 75 vs 74.7 Å³ (+0.4%) differs by more than twice "
        "the combined error.",
    ),
    # A thermal coefficient takes the thermal model's unit; exponent notation
    # follows the stored spelling's significant digits.
    "exponent": (
        record(
            "similar",
            rows=[
                row(
                    "alpha1",
                    8.7e-9,
                    5.6123e-9,
                    similar=False,
                    within_combined_2sigma=True,
                )
            ],
            **thermal("AlphaKT", {"Tr": 300.0, "alpha0": 3e-5, "alpha1": 8.7e-9}),
        ),
        "Refit vs published: alpha1 5.61e-09 vs 8.7e-09 K^-2 (-35.5%) differs by "
        "more than the agreement limit.",
    ),
    # A difference of 100% or more gives the two values, no percentage.
    "large": (
        record(
            "not_reproduced",
            rows=[
                row("theta0", 600.0, 766.2634, similar=False),
                row("gamma0", 2.609, 8.0, similar=False),
            ],
            **thermal("MieGruneisenDebye", {"theta0": 600.0, "gamma0": 2.609}),
        ),
        "Refit vs published: theta0 766.3 vs 600 K (+27.7%) and gamma0 8 vs 2.609 "
        "differ by more than the agreement limit.",
    ),
    # Too long with values: percentages. One clause per criterion, the one
    # with most coefficients first; a coefficient already named is not
    # repeated.
    "mixed": (
        record(
            "not_reproduced",
            rows=[
                row("V0", 74.7, 75.5, similar=True, within_combined_2sigma=False),
                row("K0", 160, 200, similar=False, within_combined_2sigma=False),
                row("K0_prime", 4, 0.0, similar=False, within_combined_2sigma=True),
            ],
        ),
        "Refit vs published: V0 (+1.1%) and K0 (+25.0%) differ by more than twice "
        "the combined error; K0' (0 vs 4) differs by more than the agreement "
        "limit, as does K0.",
    ),
    # A coefficient failing two criteria is named once; "also" for the second.
    "also": (
        record(
            "not_reproduced",
            rows=[
                row("V0", 74.7, 76.87, similar=True, within_combined_2sigma=False),
                row("K0", 160, 135.68, similar=False, within_combined_2sigma=False),
            ],
        ),
        "Refit vs published: V0 (+2.9%) and K0 (-15.2%) differ by more than twice "
        "the combined error; K0 also by more than the agreement limit.",
    ),
    # Four or more failing coefficients are all named.
    "four": (
        record(
            "not_reproduced",
            rows=[
                row("V0", 74.7, 75.0, similar=False),
                row("K0", 160, 200, similar=False),
                row("K0_prime", 4, 2, similar=False),
                row("gamma0", 1.5, 1.6, similar=False),
            ],
        ),
        "Refit vs published: V0 (+0.4%), K0 (+25.0%), K0' (-50.0%) and gamma0 "
        "(+6.7%) differ by more than the agreement limit.",
    ),
    # More than four failing coefficients are counted; the largest keep their
    # percentage.
    "largest": (
        record(
            "not_reproduced",
            rows=[
                row(f"theta_{m}0", 1888.0, 1888.0 * (1.1 + 0.1 * i), similar=False)
                for i, m in enumerate("abcdefg")
            ],
        ),
        "Refit vs published: 7 coefficients differ by more than the agreement "
        "limit; largest: theta_g0 (+70.0%), theta_f0 (+60.0%) and theta_e0 "
        "(+50.0%).",
    ),
    # A count and the largest, which keep their values while the line fits.
    "count": (
        record(
            "not_reproduced",
            rows=[
                row(f"coefficient_{i}", 1.25, 1.5 + i, similar=False) for i in range(12)
            ],
        ),
        "Refit vs published: 12 coefficients differ by more than the agreement "
        "limit; largest: coefficient_11 (12.5 vs 1.25), coefficient_10 and "
        "coefficient_9.",
    ),
    # Nothing failed: the coefficients without a combined-error test, and
    # which error is missing.
    "paper_gives_no_error": (
        record(
            "similar",
            rows=[
                with_errors(
                    row("K0", 160, 161, similar=True, within_combined_2sigma=None),
                    None,
                    1.2,
                ),
                with_errors(
                    row("K0_prime", 4, 4.1, similar=True, within_combined_2sigma=None),
                    None,
                    0.1,
                ),
            ],
        ),
        "All coefficients are within the agreement limit; the paper gives no "
        "error for K0 and K0', so they could not be tested against twice the "
        "combined error.",
    ),
    "refit_gives_no_error": (
        record(
            "similar",
            rows=[
                with_errors(
                    row("V0", 74.7, 74.71, similar=True, within_combined_2sigma=None),
                    0.1,
                    None,
                ),
                with_errors(
                    row("K0", 160, 161, similar=True, within_combined_2sigma=True),
                    2,
                    1,
                ),
            ],
        ),
        "All coefficients are within the agreement limit; the refit gives no "
        "error for V0, so it could not be tested against twice the combined error.",
    ),
    "neither_gives_an_error": (
        record(
            "similar",
            rows=[row("K0", 160, 161, similar=True, within_combined_2sigma=None)],
        ),
        "All coefficients are within the agreement limit; no published or refit "
        "error for K0, so it could not be tested against twice the combined error.",
    ),
    "several_reasons": (
        record(
            "similar",
            rows=[
                with_errors(
                    row("V0", 74.7, 74.71, similar=True, within_combined_2sigma=None),
                    0.1,
                    None,
                ),
                row("K0", 160, 161, similar=True, within_combined_2sigma=None),
                row("K0_prime", 4, 4.1, similar=True, within_combined_2sigma=None),
            ],
        ),
        "V0 has no refit error and K0 and K0' no published or refit error, so none "
        "could be tested against twice the combined error.",
    ),
    "two_reasons": (
        record(
            "similar",
            rows=[
                row("K0", 160, 161, similar=True, within_combined_2sigma=None),
                with_errors(
                    row("K0_prime", 4, 4.1, similar=True, within_combined_2sigma=None),
                    0.1,
                    None,
                ),
                row("V0", 74.7, 74.71, similar=True, within_combined_2sigma=True),
            ],
        ),
        # With "All coefficients are within the agreement limit; " too long.
        "K0 has no published or refit error and K0' no refit error, so neither "
        "could be tested against twice the combined error.",
    ),
    # A coefficient the record holds fixed is not a candidate for the test.
    "fixed_not_untested": (
        record(
            "similar",
            rows=[
                with_errors(
                    row("V0", 74.7, 74.71, similar=True, within_combined_2sigma=None),
                    0.1,
                    None,
                ),
                with_errors(
                    row("K0", 160, 161, similar=True, within_combined_2sigma=None),
                    2,
                    None,
                ),
                row("K0_prime", 4, 4.0, similar=True, within_combined_2sigma=None),
            ],
            parameters={
                "V0": {"value": 74.7, "error": 0.1},
                "K0": {"value": 160, "error": 2},
                "K0_prime": {"value": 4, "fixed": True},
            },
        ),
        "All coefficients are within the agreement limit; the refit gives no "
        "error for V0 and K0, so they could not be tested against twice the "
        "combined error.",
    ),
    # Names and numbers print as the record summary prints them: dK/dT, and
    # e-notation below 1e-4.
    "summary_spelling": (
        record(
            "similar",
            rows=[
                row(
                    "alpha0",
                    0.0001158,
                    9.2237e-5,
                    similar=True,
                    within_combined_2sigma=False,
                ),
                row(
                    "dK_dT", -0.008, -0.0123, similar=True, within_combined_2sigma=False
                ),
            ],
            **thermal("AlphaKT", {"Tr": 300.0, "alpha0": 0.0001158, "dK_dT": -0.008}),
        ),
        "Refit vs published: alpha0 9.224e-05 vs 0.0001158 K^-1 (-20.3%) and dK/dT "
        "-0.0123 vs -0.008 GPa/K (-53.8%) differ by more than twice the combined "
        "error.",
    ),
    # No stored errors at all, and not every row passed the agreement limit.
    "unknown": (
        record(
            "similar",
            rows=[
                {
                    "parameter": "K0",
                    "published": 160,
                    "refit": 161,
                    "within_combined_2sigma": None,
                },
                row("V0", 74.7, 74.71, similar=True, within_combined_2sigma=True),
            ],
        ),
        "No error comparison is recorded for K0, so it could not be tested against "
        "twice the combined error.",
    ),
    # A published value with more than six significant digits (a computed
    # conversion) is not repeated.
    "long_published": (
        record(
            "similar",
            rows=[
                row(
                    "V0",
                    41.480265898,
                    41.5535840287,
                    similar=True,
                    within_combined_2sigma=False,
                )
            ],
        ),
        "Refit vs published: V0 (+0.2%) differs by more than twice the combined error.",
    ),
}


@pytest.mark.parametrize("case", sorted(WHY_CASES))
def test_why_says_what_differs_in_plain_words(case):
    source, expected = WHY_CASES[case]
    assert reproduction_why(source) == expected
    assert len(expected) <= WHY_MAX_LENGTH


def test_why_states_nothing_the_table_does_not_support():
    close = WHY_CASES["values"][0]
    # Nothing failed and nothing untested: no statement is invented.
    passing = record(
        "not_reproduced",
        rows=[row("K0", 160, 161, similar=True, within_combined_2sigma=True)],
    )
    assert reproduction_why(passing) is None
    assert (
        reproduction_why(
            record(
                "parity",
                rows=close["evidence"]["reproduction"]["assessment"]["parameters"],
            )
        )
        is None
    )
    assert reproduction_why(record("similar")) is None
    # An untested coefficient explains a close record only.
    untested = WHY_CASES["neither_gives_an_error"][0]
    rows = untested["evidence"]["reproduction"]["assessment"]["parameters"]
    assert reproduction_why(record("not_reproduced", rows=rows)) is None


def printed(value):
    """A stored number as the why line prints it: as the record summary does."""
    return format_number(value)


def test_bundled_why_lines_are_short_and_rounded():
    with_why = 0
    without = set()
    for material in list_material_documents():
        full = get_material_document(material, format_version=5)
        for complete in full["eos_records"]:
            check = complete["check"]
            if check["reproduction"] not in ("close", "disagrees"):
                assert "why" not in check
                continue
            if "why" not in check:
                without.add(complete["identifier"])
                continue
            with_why += complete["identifier"] not in added_records()
            why = check["why"]
            assert len(why) <= WHY_MAX_LENGTH
            # Every coefficient is named or counted; no percentage of 100% or
            # more is printed (the two values are given instead).
            assert not re.search(r"and \d+ more", why)
            assert all(abs(float(p)) < 100 for p in re.findall(r"([-+][\d.]+)%", why))
            rows = complete["evidence"]["reproduction"]["assessment"]["parameters"]
            number = r"-?[\d.]+(?:e[-+]\d+)?"
            for match in re.finditer(rf"({number}) vs ({number})\b", why):
                refit, published = match.groups()
                # Decimals of the mantissa ("6e-10" has none, "1.0e-09" one).
                decimals = len(published.partition("e")[0].partition(".")[2])
                assert len(refit.partition("e")[0].partition(".")[2]) <= decimals + 1
                assert any(printed(r["published"]) == published for r in rows)
    assert without == RECORDS_WITHOUT_WHY
    assert with_why == 127


def test_a_fixed_coefficient_is_never_named_untested():
    # Native flags and v3 export lists, for the reference and thermal models.
    native = {
        "parameters": {"V0": 74.7, "K0_prime": {"value": 4, "fixed": True}},
        "thermal": {"parameters": {"gamma0": {"value": 1.5, "fixed": True}}},
    }
    legacy = {
        "eos": {"parameters": {"V0": 74.7, "K0_prime": 4}},
        "fixed_parameters": ["K0_prime"],
        "thermal": {"parameters": {"gamma0": 1.5}, "fixed_parameters": ["gamma0"]},
    }
    for source in (native, legacy):
        assert _is_fixed(source, "K0_prime")
        assert _is_fixed(source, "rt_eos.K0_prime")
        assert _is_fixed(source, "gamma0")
        assert not _is_fixed(source, "V0")
    for material in list_material_documents():
        full = get_material_document(material, format_version=5)
        for complete in full["eos_records"]:
            why = complete["check"].get("why") or ""
            if "could not be tested" not in why and "could be tested" not in why:
                continue
            rows = complete["evidence"]["reproduction"]["assessment"]["parameters"]
            for entry in rows:
                if _is_fixed(complete, entry["parameter"]):
                    name = _parameter_name(entry["parameter"])
                    assert not re.search(rf"(^|[ ;]){re.escape(name)}\b", why), (
                        complete["identifier"]
                    )


def test_check_fields_are_derived_without_a_check_block():
    source = record(
        "similar",
        rows=[row("K0", 256, 263.73, similar=True, within_combined_2sigma=False)],
    )
    plain = Material.from_eosmat(document(source)).eos_records[0]
    with_check = copy.deepcopy(source)
    with_check["check"] = check_block(source)
    stored = Material.from_eosmat(document(with_check)).eos_records[0]
    for name in (
        "reproduction",
        "compared",
        "why",
        "warning",
        "source_check",
        "source_note",
    ):
        assert getattr(plain, name) == getattr(stored, name), name
    assert plain.why.startswith("Refit vs published: K0 263.7 vs 256 GPa")
    assert (
        "  Check: close (refit vs published fit). Refit vs published: K0 263.7 "
        in plain.summary()
    )
    assert (
        "  Source check: verified (2026-01-02). Table 2 values match."
        in plain.summary()
    )


def test_check_must_agree_with_evidence_and_old_names_are_rejected():
    contradicted = record(
        "parity", check={"reproduction": "reproduced", "compared": "none"}
    )
    with pytest.raises(EosmatError, match="cannot accompany"):
        Material.from_eosmat(document(contradicted))
    wrong_scope = record(
        "parity", "independent_refit", check={"compared": "published_fit"}
    )
    with pytest.raises(EosmatError, match="check.compared 'published_fit' contradicts"):
        Material.from_eosmat(document(wrong_scope))
    wrong_source = record("parity", check={"source_check": "deferred"})
    with pytest.raises(EosmatError, match="check.source_check 'deferred' contradicts"):
        Material.from_eosmat(document(wrong_source))
    for old, new in (("finding", "source_note"), ("source", "source_check")):
        stale = record(None, check={old: "primary_source_validated"})
        with pytest.raises(
            EosmatError, match=f"check.{old} is no longer read.*check.{new}"
        ):
            validate_eosmat_document(document(stale))
    with pytest.raises(EosmatError, match="only given for reproduction close"):
        validate_eosmat_document(
            document(record(None, check={"reproduction": "reproduced", "why": "x"}))
        )
    with pytest.raises(EosmatError, match="at most 160"):
        validate_eosmat_document(document(record(None, check={"why": "x" * 161})))


def test_bundled_files_are_in_canonical_layout():
    for material in list_material_documents():
        path = MATERIALS / f"{material}.eosmat"
        data = path.read_bytes()
        assert format_eosmat(json.loads(data)).encode("utf-8") == data, material
        audit = AUDITS / f"{material}.audit.json"
        assert (
            format_audit(json.loads(audit.read_bytes())).encode() == audit.read_bytes()
        )


def test_record_reading_order():
    stored = json.loads((MATERIALS / "mgo.eosmat").read_bytes())
    keys = list(stored)
    assert keys.index("default_records") + 1 == keys.index("eos_records")
    order = (
        "identifier",
        "label",
        "check",
        "model",
        "temperature_ref",
        "parameters",
        "thermal",
        "parameter_error_confidence",
        "parameter_covariance",
        "validity",
        "experimental_pressure_range_gpa",
        "determination_method",
        "record_kind",
        "notes",
        "fit_datasets",
        "reference",
    )
    for material in list_material_documents():
        for entry in json.loads((MATERIALS / f"{material}.eosmat").read_bytes())[
            "eos_records"
        ]:
            present = [key for key in entry if key in order]
            assert present == [key for key in order if key in entry]
            assert list(entry["check"])[:2] == ["reproduction", "compared"]


def test_default_records_map_selects_exactly_like_record_flags():
    first = record("parity", model="BM3")
    second = {**record("parity"), "identifier": "example_vinet", "model": "Vinet"}
    mapped = document(first, second, default_records={"equilibrium": "example_vinet"})
    flagged = document(first, {**second, "default": True})
    by_map = Material.from_eosmat(mapped)
    by_flag = Material.from_eosmat(flagged)
    assert by_map.default_record().identifier == "example_vinet"
    assert [r.is_default for r in by_map.eos_records] == [
        r.is_default for r in by_flag.eos_records
    ]
    for version in (3, 4):
        exported = convert_eosmat(mapped, version=version)
        assert "default_records" not in exported
        assert exported["eos_records"][1]["default_for"] == "equilibrium"
        assert "default_for" not in exported["eos_records"][0]
    assert convert_eosmat(mapped, version=5)["default_records"] == {
        "equilibrium": "example_vinet"
    }
    # Both forms in one file are accepted only when they agree.
    validate_eosmat_document(
        document(
            first,
            {**second, "default_for": "equilibrium"},
            default_records={"equilibrium": "example_vinet"},
        )
    )
    for bad, message in (
        (
            document(
                {**first, "default": True},
                second,
                default_records={"equilibrium": "example_vinet"},
            ),
            "flagged as the default equilibrium record, but default_records names 'example_vinet'",
        ),
        (
            document(
                first,
                {**second, "default": False},
                default_records={"equilibrium": "example_vinet"},
            ),
            "flagged default: false",
        ),
        (
            document(first, second, default_records={"equilibrium": "missing"}),
            "unknown record",
        ),
        (
            document(first, second, default_records={"hugoniot": "example_vinet"}),
            "which is not a hugoniot record",
        ),
        (
            document(first, second, default_records={"best": "example_vinet"}),
            "keys must be",
        ),
    ):
        with pytest.raises(EosmatError, match=re.escape(message)):
            validate_eosmat_document(bad)
    with pytest.raises(EosmatError, match="requires EOSMAT format_version 5"):
        validate_eosmat_document({**mapped, "format_version": 4})


def test_bundled_default_maps_and_material_summary():
    flagged = 0
    for material in list_material_documents():
        stored = json.loads((MATERIALS / f"{material}.eosmat").read_bytes())
        for entry in stored["eos_records"]:
            assert not {"default", "default_for"} & set(entry)
        flagged += "default_records" in stored
    assert flagged == 80
    mgo = get_material("mgo")
    default = mgo.default_record()
    lines = mgo.summary().splitlines()
    assert f"  Default: {default.name} ({default.identifier})" in lines
    assert "default_record()" not in mgo.summary()
    header = next(
        i for i, line in enumerate(lines) if line.split()[:2] == ["identifier", "model"]
    )
    assert lines[header + 1].split()[1] == default.identifier
    iron = get_material("iron")
    assert iron.default_record() is not None and not iron.default_record().is_default
    assert (
        f"  Default: none designated; the first record is used "
        f"({iron.default_record().identifier})"
    ) in iron.summary().splitlines()
    palladium = get_material("palladium").summary()
    assert (
        "  Default: none; selection on hold. Frost experimental fit reproduction "
        "remains unresolved."
    ) in palladium
    assert sum(line.startswith("  Default") for line in palladium.splitlines()) == 1
