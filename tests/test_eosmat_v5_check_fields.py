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
    check_block,
    compared_scope,
    reproduction_why,
)
from peritheos._eosmat_format import format_audit, format_eosmat
from peritheos.eosmat import list_material_documents, validate_eosmat_document
from peritheos.errors import EosmatError
from tests.catalog_state import added_records

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
AUDITS = ROOT / "peritheos/data/audits"

#: Close/disagrees records whose stored comparison table names no failed
#: criterion (and, for close, no untested combined-error rule): no ``why``.
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


def test_why_names_failed_coefficients_from_the_stored_table():
    close = record(
        "similar",
        rows=[
            row("V0", 162.4, 162.41, similar=True, within_combined_2sigma=True),
            row("K0", 256, 263.7307, similar=True, within_combined_2sigma=False),
        ],
    )
    assert reproduction_why(close) == (
        "K0: refit 263.7 vs published 256 GPa (+3.0%), outside the combined error."
    )
    # The second spelling of the table, and two failed criteria in one row.
    disagrees = record(
        "not_reproduced",
        rows=[
            {
                "parameter": "rt_eos.K0_prime",
                "published": 4.91,
                "published_sigma": 0.12,
                "refit": 4.183456,
                "refit_sigma": None,
                "delta": -0.726544,
                "relative_delta": -0.148,
                "similar": False,
                "within_combined_2sigma": False,
                "within_reported_error": False,
            }
        ],
    )
    assert reproduction_why(disagrees) == (
        "K0': refit 4.183 vs published 4.91 (-14.8%), outside the combined "
        "error, beyond the agreement limit and outside the published error."
    )
    untested_rule = record(
        "similar",
        rows=[
            row("K0", 160, 161, similar=True, within_combined_2sigma=None),
            row("K0_prime", 4, 4.1, similar=True, within_combined_2sigma=None),
        ],
    )
    assert reproduction_why(untested_rule) == (
        "K0 and K0' have no error to test against the combined-error rule; "
        "all within the agreement limit."
    )
    single = record(
        "similar",
        rows=[row("V0", 74.7, 74.71, similar=True, within_combined_2sigma=None)],
    )
    assert reproduction_why(single) == (
        "V0 has no error to test against the combined-error rule; "
        "all within the agreement limit."
    )
    # Published values print as stored, without a trailing ".0"; the refit
    # gets at most one more decimal.
    integral = record(
        "similar",
        rows=[row("theta0", 600.0, 766.2634, similar=False)],
    )
    assert reproduction_why(integral) == (
        "theta0: refit 766.3 vs published 600 K (+27.7%), beyond the agreement limit."
    )
    # More than three failures: the three largest relative differences are
    # named and the rest counted.
    four = record(
        "not_reproduced",
        rows=[
            row("V0", 74.7, 75.0, similar=False),
            row("K0", 160, 200, similar=False),
            row("K0_prime", 4, 2, similar=False),
            row("gamma0", 1.5, 1.6, similar=False),
        ],
    )
    assert reproduction_why(four) == (
        "K0' (-50.0%), K0 (+25.0%), gamma0 (+6.7%) and 1 more: "
        "beyond the agreement limit."
    )
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
    many = record(
        "not_reproduced",
        rows=[row(f"coefficient_{i}", 1.25, 1.5 + i, similar=False) for i in range(12)],
    )
    text = reproduction_why(many)
    assert text is None or len(text) <= WHY_MAX_LENGTH


def printed(value):
    """A stored number as the why line prints it (no trailing ``.0``)."""
    integral = isinstance(value, float) and value.is_integer()
    return json.dumps(int(value) if integral else value)


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
            assert len(check["why"]) <= WHY_MAX_LENGTH
            rows = complete["evidence"]["reproduction"]["assessment"]["parameters"]
            for match in re.finditer(r"refit (\S+) vs published (\S+)", check["why"]):
                refit, published = match.groups()
                # Decimals of the mantissa ("6e-10" has none, "1.0e-09" one).
                decimals = len(published.partition("e")[0].partition(".")[2])
                assert len(refit.partition("e")[0].partition(".")[2]) <= decimals + 1
                assert any(printed(r["published"]) == published for r in rows)
    assert without == RECORDS_WITHOUT_WHY
    assert with_why == 127


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
    assert plain.why.startswith("K0: refit 263.7 vs published 256 GPa")
    assert (
        "  Check: close (refit vs published fit). K0: refit 263.7 " in plain.summary()
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
