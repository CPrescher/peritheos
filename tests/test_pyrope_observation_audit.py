"""Source-table preservation and independent adiabatic-conversion checks."""

import hashlib
import json

import numpy as np

from scripts.audit_garnet_angel_2022 import independent_pressure, source
from scripts.audit_pyrope_angel_2022_observations import (
    DATA,
    audit,
    hartwig_audit,
    observations,
    properties,
)


def test_transcriptions_preserve_source_rows_tokens_and_missing_errors():
    manifest = json.loads((DATA / "sources.json").read_text())
    for item, count in zip(manifest["sources"], (33, 40)):
        path = DATA / item["transcription"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == item["transcription_sha256"]
        )
        rows = observations(item["transcription"])
        assert [int(row["source_row"]) for row in rows] == list(range(1, count + 1))
        assert all(row["temperature_sd_k"] == "" for row in rows)
    milani = observations("milani-2015-table-a1.csv")
    assert milani[0]["volume_token"] == "1503.98(25)"
    assert milani[-1]["volume_token"] == "1538.29(16)"
    # A separate crystallographic consistency check on all transcribed rows.
    for row in milani:
        assert (
            abs(
                float(row["lattice_a_angstrom"]) ** 3
                - float(row["unit_cell_volume_a3"])
            )
            < 0.03
        )
    zou = observations("zou-2012b-table-i.csv")
    assert zou[0]["ks_gpa_token"] == "204.9(12)"
    assert zou[-1]["ks_gpa_token"] == "235.8(39)"
    assert zou[36]["density_g_cm3_token"] == "3.8086(1)"


def test_adiabatic_properties_match_official_console_and_pressure_derivatives():
    p, _ = source("pyrope")
    receipt = json.loads((DATA / "eosfit-properties/checks.json").read_text())
    for state in receipt["states"]:
        result = properties(
            p, state["pressure_gpa"], state["temperature_k"], receipt["gas_constant"]
        )
        for key, value in result.items():
            assert abs(value - state[key]) < receipt["absolute_tolerances"][key]
        # Validate analytical KT and heat-capacity-derived alpha against
        # finite differences of the independent pressure/energy equation.
        modern = properties(p, state["pressure_gpa"], state["temperature_k"])
        v, t = modern["volume_cm3_mol"], state["temperature_k"]
        h, ht = 0.0001, 0.01
        kt = (
            -v
            * (independent_pressure(p, v + h, t) - independent_pressure(p, v - h, t))
            / (2 * h)
        )
        alpha = (
            independent_pressure(p, v, t + ht) - independent_pressure(p, v, t - ht)
        ) / (2 * ht * kt)
        np.testing.assert_allclose(
            [modern["kt_gpa"], modern["alpha_per_k"]], [kt, alpha], rtol=1e-7
        )


def test_partial_scope_and_selections_stay_explicit():
    result = audit()
    assert result["author_workflow_reproduced"] is False
    assert result["milani_expansion"]["residual_summary_a3"]["count"] == 33
    rows = result["zou_moduli"]["rows"]
    assert len(rows) == 40
    assert sum(row["volume_used_by_angel_2022"] for row in rows) == 7
    assert all(
        ("observed_minus_calculated_pressure_gpa" in row)
        == (row["temperature_k"] == 300)
        for row in rows
    )
    assert (
        result["sinogeikin_2000_scale"][
            "calculated_modulus_multiplier_to_observed_scale"
        ]
        == 0.97828
    )
    assert result["milani_expansion"]["residual_summary_a3"]["rms"] > 0.4
    assert result["zou_moduli"]["residual_summary_gpa"]["rms"] < 1


def test_recovered_hartwig_values_are_not_silently_repaired():
    meta = json.loads((DATA / "hartwig-galkin-2021-source.json").read_text())
    assert (
        hashlib.sha256((DATA / meta["transcription"]).read_bytes()).hexdigest()
        == meta["transcription_sha256"]
    )
    rows = observations(meta["transcription"])
    assert len(rows) == 12
    assert [float(row["temperature_k"]) for row in rows] == [
        103,
        153,
        203,
        253,
        297.5,
        347.5,
        391,
        451,
        500,
        550,
        600,
        650,
    ]
    assert rows[4]["source_calculated_molar_volume_cm3_mol_token"] == "113.2946"
    assert rows[10]["lattice_a_angstrom_token"] == "11.4827"
    assert rows[10]["measured_molar_volume_cm3_mol_token"] == "114.121"
    assert all(
        row["temperature_error_k"]
        == row["lattice_a_error_angstrom"]
        == row["molar_volume_error_cm3_mol"]
        == ""
        for row in rows
    )
    p, _ = source("pyrope")
    checked = hartwig_audit(p)
    assert checked["temperature_range_k"] == [103, 650]
    # Independently recognizable source conflicts survive the audit.
    assert checked["rows"][4]["source_delta_discrepancy_exceeds_printed_rounding"]
    assert checked["rows"][10]["lattice_volume_discrepancy_exceeds_printed_rounding"]
    assert (
        abs(checked["rows"][10]["measured_volume_minus_lattice_volume_cm3_mol"]) > 0.15
    )
