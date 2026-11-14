"""Preserve final/preprint source equality and its effect on diagnostic inputs."""

import hashlib
import json
from unittest.mock import patch

import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA, observations
from scripts.compare_grossular_du_final_table import (
    EVIDENCE,
    FINAL,
    OUTPUT,
    PREPRINT,
    compare,
    report,
)
from scripts.refit_grossular_du_bosenick import cases, prepare


def test_final_table_tokens_widths_and_frozen_source_hashes():
    result = report()
    assert result == json.loads((DATA / OUTPUT).read_text())
    assert result["comparison"]["rows"] == 13
    assert result["comparison"]["identical"]
    assert set(result["comparison"]["maximum_absolute_differences"].values()) == {0}
    assert result["files"][PREPRINT] == (
        "d569ab2f57e87362027517a7bbfbfeb9ba76ddca0cc992eb39ae396bd1ac2b00"
    )
    assert not result["author_workflow_reproduced"]
    # The preceding fit's evidence is a dated receipt, not retroactively rewritten.
    manifest = json.loads((DATA / "du-bosenick-artifact-manifest.json").read_text())
    for item in manifest["files"]:
        assert (
            hashlib.sha256((DATA / item["file"]).read_bytes()).hexdigest()
            == (item["sha256"])
        )


def test_changed_temperature_or_printed_uncertainty_is_detected():
    old, new = observations(PREPRINT), observations(FINAL)
    alpha = json.loads((DATA / EVIDENCE).read_text())[
        "preprint_alpha_tokens_visually_transcribed"
    ]
    new[1]["temperature_k"] = "393.8"
    new[2]["cell_volume_token"] = "1668.3(2)"
    result = compare(old, new, alpha)
    assert not result["identical"]
    assert result["differing_rows"] == [
        dict(source_row="2", fields=["temperature_k"]),
        dict(source_row="3", fields=["cell_volume_token"]),
    ]


def test_final_rows_leave_all32_prepared_diagnostic_inputs_unchanged():
    # Replace just the source reader used by the latest Du controls, without
    # changing any files. Earlier background reductions are preserved because
    # the complete original nine-field source projection is exactly identical.
    final = observations(FINAL)
    assert [{key: row[key] for key in observations(PREPRINT)[0]} for row in final] == (
        observations(PREPRINT)
    )

    def final_reader(name):
        return final if name == PREPRINT else observations(name)

    for case in cases():
        before = prepare(case)
        with patch("scripts.refit_grossular_du_bosenick.observations", final_reader):
            after = prepare(case)
        assert before[0] == after[0]  # All219rows, including all3Conradpoints.
        assert before[1] == after[1]
        assert before[3] == after[3]
        assert before[2].keys() == after[2].keys()
        for key in before[2]:
            np.testing.assert_array_equal(before[2][key], after[2][key])
