"""Portable audit replay must not weaken source or conclusion preservation."""

import copy
import json
import math

import pytest

from scripts import reproduce_bykova_2018_experimental as audit
from scripts.bykova_replay_checks import check_report


def test_calculated_graph_roundoff_is_portable():
    saved = json.loads(audit.RESULT.read_text())
    current = copy.deepcopy(saved)
    current["curve_check"]["bounds"]["pressure_bound_gpa"] += 1e-13
    current["curve_check"]["checkpoints"][0]["volume_a3_conventional_cell"] += 1e-12
    check_report(saved, current)


@pytest.mark.parametrize(
    "change", ["coefficient", "fixed_pressure", "hash", "flag", "count", "curve"]
)
def test_source_and_scientific_changes_are_rejected(change):
    saved = json.loads(audit.RESULT.read_text())
    current = copy.deepcopy(saved)
    if change == "coefficient":
        current["published_parameters"]["K0"] = math.nextafter(254.0, math.inf)
    elif change == "fixed_pressure":
        current["diagnostic_refit"]["parameters"]["P0"] += 1e-13
    elif change == "hash":
        current["table_sha256"] = "changed"
    elif change == "flag":
        current["curve_check"]["all_within_graphical_bounds"] = False
    elif change == "count":
        current["observations_check"]["observations"] += 1
    else:
        current["curve_check"]["checkpoints"][0]["pressure_gpa"] += 1e-4
    with pytest.raises(ValueError):
        check_report(saved, current)
