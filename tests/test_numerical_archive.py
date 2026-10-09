"""Archive comparisons preserve source strings, holes and calculated precision."""

import pytest

from scripts.check_numerical_archive import check_csv


def test_calculated_roundoff_is_allowed_but_source_strings_are_exact():
    saved = "source,predicted\n1.00,2.0\n"
    check_csv(saved, "source,predicted\n1.00,2.00000001\n", {"predicted"})
    with pytest.raises(AssertionError):
        check_csv(saved, "source,predicted\n1.0,2.0\n", {"predicted"})
    with pytest.raises(AssertionError):
        check_csv(saved, "source,predicted\n1.00,2.01\n", {"predicted"})


def test_missing_calculated_values_and_columns_remain_exact():
    saved = "source,predicted\n1.00,\n"
    with pytest.raises(AssertionError):
        check_csv(saved, "source,predicted\n1.00,0\n", {"predicted"})
    with pytest.raises(AssertionError):
        check_csv(saved, "source,other\n1.00,\n", {"predicted"})
