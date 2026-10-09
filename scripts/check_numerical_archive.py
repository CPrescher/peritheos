"""Compare calculated CSV columns without weakening source-column checks."""

import csv
import io

import numpy as np


def check_rows(saved, current, numerical_columns, *, rtol=1e-7, atol=1e-8):
    """Keep row order, columns and source strings exact; tolerate calculation roundoff."""
    assert len(saved) == len(current)
    for index, (left, right) in enumerate(zip(saved, current)):
        assert left.keys() == right.keys(), index
        for column in left:
            if column in numerical_columns and left[column] and right[column]:
                assert np.isclose(
                    float(left[column]), float(right[column]), rtol=rtol, atol=atol
                ), (index, column, left[column], right[column])
            else:
                assert left[column] == right[column], (index, column)


def check_csv(saved, current, numerical_columns, **tolerances):
    """Parse CSV archives and compare only explicitly named calculated columns numerically."""
    check_rows(
        list(csv.DictReader(io.StringIO(saved))),
        list(csv.DictReader(io.StringIO(current))),
        numerical_columns,
        **tolerances,
    )
