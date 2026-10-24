"""Bound calculated audit roundoff while retaining exact Bykova source evidence."""

from __future__ import annotations

import math

CALCULATED_SECTIONS = {
    "fit",
    "rounding_diagnostic",
    "figure_check",
    "curve_check",
    "observations_check",
    "diagnostic_refit",
}
SOURCE_FIELDS = {
    "stored_coefficients_unchanged",
    "published_parameters",
    "reference_volume_per_formula_unit_a3",
    "reference_pressure_gpa",
    "P0",
    "volume_range_a3",
    "marker",
    "figure_x_pixel",
    "figure_y_pixel",
    "table10_counterpart",
}


def check_report(saved, current, path="audit", calculated=False):
    """Compare source facts exactly and only derived floating values tolerantly.

    Graph-coordinate transforms, bounds, residuals and rounding intervals allow
    1e-12 relative / 1e-10 absolute roundoff in their stated units. Optimized
    diagnostics allow 1e-5 relative / 1e-10 absolute drift. Source parameters,
    pixel readings, hashes, row identities, counts, flags and structure are exact.
    The archived values and source CSV bytes are never rewritten by this check.
    """
    if type(saved) is not type(current):
        raise ValueError(f"{path}: type differs")
    if isinstance(current, dict):
        if saved.keys() != current.keys():
            raise ValueError(f"{path}: fields differ")
        for key, value in current.items():
            child = f"{path}.{key}"
            if key in SOURCE_FIELDS:
                if saved[key] != value:
                    raise ValueError(f"{child}: source evidence differs")
            else:
                check_report(
                    saved[key], value, child, calculated or key in CALCULATED_SECTIONS
                )
    elif isinstance(current, list):
        if len(saved) != len(current):
            raise ValueError(f"{path}: length differs")
        for index, (left, right) in enumerate(zip(saved, current)):
            check_report(left, right, f"{path}[{index}]", calculated)
    elif isinstance(current, float) and calculated:
        optimized = path.startswith(("audit.fit.", "audit.diagnostic_refit."))
        relative = 1e-5 if optimized else 1e-12
        if not math.isclose(saved, current, rel_tol=relative, abs_tol=1e-10):
            raise ValueError(
                f"{path}: calculated value differs ({saved!r}, {current!r})"
            )
    elif saved != current:
        raise ValueError(f"{path}: source/assessment value differs")
