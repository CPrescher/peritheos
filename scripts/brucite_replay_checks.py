"""Compare brucite audit replays without changing source facts or archived fits."""

from __future__ import annotations

import math

# These fields are copied from primary evidence, not obtained by optimization.
SOURCE_FIELDS = {
    "source",
    "primary_source",
    "datasets",
    "source_json_sha256",
    "volume_basis",
    "catalog_scope",
    "published_parameters",
    "published_moduli",
    "fixed_published_moduli",
    "fixed_parameters",
    "pressure_range_gpa",
    "temperature_k",
    "source_rows",
    "observations",
    "fixed_temperature_k",
    "validation_dataset",
    "fit_dataset",
}


def check_report(actual, saved, path="report"):
    """Require exact evidence/structure and bounded calculated float agreement.

    Source fields, row identities, hashes, flags and classifications are exact.
    Diagnostic solver floats allow 1e-4 relative drift. Residuals allow 2e-7
    in their stated
    units; coefficients and covariance use a 1e-24 floor so small coefficients
    cannot disappear behind a unit-sized absolute tolerance. Software versions
    describe the archived run, while optimizer descriptions still match exactly.
    """
    if type(actual) is not type(saved):
        raise ValueError(f"{path}: type differs")
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError(f"{path}: fields differ")
        for key, value in actual.items():
            child = f"{path}.{key}"
            if key in SOURCE_FIELDS:
                if value != saved[key]:
                    raise ValueError(f"{child}: source evidence differs")
            elif path.endswith(".software") and key in {"python", "numpy", "scipy"}:
                continue
            else:
                check_report(value, saved[key], child)
    elif isinstance(actual, list):
        if len(actual) != len(saved):
            raise ValueError(f"{path}: length differs")
        for index, (left, right) in enumerate(zip(actual, saved)):
            check_report(left, right, f"{path}[{index}]")
    elif isinstance(actual, float):
        coefficient = any(
            word in path
            for word in (
                "parameters",
                "covariance",
                "standard_error",
                "alpha",
                "V300",
                "V0_a3",
            )
        )
        floor = 1e-24 if coefficient else 2e-7
        relative = 1e-4
        if not math.isclose(actual, saved, rel_tol=relative, abs_tol=floor):
            raise ValueError(
                f"{path}: calculated value differs ({actual!r}, {saved!r})"
            )
    elif actual != saved:
        raise ValueError(f"{path}: value differs")
