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

CONVERGED_LEAST_SQUARES_MESSAGES = {
    "`gtol` termination condition is satisfied.",
    "`ftol` termination condition is satisfied.",
    "`xtol` termination condition is satisfied.",
    "Both `ftol` and `xtol` termination conditions are satisfied.",
}


def check_report(actual, saved, path="report"):
    """Require exact evidence/structure and bounded calculated float agreement.

    Source fields, row identities, hashes, flags and classifications are exact.
    Diagnostic solver floats allow 1e-4 relative drift. Residuals allow 2e-7
    in their stated
    units; coefficients and covariance use a 1e-24 floor so small coefficients
    cannot disappear behind a unit-sized absolute tolerance. Software versions
    describe the archived run. Successful least-squares termination descriptions
    may differ between supported solvers; success and every numerical check still
    agree. Latent-coordinate corrections allow 1e-5 of one reported sigma.
    Successful multistart evaluation counts and gradient optimality describe
    the solver run; counts stay positive and scaled optimality stays below
    1e-3 or twice the archived value. They are not scientific classifications.
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
            elif (
                key == "solver_message"
                and actual.get("solver_success") is True
                and saved.get("solver_success") is True
                and value in CONVERGED_LEAST_SQUARES_MESSAGES
                and saved[key] in CONVERGED_LEAST_SQUARES_MESSAGES
            ):
                continue
            elif ".multistart[" in path and key in {"nfev", "optimality"}:
                if (
                    actual.get("success") is not True
                    or saved.get("success") is not True
                ):
                    check_report(value, saved[key], child)
                elif key == "nfev":
                    if (
                        type(value) is not int
                        or type(saved[key]) is not int
                        or min(value, saved[key]) < 1
                    ):
                        raise ValueError(f"{child}: invalid solver evaluation count")
                elif not math.isfinite(value) or not 0 <= value <= max(
                    1e-3, 2 * saved[key]
                ):
                    raise ValueError(
                        f"{child}: solver optimality exceeded replay bound"
                    )
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
        if (
            (".projection." in path and "_sigma[" in path)
            or ".signed_orthogonal_distance_sigma[" in path
            or path.endswith(".maximum_parameter_difference_in_native_standard_errors")
        ):
            floor = 1e-5
        relative = 1e-4
        if not math.isclose(actual, saved, rel_tol=relative, abs_tol=floor):
            raise ValueError(
                f"{path}: calculated value differs ({actual!r}, {saved!r})"
            )
    elif actual != saved:
        raise ValueError(f"{path}: value differs")
