"""Display-ready overviews of materials and EOS records, and their text form.

An overview is a plain JSON-serialisable dictionary built only from stored
record fields and facts the library derives in code (model names, default
selection, the range the reader enforces). Nothing is inferred: a missing
value stays ``None``. ``EOSRecord.summary()`` and ``Material.summary()`` are
rendered from these dictionaries, so the text and the data cannot diverge.

Building an overview never reads an audit sidecar unless ``include_audit=True``.
All display wording is defined in the fixed tables below, in one place.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from decimal import ROUND_HALF_UP, Decimal, localcontext
from statistics import NormalDist
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from peritheos.materials import EOSRecord, Material

#: Version of the overview dictionary layout. A key is never removed or given
#: a different meaning without incrementing it; keys may be added.
OVERVIEW_VERSION = 1

# --- Fixed display tables ---------------------------------------------------

#: Human names for model identifiers, taken from the names the documentation
#: uses (equation pages and model inventories). Unlisted models show their
#: identifier.
MODEL_NAMES = {
    "BM2": "Birch-Murnaghan, second order",
    "BM2Reference": "Birch-Murnaghan, second order, finite-pressure reference",
    "BM3": "Birch-Murnaghan, third order",
    "BM4": "Birch-Murnaghan, fourth order",
    "Murnaghan": "Murnaghan",
    "NaturalStrain2": "Natural strain, second order",
    "NaturalStrain3": "Natural strain, third order",
    "NaturalStrain4": "Natural strain, fourth order",
    "ModifiedTait": "Modified Tait",
    "Baonza": "Baonza pseudospinodal",
    "Vinet": "Vinet",
    "RydbergStacey": "Generalized Rydberg-Stacey",
    "Holzapfel": "Holzapfel",
    "MieGruneisenDebye": "Mie-Gruneisen-Debye",
    "MieGruneisenEinstein": "Mie-Gruneisen-Einstein",
    "Tange2009Debye": "Tange 2009 MgO thermal model",
    "AsymptoticPowerLawMieGruneisenDebyeExcess": (
        "Asymptotic-power-law Debye with T2 excess"
    ),
    "QCompromiseDebye": "Q-compromise Debye thermal pressure",
    "DebyeQuadraticThermalPressure": "Debye plus quadratic thermal pressure",
    "DebyeAnharmonicHelmholtz": "Monatomic Debye and anharmonic Helmholtz energy",
    "SoundVelocityDebyeHelmholtz": "Sound-velocity quasi-Debye Helmholtz EOS",
    "Dewaele2006": "Dewaele 2006 hcp-Fe thermal pressure scale",
    "DorogokupetsOganov2007": "Dorogokupets-Oganov 2007 four-oscillator model",
    "DoubleDebyeHelmholtz": "Double-Debye Helmholtz",
    "LinearThermalPressure": "Linear thermal pressure",
    "LogVolumeThermalPressure": "Logarithmic-volume linear thermal pressure",
    "SecondOrderTaylorThermalPressure": (
        "Second-order temperature-compression thermal pressure"
    ),
    "AlphaKT": "Temperature-dependent reference state",
    "ThermalReferenceStateEOS": "Temperature-dependent reference state",
    "HollandPowellThermalPressure": "Holland-Powell thermal pressure",
    "Sokolova2016": "Multi-oscillator Gruneisen thermal pressure",
    "MultiOscillatorGruneisenThermalEOS": "Multi-oscillator Gruneisen thermal pressure",
}

#: Short labels for coefficient names; other names are shown unchanged.
PARAMETER_DISPLAY_NAMES = {
    "K0_prime": "K0'",
    "K0_double_prime": "K0''",
    "K_infinity_prime": "K_infinity'",
    "dK_dT": "dK/dT",
    "dK_dT_V": "(dK/dT)_V",
    "dKprime_dT": "dK'/dT",
}

#: ``check.compared`` as a short phrase.
COMPARED_PHRASES = {
    "published_fit": "refit vs published fit",
    "peritheos_refit": "replay of a Peritheos refit",
    "source_equation": "source equation rebuilt",
    "none": "nothing could be refitted",
}

#: One plain sentence for each valid (reproduction, compared) pair, following
#: the definitions in docs/eos-labels.md.
CHECK_EXPLANATIONS = {
    ("reproduced", "published_fit"): (
        "A refit of the source data gave the published coefficients."
    ),
    ("reproduced", "peritheos_refit"): (
        "A stored Peritheos refit was recovered under its own inputs; this says "
        "nothing about the original paper."
    ),
    ("reproduced", "source_equation"): (
        "The rebuilt source equation met its checkpoints."
    ),
    ("close", "published_fit"): (
        "A refit of the source data gave similar coefficients but did not meet "
        "the full criterion."
    ),
    ("close", "peritheos_refit"): (
        "Replaying the stored Peritheos refit gave similar coefficients but did "
        "not meet the full criterion."
    ),
    ("close", "source_equation"): (
        "The rebuilt source equation came close to its checkpoints but did not "
        "meet the full criterion."
    ),
    ("disagrees", "published_fit"): (
        "A refit of the source data did not match the published coefficients; "
        "the cause is not established."
    ),
    ("disagrees", "peritheos_refit"): (
        "Replaying the stored Peritheos refit did not match it; the cause is not "
        "established."
    ),
    ("disagrees", "source_equation"): (
        "The rebuilt source equation did not match its checkpoints; the cause is "
        "not established."
    ),
    ("untested", "none"): (
        "Nothing could be refitted or compared; this is not a negative result."
    ),
    ("untested", "source_equation"): (
        "Only the source equation or its checkpoints were rebuilt; no data were "
        "refitted. This is not a negative result."
    ),
}

# Units of coefficients other than volumes and pressures, as the implementation
# uses the stored value (docs/units.md, the equation pages and the model
# docstrings and equations). "1" marks a dimensionless coefficient. A source
# scaling the code applies is part of the unit ("10^-6 K^-1"). Thermal
# coefficients are never converted: molar thermal models use J bar^-1 mol^-1,
# the others the record volume. Anything not listed has an unknown unit (None).
_KELVIN = "K"
_DIMENSIONLESS = "1"
_MOLAR_VOLUME = "J bar^-1 mol^-1"
#: Units every model shares: K0' is dimensionless and K0'' an inverse pressure
#: (docs/equation-reference.md); temperatures are in K (docs/units.md); the
#: Gruneisen parameter gamma0 is dimensionless (MGD equation page).
_COMMON_UNITS = {
    "K0_prime": _DIMENSIONLESS,
    "K_infinity_prime": _DIMENSIONLESS,
    "K0_double_prime": "GPa^-1",
    "Tr": _KELVIN,
    "gamma0": _DIMENSIONLESS,
}
#: ``n`` counts atoms per formula unit in Debye and Einstein models
#: (docs/models/thermal.md).
_ATOMS = {"n": _DIMENSIONLESS}
_DEBYE = {"theta0": _KELVIN, "gamma0": _DIMENSIONLESS, "n": _DIMENSIONLESS}
_MGD = {**_DEBYE, "q": _DIMENSIONLESS}
_DOUBLE_DEBYE = {
    "Vp": _MOLAR_VOLUME,
    "Ve": _MOLAR_VOLUME,
    "n": _DIMENSIONLESS,
    "phi0": "J mol^-1",
    **{name: _KELVIN for name in ("theta_a0", "theta_b0", "theta_1_0", "theta_0_0")},
    **{f"a_{mode}": f"({_MOLAR_VOLUME})^-1" for mode in ("a", "b", "1", "0")},
}
#: Dimensionless exponent and shape coefficients of the double-Debye
#: temperature laws (V/Vp)^-b and of the anharmonic term (V/Ve)^kappa.
_DOUBLE_DEBYE_EXPONENTS = {
    name: _DIMENSIONLESS for name in ("b_a", "b_b", "b_1", "b_0", "kappa")
}
#: gamma = gamma0 [1 + a((V/V0)^b - 1)] (Tange et al. 2009).
_ASYMPTOTIC_GAMMA = {
    "theta0": _KELVIN,
    **_ATOMS,
    "a": _DIMENSIONLESS,
    "b": _DIMENSIONLESS,
}
#: Sokolova et al. (2016) multi-oscillator model; the code divides a_0 and e_0
#: by 10^6 (eos/thermal/multi_oscillator.py).
_MULTI_OSCILLATOR = {
    **_ATOMS,
    "QE1o": _KELVIN,
    "QE2o": _KELVIN,
    **{name: _DIMENSIONLESS for name in ("mE1", "mE2", "delta", "t", "m", "g")},
    "a_0": "10^-6 K^-1",
    "e_0": "10^-6 K^-1",
}
#: AlphaKT: alpha0 + alpha1 T in K^-1, K0(T) = K0 + dK_dT (T - Tr) or
#: 1/K0(T) = 1/K0 + beta1 dT + beta2 d(T^2) + beta3 d(T^3), and K0' changed by
#: kprime_log_coefficient (T - Tr) ln(T/Tr) (eos/thermal/linear.py).
_THERMAL_REFERENCE_STATE = {
    "alpha0": "K^-1",
    "alpha1": "K^-2",
    "dK_dT": "GPa/K",
    "beta1": "GPa^-1 K^-1",
    "beta2": "GPa^-1 K^-2",
    "beta3": "GPa^-1 K^-3",
    "kprime_log_coefficient": "K^-1",
    "dKprime_dT": "K^-1",
}
MODEL_PARAMETER_UNITS: Mapping[str, Mapping[str, str]] = {
    "Holzapfel": {"n": _DIMENSIONLESS, "Z": _DIMENSIONLESS},
    "OddInversePower": {name: "GPa" for name in ("C3", "C5", "C7", "C9")},
    # P in GPa and rho in g/cm^3 (eos/rt/density_polynomial.py).
    "DensityPolynomial3": {
        "rho0": "g/cm^3",
        "c0": "GPa",
        "c1": "GPa cm^3/g",
        "c2": "GPa cm^6/g^2",
        "c3": "GPa cm^9/g^3",
    },
    "Vinet3": {name: _DIMENSIONLESS for name in ("eta", "beta", "psi")},
    # Us = c0 + s up with Us and up in km/s.
    "LinearUsUpHugoniot": {"rho0": "g/cm^3", "c0": "km/s", "s": _DIMENSIONLESS},
    # (r_cut/r_NN)^2, an integer (eos/experimental/maltby_2024.py).
    "Maltby2024": {"shell_cutoff_squared": _DIMENSIONLESS},
    "MieGruneisenDebye": {**_MGD, "Cvmax": "J mol^-1 K^-1"},
    "MieGruneisenEinstein": _MGD,
    "DebyeTabulatedThermalPressure": _MGD,
    "DebyeQuadraticThermalPressure": {**_MGD, "A": "GPa/K^2", "m": _DIMENSIONLESS},
    "VariableQMieGruneisenDebye": {
        **_DEBYE,
        "q0": _DIMENSIONLESS,
        "q1": _DIMENSIONLESS,
    },
    "QCompromiseDebye": _DEBYE,
    "Tange2009Debye": _ASYMPTOTIC_GAMMA,
    "AsymptoticPowerLawMieGruneisenDebye": _ASYMPTOTIC_GAMMA,
    "AsymptoticPowerLawMieGruneisenDebyeExcess": {
        **_ASYMPTOTIC_GAMMA,
        "beta0": "J mol^-1 K^-2",
        "m": _DIMENSIONLESS,
    },
    "AsymptoticDebyeTabulatedPressure": {
        **_ASYMPTOTIC_GAMMA,
        "cold_volume_ratio": _DIMENSIONLESS,
    },
    "DebyeAnharmonicHelmholtz": {
        "theta0": _KELVIN,
        **_ATOMS,
        **{name: _DIMENSIONLESS for name in ("q", "b1", "b2", "b3")},
    },
    "Dewaele2006": {
        **_ATOMS,
        "theta0": _KELVIN,
        "anharmonic_a": "K^-1",
        "electronic_e": "K^-1",
        **{
            name: _DIMENSIONLESS
            for name in ("gamma_inf", "beta", "anharmonic_m", "electronic_g")
        },
    },
    # The code multiplies anharmonic_a and electronic_e by 10^-6
    # (eos/thermal/dorogokupets_oganov.py); defect_H and defect_S are H/R and
    # S/R in exp(S/x - H/(T x^2)).
    "DorogokupetsOganov2007": {
        **_ATOMS,
        **{f"theta_{mode}": _KELVIN for mode in ("B1", "B2", "E1", "E2")},
        **{f"m_{mode}": _DIMENSIONLESS for mode in ("B1", "B2", "E1", "E2")},
        "d_B1": _DIMENSIONLESS,
        "d_B2": _DIMENSIONLESS,
        **{
            name: _DIMENSIONLESS
            for name in ("gamma_inf", "beta", "anharmonic_m", "electronic_g")
        },
        "anharmonic_a": "10^-6 K^-1",
        "electronic_e": "10^-6 K^-1",
        "defect_H": _KELVIN,
        "defect_S": _DIMENSIONLESS,
    },
    "DoubleDebyeHelmholtz": {
        **_DOUBLE_DEBYE,
        **_DOUBLE_DEBYE_EXPONENTS,
        "alpha0": "K^-1",
    },
    "DoubleDebyeLogMomentHelmholtz": {
        **_DOUBLE_DEBYE,
        **_DOUBLE_DEBYE_EXPONENTS,
        "anharmonic_a": "K^-1",
    },
    "HollandPowellThermalPressure": {**_ATOMS, "theta": _KELVIN, "alpha0": "K^-1"},
    "ThermalModifiedTait": {**_ATOMS, "theta": _KELVIN, "alpha0": "K^-1"},
    "LinearThermalPressure": {"alpha_KT": "GPa/K"},
    "LogVolumeThermalPressure": {"alpha_KT_ref": "GPa/K", "dK_dT_V": "GPa/K"},
    "SecondOrderTaylorThermalPressure": {
        "c0": "GPa",
        "c1": "GPa",
        "c3": "GPa",
        "c2": "GPa/K",
        "c5": "GPa/K",
        "c4": "GPa/K^2",
        "eta0": _DIMENSIONLESS,
    },
    "Sokolova2016": _MULTI_OSCILLATOR,
    "MultiOscillatorGruneisen": _MULTI_OSCILLATOR,
    "MultiOscillatorGruneisenThermalEOS": _MULTI_OSCILLATOR,
    "AlphaKT": _THERMAL_REFERENCE_STATE,
    "ThermalReferenceStateEOS": _THERMAL_REFERENCE_STATE,
}


def parameter_unit(
    model: str, name: str, *, volume_unit: str, pressure_unit: str
) -> str | None:
    """Documented unit of one stored coefficient, or ``None`` when not stated."""
    if name == "V0":
        return volume_unit
    if name in ("K0", "P0"):
        return pressure_unit
    if name == "K0_double_prime":
        return f"{pressure_unit}^-1"
    return MODEL_PARAMETER_UNITS.get(model, {}).get(name) or _COMMON_UNITS.get(name)


# --- Number display ---------------------------------------------------------
#
# Display only: stored values are never changed. A number that is stored with
# at most six significant digits prints as stored. A longer (computed) number
# is rounded as in the reconstruction reports (scripts/report_numbers.py): to
# the decimal place of its error, the error to two significant digits when its
# leading digit is 1 or 2 and otherwise one, or else to six significant digits.
# In the summary of a refit, derived or diagnostic record, the error of an
# adjustable coefficient stored with more than three significant digits is
# rounded by the same rule, and its value to the error's decimal place
# (``format_fitted_value_error``). Published records, fixed coefficients and
# coefficients without an error keep the rule above, so a published pair is
# never shown with fewer digits than a paper may have printed.

_STORED_DIGITS = 6


def _significant_digits(value: float) -> int:
    mantissa = repr(float(value)).lower().partition("e")[0]
    digits = mantissa.replace("-", "").replace(".", "").lstrip("0")
    return len(digits.rstrip("0")) or 1


def _is_number(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float, np.number))
        and math.isfinite(float(value))
    )


def format_number(value: Any, significant: int = _STORED_DIGITS) -> str:
    """A stored number for display; long computed tails are rounded."""
    if value is None:
        return "n/a"
    if not _is_number(value):
        return str(value)
    return format(float(value), f".{significant}g")


def _report_rounding(value: float, error: float) -> tuple[str, str]:
    """The reconstruction-report rounding of a value and its positive error."""
    with localcontext() as context:
        context.prec = 50
        uncertainty = Decimal(repr(float(error)))
        exponent = uncertainty.adjusted()
        leading = int(uncertainty.scaleb(-exponent))
        place = exponent - (1 if leading in (1, 2) else 0)
        quantum = Decimal(1).scaleb(place)
        rounded_error = uncertainty.quantize(quantum, rounding=ROUND_HALF_UP)
        rounded_value = Decimal(repr(float(value))).quantize(
            quantum, rounding=ROUND_HALF_UP
        )
        if rounded_value == 0:
            rounded_value = abs(rounded_value)
        if exponent < -4 or exponent >= 6:
            digits = max(0, exponent - place)
            return (
                f"{rounded_value.scaleb(-exponent):.{digits}f}e{exponent:+d}",
                f"{rounded_error.scaleb(-exponent):.{digits}f}e{exponent:+d}",
            )
        return format(rounded_value, "f"), format(rounded_error, "f")


def _is_long(value: Any) -> bool:
    return _is_number(value) and _significant_digits(float(value)) > _STORED_DIGITS


def format_value_error(value: Any, error: Any) -> tuple[str, str | None]:
    """Display texts for a value and its error (``None`` when no error)."""
    value_text = format_number(value)
    error_text = None if error is None else format_number(error)
    if (
        _is_number(error)
        and float(error) > 0
        and _is_number(value)
        and (_is_long(value) or _is_long(error))
    ):
        rounded_value, rounded_error = _report_rounding(float(value), float(error))
        if _is_long(value):
            value_text = rounded_value
        if _is_long(error):
            error_text = rounded_error
    return value_text, error_text


def _decimal_text(number: Decimal, place: int) -> str:
    """``number`` (already rounded to ``10**place``) as fixed or e-notation.

    The notation switches where ``format(x, ".6g")`` does, so a rounded number
    reads like a stored one: ``1.81e-05``, ``0.0131``, ``232``.
    """
    if number == 0:
        number = abs(number)
        exponent = place
    else:
        exponent = number.adjusted()
    if exponent < -4 or exponent >= 6:
        digits = max(0, exponent - place)
        return f"{number.scaleb(-exponent):.{digits}f}e{exponent:+03d}"
    return f"{number:.{max(0, -place)}f}"


def format_fitted_value_error(value: Any, error: Any) -> tuple[str, str | None]:
    """Display texts for an adjustable coefficient and its error.

    For records whose errors Peritheos computed (refit, derived, diagnostic);
    a published record uses ``format_value_error``. An error stored with more
    than three significant digits is rounded as the authoring guide prescribes
    (two significant digits when its leading digit is 1 or 2, otherwise one)
    and the value to the same decimal place; a value already stored more
    coarsely keeps its stored text. An error with at most three significant
    digits, and its value, print as ``format_value_error`` prints them.
    Display only: the stored numbers are not changed.
    """
    if not (
        _is_number(value)
        and _is_number(error)
        and float(error) > 0
        and _significant_digits(float(error)) > 3
    ):
        return format_value_error(value, error)
    with localcontext() as context:
        context.prec = 50
        stored_error = Decimal(repr(float(error)))
        stored_value = Decimal(repr(float(value)))
        exponent = stored_error.adjusted()
        leading = int(stored_error.scaleb(-exponent))
        place = exponent - (1 if leading in (1, 2) else 0)
        quantum = Decimal(1).scaleb(place)
        error_text = _decimal_text(stored_error.quantize(quantum, ROUND_HALF_UP), place)
        value_digits = _significant_digits(float(value))
        value_place = (
            stored_value.adjusted() - value_digits + 1 if stored_value else place
        )
        if value_place >= place and not _is_long(value):
            return format_number(value), error_text
        rounded = stored_value.quantize(quantum, ROUND_HALF_UP)
        return _decimal_text(rounded, place), error_text


def format_range(bounds: Any) -> str | None:
    """``0-203``, ``>=5``, ``<=40`` or ``None`` when no bound is stated."""
    if bounds is None:
        return None
    lower, upper = bounds

    def text(bound: Any) -> str:
        if _is_number(bound) and _significant_digits(float(bound)) > _STORED_DIGITS:
            return format_number(bound, 4)
        return format_number(bound)

    if lower is None and upper is None:
        return None
    if upper is None:
        return f">={text(lower)}"
    if lower is None:
        return f"<={text(upper)}"
    if lower == upper:
        return text(lower)
    separator = " to " if float(lower) < 0 else "-"
    return f"{text(lower)}{separator}{text(upper)}"


# --- Building the overview --------------------------------------------------


def _plain(value: Any) -> Any:
    """JSON-safe copy: mappings to dicts, tuples to lists, numbers to floats."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _surname(author: str) -> str:
    author = " ".join(str(author).split())
    if "," in author:
        return author.split(",", 1)[0].strip()
    words = author.split()
    if any(word.endswith(".") for word in words[:-1]):
        return words[-1]
    return author


def _authors_short(authors: list[str]) -> str:
    names = [_surname(author) for author in authors if str(author).strip()]
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return f"{names[0]} and {names[1]}"
    return f"{names[0]} et al."


def _split_authors(authors: str) -> list[str]:
    """Split a joined author string; ``Surname, I.`` pairs stay together."""
    if ";" in authors:
        return [part.strip() for part in authors.split(";") if part.strip()]
    joined: list[str] = []
    for part in authors.replace(" and ", ", ").split(","):
        words = part.split()
        if not words:
            continue
        if joined and all(word.endswith(".") for word in words):
            joined[-1] = f"{joined[-1]}, {' '.join(words)}"
        else:
            joined.append(" ".join(words))
    return joined


def _citation(record: EOSRecord, metadata: Mapping[str, Any]) -> dict[str, Any]:
    """``citation``, ``short_citation``, ``title`` and ``doi`` of a record."""
    stored = metadata.get("reference")
    if isinstance(stored, str):
        return {
            "citation": stored,
            "short_citation": stored,
            "title": None,
            "doi": None,
        }
    if isinstance(stored, Mapping):
        authors = [str(author) for author in stored.get("authors") or ()]
        year = stored.get("year")
        journal = stored.get("source")
        volume = stored.get("volume")
        locator = stored.get("locator")
        title = stored.get("title")
        doi = stored.get("doi") or None
    else:
        reference = record.reference
        authors = _split_authors(reference.authors)
        year = reference.year or None
        journal = volume = locator = None
        title = reference.title or None
        doi = reference.doi or None
    lead = _authors_short(authors)
    short = " ".join(str(part) for part in (lead, year) if part)
    citation = lead + (f" ({year})" if year else "")
    venue = " ".join(str(part) for part in (journal, volume) if part)
    if locator:
        venue = f"{venue}, {locator}" if venue else str(locator)
    if venue:
        citation = f"{citation}. {venue}" if citation else venue
    if doi:
        citation = f"{citation}. doi:{doi}" if citation else f"doi:{doi}"
    return {
        "citation": citation or None,
        "short_citation": short or None,
        "title": None if title is None else " ".join(str(title).split()),
        "doi": doi,
    }


def _model_name(identifier: str | None) -> dict[str, Any] | None:
    if identifier is None:
        return None
    return {"identifier": identifier, "name": MODEL_NAMES.get(identifier, identifier)}


def _model_identifiers(
    record: EOSRecord, metadata: Mapping[str, Any]
) -> tuple[str, str | None]:
    from peritheos.eos import ThermalEOS

    component = metadata.get("eos")
    if isinstance(component, Mapping) and component.get("type"):
        thermal = metadata.get("thermal")
        thermal_type = (
            str(thermal.get("type") or thermal.get("model"))
            if isinstance(thermal, Mapping)
            else None
        )
        return str(component["type"]), thermal_type
    eos = record.eos
    if isinstance(eos, ThermalEOS):
        return type(eos.rt_eos).__name__, type(eos).__name__
    return type(eos).__name__, None


def _parameter_rows(
    record: EOSRecord,
    metadata: Mapping[str, Any],
    model: str,
    thermal_model: str | None,
) -> list[dict[str, Any]]:
    from peritheos.eos import ThermalEOS
    from peritheos.materials import _own_parameter_values

    volume_unit = str(record.volume_unit)
    pressure_unit = "GPa"
    eos = record.eos
    reference_eos = eos.rt_eos if isinstance(eos, ThermalEOS) else eos
    rows: list[dict[str, Any]] = []

    def row(component: str, name: str, value: Any, error: Any, fixed: Any) -> None:
        owner = model if component == "model" else str(thermal_model)
        rows.append(
            {
                "component": component,
                "name": name,
                "display_name": PARAMETER_DISPLAY_NAMES.get(name, name),
                "value": _plain(value),
                "error": None if error is None else _plain(error),
                "fixed": fixed,
                "unit": parameter_unit(
                    owner,
                    name,
                    volume_unit=volume_unit,
                    pressure_unit=pressure_unit,
                ),
            }
        )

    def model_order(stored: Mapping[str, Any], executable: Any) -> list[str]:
        try:
            order = list(_own_parameter_values(executable))
        except Exception:
            order = []
        return [name for name in order if name in stored] + [
            name for name in stored if name not in order
        ]

    component = metadata.get("eos")
    if isinstance(component, Mapping) and isinstance(
        component.get("parameters"), Mapping
    ):
        values = component["parameters"]
        errors = dict(metadata.get("parameter_errors") or {})
        if (
            errors.get("V0") is not None
            and metadata.get("parameter_error_volume_basis", "model") != "public"
            and record.volume_scale != 1.0
        ):
            # Historical v3/v4 files store the V0 error in model units.
            errors["V0"] = float(errors["V0"]) / record.volume_scale
        fixed = set(metadata.get("fixed_parameters") or ())
        for name in model_order(values, reference_eos):
            row("model", str(name), values[name], errors.get(name), name in fixed)
        thermal = metadata.get("thermal")
        if isinstance(thermal, Mapping) and isinstance(
            thermal.get("parameters"), Mapping
        ):
            values = thermal["parameters"]
            thermal_errors = thermal.get("parameter_errors") or {}
            thermal_fixed = set(thermal.get("fixed_parameters") or ())
            executable = eos if isinstance(eos, ThermalEOS) else None
            for name in model_order(values, executable):
                row(
                    "thermal_model",
                    str(name),
                    values[name],
                    thermal_errors.get(name),
                    name in thermal_fixed,
                )
        return rows
    # Records constructed in Python carry no record file: report the executable
    # values with the reference volume in the public unit. Fixed flags are not
    # recorded there.
    errors = dict(record.parameter_errors or {})
    prefix = "rt_eos." if isinstance(eos, ThermalEOS) else ""
    for name, value in _own_parameter_values(reference_eos).items():
        error = errors.get(prefix + name)
        if name == "V0":
            value = value / record.volume_scale
            error = None if error is None else error / record.volume_scale
        row("model", name, float(value), error, None)
    if isinstance(eos, ThermalEOS):
        for name, value in eos.parameter_values(include_reference=False).items():
            row("thermal_model", name, float(value), errors.get(name), None)
    return rows


def _bound(value: float, *, lower: bool, temperature: bool = False) -> float | None:
    if math.isinf(value):
        return None
    if lower and temperature and value == np.nextafter(0.0, 1.0):
        return None
    return float(value)


def _effective_range(bounds: Any, *, temperature: bool = False) -> list | None:
    if bounds is None:
        return None
    lower = _bound(float(bounds[0]), lower=True, temperature=temperature)
    upper = _bound(float(bounds[1]), lower=False)
    if lower is None and upper is None:
        return None
    return [lower, upper]


def _stored_range(value: Any) -> list | None:
    return None if value is None else [_plain(item) for item in value]


def _ranges(record: EOSRecord, metadata: Mapping[str, Any]) -> dict[str, Any]:
    validity = metadata.get("validity")
    stored_validity = None
    if isinstance(validity, Mapping):
        stored_validity = {
            "pressure_gpa": _stored_range(validity.get("pressure_gpa")),
            "temperature_k": _stored_range(validity.get("temperature_k")),
            "volume_ratio": _stored_range(validity.get("volume_ratio")),
        }
    experimental_keys = (
        "experimental_pressure_range_gpa",
        "experimental_temperature_range_k",
    )
    experimental = None
    if any(metadata.get(key) is not None for key in experimental_keys):
        experimental = {
            "pressure_gpa": _stored_range(metadata.get(experimental_keys[0])),
            "temperature_k": _stored_range(metadata.get(experimental_keys[1])),
            "status": metadata.get("pressure_range_status"),
        }
    enforced = record.validity
    effective = {
        "pressure_gpa": _effective_range(enforced.pressure_gpa),
        "temperature_k": _effective_range(enforced.temperature_k, temperature=True),
        "volume_ratio": _effective_range(enforced.volume_ratio),
        "particle_velocity_km_s": None,
    }
    if not any(
        effective[key] is not None
        for key in ("pressure_gpa", "temperature_k", "volume_ratio")
    ):
        source = "none"
    elif isinstance(validity, Mapping) or experimental is None:
        source = "validity"
    else:
        source = "experimental"
    domain = getattr(record, "branch_domain", None)
    if domain is not None:
        effective["particle_velocity_km_s"] = [
            float(bound) for bound in domain.particle_velocity_km_s
        ]
    return {
        "validity": stored_validity,
        "experimental": experimental,
        "effective": {"source": source, **effective},
    }


def _pressure_scale(record: EOSRecord) -> dict[str, Any] | None:
    calibration = record.audit.get("pressure_calibration")
    if not isinstance(calibration, Mapping) or not calibration:
        return None
    methods = []
    for method in calibration.get("methods") or ():
        if not isinstance(method, Mapping):
            continue
        reference = method.get("reference")
        citation = None
        if isinstance(reference, (Mapping, str)):
            citation = _citation(record, {"reference": reference})["citation"]
        methods.append(
            {
                "kind": method.get("kind"),
                "material": method.get("material"),
                "reference_eos_record": method.get("reference_eos_record"),
                "reference_calibration_record": method.get(
                    "reference_calibration_record"
                ),
                "citation": citation,
                "scope": method.get("scope"),
            }
        )
    return {"status": calibration.get("status"), "methods": methods}


def _source_document(record: EOSRecord) -> Mapping[str, Any] | None:
    source = record.__dict__.get("_audit_source")
    return None if source is None else source._document


def _containing_material(record: EOSRecord) -> Material | None:
    """The bundled catalog material holding exactly this record object."""
    document = _source_document(record)
    if document is None:
        return None
    from peritheos.catalog import _catalog_index

    material = _catalog_index().materials.get(str(document["identifier"]))
    if material is not None and any(item is record for item in material.eos_records):
        return material
    return None


def _selection(record: EOSRecord, material: Material | None) -> tuple[bool, bool]:
    """``(is_default, selected_when_no_default)`` within ``material``."""
    from peritheos.materials import HugoniotRecord

    if material is None:
        return bool(record.is_default), False
    selected = (
        material.default_hugoniot_record()
        if isinstance(record, HugoniotRecord)
        else material.default_equilibrium_record()
    )
    if selected is not record:
        return False, False
    return bool(record.is_default), not record.is_default


def record_overview(
    record: EOSRecord,
    *,
    include_audit: bool = False,
    material: Material | None = None,
) -> dict[str, Any]:
    """Overview dictionary of one record; see ``EOSRecord.overview``."""
    from peritheos._eosmat_audit import derived_report_links
    from peritheos._eosmat_v4 import default_category

    metadata = record._record_metadata()
    headline = record._headline()
    if material is None:
        material = _containing_material(record)
    document = _source_document(record)
    if material is not None:
        material_identity = {
            "identifier": material.identifier,
            "name": material.name,
            "formula": material.formula,
            "phase": material.phase,
        }
    else:
        material_identity = {
            "identifier": None if document is None else str(document["identifier"]),
            "name": (
                None
                if document is None or document.get("name") is None
                else str(document["name"])
            ),
            "formula": record.material,
            "phase": record.phase,
        }
    model, thermal_model = _model_identifiers(record, metadata)
    is_default, fallback = _selection(record, material)
    reproduction = str(headline["reproduction"])
    compared = str(headline["compared"])
    validity_notes = metadata.get("validity")
    validity_notes = (
        validity_notes.get("notes")
        if isinstance(validity_notes, Mapping)
        else list(record.validity.notes)
    )
    thermal = metadata.get("thermal")
    thermal_notes = thermal.get("notes") if isinstance(thermal, Mapping) else None
    check = metadata.get("check")
    audit_reference = check.get("audit") if isinstance(check, Mapping) else None
    material_identifier = material_identity["identifier"]
    report_url = (
        derived_report_links(str(material_identifier), record.identifier)["report_url"]
        if document is not None
        else None
    )
    category = (
        default_category(metadata)
        if metadata.get("eos") is not None
        else ("hugoniot" if record.is_hugoniot else "equilibrium")
    )
    return {
        "overview_version": OVERVIEW_VERSION,
        "identifier": record.identifier,
        "label": record.name,
        "aliases": list(record.aliases),
        "material": material_identity,
        **_citation(record, metadata),
        "model": _model_name(model),
        "thermal_model": _model_name(thermal_model),
        "category": category,
        "equation_kind": record.equation_kind,
        "record_kind": metadata.get("record_kind"),
        "determination_method": metadata.get("determination_method"),
        "reference_temperature_k": float(record.reference_temperature),
        "is_default": is_default,
        "selected_when_no_default": fallback,
        "catalog_access": metadata.get("catalog_access"),
        "pressure_scale_variant_of": record.pressure_scale_variant_of,
        "check": {
            "reproduction": reproduction,
            "compared": compared,
            "why": headline["why"],
            "warning": headline["warning"],
            "warning_reason": headline["warning_reason"],
            "source_check": headline["source_check"],
            "source_checked_on": headline["source_checked_on"],
            "source_note": headline["source_note"],
            "display": {
                "compared": COMPARED_PHRASES.get(compared, compared),
                "explanation": CHECK_EXPLANATIONS.get((reproduction, compared)),
            },
        },
        "parameters": _parameter_rows(record, metadata, model, thermal_model),
        "error_confidence": record.parameter_error_confidence,
        "covariance_parameters": (
            None
            if record.parameter_covariance is None
            else [str(name) for name in record.covariance_parameters or ()]
        ),
        "correlation": _correlation(record),
        "ranges": _ranges(record, metadata),
        "pressure_scale": _pressure_scale(record) if include_audit else None,
        "notes": [str(note) for note in record.notes],
        "thermal_notes": None if thermal_notes is None else _plain(thermal_notes),
        "validity_notes": [str(note) for note in validity_notes or ()],
        "audit": audit_reference,
        "report_url": report_url,
    }


def _hold(material: Material) -> Mapping[str, Any] | None:
    hold = material.eosmat_metadata.get("recommendation_hold")
    if isinstance(hold, Mapping) and hold.get("status") == "on_hold":
        return hold
    return None


def material_overview(
    material: Material, *, include_audit: bool = False
) -> dict[str, Any]:
    """Overview dictionary of a material; see ``Material.overview``."""
    default = material.default_record()
    designated = {}
    for category, selected in (
        ("equilibrium", material.default_equilibrium_record()),
        ("hugoniot", material.default_hugoniot_record()),
    ):
        designated[category] = (
            {"identifier": selected.identifier, "label": selected.name}
            if selected is not None and selected.is_default
            else None
        )
    hold = _hold(material)
    listed = ([default] if default is not None else []) + [
        record for record in material.eos_records if record is not default
    ]
    notes = material.eosmat_metadata.get("notes")
    return {
        "overview_version": OVERVIEW_VERSION,
        "identifier": material.identifier,
        "name": material.name,
        "formula": material.formula,
        "phase": material.phase,
        "cell_contents": material.cell_contents,
        "volume_unit": material.volume_unit,
        "formula_units_per_cell": material.formula_units_per_cell,
        "aliases": list(material.aliases),
        "family_id": material.family_id,
        "notes": None if notes is None else _plain(notes),
        "default_record": None if default is None else default.identifier,
        "default_records": designated,
        "default_hold": (
            None if hold is None or not hold.get("reason") else str(hold["reason"])
        ),
        "records": [
            record_overview(record, include_audit=include_audit, material=material)
            for record in listed
        ],
    }


# --- Text rendering ---------------------------------------------------------


def _correlation(record: EOSRecord) -> list[list[float | None]] | None:
    """Correlation coefficients of the stored covariance, in its parameter order.

    An entry is ``None`` where a variance is not positive.
    """
    covariance = record.parameter_covariance
    if covariance is None:
        return None
    sigma = [
        math.sqrt(row[index]) if row[index] > 0 else None
        for index, row in enumerate(covariance)
    ]
    return [
        [
            None if a is None or b is None else value / (a * b)
            for value, b in zip(row, sigma)
        ]
        for row, a in zip(covariance, sigma)
    ]


def _correlation_text(
    names: list[str], correlation: list[list[float | None]] | None
) -> list[str]:
    """The correlation matrix as an aligned table (two decimals)."""
    if not correlation or len(names) != len(correlation):
        return []
    width = max(6, *(len(name) for name in names)) + 2
    lines = ["    Parameter correlations (from the covariance):"]
    lines.append("      " + " " * width + "".join(name.rjust(width) for name in names))
    for name, row in zip(names, correlation):
        cells = "".join(
            ("n/a" if value is None else f"{value:.2f}").rjust(width) for value in row
        )
        lines.append("      " + name.ljust(width) + cells)
    return lines


def _line(text: Any) -> str:
    return " ".join(str(text).split())


def _model_text(model: Mapping[str, Any] | None) -> str:
    if model is None:
        return ""
    name, identifier = model["name"], model["identifier"]
    return identifier if name == identifier else f"{name} ({identifier})"


def _ranges_text(ranges: Mapping[str, Any]) -> list[str]:
    def pair(values: Mapping[str, Any] | None, *keys: str) -> str:
        if values is None:
            return ""
        labels = {"pressure_gpa": "P", "temperature_k": "T", "volume_ratio": "V/V0"}
        parts = [
            f"{labels[key]} {text}"
            for key in keys
            if (text := format_range(values.get(key))) is not None
        ]
        return ", ".join(parts)

    effective = ranges["effective"]
    lines = []
    keys = ("pressure_gpa", "temperature_k", "volume_ratio")
    validity = pair(ranges["validity"], *keys)
    experimental = pair(ranges["experimental"], *keys[:2])
    enforced = " (enforced by check_validity=True)"
    if validity:
        lines.append(
            "Validity range"
            + (enforced if effective["source"] == "validity" else "")
            + f": {validity}"
        )
    if experimental:
        lines.append(
            "Experimental data range"
            + (enforced if effective["source"] == "experimental" else "")
            + f": {experimental}"
        )
    if not validity and not experimental:
        stated = pair(effective, *keys)
        if stated:
            lines.append(f"Validity range{enforced}: {stated}")
    velocity = format_range(effective.get("particle_velocity_km_s"))
    if velocity:
        lines.append(f"Branch domain (always enforced): up {velocity} km/s")
    if not lines:
        return []
    return ["  Ranges (P in GPa, T in K):", *(f"    {line}" for line in lines)]


_CATEGORY_WORDS = {"equilibrium": "equilibrium", "hugoniot": "Hugoniot"}

#: Two-sided normal confidence of one and two standard deviations.
_ONE_SIGMA = 0.6826894921370859
_TWO_SIGMA = 0.9544997361036416


def error_convention_text(confidence: float | None) -> str:
    """How the printed parameter errors relate to a standard deviation.

    A stored level equal to the 1-sigma or 2-sigma level to four decimals
    (0.6827, 0.682689..., 0.9545) reads "1 sigma" or "2 sigma"; any other
    level is printed with the normal factor computed from the stored value,
    the factor uncertainty propagation divides by (0.95: 1.96 sigma, 0.68:
    0.994 sigma). The percentage never rounds to 0% or 100%.
    """
    if confidence is None:
        return "confidence level not stated; used as 1 sigma"
    if round(confidence, 4) == round(_ONE_SIGMA, 4):
        return "1 sigma"
    if round(confidence, 4) == round(_TWO_SIGMA, 4):
        return "2 sigma"
    # ndtri(0.5 * (1 + c)); the upper tail keeps it finite for c near 1.
    factor = (
        NormalDist().inv_cdf(0.5 * (1.0 + confidence))
        if confidence < 0.5
        else -NormalDist().inv_cdf(0.5 * (1.0 - confidence))
    )
    percent = confidence * 100
    for digits in range(4, 18):
        text = f"{percent:.{digits}g}"
        if 0.0 < float(text) < 100.0:
            break
    else:
        text = repr(percent)
    # A level near but not at 1 or 2 sigma never prints the round factor.
    digits = 3
    while digits < 17 and float(f"{factor:.{digits}g}") in (1.0, 2.0):
        digits += 1
    return f"{text}% half-widths ({factor:.{digits}g} sigma)"


#: What the printed errors are, by stored ``record_kind``; any other kind
#: (``published`` or none stated) prints "errors as published".
_ERROR_ORIGINS = {
    "refit": "errors from the Peritheos refit",
    "derived": "errors as stored for this derived record",
    "diagnostic": "errors as stored for this diagnostic record",
}

#: Crystal systems as EOSMAT stores them in ``symmetry``.
_CRYSTAL_SYSTEMS = frozenset(
    {
        "CUBIC",
        "HEXAGONAL",
        "TETRAGONAL",
        "ORTHORHOMBIC",
        "TRIGONAL",
        "RHOMBOHEDRAL",
        "MONOCLINIC",
        "TRICLINIC",
    }
)


def phase_text(phase: str) -> str:
    """A phase for display on one line.

    A material that states no ``phase`` takes its ``symmetry``, stored in
    capitals (``CUBIC``); that crystal system prints in lower case.
    """
    phase = _line(phase)
    return phase.lower() if phase in _CRYSTAL_SYSTEMS else phase


def error_origin_text(record_kind: Any) -> str:
    """``errors as published`` or what else the record's errors are."""
    return _ERROR_ORIGINS.get(str(record_kind), "errors as published")


def record_summary(overview: Mapping[str, Any]) -> str:
    """Plain text for one record overview."""
    check = overview["check"]
    material = overview["material"]
    identity = f"{material['formula']}, {phase_text(material['phase'])}"
    if material["identifier"]:
        name = material["name"] or material["identifier"]
        identity = f"{name} ({material['identifier']}); {identity}"
    lines = [
        overview["identifier"],
        f"  {_line(overview['label'])}",
        f"  Material: {identity}",
    ]
    category = _CATEGORY_WORDS.get(overview["category"], overview["category"])
    if overview["is_default"]:
        lines.append(f"  Default: the designated {category} record")
    elif overview["selected_when_no_default"]:
        lines.append(
            f"  Default: used because the material designates no {category} default"
        )
    model = _model_text(overview["model"])
    if overview["thermal_model"]:
        model += " + " + _model_text(overview["thermal_model"])
    lines.append(
        f"  Model: {model}; reference temperature "
        f"{format_number(overview['reference_temperature_k'])} K"
    )
    explanation = check["display"]["explanation"]
    check_text = f"{check['reproduction']} ({check['display']['compared']})."
    if check["why"]:
        check_text += f" {_line(check['why'])}"
    elif explanation:
        check_text += f" {explanation}"
    lines.append(f"  Check: {check_text}")
    if check["warning"]:
        warning = check["warning"]
        if check["warning_reason"]:
            warning += f" {check['warning_reason']}"
        lines.append(f"  Warning: {_line(warning)}")
    lines.append(
        f"  Parameters, {error_origin_text(overview['record_kind'])}: "
        + error_convention_text(overview["error_confidence"])
    )
    rows = overview["parameters"]
    width = max((len(row["display_name"]) for row in rows), default=0)
    # Only errors that Peritheos computed are rounded further; a published
    # error may be a paper's printed digits.
    format_adjustable = (
        format_fitted_value_error
        if str(overview["record_kind"]) in _ERROR_ORIGINS
        else format_value_error
    )
    texts = []
    for row in rows:
        if row["fixed"]:
            value, error = format_value_error(row["value"], row["error"])
            text = f"{value} fixed"
            if error is not None:
                text += f" (published +/- {error})"
            texts.append(text)
            continue
        value, error = format_adjustable(row["value"], row["error"])
        if row["value"] is None and error is None:
            text = "not stated"
        elif error is not None:
            text = f"{value} +/- {error}"
        else:
            text = f"{value} (no error reported)"
        texts.append(text)
    value_width = max((len(text) for text in texts), default=0)
    unknown = []
    for component in ("model", "thermal_model"):
        members = [
            (row, text)
            for row, text in zip(rows, texts)
            if row["component"] == component
        ]
        if not members:
            continue
        lines.append(f"    {_model_text(overview[component])}:")
        for row, text in members:
            unit = row["unit"]
            if unit is None:
                unknown.append(row["display_name"])
            unit_text = (
                ""
                if unit is None
                else ("dimensionless" if unit == _DIMENSIONLESS else unit)
            )
            lines.append(
                f"      {row['display_name']:<{width}} = {text:<{value_width}}  "
                f"{unit_text}".rstrip()
            )
    if unknown:
        lines.append(
            "    No documented unit for: "
            + ", ".join(unknown)
            + " (see the model's equation page)."
        )
    if overview["covariance_parameters"]:
        stripped = [
            name.removeprefix("rt_eos.") for name in overview["covariance_parameters"]
        ]
        names = [PARAMETER_DISPLAY_NAMES.get(name, name) for name in stripped]
        lines.append(f"    Covariance reported for: {', '.join(names)}.")
        lines.extend(_correlation_text(names, overview["correlation"]))
    lines.extend(_ranges_text(overview["ranges"]))
    if overview["pressure_scale_variant_of"]:
        lines.append(
            f"  Pressure-scale variant of: {overview['pressure_scale_variant_of']}"
        )
    if overview["citation"]:
        lines.append(f"  Reference: {_line(overview['citation'])}")
    if overview["title"]:
        lines.append(f"    {_line(overview['title'])}")
    if check["source_check"] is not None:
        text = f"  Source check: {check['source_check']}"
        if check["source_checked_on"]:
            text += f" ({check['source_checked_on']})"
        if check["source_note"]:
            text += f". {_line(check['source_note'])}"
        lines.append(text)
    if overview["report_url"]:
        lines.append(f"  Report: {overview['report_url']}")
    return "\n".join(lines)


def material_summary(overview: Mapping[str, Any]) -> str:
    """Plain text for one material overview."""
    title = overview["name"] or overview["identifier"]
    lines = [
        f"{title} ({overview['identifier']})"
        if title != overview["identifier"]
        else title,
        f"  Formula: {overview['formula']}    Phase: {phase_text(overview['phase'])}",
        f"  Cell contents: {_line(overview['cell_contents'])}",
        f"  Volume unit: {overview['volume_unit']}",
    ]
    if overview["aliases"]:
        lines.append(f"  Aliases: {', '.join(overview['aliases'])}")
    if overview["family_id"]:
        lines.append(f"  Family: {overview['family_id']}")
    records = overview["records"]
    by_id = {record["identifier"]: record for record in records}
    default = overview["default_record"]
    designated = overview["default_records"]
    hold = overview["default_hold"]
    if default is None:
        lines.append(
            "  Default: none; selection on hold. " + _line(hold)
            if hold
            else "  Default: none"
        )
    elif by_id[default]["is_default"]:
        lines.append(f"  Default: {_line(by_id[default]['label'])} ({default})")
    else:
        lines.append(
            f"  Default: none designated; the first record is used ({default})"
        )
    if default is not None and hold:
        lines.append(f"  Default hold: {_line(hold)}")
    hugoniot = designated.get("hugoniot")
    if hugoniot is not None and hugoniot["identifier"] != default:
        lines.append(
            f"  Hugoniot default: {_line(hugoniot['label'])} ({hugoniot['identifier']})"
        )
    count = len(records)
    lines.append(
        f"  {count} EOS record{'s' if count != 1 else ''} (* default, ! warning; "
        "P in GPa and T in K are the range enforced by check_validity=True, "
        "- none stated):"
    )
    header = ("", "identifier", "model", "P (GPa)", "T (K)", "check")
    rows = [header]
    references = ["reference"]
    for record in records:
        markers = ("*" if record["is_default"] else "") + (
            "!" if record["check"]["warning"] else ""
        )
        model = record["model"]["identifier"]
        if record["thermal_model"]:
            model += " + " + record["thermal_model"]["identifier"]
        effective = record["ranges"]["effective"]
        rows.append(
            (
                markers,
                record["identifier"],
                model,
                format_range(effective["pressure_gpa"]) or "-",
                format_range(effective["temperature_k"]) or "-",
                record["check"]["reproduction"],
            )
        )
        reference = record["short_citation"] or ""
        if record["record_kind"] == "refit":
            reference = (
                f"Peritheos refit: {reference}" if reference else "Peritheos refit"
            )
        references.append(reference)
    widths = [max(len(row[column]) for row in rows) for column in range(6)]
    for row, reference in zip(rows, references):
        cells = [f"{cell:<{width}}" for cell, width in zip(row, widths)]
        lines.append(("  " + "  ".join(cells) + "  " + reference).rstrip())
    grouped: dict[str, list[str]] = {}
    for record in records:
        if record["check"]["warning"]:
            grouped.setdefault(_line(record["check"]["warning"]), []).append(
                record["identifier"]
            )
    if grouped:
        lines.append("  Warnings:")
        for warning, identifiers in grouped.items():
            lines.append(f"    ! {warning}")
            lines.append(f"      {', '.join(identifiers)}")
    return "\n".join(lines)
