"""Peritheos: thermodynamic equations of state calculations."""

from peritheos.uncertainty import (
    EOSUncertainty,
    ParameterUncertainty,
    PredictionUncertainty,
)

__version__ = "0.4.0"

__all__ = [
    "EOSUncertainty",
    "ParameterUncertainty",
    "PredictionUncertainty",
    "__version__",
]
