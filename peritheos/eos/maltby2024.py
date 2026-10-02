"""Executable published Maltby model; unresolved source reproduction, not validated."""

from __future__ import annotations

import numpy as np

from peritheos.eos import EosBase
from peritheos.eos.experimental.maltby_2024 import Maltby2024Published
from peritheos.errors import UnsupportedOperationError


class Maltby2024(EosBase):
    """Complete Helmholtz EOS with fixed published coefficients and explicit cutoff.

    Volumes use J/bar/mol; pressure uses GPa. ``reference_volume`` is the
    characteristic volume in the source laws, not a zero-pressure volume.
    Neither parameter fitting nor DAC confinement is supported by this adapter.
    """

    reference_volume = 2.256
    Tr = 300.0

    def __init__(self, shell_cutoff_squared: int):
        if (
            isinstance(shell_cutoff_squared, bool)
            or not isinstance(shell_cutoff_squared, (int, float))
            or not np.isfinite(shell_cutoff_squared)
            or shell_cutoff_squared != int(shell_cutoff_squared)
        ):
            raise ValueError("shell_cutoff_squared must be an integer")
        shell_cutoff_squared = int(shell_cutoff_squared)
        self._model = Maltby2024Published(shell_cutoff_squared)
        self.shell_cutoff_squared = shell_cutoff_squared

    def _evaluate(self, method, first, temperature):
        result = np.vectorize(getattr(self._model, method), otypes=[float])(
            first, temperature
        )
        return float(result) if result.ndim == 0 else result

    def pressure(self, volume, temperature=300.0):
        return self._evaluate("pressure", volume, temperature)

    def volume(self, pressure, temperature=300.0):
        return self._evaluate("volume", pressure, temperature)

    def bulk_modulus(self, volume, temperature=300.0):
        return self._evaluate("bulk_modulus", volume, temperature)

    def thermal_pressure_increment(self, volume, temperature):
        return self._evaluate("thermal_pressure_increment", volume, temperature)

    def dac_thermal_pressure(self, *args, **kwargs):
        raise UnsupportedOperationError(
            "DAC confinement is not supported for Maltby2024"
        )

    def volume_with_dac_confinement(self, *args, **kwargs):
        raise UnsupportedOperationError(
            "DAC confinement is not supported for Maltby2024"
        )

    def temperature_from_volumes(self, *args, **kwargs):
        raise UnsupportedOperationError(
            "Temperature inversion is not supported for Maltby2024"
        )
