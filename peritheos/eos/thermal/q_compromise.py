"""Constant-theta, constant-gamma/V Debye thermal pressure."""

from __future__ import annotations

from peritheos.eos import (
    EosBase,
    NumericType,
    ThermalEOS,
    _native_for_exact_model,
    _native_thermal_evaluate,
    validate_finite_scalar,
    validate_positive_scalar,
)

from .mie_gruneisen import MieGruneisenDebye


class QCompromiseDebye(ThermalEOS):
    """Add ``gamma0/V0 * [E_D(theta0,T)-E_D(theta0,Tr)] / 1e4``.

    Volumes are J/bar/mol of formula units; ``n`` is atoms per formula unit.
    Both theta and gamma/V are constant. There is no q parameter. This is
    the pressure approximation of Angel et al. (2022),
    doi:10.1007/s00410-022-01918-5, not a full quasiharmonic potential.
    Pressure, derivatives, inversion, fitting and uncertainty are supported;
    caloric potentials and heat capacities are deliberately not exposed.
    """

    def __init__(
        self, rt_eos: EosBase, Tr: float, theta0: float, gamma0: float, n: float
    ) -> None:
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.theta0 = validate_positive_scalar(theta0, "theta0")
        self.gamma0 = validate_finite_scalar(gamma0, "gamma0")
        self.n = validate_positive_scalar(n, "n")
        # Energy evaluated only at V0: this helper never supplies the pressure law.
        self._energy_model = MieGruneisenDebye(
            rt_eos, self.Tr, self.theta0, self.gamma0, 0.0, self.n
        )
        reference = _native_for_exact_model(rt_eos)
        if reference is not None and type(self) is QCompromiseDebye:
            from peritheos import _rust

            self._native = _rust.ThermalEos.q_compromise_debye(
                reference, self.Tr, self.theta0, self.gamma0, self.n
            )

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        volumes, temperatures = self._broadcast_state(V, T)
        if hasattr(self, "_native"):
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        energy = self._energy_model.thermal_energy(self.rt_eos.V0, temperatures)
        reference = self._energy_model.thermal_energy(self.rt_eos.V0, self.Tr)
        return self._scalar_or_array(
            self.gamma0 * (energy - reference) / self.rt_eos.V0 / 1e4
        )
