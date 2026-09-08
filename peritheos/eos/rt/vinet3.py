"""Third-order Vinet pressure form of Fratanduono et al. (2020)."""

import numpy as np

from peritheos.eos import (
    EosBase,
    NumericType,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class Vinet3(EosBase):
    """Vinet form with a cubic polynomial in the exponential.

    The "3" denotes the cubic polynomial in the exponent.

    Supplemental Section S4 Eq. (2), DOI 10.1103/PhysRevLett.124.015701:
    ``P = 3*K0*(1-x)/x**2 * exp(eta*(1-x)+beta*(1-x)**2+psi*(1-x)**3)``,
    where ``x = (V/V0)**(1/3) = (rho0/rho)**(1/3)``. Volumes may use any
    consistent unit; pressure uses the unit of K0. The coefficients eta,
    beta and psi are independent; eta is not K0_prime. Ordinary Vinet is
    recovered for beta=psi=0 and eta=3*(K0_prime-1)/2.
    """

    def __init__(
        self, V0: float, K0: float, eta: float, beta: float, psi: float
    ) -> None:
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.eta = validate_finite_scalar(eta, "eta")
        self.beta = validate_finite_scalar(beta, "beta")
        self.psi = validate_finite_scalar(psi, "psi")

    def pressure(self, V: NumericType) -> NumericType:
        """Return the published pressure at volume V."""
        x = np.cbrt(np.asarray(validate_volume(V), dtype=float) / self.V0)
        y = 1.0 - x
        result = (
            3.0
            * self.K0
            * y
            / x**2
            * np.exp(y * (self.eta + y * (self.beta + y * self.psi)))
        )
        return float(result) if result.ndim == 0 else result

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """Return analytic -V*dP/dV, regular also at V0."""
        x = np.cbrt(np.asarray(validate_volume(V), dtype=float) / self.V0)
        y = 1.0 - x
        slope = self.eta + 2.0 * self.beta * y + 3.0 * self.psi * y**2
        result = (
            self.K0
            / x**2
            * np.exp(y * (self.eta + y * (self.beta + y * self.psi)))
            * (x + 2.0 * y + x * y * slope)
        )
        return float(result) if result.ndim == 0 else result
