import numpy as np

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class Holzapfel(EosBase):
    """
    Holzapfel equation of state. Introduced in Holzapfel et al. 2001. And here implemented as described in
    Sokolova et al. 2016.
    """

    def __init__(
        self, V0: float, K0: float, K0_prime: float, n: float, Z: float
    ) -> None:
        """
        Initialize the Holzapfel equation of state.

        **Note:**
        Contrary to other EOS implementations, the reference values need to be provided in the correct units.
        For consistency with the equation of states, the pressure unit is GPa - contrary to excel spreadsheets from
        Sokolova et al. 2016 where the bulk modulus is given in kbar.

        Parameters
        ----------
        V0 : float
            Reference volume [JBar^-1] (same as [cm^3/mol]/10)
        K0 : float
            Bulk modulus at reference volume in [GPa]
        K0_prime : float
            Pressure derivative of bulk modulus at reference volume
        n : float
            Number of atoms in a chemical formula
        Z : float
            Atomic number for an element or effective atomic number for a
            compound, as defined by the selected Holzapfel parameterization
        """
        super().__init__()
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self.n = validate_positive_scalar(n, "n")
        self.Z = validate_positive_scalar(Z, "Z")

        self._native = _rust.RtEos.holzapfel(
            self.V0, self.K0, self.K0_prime, self.n, self.Z
        )

    def pressure(self, V: NumericType) -> NumericType:
        """
        Calculate pressure using the Holzapfel equation of state.

        Parameters
        ----------
        V : float or np.ndarray
            Volume in [JBar^-1] (same as [cm^3/mol]/10)

        Returns
        -------
        float or np.ndarray
            Pressure in [GPa]
        """
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "pressure", V)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Calculate bulk modulus using the Holzapfel equation of state.

        Parameters
        ----------
        V : float or np.ndarray
            Volume in [JBar^-1] (same as [cm^3/mol]/10)

        Returns
        -------
        float or np.ndarray
            Bulk modulus in [GPa]
        """
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "bulk_modulus", V)

    def bulk_modulus_derivative(self, V: NumericType, eps: float = 1e-6) -> NumericType:
        """
        Compute the pressure derivative of the bulk modulus using the native
        centered-volume convention.

        Parameters
        ----------
        V : float or np.ndarray
            Volume in [JBar^-1] (same as [cm^3/mol]/10)
        eps : float
            Step size for the numerical derivative, default is 1e-6. This value is used as the
            relative difference between the two perturbed volumes. The accuracy of the numerical
            derivative is approximately `eps * 100%`.

        Returns
        -------
        float or np.ndarray
            Pressure derivative of bulk modulus (unitless)
        """
        V = validate_volume(V)
        eps = validate_positive_scalar(eps, "eps")
        if eps >= 1:
            raise ValueError("eps must be smaller than one")
        array = np.asarray(V, dtype=float)
        if array.ndim == 0:
            return float(self._native.bulk_modulus_derivative_scalar(float(array), eps))
        return np.asarray(
            self._native.bulk_modulus_derivative_array(array, eps), dtype=float
        )


def bulk_modulus_derivative_analytical(V0, V, KT, K0, c0, c2):
    """Evaluate the historical coefficient-level derivative through Rust.

    This compatibility helper mirrors the original array broadcasting. New
    code should use :meth:`Holzapfel.bulk_modulus_derivative`.
    """
    arrays = np.broadcast_arrays(
        *[np.asarray(value, dtype=float) for value in (V0, V, KT, K0, c0, c2)]
    )
    result = np.fromiter(
        (
            _rust.holzapfel_derivative_analytical(*values)
            for values in zip(*(array.flat for array in arrays))
        ),
        dtype=float,
        count=arrays[0].size,
    ).reshape(arrays[0].shape)
    if result.ndim == 0:
        return float(result)
    return result
