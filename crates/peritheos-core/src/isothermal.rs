//! Built-in isothermal equations of state.

use crate::validation::{finite_parameter, finite_result, positive_parameter, positive_state};
use crate::{EosError, EosResult, IsothermalEos};

fn eulerian_terms(v0: f64, volume: f64) -> (f64, f64, f64) {
    let eta = (v0 / volume).cbrt();
    let eta2 = eta * eta;
    let eta5 = eta2 * eta2 * eta;
    (0.5 * (eta2 - 1.0), eta5, eta5 * eta2)
}

/// Second-order Birch--Murnaghan equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BM2 {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
}

impl BM2 {
    /// Construct a second-order Birch--Murnaghan model.
    ///
    /// # Errors
    ///
    /// Returns an error unless `v0` and `k0` are positive and finite.
    pub fn new(v0: f64, k0: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
        })
    }
}

impl IsothermalEos for BM2 {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let (strain, eta5, _) = eulerian_terms(self.v0, volume);
        finite_result(3.0 * self.k0 * strain * eta5)
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let (strain, eta5, _) = eulerian_terms(self.v0, volume);
        finite_result(self.k0 * (1.0 + 7.0 * strain) * eta5)
    }
}

/// Third-order Birch--Murnaghan equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BM3 {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
}

impl BM3 {
    /// Construct a third-order Birch--Murnaghan model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite constructor parameters.
    pub fn new(v0: f64, k0: f64, k0_prime: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
            k0_prime: finite_parameter(k0_prime, "K0_prime")?,
        })
    }
}

impl IsothermalEos for BM3 {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let (strain, eta5, _) = eulerian_terms(self.v0, volume);
        finite_result(3.0 * self.k0 * strain * eta5 * (1.0 + 1.5 * (self.k0_prime - 4.0) * strain))
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let (strain, eta5, _) = eulerian_terms(self.v0, volume);
        finite_result(
            eta5 * (self.k0
                + (3.0 * self.k0 * self.k0_prime - 5.0 * self.k0) * strain
                + 13.5 * (self.k0 * self.k0_prime - 4.0 * self.k0) * strain * strain),
        )
    }
}

/// Fourth-order Birch--Murnaghan equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BM4 {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
    /// Reference second pressure derivative of the bulk modulus.
    pub k0_double_prime: f64,
}

impl BM4 {
    /// Construct a fourth-order Birch--Murnaghan model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite constructor parameters.
    pub fn new(v0: f64, k0: f64, k0_prime: f64, k0_double_prime: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
            k0_prime: finite_parameter(k0_prime, "K0_prime")?,
            k0_double_prime: finite_parameter(k0_double_prime, "K0_double_prime")?,
        })
    }

    fn coefficients(&self) -> (f64, f64) {
        let zeta = 0.75 * (4.0 - self.k0_prime);
        let xi = 0.375
            * (self.k0 * self.k0_double_prime
                + (self.k0_prime - 4.0) * (self.k0_prime - 3.0)
                + 35.0 / 9.0);
        (zeta, xi)
    }
}

impl IsothermalEos for BM4 {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let (strain, eta5, _) = eulerian_terms(self.v0, volume);
        let (zeta, xi) = self.coefficients();
        finite_result(
            3.0 * self.k0
                * strain
                * eta5
                * (1.0 - 2.0 * zeta * strain + 4.0 * xi * strain * strain),
        )
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let (strain, eta5, eta7) = eulerian_terms(self.v0, volume);
        let (zeta, xi) = self.coefficients();
        let correction = 1.0 - 2.0 * zeta * strain + 4.0 * xi * strain * strain;
        let derivative = 1.0 - 4.0 * zeta * strain + 12.0 * xi * strain * strain;
        finite_result(5.0 * strain * self.k0 * eta5 * correction + self.k0 * eta7 * derivative)
    }
}

/// Murnaghan equation of state, including its continuous zero-derivative limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Murnaghan {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
}

impl Murnaghan {
    /// Construct a Murnaghan model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite constructor parameters.
    pub fn new(v0: f64, k0: f64, k0_prime: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
            k0_prime: finite_parameter(k0_prime, "K0_prime")?,
        })
    }
}

impl IsothermalEos for Murnaghan {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    #[allow(clippy::float_cmp)]
    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let compression = (self.v0 / volume).ln();
        let pressure = if self.k0_prime == 0.0 {
            self.k0 * compression
        } else {
            self.k0 * (self.k0_prime * compression).exp_m1() / self.k0_prime
        };
        finite_result(pressure)
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        finite_result(self.k0 * (self.k0_prime * (self.v0 / volume).ln()).exp())
    }
}

fn natural_strain_pressure(v0: f64, k0: f64, a: f64, b: f64, volume: f64) -> EosResult<f64> {
    let volume = positive_state(volume, "volume")?;
    let strain = (v0 / volume).ln() / 3.0;
    finite_result(3.0 * k0 * (v0 / volume) * (strain + a * strain.powi(2) + b * strain.powi(3)))
}

fn natural_strain_bulk_modulus(v0: f64, k0: f64, a: f64, b: f64, volume: f64) -> EosResult<f64> {
    let volume = positive_state(volume, "volume")?;
    let strain = (v0 / volume).ln() / 3.0;
    finite_result(
        k0 * (v0 / volume)
            * (1.0
                + (3.0 + 2.0 * a) * strain
                + (3.0 * a + 3.0 * b) * strain.powi(2)
                + 3.0 * b * strain.powi(3)),
    )
}

/// Second-order natural-strain (Poirier--Tarantola) equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NaturalStrain2 {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
}

impl NaturalStrain2 {
    /// Construct a second-order natural-strain model.
    ///
    /// # Errors
    ///
    /// Returns an error unless `v0` and `k0` are positive and finite.
    pub fn new(v0: f64, k0: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
        })
    }
}

impl IsothermalEos for NaturalStrain2 {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        natural_strain_pressure(self.v0, self.k0, 0.0, 0.0, volume)
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        natural_strain_bulk_modulus(self.v0, self.k0, 0.0, 0.0, volume)
    }
}

/// Third-order natural-strain (Poirier--Tarantola) equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NaturalStrain3 {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
}

impl NaturalStrain3 {
    /// Construct a third-order natural-strain model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite constructor parameters.
    pub fn new(v0: f64, k0: f64, k0_prime: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
            k0_prime: finite_parameter(k0_prime, "K0_prime")?,
        })
    }

    fn a(&self) -> f64 {
        1.5 * (self.k0_prime - 2.0)
    }
}

impl IsothermalEos for NaturalStrain3 {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        natural_strain_pressure(self.v0, self.k0, self.a(), 0.0, volume)
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        natural_strain_bulk_modulus(self.v0, self.k0, self.a(), 0.0, volume)
    }
}

/// Fourth-order natural-strain (Poirier--Tarantola) equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NaturalStrain4 {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
    /// Reference second pressure derivative of the bulk modulus.
    pub k0_double_prime: f64,
}

impl NaturalStrain4 {
    /// Construct a fourth-order natural-strain model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite constructor parameters.
    pub fn new(v0: f64, k0: f64, k0_prime: f64, k0_double_prime: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
            k0_prime: finite_parameter(k0_prime, "K0_prime")?,
            k0_double_prime: finite_parameter(k0_double_prime, "K0_double_prime")?,
        })
    }

    fn coefficients(&self) -> (f64, f64) {
        let difference = self.k0_prime - 2.0;
        (
            1.5 * difference,
            1.5 * (self.k0 * self.k0_double_prime + 1.0 + difference + difference * difference),
        )
    }
}

impl IsothermalEos for NaturalStrain4 {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let (a, b) = self.coefficients();
        natural_strain_pressure(self.v0, self.k0, a, b, volume)
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let (a, b) = self.coefficients();
        natural_strain_bulk_modulus(self.v0, self.k0, a, b, volume)
    }
}

/// Modified Tait equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModifiedTait {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
    /// Reference second pressure derivative of the bulk modulus.
    pub k0_double_prime: f64,
    a: f64,
    b: f64,
    c: f64,
}

impl ModifiedTait {
    /// Construct a modified Tait model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid parameters or a singular coefficient set.
    #[allow(clippy::float_cmp)]
    pub fn new(v0: f64, k0: f64, k0_prime: f64, k0_double_prime: f64) -> EosResult<Self> {
        let v0 = positive_parameter(v0, "V0")?;
        let k0 = positive_parameter(k0, "K0")?;
        let k0_prime = finite_parameter(k0_prime, "K0_prime")?;
        let k0_double_prime = finite_parameter(k0_double_prime, "K0_double_prime")?;
        let one_plus_prime = 1.0 + k0_prime;
        let numerator_c = one_plus_prime + k0 * k0_double_prime;
        let denominator_c = k0_prime * k0_prime + k0_prime - k0 * k0_double_prime;
        if one_plus_prime == 0.0 || numerator_c == 0.0 || denominator_c == 0.0 {
            return Err(EosError::InvalidParameter {
                name: "K0_prime/K0_double_prime",
                reason: "produce a singular modified Tait EOS",
            });
        }

        Ok(Self {
            v0,
            k0,
            k0_prime,
            k0_double_prime,
            a: one_plus_prime / numerator_c,
            b: k0_prime / k0 - k0_double_prime / one_plus_prime,
            c: numerator_c / denominator_c,
        })
    }

    fn compression_base(&self, volume: f64) -> EosResult<(f64, f64)> {
        let volume = positive_state(volume, "volume")?;
        let relative_volume = volume / self.v0;
        let base = (relative_volume + self.a - 1.0) / self.a;
        if !base.is_finite() || base <= 0.0 {
            return Err(EosError::InvalidState {
                name: "volume",
                reason: "is outside the modified Tait EOS domain",
            });
        }
        Ok((relative_volume, base))
    }
}

impl IsothermalEos for ModifiedTait {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let (_, base) = self.compression_base(volume)?;
        finite_result((-base.ln() / self.c).exp_m1() / self.b)
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let (relative_volume, base) = self.compression_base(volume)?;
        finite_result(self.k0 * relative_volume * ((-1.0 / self.c - 1.0) * base.ln()).exp())
    }
}

/// Vinet equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vinet {
    /// Reference volume.
    pub v0: f64,
    /// Reference bulk modulus.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
}

impl Vinet {
    /// Construct a Vinet model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite constructor parameters.
    pub fn new(v0: f64, k0: f64, k0_prime: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            k0: positive_parameter(k0, "K0")?,
            k0_prime: finite_parameter(k0_prime, "K0_prime")?,
        })
    }
}

impl IsothermalEos for Vinet {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let compression = (volume / self.v0).cbrt();
        let eta = 1.5 * (self.k0_prime - 1.0);
        finite_result(
            3.0 * self.k0 * (1.0 - compression) / compression.powi(2)
                * (eta * (1.0 - compression)).exp(),
        )
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let compression = (volume / self.v0).cbrt();
        let eta = 1.5 * (self.k0_prime - 1.0);
        finite_result(
            self.k0 / compression.powi(2)
                * (1.0 + (eta * compression + 1.0) * (1.0 - compression))
                * (eta * (1.0 - compression)).exp(),
        )
    }
}

/// Holzapfel equation of state in the Sokolova et al. convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Holzapfel {
    /// Reference molar volume in J bar^-1 mol^-1.
    pub v0: f64,
    /// Reference bulk modulus in `GPa`.
    pub k0: f64,
    /// Reference pressure derivative of the bulk modulus.
    pub k0_prime: f64,
    /// Number of atoms in the formula unit.
    pub n: f64,
    /// Atomic number of the formula unit.
    pub z: f64,
    c0: f64,
    c2: f64,
}

impl Holzapfel {
    /// Construct a Holzapfel model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite constructor parameters.
    pub fn new(v0: f64, k0: f64, k0_prime: f64, n: f64, z: f64) -> EosResult<Self> {
        let v0 = positive_parameter(v0, "V0")?;
        let k0 = positive_parameter(k0, "K0")?;
        let k0_prime = finite_parameter(k0_prime, "K0_prime")?;
        let n = positive_parameter(n, "n")?;
        let z = positive_parameter(z, "Z")?;
        let fermigas_pressure = 1003.6 * (z * n / (v0 * 10.0)).powf(5.0 / 3.0);
        let c0 = -(3.0 * k0 / fermigas_pressure).ln();
        let c2 = 1.5 * (k0_prime - 3.0) - c0;
        Ok(Self {
            v0,
            k0,
            k0_prime,
            n,
            z,
            c0,
            c2,
        })
    }

    /// Pressure derivative of the bulk modulus using the public numerical convention.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or unless `epsilon` lies in `(0, 1)`.
    pub fn bulk_modulus_derivative(&self, volume: f64, epsilon: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let epsilon = positive_parameter(epsilon, "eps")?;
        if epsilon >= 1.0 {
            return Err(EosError::InvalidParameter {
                name: "eps",
                reason: "must be smaller than one",
            });
        }
        let upper = volume * (1.0 + epsilon);
        let lower = volume * (1.0 - epsilon);
        finite_result(
            (self.bulk_modulus(upper)? - self.bulk_modulus(lower)?)
                / (self.pressure(upper)? - self.pressure(lower)?),
        )
    }
}

impl IsothermalEos for Holzapfel {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let x = (volume / self.v0).cbrt();
        finite_result(
            3.0 * self.k0
                * (self.c0 * (1.0 - x)).exp()
                * (x.powi(-5) - x.powi(-4))
                * (1.0 + self.c2 * x - self.c2 * x * x),
        )
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let x = (volume / self.v0).cbrt();
        let correction = 1.0 + self.c2 * x * (1.0 - x);
        let bracket_1 = (5.0 - 4.0 * x) * correction;
        let bracket_2 = self.c0 * x * (1.0 - x) * correction;
        let bracket_3 = -(1.0 - x) * (self.c2 * x - 2.0 * self.c2 * x * x);
        finite_result(
            self.k0
                * x.powi(-5)
                * (self.c0 * (1.0 - x)).exp()
                * (bracket_1 + bracket_2 + bracket_3),
        )
    }

    fn bulk_modulus_derivative(&self, volume: f64, relative_step: f64) -> EosResult<f64> {
        Holzapfel::bulk_modulus_derivative(self, volume, relative_step)
    }
}

/// Evaluate the legacy Sokolova-workbook analytical Holzapfel derivative.
///
/// This coefficient-level function exists to preserve the historical Python
/// helper. New code should construct [`Holzapfel`] and call
/// [`Holzapfel::bulk_modulus_derivative`] instead.
///
/// # Errors
///
/// Returns an error for invalid inputs or a non-finite result.
pub fn holzapfel_bulk_modulus_derivative_analytical(
    v0: f64,
    volume: f64,
    bulk_modulus: f64,
    k0: f64,
    c0: f64,
    c2: f64,
) -> EosResult<f64> {
    let v0 = positive_parameter(v0, "V0")?;
    let volume = positive_state(volume, "volume")?;
    let bulk_modulus = finite_parameter(bulk_modulus, "KT")?;
    let k0 = positive_parameter(k0, "K0")?;
    let c0 = finite_parameter(c0, "c0")?;
    let c2 = finite_parameter(c2, "c2")?;
    let x = (volume / v0).cbrt();
    let correction = 1.0 + c2 * x - c2 * x * x;
    let common = (-5.0 / x.powi(2) + 4.0 / x) * correction - (1.0 / x - 1.0) * correction * c0
        + (1.0 / x - 1.0) * (c2 - 2.0 * c2 * x);
    let exponential = (c0 * (1.0 - x)).exp();
    let term_1 = 3.0 / x.powi(4) * k0 * exponential * common;
    let term_2 = k0 * exponential * c0 * common / x.powi(3);
    let term_3 = k0
        * exponential
        * ((10.0 / x.powi(3) - 4.0 / x.powi(2)) * correction
            + (-5.0 / x.powi(2) + 4.0 / x) * (c2 - 2.0 * c2 * x)
            + c0 / x.powi(2) * correction
            - (1.0 / x - 1.0) * (c2 - 2.0 * c2 * x) * c0
            - (c2 - 2.0 * c2 * x) / x.powi(2)
            - 2.0 * c2 * (1.0 / x - 1.0))
        / x.powi(3);
    finite_result((term_1 + term_2 - term_3) / (-bulk_modulus / x) / 3.0)
}
