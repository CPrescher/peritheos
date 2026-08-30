//! Native numerical core for Peritheos.
//!
//! The core API is scalar-first so it can be used without an array framework.
//! Batch and language-binding layers build on the same checked scalar model
//! methods. Working units follow the public Peritheos conventions.

mod error;

pub use error::EosError;

/// Convenient result alias for EOS construction and evaluation.
pub type EosResult<T> = Result<T, EosError>;

/// Common behavior of an isothermal equation of state.
pub trait IsothermalEos {
    /// Reference volume of the model.
    fn reference_volume(&self) -> f64;

    /// Pressure at a positive finite volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is invalid, outside the model domain,
    /// or produces a non-finite result.
    fn pressure(&self, volume: f64) -> EosResult<f64>;

    /// Isothermal bulk modulus at a positive finite volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is invalid, outside the model domain,
    /// or produces a non-finite result.
    fn bulk_modulus(&self, volume: f64) -> EosResult<f64>;

    /// Volume on the invertible branch nearest the reference volume.
    ///
    /// # Errors
    ///
    /// Returns an error when pressure is invalid or the requested root cannot
    /// be bracketed and converged on the supported branch.
    fn volume(&self, pressure: f64) -> EosResult<f64>;
}

/// Common mechanical behavior of a thermal equation of state.
pub trait ThermalEos {
    /// Isothermal reference model used by the thermal model.
    type Reference: IsothermalEos;

    /// Reference EOS instance.
    fn reference_eos(&self) -> &Self::Reference;

    /// Reference temperature in kelvin.
    fn reference_temperature(&self) -> f64;

    /// Thermal pressure relative to the reference temperature.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid states, model-domain failures, or
    /// non-finite results.
    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64>;

    /// Total pressure at a volume and temperature.
    ///
    /// # Errors
    ///
    /// Returns an error from either component EOS or when their sum is not
    /// finite.
    fn pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let thermal = self.thermal_pressure(volume, temperature)?;
        let reference = self.reference_eos().pressure(volume)?;
        let result = reference + thermal;
        if result.is_finite() {
            Ok(result)
        } else {
            Err(EosError::NonFiniteResult)
        }
    }

    /// Volume on the invertible branch nearest the reference volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the requested root cannot be bracketed and
    /// converged on the supported branch.
    fn volume(&self, pressure: f64, temperature: f64) -> EosResult<f64>;

    /// Positive-temperature root nearest the reference temperature.
    ///
    /// # Errors
    ///
    /// Returns an error when the state is invalid or a positive-temperature
    /// root cannot be bracketed and converged.
    fn temperature(&self, pressure: f64, volume: f64) -> EosResult<f64>;
}

/// Caloric quantities supplied by thermal models with a defined potential.
pub trait CaloricEos: ThermalEos {
    /// Constant-volume molar heat capacity in J mol^-1 K^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid states, unavailable model domains, or
    /// non-finite results.
    fn molar_heat_capacity_v(&self, volume: f64, temperature: f64) -> EosResult<f64>;
}

#[cfg(test)]
mod tests {
    use super::EosError;

    #[test]
    fn errors_have_descriptive_categories() {
        let error = EosError::InvalidParameter {
            name: "V0",
            reason: "must be positive and finite",
        };
        assert_eq!(
            error.to_string(),
            "invalid parameter V0: must be positive and finite"
        );
    }
}
