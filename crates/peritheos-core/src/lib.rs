//! Native numerical core for Peritheos.
//!
//! The core API is scalar-first so it can be used without an array framework.
//! Batch and language-binding layers build on the same checked scalar model
//! methods. Working units follow the public Peritheos conventions.

mod error;
mod quadrature;
mod root;
mod validation;

pub mod isothermal;
pub mod thermal;

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
    fn volume(&self, pressure: f64) -> EosResult<f64> {
        root::solve_volume(self, pressure)
    }
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

    /// Isothermal bulk modulus from the same centered relative-volume
    /// derivative convention as the Python API.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid states, step sizes, or model evaluation.
    fn bulk_modulus(&self, volume: f64, temperature: f64, relative_step: f64) -> EosResult<f64> {
        let volume = validation::positive_state(volume, "volume")?;
        let temperature = validation::positive_state(temperature, "temperature")?;
        let relative_step = validation::positive_state(relative_step, "relative_step")?;
        let step = relative_step * volume;
        let derivative = (self.pressure(volume + step, temperature)?
            - self.pressure(volume - step, temperature)?)
            / (2.0 * step);
        validation::finite_result(-volume * derivative)
    }

    /// Isothermal compressibility in `GPa^-1`.
    ///
    /// # Errors
    ///
    /// Returns an error when bulk-modulus evaluation fails or is non-finite.
    fn isothermal_compressibility(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        validation::finite_result(1.0 / self.bulk_modulus(volume, temperature, 1.0e-6)?)
    }

    /// Volumetric thermal expansivity in `K^-1`.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid states, step sizes, or model evaluation.
    fn thermal_expansivity(
        &self,
        volume: f64,
        temperature: f64,
        relative_step: f64,
    ) -> EosResult<f64> {
        let volume = validation::positive_state(volume, "volume")?;
        let temperature = validation::positive_state(temperature, "temperature")?;
        let relative_step = validation::positive_state(relative_step, "relative_step")?;
        let step = (relative_step * temperature).min(0.49 * temperature);
        let pressure_derivative = (self.pressure(volume, temperature + step)?
            - self.pressure(volume, temperature - step)?)
            / (2.0 * step);
        validation::finite_result(
            pressure_derivative / self.bulk_modulus(volume, temperature, 1.0e-6)?,
        )
    }

    /// Volume on the invertible branch nearest the reference volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the requested root cannot be bracketed and
    /// converged on the supported branch.
    fn volume(&self, pressure: f64, temperature: f64) -> EosResult<f64> {
        let temperature = validation::positive_state(temperature, "temperature")?;
        root::solve_volume_function(
            |volume| self.pressure(volume, temperature),
            pressure,
            self.reference_eos().reference_volume(),
        )
    }

    /// Positive-temperature root nearest the reference temperature.
    ///
    /// # Errors
    ///
    /// Returns an error when the state is invalid or a positive-temperature
    /// root cannot be bracketed and converged.
    fn temperature(&self, pressure: f64, volume: f64) -> EosResult<f64> {
        let pressure = validation::finite_state(pressure, "pressure")?;
        let volume = validation::positive_state(volume, "volume")?;
        let target_thermal_pressure = pressure - self.reference_eos().pressure(volume)?;
        root::solve_temperature_function(
            |temperature| self.thermal_pressure(volume, temperature),
            target_thermal_pressure,
            self.reference_temperature(),
        )
    }

    /// Effective pressure contribution retained by a DAC confinement fraction.
    ///
    /// # Errors
    ///
    /// Returns an error unless `f_dac` is finite and lies in `[0, 1)`.
    fn dac_thermal_pressure(&self, volume: f64, temperature: f64, f_dac: f64) -> EosResult<f64> {
        let f_dac = validation::finite_state(f_dac, "f_dac")?;
        if !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        validation::finite_result(f_dac * self.thermal_pressure(volume, temperature)?)
    }

    /// Infer a heated-state temperature from ambient and heated volumes.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volumes, confinement fraction, an
    /// unbracketed root, or a root below the reference temperature.
    fn temperature_from_volumes(
        &self,
        ambient_volume: f64,
        heated_volume: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        let ambient_volume = validation::positive_state(ambient_volume, "volume")?;
        let heated_volume = validation::positive_state(heated_volume, "volume")?;
        let f_dac = validation::finite_state(f_dac, "f_dac")?;
        if !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        let target = (self.reference_eos().pressure(ambient_volume)?
            - self.reference_eos().pressure(heated_volume)?)
            / (1.0 - f_dac);
        if target < 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        let temperature = root::solve_temperature_function(
            |value| self.thermal_pressure(heated_volume, value),
            target,
            self.reference_temperature(),
        )?;
        let tolerance = 1.0e-10 * self.reference_temperature().max(1.0);
        if temperature < self.reference_temperature() - tolerance {
            Err(EosError::OutsideInvertibleRange)
        } else {
            Ok(temperature)
        }
    }
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

    /// Constant-pressure molar heat capacity in J mol^-1 K^-1.
    ///
    /// # Errors
    ///
    /// Returns an error when any required thermoelastic quantity fails.
    fn molar_heat_capacity_p(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let cv = self.molar_heat_capacity_v(volume, temperature)?;
        let alpha = self.thermal_expansivity(volume, temperature, 1.0e-5)?;
        let bulk_modulus = self.bulk_modulus(volume, temperature, 1.0e-6)?;
        validation::finite_result(cv + alpha * alpha * bulk_modulus * volume * temperature * 1.0e4)
    }

    /// Thermodynamic Gruneisen parameter.
    ///
    /// # Errors
    ///
    /// Returns an error when any required thermoelastic quantity fails.
    fn gruneisen_parameter(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let alpha = self.thermal_expansivity(volume, temperature, 1.0e-5)?;
        let bulk_modulus = self.bulk_modulus(volume, temperature, 1.0e-6)?;
        let cv = self.molar_heat_capacity_v(volume, temperature)?;
        validation::finite_result(alpha * bulk_modulus * volume * 1.0e4 / cv)
    }

    /// Adiabatic bulk modulus in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error when any required thermoelastic quantity fails.
    fn adiabatic_bulk_modulus(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let kt = self.bulk_modulus(volume, temperature, 1.0e-6)?;
        let cv = self.molar_heat_capacity_v(volume, temperature)?;
        let cp = self.molar_heat_capacity_p(volume, temperature)?;
        validation::finite_result(kt * cp / cv)
    }
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
