//! Equations of state, thermodynamic properties, inversion, and material files.
//!
//! `peritheos` is the Rust API behind Peritheos. It combines checked EOS
//! calculations, ordered batch evaluation, executable `.eosmat` models,
//! fitting, and uncertainty propagation without a Python runtime.
//!
//! # Quick start: construct, evaluate, and invert
//!
//! Model constructors validate their parameters. Evaluation and inversion use
//! the common [`IsothermalEos`] trait, so changing equation families does not
//! change the calling pattern.
//!
//! ```
//! use peritheos::{isothermal::BM3, IsothermalEos};
//!
//! let eos = BM3::new(10.0, 160.0, 4.0)?;
//! let pressure = eos.pressure(9.0)?;
//! let recovered_volume = eos.volume(pressure)?;
//!
//! assert!(pressure > 0.0);
//! assert!((recovered_volume - 9.0).abs() < 1.0e-10);
//! # Ok::<(), peritheos::EosError>(())
//! ```
//!
//! # Thermal and caloric calculations
//!
//! Thermal models wrap a reference isotherm. The [`ThermalEos`] trait provides
//! total pressure and both P-T-to-V and P-V-to-T inversion; [`CaloricEos`]
//! adds heat capacities and related quantities when the model defines them.
//!
//! ```
//! use peritheos::{
//!     isothermal::BM3,
//!     thermal::MieGruneisenDebye,
//!     CaloricEos, ThermalEos,
//! };
//!
//! let reference = BM3::new(1.02, 165.0, 5.0)?;
//! let eos = MieGruneisenDebye::new(reference, 300.0, 170.0, 2.9, 1.0, 1.0)?;
//! let pressure = eos.pressure(0.95, 1_500.0)?;
//! let recovered_volume = eos.volume(pressure, 1_500.0)?;
//! let heat_capacity = eos.molar_heat_capacity_v(0.95, 1_500.0)?;
//!
//! assert!((recovered_volume - 0.95).abs() < 1.0e-9);
//! assert!(heat_capacity > 0.0);
//! # Ok::<(), peritheos::EosError>(())
//! ```
//!
//! # Material records instead of handwritten parameters
//!
//! [`load_eosmat`] and [`load_eosmat_str`] parse and validate canonical
//! Peritheos format-3 and legacy Dioptas format-2 files. A loaded [`EosRecord`]
//! is runtime dispatched but exposes the same pressure and inversion workflow:
//!
//! ```no_run
//! use peritheos::load_eosmat;
//!
//! let material = load_eosmat("gold.eosmat")?;
//! let record = material.default_record().expect("material has a default EOS");
//! let pressure = record.pressure(60.0, 300.0)?;
//! material.save("gold-copy.eosmat")?;
//! println!("{}: {pressure:.3} GPa", material.name);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Ordered batch evaluation
//!
//! Import the extension traits from [`batch`] when inputs already live in
//! slices. Batch methods preserve input order and return the first scalar
//! error; they intentionally do not impose an array or threading framework.
//!
//! ```
//! use peritheos::{batch::IsothermalEosBatch, isothermal::Vinet};
//!
//! let eos = Vinet::new(10.0, 160.0, 4.0)?;
//! let pressures = eos.pressures(&[10.0, 9.5, 9.0])?;
//! assert_eq!(pressures.len(), 3);
//! assert!(pressures.windows(2).all(|pair| pair[0] < pair[1]));
//! # Ok::<(), peritheos::EosError>(())
//! ```
//!
//! # Fitting and uncertainty
//!
//! The [`fit`] module provides bounded robust least squares, EOS-specific
//! observation types, covariance estimation, and local or Monte Carlo
//! uncertainty propagation. Its high-level routines accept model factories,
//! so the same fitting workflow works with any model implementing these
//! traits. See [`fit`] for complete examples.
//!
//! # Units and conventions
//!
//! - Pressure and bulk modulus use `GPa`.
//! - Temperature uses kelvin and must be positive.
//! - Volumes must be positive. Isothermal models accept any volume unit that
//!   is used consistently with `V0`.
//! - Thermal energy models use molar volume in `J bar^-1 mol^-1`; heat
//!   capacities and energies use molar SI units.
//! - Inversion returns the supported branch nearest the reference state.
//!
//! Browse [`isothermal`], [`thermal`], and [`fit`] for the main API families,
//! or start with the complete examples shipped with the crate.

mod error;
mod quadrature;
mod root;
mod validation;

pub mod batch;
pub mod eosmat;
pub mod fit;
pub mod isothermal;
pub mod pressure_calibration;
pub mod thermal;

pub use batch::{CaloricEosBatch, IsothermalEosBatch, ThermalEosBatch};
pub use eosmat::{
    load_eosmat, load_eosmat_reader, load_eosmat_str, material_from_value, save_eosmat,
    serialize_eosmat, validate_eosmat_document, EosRecord, EosmatError, EosmatErrorKind,
    IsothermalModel, LoadedEos, Material, ThermalModel, EOSMAT_FORMAT, EOSMAT_FORMAT_VERSION,
};
pub use error::{EosError, EosErrorKind};
pub use pressure_calibration::{
    diamond_raman_calibration, recalculate_diamond_raman_pressure, recalculate_ruby_pressure,
    ruby_calibration, DiamondRamanCalibration, DiamondRamanCalibrationModel, RubyCalibrationModel,
    RubyFluorescenceCalibration, DIAMOND_RAMAN_CALIBRATIONS, RUBY_CALIBRATIONS,
};

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

    /// Pressure derivative of the bulk modulus, `dK/dP`.
    ///
    /// The default uses the same centered relative-volume convention as the
    /// public Python fallback. Models with a stable analytical expression may
    /// override it.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid states, step sizes, or derivatives.
    fn bulk_modulus_derivative(&self, volume: f64, relative_step: f64) -> EosResult<f64> {
        let volume = validation::positive_state(volume, "volume")?;
        let relative_step = validation::positive_state(relative_step, "relative_step")?;
        if relative_step >= 1.0 {
            return Err(EosError::InvalidState {
                name: "relative_step",
                reason: "must be smaller than one",
            });
        }
        let step = relative_step * volume;
        let pressure_difference = self.pressure(volume + step)? - self.pressure(volume - step)?;
        if pressure_difference == 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        validation::finite_result(
            (self.bulk_modulus(volume + step)? - self.bulk_modulus(volume - step)?)
                / pressure_difference,
        )
    }

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

    /// Thermal-pressure increase above the reference-temperature isotherm.
    ///
    /// For reference-isotherm models this is identical to
    /// [`Self::thermal_pressure`]. Absolute free-energy models override it to
    /// remove their non-zero reference-temperature pressure contribution.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid states, model-domain failures, or
    /// non-finite results.
    fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        self.thermal_pressure(volume, temperature)
    }

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

    /// Heated volume for a cold pressure and retained DAC thermal pressure.
    ///
    /// Solves `P(V,T) = P_cold + f_dac * Delta P_thermal(V,T)` on the
    /// invertible branch nearest the reference volume. The total hot pressure
    /// is available by evaluating [`Self::pressure`] at the returned state.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid pressure, temperature, or confinement
    /// fraction, or when the volume root cannot be bracketed and converged.
    fn volume_with_dac_confinement(
        &self,
        cold_pressure: f64,
        temperature: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        let temperature = validation::positive_state(temperature, "temperature")?;
        let f_dac = validation::finite_state(f_dac, "f_dac")?;
        if !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        root::solve_volume_function(
            |volume| {
                Ok(self.pressure(volume, temperature)?
                    - self.dac_thermal_pressure(volume, temperature, f_dac)?)
            },
            cold_pressure,
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
        validation::finite_result(f_dac * self.thermal_pressure_increment(volume, temperature)?)
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
    use super::{EosError, EosErrorKind};

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
        assert_eq!(error.kind(), EosErrorKind::InvalidParameter);
        assert_eq!(error.code(), "eos.invalid_parameter");
        assert_eq!(error.field(), Some("V0"));
        assert!(error.is_validation());

        let cases = [
            (
                EosError::InvalidState {
                    name: "temperature",
                    reason: "must be positive and finite",
                },
                EosErrorKind::InvalidState,
                "eos.invalid_state",
                Some("temperature"),
                true,
            ),
            (
                EosError::OutsideInvertibleRange,
                EosErrorKind::OutsideInvertibleRange,
                "eos.outside_invertible_range",
                None,
                true,
            ),
            (
                EosError::BracketingFailed,
                EosErrorKind::BracketingFailed,
                "eos.bracketing_failed",
                None,
                true,
            ),
            (
                EosError::ConvergenceFailed,
                EosErrorKind::ConvergenceFailed,
                "eos.convergence_failed",
                None,
                false,
            ),
            (
                EosError::NonFiniteResult,
                EosErrorKind::NonFiniteResult,
                "eos.non_finite_result",
                None,
                false,
            ),
        ];
        for (error, kind, code, field, validation) in cases {
            assert_eq!(error.kind(), kind);
            assert_eq!(error.code(), code);
            assert_eq!(error.field(), field);
            assert_eq!(error.is_validation(), validation);
            assert!(!error.to_string().is_empty());
        }
    }
}
