//! Executable pressure calibrations and pressure-scale conversion.

use crate::{validation, EosError, EosResult};

/// Functional form of a ruby R1 fluorescence pressure calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RubyCalibrationModel {
    /// `P = (A/B) * (ratio^B - 1)`.
    PowerLaw { a_gpa: f64, b: f64 },
    /// `P = A*x*(1+m*x)`, where `x = ratio - 1`.
    QuadraticShift { a_gpa: f64, m: f64 },
    /// Holzapfel's modified Freund-Ingalls form.
    HolzapfelFreundIngalls { a_gpa: f64, b: f64, c: f64 },
}

/// A published room-temperature ruby R1 fluorescence pressure scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RubyFluorescenceCalibration {
    pub identifier: &'static str,
    pub label: &'static str,
    pub doi: &'static str,
    pub reference_wavelength_nm: f64,
    pub model: RubyCalibrationModel,
}

impl RubyFluorescenceCalibration {
    /// Calculate pressure in `GPa` from the corrected wavelength ratio `lambda/lambda0`.
    ///
    /// # Errors
    ///
    /// Returns an error for a non-finite ratio, a ratio below one, or a
    /// non-finite result.
    pub fn pressure_from_ratio(self, wavelength_ratio: f64) -> EosResult<f64> {
        let ratio = validation::finite_state(wavelength_ratio, "wavelength_ratio")?;
        if ratio < 1.0 {
            return Err(EosError::InvalidState {
                name: "wavelength_ratio",
                reason: "must be at least one for non-negative pressure",
            });
        }
        let pressure = match self.model {
            RubyCalibrationModel::PowerLaw { a_gpa, b } => {
                (a_gpa / b) * ratio.powf(b).mul_add(1.0, -1.0)
            }
            RubyCalibrationModel::QuadraticShift { a_gpa, m } => {
                let x = ratio - 1.0;
                a_gpa * x * x.mul_add(m, 1.0)
            }
            RubyCalibrationModel::HolzapfelFreundIngalls { a_gpa, b, c } => {
                let exponent = ((b + c) / c) * (1.0 - ratio.powf(-c));
                (a_gpa / (b + c)) * exponent.exp_m1()
            }
        };
        validation::finite_result(pressure)
    }

    /// Invert pressure in `GPa` to the corrected wavelength ratio `lambda/lambda0`.
    ///
    /// # Errors
    ///
    /// Returns an error for negative or non-finite pressure, a state outside
    /// the invertible model domain, or a non-finite result.
    pub fn wavelength_ratio(self, pressure_gpa: f64) -> EosResult<f64> {
        let pressure = validation::finite_state(pressure_gpa, "pressure_gpa")?;
        if pressure < 0.0 {
            return Err(EosError::InvalidState {
                name: "pressure_gpa",
                reason: "must be non-negative",
            });
        }
        let ratio = match self.model {
            RubyCalibrationModel::PowerLaw { a_gpa, b } => {
                b.mul_add(pressure / a_gpa, 1.0).powf(1.0 / b)
            }
            RubyCalibrationModel::QuadraticShift { a_gpa, m } => {
                let x = (0.5 * (4.0 * m * pressure / a_gpa).ln_1p()).exp_m1() / (2.0 * m);
                1.0 + x
            }
            RubyCalibrationModel::HolzapfelFreundIngalls { a_gpa, b, c } => {
                let base = 1.0 - (c / (b + c)) * ((b + c) * pressure / a_gpa).ln_1p();
                if base <= 0.0 {
                    return Err(EosError::OutsideInvertibleRange);
                }
                base.powf(-1.0 / c)
            }
        };
        validation::finite_result(ratio)
    }

    /// Calculate pressure from a temperature-corrected R1 wavelength in nm.
    ///
    /// # Errors
    ///
    /// Returns an error for a non-positive wavelength or invalid ratio.
    pub fn pressure_from_wavelength(self, wavelength_nm: f64) -> EosResult<f64> {
        let wavelength = validation::positive_state(wavelength_nm, "wavelength_nm")?;
        self.pressure_from_ratio(wavelength / self.reference_wavelength_nm)
    }

    /// Calculate the temperature-corrected R1 wavelength in nm implied by pressure.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid pressure or a non-finite result.
    pub fn wavelength_from_pressure(self, pressure_gpa: f64) -> EosResult<f64> {
        validation::finite_result(
            self.wavelength_ratio(pressure_gpa)? * self.reference_wavelength_nm,
        )
    }
}

/// Bundled executable ruby calibrations, in chronological order.
pub const RUBY_CALIBRATIONS: [RubyFluorescenceCalibration; 5] = [
    RubyFluorescenceCalibration {
        identifier: "ruby_mao_1978",
        label: "Mao et al. (1978) ruby scale",
        doi: "10.1063/1.325277",
        reference_wavelength_nm: 694.2,
        model: RubyCalibrationModel::PowerLaw {
            a_gpa: 1904.0,
            b: 5.0,
        },
    },
    RubyFluorescenceCalibration {
        identifier: "ruby_mao_1986",
        label: "Mao, Xu, and Bell (1986) quasihydrostatic ruby scale",
        doi: "10.1029/JB091iB05p04673",
        reference_wavelength_nm: 694.24,
        model: RubyCalibrationModel::PowerLaw {
            a_gpa: 1904.0,
            b: 7.665,
        },
    },
    RubyFluorescenceCalibration {
        identifier: "ruby_dewaele_2004",
        label: "Dewaele, Loubeyre, and Mezouar (2004) revised ruby scale",
        doi: "10.1103/PhysRevB.70.094112",
        reference_wavelength_nm: 694.24,
        model: RubyCalibrationModel::PowerLaw {
            a_gpa: 1904.0,
            b: 9.5,
        },
    },
    RubyFluorescenceCalibration {
        identifier: "ruby_holzapfel_2005",
        label: "Holzapfel (2005) practical ruby scale",
        doi: "10.1080/09511920500147501",
        reference_wavelength_nm: 694.24,
        model: RubyCalibrationModel::HolzapfelFreundIngalls {
            a_gpa: 1845.0,
            b: 14.7,
            c: 7.5,
        },
    },
    RubyFluorescenceCalibration {
        identifier: "ruby_dorogokupets_oganov_2007",
        label: "Dorogokupets and Oganov (2007) ruby scale",
        doi: "10.1103/PhysRevB.75.024115",
        reference_wavelength_nm: 694.24,
        model: RubyCalibrationModel::QuadraticShift {
            a_gpa: 1884.0,
            m: 5.5,
        },
    },
];

/// Find a bundled ruby calibration by stable identifier.
#[must_use]
pub fn ruby_calibration(identifier: &str) -> Option<RubyFluorescenceCalibration> {
    RUBY_CALIBRATIONS
        .iter()
        .copied()
        .find(|calibration| calibration.identifier == identifier)
}

/// Convert ruby-derived pressure between scales through the corrected R1 ratio.
///
/// # Errors
///
/// Returns an error when source-scale inversion or target-scale evaluation
/// rejects the state.
pub fn recalculate_ruby_pressure(
    pressure_gpa: f64,
    source: RubyFluorescenceCalibration,
    target: RubyFluorescenceCalibration,
) -> EosResult<f64> {
    target.pressure_from_ratio(source.wavelength_ratio(pressure_gpa)?)
}
