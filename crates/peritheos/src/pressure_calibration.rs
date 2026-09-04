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
pub const RUBY_CALIBRATIONS: [RubyFluorescenceCalibration; 7] = [
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
        identifier: "ruby_dewaele_2008",
        label: "Dewaele et al. (2008) hydrostatic ruby scale",
        doi: "10.1103/PhysRevB.78.104102",
        reference_wavelength_nm: 694.24,
        model: RubyCalibrationModel::PowerLaw {
            a_gpa: 1920.0,
            b: 9.61,
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
    RubyFluorescenceCalibration {
        identifier: "ruby_shen_2020",
        label: "IPPS-Ruby2020 (Shen et al., 2020)",
        doi: "10.1080/08957959.2020.1791107",
        reference_wavelength_nm: 694.25,
        model: RubyCalibrationModel::QuadraticShift {
            a_gpa: 1870.0,
            m: 5.63,
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

/// Functional form of a diamond-anvil high-frequency Raman-edge calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiamondRamanCalibrationModel {
    /// `P = A*x + B*x^2`, where `x = omega/omega0 - 1`.
    NormalizedQuadratic { a_gpa: f64, b_gpa: f64 },
    /// `P = K0*x*(1 + 0.5*(K0' - 1)*x)`.
    AkahamaQuadratic { k0_gpa: f64, k0_prime: f64 },
}

/// A published room-temperature diamond-anvil Raman-edge pressure scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiamondRamanCalibration {
    pub identifier: &'static str,
    pub label: &'static str,
    pub doi: &'static str,
    pub reference_wavenumber_cm1: f64,
    pub model: DiamondRamanCalibrationModel,
}

impl DiamondRamanCalibration {
    /// Calculate pressure in `GPa` from the edge ratio `omega/omega0`.
    ///
    /// # Errors
    ///
    /// Returns an error for a non-finite ratio, a ratio below one, or a
    /// non-finite result.
    pub fn pressure_from_ratio(self, wavenumber_ratio: f64) -> EosResult<f64> {
        let ratio = validation::finite_state(wavenumber_ratio, "wavenumber_ratio")?;
        if ratio < 1.0 {
            return Err(EosError::InvalidState {
                name: "wavenumber_ratio",
                reason: "must be at least one for non-negative pressure",
            });
        }
        let shift = ratio - 1.0;
        let pressure = match self.model {
            DiamondRamanCalibrationModel::NormalizedQuadratic { a_gpa, b_gpa } => {
                shift.mul_add(b_gpa * shift, a_gpa * shift)
            }
            DiamondRamanCalibrationModel::AkahamaQuadratic { k0_gpa, k0_prime } => {
                k0_gpa * shift * (0.5 * (k0_prime - 1.0) * shift + 1.0)
            }
        };
        validation::finite_result(pressure)
    }

    /// Invert pressure in `GPa` to the edge ratio `omega/omega0`.
    ///
    /// # Errors
    ///
    /// Returns an error for negative or non-finite pressure, or a non-finite
    /// result.
    pub fn wavenumber_ratio(self, pressure_gpa: f64) -> EosResult<f64> {
        let pressure = validation::finite_state(pressure_gpa, "pressure_gpa")?;
        if pressure < 0.0 {
            return Err(EosError::InvalidState {
                name: "pressure_gpa",
                reason: "must be non-negative",
            });
        }
        let (linear, quadratic) = match self.model {
            DiamondRamanCalibrationModel::NormalizedQuadratic { a_gpa, b_gpa } => (a_gpa, b_gpa),
            DiamondRamanCalibrationModel::AkahamaQuadratic { k0_gpa, k0_prime } => {
                (k0_gpa, 0.5 * k0_gpa * (k0_prime - 1.0))
            }
        };
        let shift =
            ((linear * linear + 4.0 * quadratic * pressure).sqrt() - linear) / (2.0 * quadratic);
        validation::finite_result(1.0 + shift)
    }

    /// Calculate pressure from the diamond-edge wavenumber in `cm^-1`.
    ///
    /// # Errors
    ///
    /// Returns an error for a non-positive wavenumber or invalid ratio.
    pub fn pressure_from_wavenumber(self, wavenumber_cm1: f64) -> EosResult<f64> {
        let wavenumber = validation::positive_state(wavenumber_cm1, "wavenumber_cm1")?;
        self.pressure_from_ratio(wavenumber / self.reference_wavenumber_cm1)
    }

    /// Calculate the diamond-edge wavenumber in `cm^-1` implied by pressure.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid pressure or a non-finite result.
    pub fn wavenumber_from_pressure(self, pressure_gpa: f64) -> EosResult<f64> {
        validation::finite_result(
            self.wavenumber_ratio(pressure_gpa)? * self.reference_wavenumber_cm1,
        )
    }
}

/// Bundled executable diamond Raman calibrations, in chronological order.
pub const DIAMOND_RAMAN_CALIBRATIONS: [DiamondRamanCalibration; 2] = [
    DiamondRamanCalibration {
        identifier: "diamond_raman_akahama_2006",
        label: "Akahama and Kawamura (2006) diamond-anvil Raman-edge scale",
        doi: "10.1063/1.2335683",
        reference_wavenumber_cm1: 1334.0,
        model: DiamondRamanCalibrationModel::AkahamaQuadratic {
            k0_gpa: 547.0,
            k0_prime: 3.75,
        },
    },
    DiamondRamanCalibration {
        identifier: "diamond_raman_eremets_2023",
        label: "Eremets et al. (2023) universal diamond-edge Raman scale",
        doi: "10.1038/s41467-023-36429-9",
        reference_wavenumber_cm1: 1332.5,
        model: DiamondRamanCalibrationModel::NormalizedQuadratic {
            a_gpa: 517.0,
            b_gpa: 764.0,
        },
    },
];

/// Find a bundled diamond Raman calibration by stable identifier.
#[must_use]
pub fn diamond_raman_calibration(identifier: &str) -> Option<DiamondRamanCalibration> {
    DIAMOND_RAMAN_CALIBRATIONS
        .iter()
        .copied()
        .find(|calibration| calibration.identifier == identifier)
}

/// Convert diamond-edge pressure between scales through `omega/omega0`.
///
/// # Errors
///
/// Returns an error when source-scale inversion or target-scale evaluation
/// rejects the state.
pub fn recalculate_diamond_raman_pressure(
    pressure_gpa: f64,
    source: DiamondRamanCalibration,
    target: DiamondRamanCalibration,
) -> EosResult<f64> {
    target.pressure_from_ratio(source.wavenumber_ratio(pressure_gpa)?)
}
