//! Validation, round-trip serialization, and executable evaluation of
//! Peritheos `.eosmat` material files.
//!
//! The loader accepts canonical Peritheos format 3 and legacy Dioptas format
//! 2 documents. JSON extensions are retained in [`Material::document`] and
//! [`EosRecord::document`], while equation construction is restricted to the
//! built-in model registry.

use std::collections::HashMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;

use serde::Deserialize;
use serde_json::Value;

use crate::hugoniot::{Hugoniot, LinearUsUpHugoniot};
use crate::isothermal::{
    Holzapfel, ModifiedTait, Murnaghan, NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2,
    BM3, BM4,
};
use crate::thermal::{
    AsymptoticPowerLawMieGruneisenDebye, DebyeTemperatureLaw, DorogokupetsOganov2007,
    DorogokupetsOganov2007Parameters, DoubleDebyeHelmholtz, DoubleDebyeLogMomentHelmholtz,
    LinearThermalPressure, LogVolumeThermalPressure, MieGruneisenDebye, MieGruneisenEinstein,
    MultiOscillatorGruneisen, ReferenceStateEos, ReferenceVolumeLaw, SokolovaParameters,
    ThermalExpansionLaw, ThermalModifiedTait, ThermalReferenceState,
};
use crate::{EosError, EosResult, IsothermalEos, ThermalEos};

/// Canonical `.eosmat` format discriminator.
pub const EOSMAT_FORMAT: &str = "peritheos.material";
/// Current canonical `.eosmat` format version.
pub const EOSMAT_FORMAT_VERSION: u64 = 3;
const LEGACY_FORMAT_VERSION: u64 = 2;
const CELL_ANGSTROM3_TO_FORMULA_MOLAR_J_PER_BAR: f64 = 0.060_221_407_6;
const HUGONIOT_MASS_BASIS_RELATIVE_TOLERANCE: f64 = 1.0e-3;

/// Machine-readable category for an [`EosmatError`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum EosmatErrorKind {
    /// The underlying filesystem operation failed.
    Io,
    /// JSON decoding or typed deserialization failed.
    Json,
    /// The document structure, format, or units are invalid.
    InvalidDocument,
    /// An individual EOS record cannot be constructed.
    InvalidRecord,
}

/// Errors encountered while decoding or constructing an `.eosmat` material.
#[derive(Debug)]
#[non_exhaustive]
pub enum EosmatError {
    /// The file could not be opened, read, or written.
    Io(std::io::Error),
    /// The file is not valid JSON or does not match the expected JSON types.
    Json(serde_json::Error),
    /// The material-level document is invalid or uses an unsupported format.
    InvalidDocument(String),
    /// An EOS record could not be converted to an executable model.
    InvalidRecord {
        /// Stable identifier, or a generated identifier for legacy records.
        identifier: String,
        /// Human-readable validation or construction failure.
        reason: String,
    },
}

impl EosmatError {
    /// Return the machine-readable error category.
    #[must_use]
    pub const fn kind(&self) -> EosmatErrorKind {
        match self {
            Self::Io(_) => EosmatErrorKind::Io,
            Self::Json(_) => EosmatErrorKind::Json,
            Self::InvalidDocument(_) => EosmatErrorKind::InvalidDocument,
            Self::InvalidRecord { .. } => EosmatErrorKind::InvalidRecord,
        }
    }

    /// Return a stable, language-independent error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Io(_) => "eosmat.io",
            Self::Json(_) => "eosmat.json",
            Self::InvalidDocument(_) => "eosmat.invalid_document",
            Self::InvalidRecord { .. } => "eosmat.invalid_record",
        }
    }

    /// Return the invalid EOS record identifier, when the error is record-local.
    #[must_use]
    pub fn record_identifier(&self) -> Option<&str> {
        match self {
            Self::InvalidRecord { identifier, .. } => Some(identifier),
            Self::Io(_) | Self::Json(_) | Self::InvalidDocument(_) => None,
        }
    }
}

impl Display for EosmatError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "could not access eosmat file: {error}"),
            Self::Json(error) => write!(formatter, "invalid eosmat JSON: {error}"),
            Self::InvalidDocument(reason) => write!(formatter, "invalid eosmat document: {reason}"),
            Self::InvalidRecord { identifier, reason } => {
                write!(formatter, "invalid EOS record {identifier:?}: {reason}")
            }
        }
    }
}

impl Error for EosmatError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Json(error) => Some(error),
            Self::InvalidDocument(_) | Self::InvalidRecord { .. } => None,
        }
    }
}

impl From<std::io::Error> for EosmatError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<serde_json::Error> for EosmatError {
    fn from(error: serde_json::Error) -> Self {
        Self::Json(error)
    }
}

/// Runtime-dispatched built-in isothermal model.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum IsothermalModel {
    /// Second-order Birch--Murnaghan EOS.
    BM2(BM2),
    /// Third-order Birch--Murnaghan EOS.
    BM3(BM3),
    /// Fourth-order Birch--Murnaghan EOS.
    BM4(BM4),
    /// Holzapfel EOS.
    Holzapfel(Holzapfel),
    /// Modified Tait EOS.
    ModifiedTait(ModifiedTait),
    /// Murnaghan EOS.
    Murnaghan(Murnaghan),
    /// Second-order natural-strain EOS.
    NaturalStrain2(NaturalStrain2),
    /// Third-order natural-strain EOS.
    NaturalStrain3(NaturalStrain3),
    /// Fourth-order natural-strain EOS.
    NaturalStrain4(NaturalStrain4),
    /// Vinet EOS.
    Vinet(Vinet),
}

impl IsothermalModel {
    /// Stable mechanism-oriented `.eosmat` model identifier.
    #[must_use]
    pub const fn model_identifier(&self) -> &'static str {
        match self {
            Self::BM2(_) => "birch_murnaghan_2",
            Self::BM3(_) => "birch_murnaghan_3",
            Self::BM4(_) => "birch_murnaghan_4",
            Self::Holzapfel(_) => "holzapfel",
            Self::ModifiedTait(_) => "modified_tait",
            Self::Murnaghan(_) => "murnaghan",
            Self::NaturalStrain2(_) => "natural_strain_2",
            Self::NaturalStrain3(_) => "natural_strain_3",
            Self::NaturalStrain4(_) => "natural_strain_4",
            Self::Vinet(_) => "vinet",
        }
    }

    /// Application-facing model name used by `.eosmat`.
    #[must_use]
    pub const fn model_name(&self) -> &'static str {
        match self {
            Self::BM2(_) => "BM2",
            Self::BM3(_) => "BM3",
            Self::BM4(_) => "BM4",
            Self::Holzapfel(_) => "Holzapfel",
            Self::ModifiedTait(_) => "ModifiedTait",
            Self::Murnaghan(_) => "Murnaghan",
            Self::NaturalStrain2(_) => "NaturalStrain2",
            Self::NaturalStrain3(_) => "NaturalStrain3",
            Self::NaturalStrain4(_) => "NaturalStrain4",
            Self::Vinet(_) => "Vinet",
        }
    }
}

macro_rules! dispatch_isothermal {
    ($self:expr, $model:ident => $expression:expr) => {
        match $self {
            IsothermalModel::BM2($model) => $expression,
            IsothermalModel::BM3($model) => $expression,
            IsothermalModel::BM4($model) => $expression,
            IsothermalModel::Holzapfel($model) => $expression,
            IsothermalModel::ModifiedTait($model) => $expression,
            IsothermalModel::Murnaghan($model) => $expression,
            IsothermalModel::NaturalStrain2($model) => $expression,
            IsothermalModel::NaturalStrain3($model) => $expression,
            IsothermalModel::NaturalStrain4($model) => $expression,
            IsothermalModel::Vinet($model) => $expression,
        }
    };
}

impl IsothermalEos for IsothermalModel {
    fn reference_volume(&self) -> f64 {
        dispatch_isothermal!(self, model => model.reference_volume())
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        dispatch_isothermal!(self, model => model.pressure(volume))
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        dispatch_isothermal!(self, model => model.bulk_modulus(volume))
    }
}

impl ReferenceStateEos for IsothermalModel {
    fn reference_bulk_modulus(&self) -> f64 {
        dispatch_isothermal!(self, model => model.reference_bulk_modulus())
    }

    fn with_reference_state(&self, volume: f64, bulk_modulus: f64) -> EosResult<Self> {
        match self {
            Self::BM2(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::BM2),
            Self::BM3(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::BM3),
            Self::BM4(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::BM4),
            Self::Holzapfel(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::Holzapfel),
            Self::ModifiedTait(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::ModifiedTait),
            Self::Murnaghan(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::Murnaghan),
            Self::NaturalStrain2(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::NaturalStrain2),
            Self::NaturalStrain3(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::NaturalStrain3),
            Self::NaturalStrain4(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::NaturalStrain4),
            Self::Vinet(model) => model
                .with_reference_state(volume, bulk_modulus)
                .map(Self::Vinet),
        }
    }
}

/// Runtime-dispatched built-in shock Hugoniot model.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum HugoniotModel {
    /// Linear shock-velocity--particle-velocity relation.
    LinearUsUp(LinearUsUpHugoniot),
}

/// Coupled mechanical state on a shock Hugoniot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HugoniotState {
    /// Volume in the record's public conventional-cell unit.
    pub volume: f64,
    /// Density in g cm^-3.
    pub density: f64,
    /// Pressure in `GPa`.
    pub pressure: f64,
    /// Particle velocity in km s^-1.
    pub particle_velocity: f64,
    /// Shock velocity in km s^-1.
    pub shock_velocity: f64,
    /// Specific internal-energy increase in MJ kg^-1.
    pub specific_internal_energy_change: f64,
}

/// Initial loading history for a Hugoniot record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HugoniotLoadingPath {
    /// Single shock from the ordinary initial state.
    Principal,
    /// Single shock from a statically precompressed state.
    Precompressed,
}

/// Whether the represented states retain or transform the precursor phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HugoniotBranchKind {
    /// The represented and precursor phases are the same.
    Untransformed,
    /// The represented phase differs from the precursor phase.
    Transformed,
}

/// Scientific meaning attached to a declared branch domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HugoniotDomainKind {
    /// A phase-stability interval.
    PhaseStability,
    /// The interval covered by observations or a fit.
    ExperimentalCoverage,
    /// A source-recommended usage interval.
    Recommended,
}

/// Provenance strength of a branch-domain boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HugoniotBoundaryStatus {
    /// Numerical endpoints were reported by the source.
    ReportedExactly,
    /// The source described the interval qualitatively.
    ReportedQualitatively,
    /// The endpoints were inferred during curation.
    Inferred,
}

/// Typed precursor state retained from an `.eosmat` Hugoniot record.
#[derive(Clone, Debug, PartialEq)]
pub struct HugoniotInitialState {
    /// Precursor phase label.
    pub phase: String,
    /// Stable identifier of the precursor material document.
    pub material_identifier: String,
    /// Optional precursor equilibrium-EOS record identifier.
    pub eos_record_identifier: Option<String>,
    /// Initial temperature in kelvin.
    pub temperature_k: f64,
    /// Initial pressure in `GPa`.
    pub pressure_gpa: f64,
    /// Initial density in g cm^-3.
    pub density_g_cm3: f64,
}

/// Operational mass basis shared by `V0` and every evaluated volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HugoniotVolumeBasis {
    /// Number of chemical formula units represented by one public volume.
    pub formula_units: f64,
    /// Molar mass of one chemical formula in g mol^-1.
    pub molar_mass_g_mol: f64,
}

/// Declared particle-velocity interval for one phase-specific branch.
#[derive(Clone, Debug, PartialEq)]
pub struct HugoniotBranchDomain {
    /// Inclusive particle-velocity interval in km s^-1.
    pub particle_velocity_km_s: [f64; 2],
    /// Scientific meaning of the interval.
    pub kind: HugoniotDomainKind,
    /// Provenance strength of its endpoints.
    pub boundary_status: HugoniotBoundaryStatus,
    /// Qualifications retained from the record.
    pub notes: Vec<String>,
}

/// Typed metadata required by a phase-specific Hugoniot record.
#[derive(Clone, Debug, PartialEq)]
pub struct HugoniotRecordMetadata {
    /// Initial loading history.
    pub loading_path: HugoniotLoadingPath,
    /// Whether this branch transforms the precursor phase.
    pub branch_kind: HugoniotBranchKind,
    /// Resolvable precursor state.
    pub initial_state: HugoniotInitialState,
    /// Enforced mass basis.
    pub volume_basis: HugoniotVolumeBasis,
    /// Declared branch domain.
    pub branch_domain: HugoniotBranchDomain,
}

impl HugoniotModel {
    /// Stable mechanism-oriented `.eosmat` model identifier.
    #[must_use]
    pub const fn model_identifier(&self) -> &'static str {
        match self {
            Self::LinearUsUp(_) => "linear_us_up_hugoniot",
        }
    }

    /// Application-facing model name used by `.eosmat`.
    #[must_use]
    pub const fn model_name(&self) -> &'static str {
        match self {
            Self::LinearUsUp(_) => "LinearUsUpHugoniot",
        }
    }

    fn shock_velocity_from_particle_velocity(&self, particle_velocity: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => {
                model.shock_velocity_from_particle_velocity(particle_velocity)
            }
        }
    }

    fn pressure_from_particle_velocity(&self, particle_velocity: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.pressure_from_particle_velocity(particle_velocity),
        }
    }

    fn volume_from_particle_velocity(&self, particle_velocity: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.volume_from_particle_velocity(particle_velocity),
        }
    }
}

impl Hugoniot for HugoniotModel {
    fn reference_volume(&self) -> f64 {
        match self {
            Self::LinearUsUp(model) => model.reference_volume(),
        }
    }

    fn initial_density(&self) -> f64 {
        match self {
            Self::LinearUsUp(model) => model.initial_density(),
        }
    }

    fn initial_pressure(&self) -> f64 {
        match self {
            Self::LinearUsUp(model) => model.initial_pressure(),
        }
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.pressure(volume),
        }
    }

    fn volume(&self, pressure: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.volume(pressure),
        }
    }

    fn particle_velocity(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.particle_velocity(volume),
        }
    }

    fn shock_velocity(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.shock_velocity(volume),
        }
    }

    fn density(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.density(volume),
        }
    }

    fn specific_internal_energy_change(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.specific_internal_energy_change(volume),
        }
    }

    fn tangent_modulus(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::LinearUsUp(model) => model.tangent_modulus(volume),
        }
    }
}

/// Runtime-dispatched built-in thermal model.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum ThermalModel {
    /// Tange-type asymptotic-power-law Mie--Gruneisen--Debye EOS.
    AsymptoticPowerLawMieGruneisenDebye(AsymptoticPowerLawMieGruneisenDebye<IsothermalModel>),
    /// Vinet cold curve plus absolute double-Debye Helmholtz contribution.
    DoubleDebyeHelmholtz(DoubleDebyeHelmholtz),
    /// Vinet double-Debye Helmholtz model constrained by the logarithmic moment.
    DoubleDebyeLogMomentHelmholtz(DoubleDebyeLogMomentHelmholtz),
    /// Dorogokupets--Oganov (2007) four-oscillator Helmholtz model.
    DorogokupetsOganov2007(DorogokupetsOganov2007<IsothermalModel>),
    /// Constant-slope thermal pressure.
    LinearThermalPressure(LinearThermalPressure<IsothermalModel>),
    /// Logarithmic-volume thermal pressure.
    LogVolumeThermalPressure(LogVolumeThermalPressure<IsothermalModel>),
    /// Mie--Gruneisen--Debye EOS.
    MieGruneisenDebye(MieGruneisenDebye<IsothermalModel>),
    /// Mie--Gruneisen--Einstein EOS.
    MieGruneisenEinstein(MieGruneisenEinstein<IsothermalModel>),
    /// Generic multi-oscillator Gruneisen thermal EOS.
    MultiOscillatorGruneisen(MultiOscillatorGruneisen<IsothermalModel>),
    /// Holland--Powell thermal modified Tait EOS.
    ThermalModifiedTait(ThermalModifiedTait),
    /// Temperature-dependent reference-state EOS.
    ThermalReferenceState(ThermalReferenceState<IsothermalModel>),
}

macro_rules! dispatch_thermal {
    ($self:expr, $model:ident => $expression:expr) => {
        match $self {
            ThermalModel::AsymptoticPowerLawMieGruneisenDebye($model) => $expression,
            ThermalModel::DoubleDebyeHelmholtz($model) => $expression,
            ThermalModel::DoubleDebyeLogMomentHelmholtz($model) => $expression,
            ThermalModel::DorogokupetsOganov2007($model) => $expression,
            ThermalModel::LinearThermalPressure($model) => $expression,
            ThermalModel::LogVolumeThermalPressure($model) => $expression,
            ThermalModel::MieGruneisenDebye($model) => $expression,
            ThermalModel::MieGruneisenEinstein($model) => $expression,
            ThermalModel::MultiOscillatorGruneisen($model) => $expression,
            ThermalModel::ThermalModifiedTait($model) => $expression,
            ThermalModel::ThermalReferenceState($model) => $expression,
        }
    };
}

impl ThermalModel {
    /// Stable mechanism-oriented `.eosmat` model identifier.
    #[must_use]
    pub const fn model_identifier(&self) -> &'static str {
        match self {
            Self::AsymptoticPowerLawMieGruneisenDebye(_) => {
                "asymptotic_power_law_mie_gruneisen_debye"
            }
            Self::DoubleDebyeHelmholtz(_) => "double_debye_helmholtz",
            Self::DoubleDebyeLogMomentHelmholtz(_) => "double_debye_log_moment_helmholtz",
            Self::DorogokupetsOganov2007(_) => "dorogokupets_oganov_2007",
            Self::LinearThermalPressure(_) => "linear_thermal_pressure",
            Self::LogVolumeThermalPressure(_) => "log_volume_thermal_pressure",
            Self::MieGruneisenDebye(_) => "mie_gruneisen_debye",
            Self::MieGruneisenEinstein(_) => "mie_gruneisen_einstein",
            Self::MultiOscillatorGruneisen(_) => "multi_oscillator_gruneisen_thermal_pressure",
            Self::ThermalModifiedTait(_) => "thermal_modified_tait",
            Self::ThermalReferenceState(_) => "thermal_reference_state",
        }
    }

    /// Reference temperature in kelvin.
    #[must_use]
    pub fn reference_temperature(&self) -> f64 {
        dispatch_thermal!(self, model => model.reference_temperature())
    }

    fn reference_volume(&self) -> f64 {
        dispatch_thermal!(self, model => model.reference_eos().reference_volume())
    }

    fn pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.pressure(volume, temperature))
    }

    fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.thermal_pressure_increment(volume, temperature))
    }

    fn dac_thermal_pressure(&self, volume: f64, temperature: f64, f_dac: f64) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.dac_thermal_pressure(volume, temperature, f_dac))
    }

    fn bulk_modulus(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.bulk_modulus(volume, temperature, 1.0e-6))
    }

    fn volume(&self, pressure: f64, temperature: f64) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.volume(pressure, temperature))
    }

    fn volume_with_dac_confinement(
        &self,
        cold_pressure: f64,
        temperature: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.volume_with_dac_confinement(cold_pressure, temperature, f_dac))
    }
}

/// An executable equation selected from an `.eosmat` record.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum LoadedEos {
    /// Isothermal equation with no thermal correction.
    Isothermal(IsothermalModel),
    /// Composed reference isotherm and thermal correction.
    Thermal(ThermalModel),
    /// One-dimensional pressure-volume shock path.
    Hugoniot(HugoniotModel),
}

impl LoadedEos {
    /// Whether the record has a thermal correction.
    #[must_use]
    pub const fn is_thermal(&self) -> bool {
        matches!(self, Self::Thermal(_))
    }

    /// Whether the record describes a shock Hugoniot path.
    #[must_use]
    pub const fn is_hugoniot(&self) -> bool {
        matches!(self, Self::Hugoniot(_))
    }

    /// Stable model identifier for the record's primary equation.
    #[must_use]
    pub const fn model_identifier(&self) -> &'static str {
        match self {
            Self::Isothermal(model) => model.model_identifier(),
            Self::Thermal(model) => model.model_identifier(),
            Self::Hugoniot(model) => model.model_identifier(),
        }
    }

    /// Reference-isotherm model identifier.
    ///
    /// For a Hugoniot this returns its path-model identifier for backward
    /// compatibility. Prefer [`Self::model_identifier`] when the equation kind
    /// is not already known.
    #[must_use]
    pub const fn isothermal_model_identifier(&self) -> &'static str {
        match self {
            Self::Isothermal(model) => model.model_identifier(),
            Self::Thermal(model) => match model {
                ThermalModel::AsymptoticPowerLawMieGruneisenDebye(value) => {
                    value.rt_eos.model_identifier()
                }
                ThermalModel::DoubleDebyeHelmholtz(_)
                | ThermalModel::DoubleDebyeLogMomentHelmholtz(_) => "vinet",
                ThermalModel::DorogokupetsOganov2007(value) => value.rt_eos.model_identifier(),
                ThermalModel::LinearThermalPressure(value) => value.rt_eos.model_identifier(),
                ThermalModel::LogVolumeThermalPressure(value) => value.rt_eos.model_identifier(),
                ThermalModel::MieGruneisenDebye(value) => value.rt_eos.model_identifier(),
                ThermalModel::MieGruneisenEinstein(value) => value.rt_eos.model_identifier(),
                ThermalModel::MultiOscillatorGruneisen(value) => value.rt_eos.model_identifier(),
                ThermalModel::ThermalModifiedTait(_) => "modified_tait",
                ThermalModel::ThermalReferenceState(value) => value.rt_eos.model_identifier(),
            },
            Self::Hugoniot(model) => model.model_identifier(),
        }
    }

    /// Thermal model identifier, if present.
    #[must_use]
    pub const fn thermal_model_identifier(&self) -> Option<&'static str> {
        match self {
            Self::Thermal(model) => Some(model.model_identifier()),
            Self::Isothermal(_) | Self::Hugoniot(_) => None,
        }
    }
}

/// One executable EOS record from a material document.
#[derive(Clone, Debug, PartialEq)]
pub struct EosRecord {
    /// Stable record identifier.
    pub identifier: String,
    /// Human-readable record label.
    pub label: String,
    /// Whether this is the document's preferred record.
    pub is_default: bool,
    /// Executable built-in EOS.
    pub eos: LoadedEos,
    /// Typed shock-path metadata, present exactly for Hugoniot records.
    pub hugoniot: Option<HugoniotRecordMetadata>,
    /// Reference temperature in kelvin.
    pub reference_temperature: f64,
    /// Original record including all extension fields.
    pub document: Value,
    volume_scale: f64,
}

/// Typed view of one phase-specific Hugoniot record.
#[derive(Clone, Copy, Debug)]
pub struct HugoniotRecord<'a> {
    /// Common EOS-record fields and executable model.
    pub record: &'a EosRecord,
    /// Required path, precursor, mass-basis, and domain metadata.
    pub metadata: &'a HugoniotRecordMetadata,
}

impl HugoniotRecord<'_> {
    fn model(&self) -> HugoniotModel {
        match self.record.eos {
            LoadedEos::Hugoniot(model) => model,
            LoadedEos::Isothermal(_) | LoadedEos::Thermal(_) => {
                unreachable!("typed Hugoniot record must contain a Hugoniot model")
            }
        }
    }

    fn validate_domain(&self, particle_velocity: f64) -> EosResult<()> {
        let [lower, upper] = self.metadata.branch_domain.particle_velocity_km_s;
        let tolerance = 16.0 * f64::EPSILON * lower.abs().max(upper.abs()).max(1.0);
        if particle_velocity.is_finite()
            && particle_velocity >= lower - tolerance
            && particle_velocity <= upper + tolerance
        {
            Ok(())
        } else {
            Err(EosError::InvalidState {
                name: "particle_velocity",
                reason: "is outside the declared Hugoniot branch domain",
            })
        }
    }

    fn model_volume_and_particle_velocity(&self, volume: f64) -> EosResult<(f64, f64)> {
        let model_volume = volume * self.record.volume_scale;
        let particle_velocity = self.model().particle_velocity(model_volume)?;
        self.validate_domain(particle_velocity)?;
        Ok((model_volume, particle_velocity))
    }

    /// Return particle velocity in km s^-1, enforcing the branch domain.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or a state outside the branch domain.
    pub fn particle_velocity(&self, volume: f64) -> EosResult<f64> {
        self.model_volume_and_particle_velocity(volume)
            .map(|(_, particle_velocity)| particle_velocity)
    }

    /// Return shock velocity in km s^-1, enforcing the branch domain.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or a state outside the branch domain.
    pub fn shock_velocity(&self, volume: f64) -> EosResult<f64> {
        let (model_volume, _) = self.model_volume_and_particle_velocity(volume)?;
        self.model().shock_velocity(model_volume)
    }

    /// Return density in g cm^-3, enforcing the branch domain.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or a state outside the branch domain.
    pub fn density(&self, volume: f64) -> EosResult<f64> {
        let (model_volume, _) = self.model_volume_and_particle_velocity(volume)?;
        self.model().density(model_volume)
    }

    /// Return specific internal-energy increase in MJ kg^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or a state outside the branch domain.
    pub fn specific_internal_energy_change(&self, volume: f64) -> EosResult<f64> {
        let (model_volume, _) = self.model_volume_and_particle_velocity(volume)?;
        self.model().specific_internal_energy_change(model_volume)
    }

    /// Return tangent stiffness in `GPa`, enforcing the branch domain.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or a state outside the branch domain.
    pub fn tangent_modulus(&self, volume: f64) -> EosResult<f64> {
        let (model_volume, _) = self.model_volume_and_particle_velocity(volume)?;
        self.model().tangent_modulus(model_volume)
    }

    /// Return a coupled shock state, enforcing the declared branch domain.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid particle velocity or a state outside the domain.
    pub fn state_from_particle_velocity(&self, particle_velocity: f64) -> EosResult<HugoniotState> {
        self.validate_domain(particle_velocity)?;
        let model = self.model();
        let model_volume = model.volume_from_particle_velocity(particle_velocity)?;
        Ok(HugoniotState {
            volume: model_volume / self.record.volume_scale,
            density: model.density(model_volume)?,
            pressure: model.pressure_from_particle_velocity(particle_velocity)?,
            particle_velocity,
            shock_velocity: model.shock_velocity_from_particle_velocity(particle_velocity)?,
            specific_internal_energy_change: model.specific_internal_energy_change(model_volume)?,
        })
    }
}

impl EosRecord {
    /// Return a typed Hugoniot view when this record represents a shock path.
    #[must_use]
    pub fn as_hugoniot(&self) -> Option<HugoniotRecord<'_>> {
        self.hugoniot.as_ref().map(|metadata| HugoniotRecord {
            record: self,
            metadata,
        })
    }
    /// Reference volume in the file's conventional-cell volume unit.
    #[must_use]
    pub fn reference_volume(&self) -> f64 {
        let model_volume = match self.eos {
            LoadedEos::Isothermal(model) => model.reference_volume(),
            LoadedEos::Thermal(model) => model.reference_volume(),
            LoadedEos::Hugoniot(model) => model.reference_volume(),
        };
        model_volume / self.volume_scale
    }

    /// Pressure in `GPa` at conventional-cell `volume` and `temperature`.
    ///
    /// Isothermal records ignore `temperature`. Hugoniot records require the
    /// coordinate to equal their precursor-state temperature.
    ///
    /// # Errors
    ///
    /// Returns an error when the state is invalid or model evaluation fails.
    pub fn pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = volume * self.volume_scale;
        match self.eos {
            LoadedEos::Isothermal(model) => model.pressure(volume),
            LoadedEos::Thermal(model) => model.pressure(volume, temperature),
            LoadedEos::Hugoniot(model) => {
                self.validate_hugoniot_temperature(temperature)?;
                self.as_hugoniot()
                    .ok_or(EosError::InvalidState {
                        name: "eos",
                        reason: "Hugoniot model is missing typed metadata",
                    })?
                    .validate_domain(model.particle_velocity(volume)?)?;
                model.pressure(volume)
            }
        }
    }

    /// Isothermal bulk modulus in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error when the state is invalid or model evaluation fails.
    pub fn bulk_modulus(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = volume * self.volume_scale;
        match self.eos {
            LoadedEos::Isothermal(model) => model.bulk_modulus(volume),
            LoadedEos::Thermal(model) => model.bulk_modulus(volume, temperature),
            LoadedEos::Hugoniot(_) => Err(EosError::InvalidState {
                name: "eos",
                reason: "must be an equilibrium isothermal or thermal EOS",
            }),
        }
    }

    /// Conventional-cell volume at a requested pressure and temperature.
    ///
    /// # Errors
    ///
    /// Returns an error when inversion fails on the supported branch.
    pub fn volume(&self, pressure: f64, temperature: f64) -> EosResult<f64> {
        let model_volume = match self.eos {
            LoadedEos::Isothermal(model) => model.volume(pressure)?,
            LoadedEos::Thermal(model) => model.volume(pressure, temperature)?,
            LoadedEos::Hugoniot(model) => {
                self.validate_hugoniot_temperature(temperature)?;
                let volume = model.volume(pressure)?;
                self.as_hugoniot()
                    .ok_or(EosError::InvalidState {
                        name: "eos",
                        reason: "Hugoniot model is missing typed metadata",
                    })?
                    .validate_domain(model.particle_velocity(volume)?)?;
                volume
            }
        };
        Ok(model_volume / self.volume_scale)
    }

    fn validate_hugoniot_temperature(&self, temperature: f64) -> EosResult<()> {
        if nearly_equal(temperature, self.reference_temperature) {
            Ok(())
        } else {
            Err(EosError::InvalidState {
                name: "temperature",
                reason: "must equal the Hugoniot initial-state temperature",
            })
        }
    }

    /// Thermal-pressure increase above the reference-temperature isotherm.
    ///
    /// The input volume uses the file's conventional-cell unit and the result
    /// is in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error for an isothermal record or invalid state.
    pub fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = volume * self.volume_scale;
        match self.eos {
            LoadedEos::Thermal(model) => model.thermal_pressure_increment(volume, temperature),
            LoadedEos::Isothermal(_) | LoadedEos::Hugoniot(_) => Err(EosError::InvalidState {
                name: "eos",
                reason: "must be thermal",
            }),
        }
    }

    /// Retained DAC pressure increment in `GPa` at a heated state.
    ///
    /// # Errors
    ///
    /// Returns an error for an isothermal record, invalid state, or invalid
    /// confinement fraction.
    pub fn dac_thermal_pressure(
        &self,
        volume: f64,
        temperature: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        let volume = volume * self.volume_scale;
        match self.eos {
            LoadedEos::Thermal(model) => model.dac_thermal_pressure(volume, temperature, f_dac),
            LoadedEos::Isothermal(_) | LoadedEos::Hugoniot(_) => Err(EosError::InvalidState {
                name: "eos",
                reason: "must be thermal",
            }),
        }
    }

    /// Heated conventional-cell volume from cold pressure and DAC confinement.
    ///
    /// # Errors
    ///
    /// Returns an error for an isothermal record, invalid state, invalid
    /// confinement fraction, or failed volume inversion.
    pub fn volume_with_dac_confinement(
        &self,
        cold_pressure: f64,
        temperature: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        let model_volume = match self.eos {
            LoadedEos::Thermal(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)?
            }
            LoadedEos::Isothermal(_) | LoadedEos::Hugoniot(_) => {
                return Err(EosError::InvalidState {
                    name: "eos",
                    reason: "must be thermal",
                });
            }
        };
        Ok(model_volume / self.volume_scale)
    }
}

/// A decoded `.eosmat` material and its executable records.
#[derive(Clone, Debug, PartialEq)]
pub struct Material {
    /// Stable material identifier.
    pub identifier: String,
    /// Human-readable material name.
    pub name: String,
    /// Chemical formula.
    pub formula: String,
    /// Represented material phase.
    pub phase: String,
    /// EOS records in source order.
    pub eos_records: Vec<EosRecord>,
    /// Original material document including all extension fields.
    pub document: Value,
}

impl Material {
    /// Find a record by its stable identifier.
    #[must_use]
    pub fn record(&self, identifier: &str) -> Option<&EosRecord> {
        self.eos_records
            .iter()
            .find(|record| record.identifier == identifier)
    }

    /// Return the preferred equilibrium record, or its first record as fallback.
    #[must_use]
    pub fn default_equilibrium_record(&self) -> Option<&EosRecord> {
        self.eos_records
            .iter()
            .filter(|record| !record.eos.is_hugoniot())
            .find(|record| record.is_default)
            .or_else(|| self.equilibrium_records().next())
    }

    /// Return the preferred Hugoniot record, or its first record as fallback.
    #[must_use]
    pub fn default_hugoniot_record(&self) -> Option<HugoniotRecord<'_>> {
        self.hugoniot_records()
            .find(|record| record.record.is_default)
            .or_else(|| self.hugoniot_records().next())
    }

    /// Return an equilibrium default first, preserving safe legacy behavior.
    #[must_use]
    pub fn default_record(&self) -> Option<&EosRecord> {
        self.default_equilibrium_record()
            .or_else(|| self.default_hugoniot_record().map(|record| record.record))
    }

    /// Iterate over equilibrium isothermal and thermal EOS records.
    pub fn equilibrium_records(&self) -> impl Iterator<Item = &EosRecord> {
        self.eos_records
            .iter()
            .filter(|record| !record.eos.is_hugoniot())
    }

    /// Iterate over phase-specific shock Hugoniot records.
    pub fn hugoniot_records(&self) -> impl Iterator<Item = HugoniotRecord<'_>> {
        self.eos_records.iter().filter_map(EosRecord::as_hugoniot)
    }

    /// Validate the retained document, including every executable EOS record.
    ///
    /// This is useful after editing [`Self::document`] directly.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid structure, metadata, or model parameters.
    pub fn validate(&self) -> Result<(), EosmatError> {
        validate_eosmat_document(&self.document)
    }

    /// Serialize the retained document as validated, pretty UTF-8 JSON.
    ///
    /// The output ends with a newline and preserves unknown extension fields.
    ///
    /// # Errors
    ///
    /// Returns an error when the retained document is no longer valid.
    pub fn to_json(&self) -> Result<String, EosmatError> {
        serialize_eosmat(&self.document)
    }

    /// Validate and save the retained document to an `.eosmat` file.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid content or file-writing failures.
    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), EosmatError> {
        save_eosmat(path, &self.document)
    }
}

#[derive(Debug, Deserialize)]
struct RawMaterial {
    format: Option<String>,
    format_version: Option<u64>,
    identifier: Option<String>,
    name: Option<String>,
    formula: Option<String>,
    phase: Option<String>,
    units: Option<RawUnits>,
    formula_units_per_cell: Option<f64>,
    #[serde(default)]
    eos_records: Vec<Value>,
}

#[derive(Debug, Deserialize)]
struct RawUnits {
    pressure: Option<String>,
    temperature: Option<String>,
    volume: Option<String>,
}

#[derive(Debug, Deserialize)]
struct RawRecord {
    identifier: Option<String>,
    label: Option<String>,
    #[serde(default)]
    default: bool,
    default_for: Option<String>,
    eos: Option<RawComponent>,
    thermal: Option<RawComponent>,
    equation_kind: Option<String>,
    loading_path: Option<String>,
    branch_kind: Option<String>,
    temperature_ref: Option<f64>,
    initial_state: Option<RawInitialState>,
    volume_basis: Option<RawVolumeBasis>,
    branch_domain: Option<RawBranchDomain>,
    volume: Option<RawVolume>,
}

#[derive(Debug, Deserialize)]
struct RawInitialState {
    phase: Option<String>,
    material_identifier: Option<String>,
    eos_record_identifier: Option<String>,
    temperature_k: Option<f64>,
    pressure_gpa: Option<f64>,
    density_g_cm3: Option<f64>,
}

#[derive(Debug, Deserialize)]
struct RawVolumeBasis {
    kind: Option<String>,
    formula_units: Option<f64>,
    molar_mass_g_mol: Option<f64>,
}

#[derive(Debug, Deserialize)]
struct RawBranchDomain {
    particle_velocity_km_s: Option<[f64; 2]>,
    kind: Option<String>,
    boundary_status: Option<String>,
    #[serde(default)]
    notes: Vec<String>,
}

#[derive(Debug, Deserialize)]
struct RawVolume {
    public_to_model_scale: Option<f64>,
}

#[derive(Debug, Deserialize)]
struct RawComponent {
    #[serde(rename = "type")]
    model_type: Option<String>,
    model: Option<String>,
    #[serde(default)]
    parameters: HashMap<String, Option<f64>>,
    debye_temperature_law: Option<String>,
    thermal_expansion_law: Option<String>,
    reference_volume_law: Option<String>,
    #[serde(default)]
    configuration: HashMap<String, Value>,
}

/// Load an `.eosmat` file and construct every EOS record through the built-in registry.
///
/// # Errors
///
/// Returns an error for I/O failures, malformed JSON, unsupported formats or
/// models, missing parameters, and invalid model parameters.
pub fn load_eosmat(path: impl AsRef<Path>) -> Result<Material, EosmatError> {
    let reader = BufReader::new(File::open(path)?);
    load_eosmat_reader(reader)
}

/// Decode an `.eosmat` document from a UTF-8 reader.
///
/// # Errors
///
/// Returns the same errors as [`load_eosmat`], except file-opening errors.
pub fn load_eosmat_reader(reader: impl Read) -> Result<Material, EosmatError> {
    let document = serde_json::from_reader(reader)?;
    material_from_value(document)
}

/// Decode an `.eosmat` document from a JSON string.
///
/// # Errors
///
/// Returns an error for malformed JSON or an invalid material document.
pub fn load_eosmat_str(source: &str) -> Result<Material, EosmatError> {
    let document = serde_json::from_str(source)?;
    material_from_value(document)
}

/// Validate a decoded `.eosmat` document without retaining an executable model.
///
/// Validation covers the canonical or legacy envelope, metadata used by the
/// Python API, unique record identifiers, covariance dimensions, ranges, and
/// construction of every referenced built-in EOS. Unknown extension fields
/// are accepted and preserved by subsequent serialization.
///
/// # Errors
///
/// Returns an error for invalid structure, unsupported models, or invalid
/// model parameters.
pub fn validate_eosmat_document(document: &Value) -> Result<(), EosmatError> {
    validate_document_structure(document)?;
    construct_material(document.clone()).map(|_| ())
}

/// Serialize a decoded `.eosmat` document as validated, pretty UTF-8 JSON.
///
/// The output ends with a newline. Unknown JSON extension fields are retained.
///
/// # Errors
///
/// Returns an error when the document is invalid or cannot be serialized.
pub fn serialize_eosmat(document: &Value) -> Result<String, EosmatError> {
    validate_eosmat_document(document)?;
    let mut serialized = serde_json::to_string_pretty(document)?;
    serialized.push('\n');
    Ok(serialized)
}

/// Validate and save a decoded `.eosmat` document.
///
/// Serialization completes before the destination is created, so an invalid
/// document cannot truncate an existing file.
///
/// # Errors
///
/// Returns an error for invalid content or file-writing failures.
pub fn save_eosmat(path: impl AsRef<Path>, document: &Value) -> Result<(), EosmatError> {
    let serialized = serialize_eosmat(document)?;
    let mut file = File::create(path)?;
    file.write_all(serialized.as_bytes())?;
    Ok(())
}

/// Construct an executable material from a decoded JSON value.
///
/// # Errors
///
/// Returns an error for an unsupported format or an invalid EOS record.
pub fn material_from_value(document: Value) -> Result<Material, EosmatError> {
    validate_document_structure(&document)?;
    construct_material(document)
}

fn construct_material(document: Value) -> Result<Material, EosmatError> {
    let raw: RawMaterial = serde_json::from_value(document.clone())?;
    let version = raw
        .format_version
        .ok_or_else(|| EosmatError::InvalidDocument("missing format_version".to_owned()))?;
    let canonical =
        raw.format.as_deref() == Some(EOSMAT_FORMAT) && version == EOSMAT_FORMAT_VERSION;
    let legacy = raw.format.is_none() && version == LEGACY_FORMAT_VERSION;
    if !canonical && !legacy {
        return Err(EosmatError::InvalidDocument(format!(
            "supported formats are {EOSMAT_FORMAT} version {EOSMAT_FORMAT_VERSION} and legacy Dioptas version {LEGACY_FORMAT_VERSION}"
        )));
    }

    let name = required_string(raw.name, "name")?;
    let formula = required_string(raw.formula, "formula")?;
    let phase = raw.phase.unwrap_or_else(|| "unspecified".to_owned());
    let identifier = match raw.identifier {
        Some(value) if !value.is_empty() => value,
        Some(_) if canonical => {
            return Err(EosmatError::InvalidDocument(
                "identifier must not be empty".to_owned(),
            ));
        }
        None if canonical => {
            return Err(EosmatError::InvalidDocument(
                "canonical format 3 requires identifier".to_owned(),
            ));
        }
        Some(_) | None => slug(&name),
    };
    if canonical {
        validate_units(raw.units.as_ref())?;
    }

    let mut eos_records = Vec::with_capacity(raw.eos_records.len());
    for (index, record_document) in raw.eos_records.into_iter().enumerate() {
        eos_records.push(record_from_value(
            record_document,
            index,
            raw.formula_units_per_cell,
            &identifier,
            &phase,
        )?);
    }

    Ok(Material {
        identifier,
        name,
        formula,
        phase,
        eos_records,
        document,
    })
}

fn invalid_document(reason: impl Into<String>) -> EosmatError {
    EosmatError::InvalidDocument(reason.into())
}

fn object_at<'a>(
    value: &'a Value,
    location: &str,
) -> Result<&'a serde_json::Map<String, Value>, EosmatError> {
    value
        .as_object()
        .ok_or_else(|| invalid_document(format!("{location} must be a JSON object")))
}

fn array_at<'a>(value: &'a Value, location: &str) -> Result<&'a [Value], EosmatError> {
    value
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| invalid_document(format!("{location} must be a JSON array")))
}

fn finite_at(value: &Value, location: &str, positive: bool) -> Result<f64, EosmatError> {
    let number = value
        .as_f64()
        .filter(|number| number.is_finite())
        .ok_or_else(|| invalid_document(format!("{location} must be a finite number")))?;
    if positive && number <= 0.0 {
        return Err(invalid_document(format!(
            "{location} must be greater than zero"
        )));
    }
    Ok(number)
}

fn validate_optional_string(
    document: &serde_json::Map<String, Value>,
    key: &str,
) -> Result<(), EosmatError> {
    if document.get(key).is_some_and(|value| !value.is_string()) {
        return Err(invalid_document(format!("{key} must be a string")));
    }
    Ok(())
}

fn validate_range(value: &Value, location: &str) -> Result<(), EosmatError> {
    let range = array_at(value, location)?;
    if range.len() != 2 {
        return Err(invalid_document(format!(
            "{location} must contain two values"
        )));
    }
    let low = finite_at(&range[0], &format!("{location}[0]"), false)?;
    let high = finite_at(&range[1], &format!("{location}[1]"), false)?;
    if low > high {
        return Err(invalid_document(format!("{location} must be ordered")));
    }
    Ok(())
}

fn validate_parameter_metadata(
    record: &serde_json::Map<String, Value>,
    location: &str,
    canonical: bool,
) -> Result<(), EosmatError> {
    match record.get("parameter_errors") {
        Some(value) => {
            for (name, error) in object_at(value, &format!("{location}.parameter_errors"))? {
                if !error.is_null() {
                    finite_at(error, &format!("{location}.parameter_errors.{name}"), false)?;
                }
            }
        }
        None if canonical => {
            return Err(invalid_document(format!(
                "{location}.parameter_errors is required"
            )));
        }
        None => {}
    }
    match record.get("fixed_parameters") {
        Some(value) => {
            if array_at(value, &format!("{location}.fixed_parameters"))?
                .iter()
                .any(|name| !name.is_string())
            {
                return Err(invalid_document(format!(
                    "{location}.fixed_parameters must contain strings"
                )));
            }
        }
        None if canonical => {
            return Err(invalid_document(format!(
                "{location}.fixed_parameters is required"
            )));
        }
        None => {}
    }
    if let Some(value) = record
        .get("parameter_error_confidence")
        .filter(|value| !value.is_null())
    {
        let confidence = finite_at(
            value,
            &format!("{location}.parameter_error_confidence"),
            false,
        )?;
        if !(0.0..1.0).contains(&confidence) || confidence == 0.0 {
            return Err(invalid_document(format!(
                "{location}.parameter_error_confidence must lie between zero and one"
            )));
        }
    }
    Ok(())
}

fn validate_covariance(value: &Value, location: &str) -> Result<(), EosmatError> {
    let covariance = object_at(value, location)?;
    let matrix = array_at(
        covariance
            .get("matrix")
            .ok_or_else(|| invalid_document(format!("{location}.matrix is required")))?,
        &format!("{location}.matrix"),
    )?;
    let order = array_at(
        covariance
            .get("parameter_order")
            .ok_or_else(|| invalid_document(format!("{location}.parameter_order is required")))?,
        &format!("{location}.parameter_order"),
    )?;
    if order.is_empty()
        || order
            .iter()
            .any(|name| name.as_str().is_none_or(str::is_empty))
        || matrix.len() != order.len()
    {
        return Err(invalid_document(format!(
            "{location} must contain a non-empty square matrix matching parameter_order"
        )));
    }
    for (row_index, row) in matrix.iter().enumerate() {
        let row = array_at(row, &format!("{location}.matrix[{row_index}]"))?;
        if row.len() != order.len() {
            return Err(invalid_document(format!(
                "{location}.matrix must be square and match parameter_order"
            )));
        }
        for (column_index, value) in row.iter().enumerate() {
            finite_at(
                value,
                &format!("{location}.matrix[{row_index}][{column_index}]"),
                false,
            )?;
        }
    }
    Ok(())
}

fn validate_pressure_reference_links(
    method: &serde_json::Map<String, Value>,
    kind: Option<&str>,
    location: &str,
) -> Result<(), EosmatError> {
    if let Some(identifier) = method.get("reference_eos_record") {
        if identifier.as_str().is_none_or(str::is_empty) {
            return Err(invalid_document(format!(
                "{location}.reference_eos_record must be a non-empty string"
            )));
        }
        if kind != Some("equation_of_state") {
            return Err(invalid_document(format!(
                "{location}.reference_eos_record requires an equation_of_state method"
            )));
        }
    }
    if let Some(identifier) = method.get("reference_calibration_record") {
        if identifier.as_str().is_none_or(str::is_empty) {
            return Err(invalid_document(format!(
                "{location}.reference_calibration_record must be a non-empty string"
            )));
        }
        if kind != Some("ruby_fluorescence") {
            return Err(invalid_document(format!(
                "{location}.reference_calibration_record requires a ruby_fluorescence method"
            )));
        }
    }
    Ok(())
}

fn validate_pressure_calibration(value: &Value, location: &str) -> Result<(), EosmatError> {
    let calibration = object_at(value, location)?;
    let status = calibration.get("status").and_then(Value::as_str);
    if !matches!(
        status,
        Some("resolved" | "partially_resolved" | "not_applicable" | "unresolved")
    ) {
        return Err(invalid_document(format!("{location}.status is invalid")));
    }
    let methods = array_at(
        calibration
            .get("methods")
            .ok_or_else(|| invalid_document(format!("{location}.methods is required")))?,
        &format!("{location}.methods"),
    )?;
    if matches!(
        status,
        Some("resolved" | "partially_resolved" | "not_applicable")
    ) && methods.is_empty()
    {
        return Err(invalid_document(format!(
            "{location}.methods must not be empty for status {status:?}"
        )));
    }
    for (index, value) in methods.iter().enumerate() {
        let method_location = format!("{location}.methods[{index}]");
        let method = object_at(value, &method_location)?;
        let kind = method.get("kind").and_then(Value::as_str);
        if !matches!(
            kind,
            Some(
                "equation_of_state"
                    | "ruby_fluorescence"
                    | "other_optical_gauge"
                    | "shock_wave"
                    | "ultrasonic"
                    | "ab_initio"
                    | "self_consistent"
                    | "ambient_pressure"
                    | "other"
            )
        ) {
            return Err(invalid_document(format!(
                "{method_location}.kind is invalid"
            )));
        }
        if method
            .get("source_location")
            .and_then(Value::as_str)
            .is_none_or(str::is_empty)
        {
            return Err(invalid_document(format!(
                "{method_location}.source_location must be a non-empty string"
            )));
        }
        if let Some(reference) = method.get("reference") {
            if !reference.is_string() && !reference.is_object() {
                return Err(invalid_document(format!(
                    "{method_location}.reference must be a string or object"
                )));
            }
        } else if kind == Some("equation_of_state") {
            return Err(invalid_document(format!(
                "{method_location}.reference is required for an EOS method"
            )));
        }
        validate_pressure_reference_links(method, kind, &method_location)?;
    }
    let recalculation = object_at(
        calibration
            .get("recalculation")
            .ok_or_else(|| invalid_document(format!("{location}.recalculation is required")))?,
        &format!("{location}.recalculation"),
    )?;
    if !matches!(
        recalculation.get("status").and_then(Value::as_str),
        Some(
            "ready"
                | "missing_calibrant_observations"
                | "reference_eos_not_bundled"
                | "reference_model_not_supported"
                | "not_applicable"
                | "not_possible"
        )
    ) {
        return Err(invalid_document(format!(
            "{location}.recalculation.status is invalid"
        )));
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn validate_document_structure(document: &Value) -> Result<(), EosmatError> {
    let document = object_at(document, "document")?;
    let format = document.get("format").and_then(Value::as_str);
    let version = document.get("format_version").and_then(Value::as_u64);
    let canonical = format == Some(EOSMAT_FORMAT) && version == Some(EOSMAT_FORMAT_VERSION);
    let legacy = document.get("format").is_none() && version == Some(LEGACY_FORMAT_VERSION);
    if !canonical && !legacy {
        return Err(invalid_document(format!(
            "supported formats are {EOSMAT_FORMAT} version {EOSMAT_FORMAT_VERSION} and legacy Dioptas version {LEGACY_FORMAT_VERSION}"
        )));
    }

    for key in ["name", "formula"] {
        if document.get(key).and_then(Value::as_str).is_none() {
            return Err(invalid_document(format!("{key} must be a string")));
        }
    }
    for key in ["identifier", "phase", "symmetry", "notes"] {
        validate_optional_string(document, key)?;
    }
    if document
        .get("cell_contents")
        .is_some_and(|value| !value.is_null() && !value.is_string())
    {
        return Err(invalid_document("cell_contents must be a string"));
    }
    if canonical
        && document
            .get("identifier")
            .and_then(Value::as_str)
            .is_none_or(str::is_empty)
    {
        return Err(invalid_document(
            "canonical format 3 requires a non-empty identifier",
        ));
    }
    if let Some(value) = document
        .get("formula_units_per_cell")
        .filter(|value| !value.is_null())
    {
        finite_at(value, "formula_units_per_cell", true)?;
    }
    for key in ["aliases", "atom_sites", "peaks"] {
        if let Some(value) = document.get(key) {
            array_at(value, key)?;
        }
    }
    let mut dataset_identifiers = std::collections::HashSet::new();
    let mut dataset_record_links = Vec::new();
    if let Some(value) = document.get("datasets") {
        for (index, value) in array_at(value, "datasets")?.iter().enumerate() {
            let location = format!("datasets[{index}]");
            let dataset = object_at(value, &location)?;
            let identifier = dataset
                .get("identifier")
                .and_then(Value::as_str)
                .filter(|value| !value.is_empty())
                .ok_or_else(|| {
                    invalid_document(format!("{location}.identifier must be a non-empty string"))
                })?;
            if !dataset_identifiers.insert(identifier) {
                return Err(invalid_document(format!(
                    "duplicate dataset identifier {identifier:?}"
                )));
            }
            for key in ["kind", "source_location"] {
                if dataset
                    .get(key)
                    .and_then(Value::as_str)
                    .is_none_or(str::is_empty)
                {
                    return Err(invalid_document(format!(
                        "{location}.{key} must be a non-empty string"
                    )));
                }
            }
            if !dataset
                .get("reference")
                .is_some_and(|value| value.is_string() || value.is_object())
            {
                return Err(invalid_document(format!(
                    "{location}.reference must be a string or object"
                )));
            }
            let columns = array_at(
                dataset
                    .get("columns")
                    .ok_or_else(|| invalid_document(format!("{location}.columns is required")))?,
                &format!("{location}.columns"),
            )?;
            if columns.is_empty() {
                return Err(invalid_document(format!(
                    "{location}.columns must be non-empty"
                )));
            }
            let mut column_names = std::collections::HashSet::new();
            for (column_index, value) in columns.iter().enumerate() {
                let column_location = format!("{location}.columns[{column_index}]");
                let column = object_at(value, &column_location)?;
                let name = column
                    .get("name")
                    .and_then(Value::as_str)
                    .filter(|value| !value.is_empty())
                    .ok_or_else(|| {
                        invalid_document(format!(
                            "{column_location}.name must be a non-empty string"
                        ))
                    })?;
                if !column_names.insert(name) {
                    return Err(invalid_document(format!(
                        "duplicate column name {name:?} in {location}"
                    )));
                }
                for key in ["quantity", "unit", "role"] {
                    if column
                        .get(key)
                        .and_then(Value::as_str)
                        .is_none_or(str::is_empty)
                    {
                        return Err(invalid_document(format!(
                            "{column_location}.{key} must be a non-empty string"
                        )));
                    }
                }
                if !matches!(
                    column.get("role").and_then(Value::as_str),
                    Some(
                        "value"
                            | "uncertainty"
                            | "standard_deviation"
                            | "standard_error"
                            | "bound"
                            | "flag"
                    )
                ) {
                    return Err(invalid_document(format!(
                        "{column_location}.role is invalid"
                    )));
                }
                if column
                    .get("of")
                    .and_then(Value::as_str)
                    .is_some_and(|name| !column_names.contains(name))
                {
                    return Err(invalid_document(format!(
                        "{column_location}.of must reference an earlier column"
                    )));
                }
            }
            match (dataset.get("rows"), dataset.get("resource")) {
                (Some(rows), None) => {
                    for (row_index, row) in array_at(rows, &format!("{location}.rows"))?
                        .iter()
                        .enumerate()
                    {
                        let row = array_at(row, &format!("{location}.rows[{row_index}]"))?;
                        if row.len() != columns.len() {
                            return Err(invalid_document(format!(
                                "{location}.rows[{row_index}] must match the column count"
                            )));
                        }
                        for (column_index, value) in row.iter().enumerate() {
                            if !value.is_null() {
                                finite_at(
                                    value,
                                    &format!("{location}.rows[{row_index}][{column_index}]"),
                                    false,
                                )?;
                            }
                        }
                    }
                }
                (None, Some(resource)) => {
                    let resource = object_at(resource, &format!("{location}.resource"))?;
                    for key in ["path", "sha256", "media_type"] {
                        if resource
                            .get(key)
                            .and_then(Value::as_str)
                            .is_none_or(str::is_empty)
                        {
                            return Err(invalid_document(format!(
                                "{location}.resource.{key} must be a non-empty string"
                            )));
                        }
                    }
                    let path = resource["path"].as_str().unwrap();
                    if Path::new(path).is_absolute()
                        || Path::new(path)
                            .components()
                            .any(|component| matches!(component, std::path::Component::ParentDir))
                    {
                        return Err(invalid_document(format!(
                            "{location}.resource.path must be relative and local"
                        )));
                    }
                    let sha256 = resource["sha256"].as_str().unwrap();
                    if sha256.len() != 64
                        || !sha256
                            .bytes()
                            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
                    {
                        return Err(invalid_document(format!(
                            "{location}.resource.sha256 must be 64 lowercase hex characters"
                        )));
                    }
                }
                _ => {
                    return Err(invalid_document(format!(
                        "{location} must contain exactly one of rows or resource"
                    )));
                }
            }
            let used_by = array_at(
                dataset.get("used_by_eos_records").ok_or_else(|| {
                    invalid_document(format!("{location}.used_by_eos_records is required"))
                })?,
                &format!("{location}.used_by_eos_records"),
            )?;
            if used_by.is_empty() {
                return Err(invalid_document(format!(
                    "{location}.used_by_eos_records must be non-empty"
                )));
            }
            for value in used_by {
                let record_identifier = value
                    .as_str()
                    .filter(|value| !value.is_empty())
                    .ok_or_else(|| {
                        invalid_document(format!(
                            "{location}.used_by_eos_records must contain non-empty strings"
                        ))
                    })?;
                dataset_record_links.push((identifier, record_identifier));
            }
        }
    }
    if let Some(value) = document.get("lattice") {
        let lattice = object_at(value, "lattice")?;
        for key in ["a", "alpha", "beta", "gamma"] {
            finite_at(
                lattice
                    .get(key)
                    .ok_or_else(|| invalid_document(format!("lattice.{key} is required")))?,
                &format!("lattice.{key}"),
                false,
            )?;
        }
        for key in ["b", "c"] {
            if let Some(value) = lattice.get(key).filter(|value| !value.is_null()) {
                finite_at(value, &format!("lattice.{key}"), false)?;
            }
        }
    }

    let empty_records = Vec::new();
    let records = match document.get("eos_records") {
        Some(value) => array_at(value, "eos_records")?,
        None => empty_records.as_slice(),
    };
    let mut identifiers = std::collections::HashSet::new();
    let mut derived_record_links = Vec::new();
    let mut fit_dataset_links = Vec::new();
    let mut default_counts = [0_u32; 2];
    for (index, value) in records.iter().enumerate() {
        let location = format!("eos_records[{index}]");
        let record = object_at(value, &location)?;
        if canonical && record.get("label").and_then(Value::as_str).is_none() {
            return Err(invalid_document(format!(
                "{location}.label must be a string"
            )));
        }
        if let Some(identifier) = record.get("identifier") {
            let identifier = identifier
                .as_str()
                .filter(|value| !value.is_empty())
                .ok_or_else(|| {
                    invalid_document(format!("{location}.identifier must be a non-empty string"))
                })?;
            if !identifiers.insert(identifier) {
                return Err(invalid_document(format!(
                    "duplicate EOS record identifier {identifier:?}"
                )));
            }
        } else if canonical {
            return Err(invalid_document(format!(
                "{location}.identifier is required"
            )));
        }
        let is_hugoniot = record
            .get("eos")
            .and_then(Value::as_object)
            .and_then(|eos| eos.get("type"))
            .and_then(Value::as_str)
            == Some("LinearUsUpHugoniot");
        let default_category = if is_hugoniot {
            "hugoniot"
        } else {
            "equilibrium"
        };
        if let Some(default_for) = record.get("default_for") {
            if default_for.as_str() != Some(default_category) {
                return Err(invalid_document(format!(
                    "{location}.default_for does not match equation_kind"
                )));
            }
        }
        if record.get("default") == Some(&Value::Bool(true))
            || record.get("default_for").and_then(Value::as_str) == Some(default_category)
        {
            default_counts[usize::from(is_hugoniot)] += 1;
        }

        let record_identifier = record.get("identifier").and_then(Value::as_str);
        let record_kind = record
            .get("record_kind")
            .and_then(Value::as_str)
            .unwrap_or("published");
        if !matches!(
            record_kind,
            "published" | "refit" | "derived" | "diagnostic"
        ) {
            return Err(invalid_document(format!(
                "{location}.record_kind is invalid"
            )));
        }
        let derived_from = record.get("derived_from_record");
        if let Some(value) = derived_from {
            let parent_identifier = value
                .as_str()
                .filter(|value| !value.is_empty())
                .ok_or_else(|| {
                    invalid_document(format!(
                        "{location}.derived_from_record must be a non-empty string"
                    ))
                })?;
            if let Some(record_identifier) = record_identifier {
                derived_record_links.push((record_identifier, parent_identifier));
            }
        }
        let fit_provenance = record.get("fit_provenance");
        if record_kind == "refit"
            && (derived_from.is_none() || !fit_provenance.is_some_and(Value::is_object))
        {
            return Err(invalid_document(format!(
                "{location} refit records require derived_from_record and fit_provenance"
            )));
        }
        if record_kind == "derived" {
            let derivation = object_at(
                record.get("derivation").ok_or_else(|| {
                    invalid_document(format!("{location}.derivation is required"))
                })?,
                &format!("{location}.derivation"),
            )?;
            if !matches!(
                derivation.get("source_kind").and_then(Value::as_str),
                Some("sesame_table" | "published_table" | "calculation")
            ) {
                return Err(invalid_document(format!(
                    "{location}.derivation.source_kind is invalid"
                )));
            }
            for key in ["source_identifier", "method"] {
                if derivation
                    .get(key)
                    .and_then(Value::as_str)
                    .is_none_or(str::is_empty)
                {
                    return Err(invalid_document(format!(
                        "{location}.derivation.{key} must be a non-empty string"
                    )));
                }
            }
        }
        if let Some(value) = fit_provenance {
            let provenance = object_at(value, &format!("{location}.fit_provenance"))?;
            let dataset_identifier = provenance
                .get("dataset")
                .and_then(Value::as_str)
                .filter(|value| !value.is_empty())
                .ok_or_else(|| {
                    invalid_document(format!(
                        "{location}.fit_provenance.dataset must be a non-empty string"
                    ))
                })?;
            if let Some(record_identifier) = record_identifier {
                fit_dataset_links.push((record_identifier, dataset_identifier));
            }
        }

        if let Some(reference) = record.get("reference") {
            if let Some(reference) = reference.as_object() {
                if !reference.get("authors").is_some_and(Value::is_array)
                    || !reference.get("year").is_some_and(Value::is_i64)
                {
                    return Err(invalid_document(format!(
                        "{location}.reference requires authors and an integer year"
                    )));
                }
            } else if !reference.is_string() {
                return Err(invalid_document(format!(
                    "{location}.reference must be a string or object"
                )));
            }
        } else if canonical {
            return Err(invalid_document(format!(
                "{location}.reference is required"
            )));
        }

        validate_parameter_metadata(record, &location, canonical)?;
        if let Some(value) = record.get("volume").filter(|value| !value.is_null()) {
            let volume = object_at(value, &format!("{location}.volume"))?;
            for key in ["reference_value", "public_to_model_scale"] {
                if let Some(value) = volume.get(key).filter(|value| !value.is_null()) {
                    finite_at(value, &format!("{location}.volume.{key}"), true)?;
                }
            }
            if volume
                .get("model_unit")
                .is_some_and(|value| !value.is_null() && !value.is_string())
            {
                return Err(invalid_document(format!(
                    "{location}.volume.model_unit must be a string"
                )));
            }
        }
        if let Some(value) = record
            .get("parameter_covariance")
            .filter(|value| !value.is_null())
        {
            validate_covariance(value, &format!("{location}.parameter_covariance"))?;
        }
        if let Some(thermal) = record.get("thermal").filter(|value| !value.is_null()) {
            validate_parameter_metadata(
                object_at(thermal, &format!("{location}.thermal"))?,
                &format!("{location}.thermal"),
                false,
            )?;
        }
        for key in [
            "experimental_pressure_range_gpa",
            "experimental_temperature_range_k",
        ] {
            if let Some(value) = record.get(key).filter(|value| !value.is_null()) {
                validate_range(value, &format!("{location}.{key}"))?;
            }
        }
        if let Some(value) = record.get("validity").filter(|value| !value.is_null()) {
            let validity = object_at(value, &format!("{location}.validity"))?;
            for key in ["pressure_gpa", "temperature_k", "volume_ratio"] {
                if let Some(value) = validity.get(key) {
                    validate_range(value, &format!("{location}.validity.{key}"))?;
                }
            }
            if validity.get("notes").is_some_and(|value| !value.is_array()) {
                return Err(invalid_document(format!(
                    "{location}.validity.notes must be an array"
                )));
            }
        }
        if let Some(value) = record
            .get("pressure_calibration")
            .filter(|value| !value.is_null())
        {
            validate_pressure_calibration(value, &format!("{location}.pressure_calibration"))?;
        }
        if canonical {
            let validation = object_at(
                record.get("scientific_validation").ok_or_else(|| {
                    invalid_document(format!("{location}.scientific_validation is required"))
                })?,
                &format!("{location}.scientific_validation"),
            )?;
            if !matches!(
                validation.get("status").and_then(Value::as_str),
                Some("primary_source_validated" | "pending_primary_source_check" | "deferred")
            ) {
                return Err(invalid_document(format!(
                    "{location}.scientific_validation.status is invalid"
                )));
            }
        }
    }
    for (category, count) in ["equilibrium", "hugoniot"].into_iter().zip(default_counts) {
        if count > 1 {
            return Err(invalid_document(format!(
                "a material may have at most one default {category} EOS record"
            )));
        }
    }
    for (dataset_identifier, record_identifier) in dataset_record_links {
        if !identifiers.contains(record_identifier) {
            return Err(invalid_document(format!(
                "dataset {dataset_identifier:?} references unknown EOS record {record_identifier:?}"
            )));
        }
    }
    for (record_identifier, parent_identifier) in derived_record_links {
        if parent_identifier == record_identifier {
            return Err(invalid_document(format!(
                "EOS record {record_identifier:?} cannot derive from itself"
            )));
        }
        if !identifiers.contains(parent_identifier) {
            return Err(invalid_document(format!(
                "EOS record {record_identifier:?} derives from unknown EOS record {parent_identifier:?}"
            )));
        }
    }
    for (record_identifier, dataset_identifier) in fit_dataset_links {
        if !dataset_identifiers.contains(dataset_identifier) {
            return Err(invalid_document(format!(
                "EOS record {record_identifier:?} fit provenance references unknown dataset {dataset_identifier:?}"
            )));
        }
    }
    Ok(())
}

fn required_string(value: Option<String>, field: &str) -> Result<String, EosmatError> {
    match value {
        Some(value) => Ok(value),
        None => Err(EosmatError::InvalidDocument(format!("missing {field}"))),
    }
}

fn validate_units(units: Option<&RawUnits>) -> Result<(), EosmatError> {
    let valid = units.is_some_and(|units| {
        units.pressure.as_deref() == Some("GPa")
            && units.temperature.as_deref() == Some("K")
            && units.volume.as_deref() == Some("angstrom^3/conventional_unit_cell")
    });
    if valid {
        Ok(())
    } else {
        Err(EosmatError::InvalidDocument(
            "canonical units must be GPa, K, and angstrom^3/conventional_unit_cell".to_owned(),
        ))
    }
}

fn slug(value: &str) -> String {
    let mut result = String::new();
    let mut separator = false;
    for character in value.chars().flat_map(char::to_lowercase) {
        if character.is_ascii_alphanumeric() {
            if separator && !result.is_empty() {
                result.push('_');
            }
            result.push(character);
            separator = false;
        } else {
            separator = true;
        }
    }
    if result.is_empty() {
        "unnamed".to_owned()
    } else {
        result
    }
}

fn record_from_value(
    document: Value,
    index: usize,
    formula_units_per_cell: Option<f64>,
    material_identifier: &str,
    represented_phase: &str,
) -> Result<EosRecord, EosmatError> {
    let raw: RawRecord =
        serde_json::from_value(document.clone()).map_err(|error| EosmatError::InvalidRecord {
            identifier: format!("record_{}", index + 1),
            reason: error.to_string(),
        })?;
    let label = raw
        .label
        .clone()
        .unwrap_or_else(|| format!("EOS record {}", index + 1));
    let identifier = raw.identifier.clone().unwrap_or_else(|| slug(&label));
    let result = build_record(
        &raw,
        identifier.clone(),
        label,
        document,
        formula_units_per_cell,
        material_identifier,
        represented_phase,
    );
    result.map_err(|reason| EosmatError::InvalidRecord { identifier, reason })
}

#[allow(clippy::too_many_lines)]
fn build_record(
    raw: &RawRecord,
    identifier: String,
    label: String,
    document: Value,
    formula_units_per_cell: Option<f64>,
    material_identifier: &str,
    represented_phase: &str,
) -> Result<EosRecord, String> {
    let reference_component = raw.eos.as_ref().ok_or("missing eos component")?;
    let reference_identifier = component_model_identifier(reference_component, false)?;
    let is_hugoniot = reference_identifier == "linear_us_up_hugoniot";
    let thermal_identifier = raw
        .thermal
        .as_ref()
        .map(|component| component_model_identifier(component, true))
        .transpose()?;
    let volume_scale = if let Some(scale) = raw
        .volume
        .as_ref()
        .and_then(|volume| volume.public_to_model_scale)
    {
        if !scale.is_finite() || scale <= 0.0 {
            return Err("volume.public_to_model_scale must be positive and finite".to_owned());
        }
        scale
    } else if thermal_identifier.is_some_and(is_molar_volume_model) {
        let formula_units = formula_units_per_cell.ok_or_else(|| {
            "formula_units_per_cell or volume.public_to_model_scale is required for molar-volume thermal EOS".to_owned()
        })?;
        if !formula_units.is_finite() || formula_units <= 0.0 {
            return Err("formula_units_per_cell must be positive and finite".to_owned());
        }
        CELL_ANGSTROM3_TO_FORMULA_MOLAR_J_PER_BAR / formula_units
    } else {
        1.0
    };

    let mut hugoniot_metadata = None;
    let eos = if is_hugoniot {
        if raw.equation_kind.as_deref() != Some("hugoniot") {
            return Err("a Hugoniot model requires equation_kind \"hugoniot\"".to_owned());
        }
        if raw.thermal.is_some() {
            return Err("a Hugoniot record cannot contain a thermal component".to_owned());
        }
        let loading_path = match raw.loading_path.as_deref() {
            Some("principal") => HugoniotLoadingPath::Principal,
            Some("precompressed") => HugoniotLoadingPath::Precompressed,
            _ => return Err("a Hugoniot record requires a supported loading_path".to_owned()),
        };
        let branch_kind = match raw.branch_kind.as_deref() {
            Some("untransformed") => HugoniotBranchKind::Untransformed,
            Some("transformed") => HugoniotBranchKind::Transformed,
            _ => return Err("a Hugoniot record requires a supported branch_kind".to_owned()),
        };
        let initial_state = raw
            .initial_state
            .as_ref()
            .ok_or("a Hugoniot record requires initial_state")?;
        let initial_phase = initial_state
            .phase
            .as_deref()
            .filter(|phase| !phase.trim().is_empty())
            .ok_or("initial_state.phase must be a non-empty string")?;
        let initial_material = initial_state
            .material_identifier
            .as_deref()
            .filter(|identifier| !identifier.is_empty())
            .ok_or("initial_state.material_identifier must be a non-empty string")?;
        if initial_state
            .eos_record_identifier
            .as_deref()
            .is_some_and(str::is_empty)
        {
            return Err("initial_state.eos_record_identifier must not be empty".to_owned());
        }
        if branch_kind == HugoniotBranchKind::Untransformed
            && (initial_phase != represented_phase || initial_material != material_identifier)
        {
            return Err(
                "an untransformed branch must reference the represented material and phase"
                    .to_owned(),
            );
        }
        if branch_kind == HugoniotBranchKind::Transformed
            && (initial_phase == represented_phase || initial_material == material_identifier)
        {
            return Err("a transformed branch requires a distinct precursor".to_owned());
        }
        let temperature = initial_state
            .temperature_k
            .ok_or("initial_state.temperature_k is required")?;
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err("initial_state.temperature_k must be positive and finite".to_owned());
        }
        let initial_pressure = initial_state
            .pressure_gpa
            .ok_or("initial_state.pressure_gpa is required")?;
        if !initial_pressure.is_finite() {
            return Err("initial_state.pressure_gpa must be finite".to_owned());
        }
        let initial_density = initial_state
            .density_g_cm3
            .ok_or("initial_state.density_g_cm3 is required")?;
        if !initial_density.is_finite() || initial_density <= 0.0 {
            return Err("initial_state.density_g_cm3 must be positive and finite".to_owned());
        }
        let volume_basis = raw
            .volume_basis
            .as_ref()
            .ok_or("a Hugoniot record requires volume_basis")?;
        if volume_basis.kind.as_deref() != Some("formula_units") {
            return Err("volume_basis.kind must be \"formula_units\"".to_owned());
        }
        let basis_formula_units = volume_basis
            .formula_units
            .ok_or("volume_basis.formula_units is required")?;
        if !basis_formula_units.is_finite() || basis_formula_units <= 0.0 {
            return Err("volume_basis.formula_units must be positive and finite".to_owned());
        }
        let material_formula_units = formula_units_per_cell
            .ok_or("formula_units_per_cell is required when a material has a Hugoniot record")?;
        if !nearly_equal(basis_formula_units, material_formula_units) {
            return Err("volume_basis.formula_units must match formula_units_per_cell".to_owned());
        }
        let molar_mass = volume_basis
            .molar_mass_g_mol
            .ok_or("volume_basis.molar_mass_g_mol is required")?;
        if !molar_mass.is_finite() || molar_mass <= 0.0 {
            return Err("volume_basis.molar_mass_g_mol must be positive and finite".to_owned());
        }
        if raw
            .temperature_ref
            .is_some_and(|reference| !nearly_equal(reference, temperature))
        {
            return Err(
                "temperature_ref must match initial_state.temperature_k for a Hugoniot record"
                    .to_owned(),
            );
        }
        let model = build_hugoniot(reference_component, volume_scale)?;
        if !nearly_equal(initial_pressure, model.initial_pressure()) {
            return Err("initial_state.pressure_gpa must match eos.parameters.P0".to_owned());
        }
        if !nearly_equal(initial_density, model.initial_density()) {
            return Err("initial_state.density_g_cm3 must match eos.parameters.rho0".to_owned());
        }
        let public_v0 = model.reference_volume() / volume_scale;
        let expected_density =
            basis_formula_units * molar_mass / (6.022_140_76e23 * public_v0 * 1.0e-24);
        if (initial_density - expected_density).abs()
            > HUGONIOT_MASS_BASIS_RELATIVE_TOLERANCE
                * initial_density.abs().max(expected_density.abs())
        {
            return Err(
                "V0, rho0, formula_units, and molar_mass_g_mol must share a mass basis".to_owned(),
            );
        }
        if loading_path == HugoniotLoadingPath::Precompressed && initial_pressure <= 0.0 {
            return Err("a precompressed Hugoniot requires positive P0".to_owned());
        }
        let domain = raw
            .branch_domain
            .as_ref()
            .ok_or("a Hugoniot record requires branch_domain")?;
        let velocity_range = domain
            .particle_velocity_km_s
            .ok_or("branch_domain.particle_velocity_km_s is required")?;
        if velocity_range[0] < 0.0
            || !velocity_range[0].is_finite()
            || !velocity_range[1].is_finite()
            || velocity_range[0] > velocity_range[1]
        {
            return Err(
                "branch-domain particle velocities must be finite, non-negative, and ordered"
                    .to_owned(),
            );
        }
        let domain_kind = match domain.kind.as_deref() {
            Some("phase_stability") => HugoniotDomainKind::PhaseStability,
            Some("experimental_coverage") => HugoniotDomainKind::ExperimentalCoverage,
            Some("recommended") => HugoniotDomainKind::Recommended,
            _ => return Err("branch_domain.kind is invalid".to_owned()),
        };
        let boundary_status = match domain.boundary_status.as_deref() {
            Some("reported_exactly") => HugoniotBoundaryStatus::ReportedExactly,
            Some("reported_qualitatively") => HugoniotBoundaryStatus::ReportedQualitatively,
            Some("inferred") => HugoniotBoundaryStatus::Inferred,
            _ => return Err("branch_domain.boundary_status is invalid".to_owned()),
        };
        hugoniot_metadata = Some(HugoniotRecordMetadata {
            loading_path,
            branch_kind,
            initial_state: HugoniotInitialState {
                phase: initial_phase.to_owned(),
                material_identifier: initial_material.to_owned(),
                eos_record_identifier: initial_state.eos_record_identifier.clone(),
                temperature_k: temperature,
                pressure_gpa: initial_pressure,
                density_g_cm3: initial_density,
            },
            volume_basis: HugoniotVolumeBasis {
                formula_units: basis_formula_units,
                molar_mass_g_mol: molar_mass,
            },
            branch_domain: HugoniotBranchDomain {
                particle_velocity_km_s: velocity_range,
                kind: domain_kind,
                boundary_status,
                notes: domain.notes.clone(),
            },
        });
        LoadedEos::Hugoniot(model)
    } else {
        let expected_kind = if raw.thermal.is_some() {
            "thermal"
        } else {
            "isothermal"
        };
        if raw
            .equation_kind
            .as_deref()
            .is_some_and(|kind| kind != expected_kind)
        {
            return Err(format!(
                "equation_kind must be {expected_kind:?} for this record"
            ));
        }
        let reference = build_isothermal(reference_component, volume_scale)?;
        match raw.thermal.as_ref() {
            None => LoadedEos::Isothermal(reference),
            Some(thermal) => LoadedEos::Thermal(build_thermal(thermal, reference)?),
        }
    };
    let reference_temperature = raw.temperature_ref.unwrap_or_else(|| match eos {
        LoadedEos::Isothermal(_) => 300.0,
        LoadedEos::Thermal(model) => model.reference_temperature(),
        LoadedEos::Hugoniot(_) => raw
            .initial_state
            .as_ref()
            .and_then(|state| state.temperature_k)
            .unwrap_or(300.0),
    });
    if !reference_temperature.is_finite() || reference_temperature <= 0.0 {
        return Err("temperature_ref must be positive and finite".to_owned());
    }
    let default_category = if is_hugoniot {
        "hugoniot"
    } else {
        "equilibrium"
    };
    if raw
        .default_for
        .as_deref()
        .is_some_and(|category| category != default_category)
    {
        return Err("default_for does not match equation_kind".to_owned());
    }

    Ok(EosRecord {
        identifier,
        label,
        is_default: raw.default || raw.default_for.as_deref() == Some(default_category),
        eos,
        hugoniot: hugoniot_metadata,
        reference_temperature,
        document,
        volume_scale,
    })
}

fn nearly_equal(left: f64, right: f64) -> bool {
    (left - right).abs() <= 1.0e-12 * left.abs().max(right.abs()).max(1.0)
}

fn is_molar_volume_model(model: &str) -> bool {
    matches!(
        model,
        "mie_gruneisen_debye"
            | "mie_gruneisen_einstein"
            | "asymptotic_power_law_mie_gruneisen_debye"
            | "double_debye_helmholtz"
            | "double_debye_log_moment_helmholtz"
            | "dorogokupets_oganov_2007"
            | "multi_oscillator_gruneisen_thermal_pressure"
            | "thermal_modified_tait"
    )
}

fn parameter(component: &RawComponent, name: &str) -> Result<f64, String> {
    component
        .parameters
        .get(name)
        .copied()
        .flatten()
        .ok_or_else(|| format!("model {:?} is missing parameter {name}", component.model))
}

fn component_model_identifier(component: &RawComponent, thermal: bool) -> Result<&str, String> {
    if let Some(model) = component.model.as_deref() {
        return Ok(model);
    }
    let model_type = component
        .model_type
        .as_deref()
        .ok_or("equation component is missing type and model")?;
    let identifier = if thermal {
        match model_type {
            "AlphaKT" => "thermal_reference_state",
            "AsymptoticPowerLawMieGruneisenDebye" => "asymptotic_power_law_mie_gruneisen_debye",
            "DoubleDebyeHelmholtz" => "double_debye_helmholtz",
            "DoubleDebyeLogMomentHelmholtz" => "double_debye_log_moment_helmholtz",
            "DorogokupetsOganov2007" => "dorogokupets_oganov_2007",
            "LinearThermalPressure" => "linear_thermal_pressure",
            "LogVolumeThermalPressure" => "log_volume_thermal_pressure",
            "MieGruneisenDebye" => "mie_gruneisen_debye",
            "MieGruneisenEinstein" => "mie_gruneisen_einstein",
            "MultiOscillatorGruneisen" | "Sokolova2016" => {
                "multi_oscillator_gruneisen_thermal_pressure"
            }
            "ThermalModifiedTait" => "thermal_modified_tait",
            _ => return Err(format!("unknown thermal type {model_type:?}")),
        }
    } else {
        match model_type {
            "BM2" => "birch_murnaghan_2",
            "BM3" => "birch_murnaghan_3",
            "BM4" => "birch_murnaghan_4",
            "Holzapfel" => "holzapfel",
            "ModifiedTait" => "modified_tait",
            "Murnaghan" => "murnaghan",
            "NaturalStrain2" => "natural_strain_2",
            "NaturalStrain3" => "natural_strain_3",
            "NaturalStrain4" => "natural_strain_4",
            "Vinet" => "vinet",
            "LinearUsUpHugoniot" => "linear_us_up_hugoniot",
            _ => return Err(format!("unknown isothermal type {model_type:?}")),
        }
    };
    Ok(identifier)
}

fn build_hugoniot(component: &RawComponent, volume_scale: f64) -> Result<HugoniotModel, String> {
    let model = component_model_identifier(component, false)?;
    let p = |name| parameter(component, name);
    match model {
        "linear_us_up_hugoniot" => {
            check_type(component, "LinearUsUpHugoniot")?;
            LinearUsUpHugoniot::new(
                p("V0")? * volume_scale,
                p("rho0")?,
                p("c0")?,
                p("s")?,
                p("P0")?,
            )
            .map(HugoniotModel::LinearUsUp)
            .map_err(|error| error.to_string())
        }
        _ => Err(format!("unknown Hugoniot model {model:?}")),
    }
}

fn build_isothermal(
    component: &RawComponent,
    volume_scale: f64,
) -> Result<IsothermalModel, String> {
    let model = component_model_identifier(component, false)?;
    let p = |name| parameter(component, name);
    let v0 = p("V0")? * volume_scale;
    let built = match model {
        "birch_murnaghan_2" => {
            check_type(component, "BM2")?;
            BM2::new(v0, p("K0")?).map(IsothermalModel::BM2)
        }
        "birch_murnaghan_3" => {
            check_type(component, "BM3")?;
            BM3::new(v0, p("K0")?, p("K0_prime")?).map(IsothermalModel::BM3)
        }
        "birch_murnaghan_4" => {
            check_type(component, "BM4")?;
            BM4::new(v0, p("K0")?, p("K0_prime")?, p("K0_double_prime")?).map(IsothermalModel::BM4)
        }
        "holzapfel" => {
            check_type(component, "Holzapfel")?;
            Holzapfel::new(v0, p("K0")?, p("K0_prime")?, p("n")?, p("Z")?)
                .map(IsothermalModel::Holzapfel)
        }
        "modified_tait" => {
            check_type(component, "ModifiedTait")?;
            ModifiedTait::new(v0, p("K0")?, p("K0_prime")?, p("K0_double_prime")?)
                .map(IsothermalModel::ModifiedTait)
        }
        "murnaghan" => {
            check_type(component, "Murnaghan")?;
            Murnaghan::new(v0, p("K0")?, p("K0_prime")?).map(IsothermalModel::Murnaghan)
        }
        "natural_strain_2" => {
            check_type(component, "NaturalStrain2")?;
            NaturalStrain2::new(v0, p("K0")?).map(IsothermalModel::NaturalStrain2)
        }
        "natural_strain_3" => {
            check_type(component, "NaturalStrain3")?;
            NaturalStrain3::new(v0, p("K0")?, p("K0_prime")?).map(IsothermalModel::NaturalStrain3)
        }
        "natural_strain_4" => {
            check_type(component, "NaturalStrain4")?;
            NaturalStrain4::new(v0, p("K0")?, p("K0_prime")?, p("K0_double_prime")?)
                .map(IsothermalModel::NaturalStrain4)
        }
        "vinet" => {
            check_type(component, "Vinet")?;
            Vinet::new(v0, p("K0")?, p("K0_prime")?).map(IsothermalModel::Vinet)
        }
        _ => return Err(format!("unknown isothermal model {model:?}")),
    };
    built.map_err(|error| error.to_string())
}

fn check_type(component: &RawComponent, expected: &str) -> Result<(), String> {
    match component.model_type.as_deref() {
        Some(actual) if actual == expected => Ok(()),
        Some(actual) => Err(format!(
            "model {:?} requires type {expected:?}, found {actual:?}",
            component.model
        )),
        None => Err(format!("model {:?} is missing type", component.model)),
    }
}

fn configuration<'a>(component: &'a RawComponent, name: &str) -> Option<&'a str> {
    match name {
        "debye_temperature_law" => component.debye_temperature_law.as_deref(),
        "thermal_expansion_law" => component.thermal_expansion_law.as_deref(),
        "reference_volume_law" => component.reference_volume_law.as_deref(),
        _ => None,
    }
    .or_else(|| component.configuration.get(name).and_then(Value::as_str))
}

#[allow(clippy::too_many_lines)]
fn build_thermal(
    component: &RawComponent,
    reference: IsothermalModel,
) -> Result<ThermalModel, String> {
    if let Some((name, _)) = component
        .parameters
        .iter()
        .find(|(name, value)| value.is_none() && name.as_str() != "Tr")
    {
        return Err(format!(
            "model {:?} parameter {name} must be a finite number",
            component.model
        ));
    }
    let model = component_model_identifier(component, true)?;
    let p = |name| parameter(component, name);
    let built = match model {
        "double_debye_helmholtz" => {
            check_type(component, "DoubleDebyeHelmholtz")?;
            let IsothermalModel::Vinet(reference) = reference else {
                return Err("double_debye_helmholtz requires a Vinet reference EOS".to_owned());
            };
            let mut model = DoubleDebyeHelmholtz::new(
                reference,
                p("Vp")?,
                p("theta_a0")?,
                p("a_a")?,
                p("b_a")?,
                p("theta_b0")?,
                p("a_b")?,
                p("b_b")?,
                p("theta_1_0")?,
                p("a_1")?,
                p("b_1")?,
                p("n")?,
                p("alpha0")?,
                p("Ve")?,
                p("kappa")?,
                p("phi0")?,
            )
            .map_err(|error| error.to_string())?;
            if let Some(tr) = component.parameters.get("Tr").copied().flatten() {
                model = model
                    .with_reference_temperature(tr)
                    .map_err(|error| error.to_string())?;
            }
            Ok(ThermalModel::DoubleDebyeHelmholtz(model))
        }
        "double_debye_log_moment_helmholtz" => {
            check_type(component, "DoubleDebyeLogMomentHelmholtz")?;
            let IsothermalModel::Vinet(reference) = reference else {
                return Err(
                    "double_debye_log_moment_helmholtz requires a Vinet reference EOS".to_owned(),
                );
            };
            let mut model = DoubleDebyeLogMomentHelmholtz::new(
                reference,
                p("Vp")?,
                p("theta_a0")?,
                p("a_a")?,
                p("b_a")?,
                p("theta_b0")?,
                p("a_b")?,
                p("b_b")?,
                p("theta_0_0")?,
                p("a_0")?,
                p("b_0")?,
                p("n")?,
                p("anharmonic_a")?,
                p("phi0")?,
            )
            .map_err(|error| error.to_string())?;
            if let Some(tr) = component.parameters.get("Tr").copied().flatten() {
                model = model
                    .with_reference_temperature(tr)
                    .map_err(|error| error.to_string())?;
            }
            Ok(ThermalModel::DoubleDebyeLogMomentHelmholtz(model))
        }
        "dorogokupets_oganov_2007" => {
            check_type(component, "DorogokupetsOganov2007")?;
            let parameters = DorogokupetsOganov2007Parameters {
                tr: p("Tr")?,
                theta_b1: p("theta_B1")?,
                d_b1: p("d_B1")?,
                m_b1: p("m_B1")?,
                theta_b2: p("theta_B2")?,
                d_b2: p("d_B2")?,
                m_b2: p("m_B2")?,
                theta_e1: p("theta_E1")?,
                m_e1: p("m_E1")?,
                theta_e2: p("theta_E2")?,
                m_e2: p("m_E2")?,
                gamma0: p("gamma0")?,
                gamma_inf: p("gamma_inf")?,
                beta: p("beta")?,
                anharmonic_a: p("anharmonic_a")?,
                anharmonic_m: p("anharmonic_m")?,
                electronic_e: p("electronic_e")?,
                electronic_g: p("electronic_g")?,
                defect_h: p("defect_H")?,
                defect_s: p("defect_S")?,
            };
            DorogokupetsOganov2007::new(reference, parameters, p("n")?)
                .map(ThermalModel::DorogokupetsOganov2007)
        }
        "linear_thermal_pressure" => {
            check_type(component, "LinearThermalPressure")?;
            LinearThermalPressure::new(reference, p("Tr")?, p("alpha_KT")?)
                .map(ThermalModel::LinearThermalPressure)
        }
        "log_volume_thermal_pressure" => {
            check_type(component, "LogVolumeThermalPressure")?;
            LogVolumeThermalPressure::new(reference, p("Tr")?, p("alpha_KT_ref")?, p("dK_dT_V")?)
                .map(ThermalModel::LogVolumeThermalPressure)
        }
        "thermal_reference_state" => {
            check_type(component, "AlphaKT")?;
            let expansion_law =
                match configuration(component, "thermal_expansion_law").unwrap_or("constant") {
                    "constant" => ThermalExpansionLaw::Constant,
                    "linear_temperature" => ThermalExpansionLaw::LinearTemperature,
                    value => return Err(format!("unknown thermal_expansion_law {value:?}")),
                };
            let volume_law = match configuration(component, "reference_volume_law")
                .unwrap_or("integrated_expansivity")
            {
                "integrated_expansivity" => ReferenceVolumeLaw::IntegratedExpansivity,
                "linear_temperature" => ReferenceVolumeLaw::LinearTemperature,
                "berman" => ReferenceVolumeLaw::Berman,
                value => return Err(format!("unknown reference_volume_law {value:?}")),
            };
            ThermalReferenceState::new(
                reference,
                p("Tr")?,
                p("alpha0")?,
                p("dK_dT")?,
                component
                    .parameters
                    .get("alpha1")
                    .copied()
                    .flatten()
                    .unwrap_or(0.0),
                expansion_law,
                volume_law,
            )
            .map(ThermalModel::ThermalReferenceState)
        }
        "mie_gruneisen_debye" => {
            check_type(component, "MieGruneisenDebye")?;
            let law = match configuration(component, "debye_temperature_law")
                .unwrap_or("integrated_gruneisen")
            {
                "integrated_gruneisen" => DebyeTemperatureLaw::IntegratedGruneisen,
                "variable_exponent" => DebyeTemperatureLaw::VariableExponent,
                value => return Err(format!("unknown debye_temperature_law {value:?}")),
            };
            MieGruneisenDebye::new_with_temperature_law(
                reference,
                p("Tr")?,
                p("theta0")?,
                p("gamma0")?,
                p("q")?,
                p("n")?,
                law,
            )
            .map(ThermalModel::MieGruneisenDebye)
        }
        "mie_gruneisen_einstein" => {
            check_type(component, "MieGruneisenEinstein")?;
            MieGruneisenEinstein::new(
                reference,
                p("Tr")?,
                p("theta0")?,
                p("gamma0")?,
                p("q")?,
                p("n")?,
            )
            .map(ThermalModel::MieGruneisenEinstein)
        }
        "asymptotic_power_law_mie_gruneisen_debye" => {
            check_type(component, "AsymptoticPowerLawMieGruneisenDebye")?;
            AsymptoticPowerLawMieGruneisenDebye::new(
                reference,
                p("Tr")?,
                p("theta0")?,
                p("gamma0")?,
                p("a")?,
                p("b")?,
                p("n")?,
            )
            .map(ThermalModel::AsymptoticPowerLawMieGruneisenDebye)
        }
        "multi_oscillator_gruneisen_thermal_pressure" => {
            let model_type = component
                .model_type
                .as_deref()
                .ok_or("thermal model is missing type")?;
            if !matches!(model_type, "MultiOscillatorGruneisen" | "Sokolova2016") {
                return Err(format!(
                    "model {model:?} requires type \"MultiOscillatorGruneisen\" or \"Sokolova2016\", found {model_type:?}"
                ));
            }
            let parameters = SokolovaParameters {
                tr: p("Tr")?,
                qe1o: p("QE1o")?,
                me1: p("mE1")?,
                qe2o: p("QE2o")?,
                me2: p("mE2")?,
                delta: p("delta")?,
                t: p("t")?,
                a_0: p("a_0")?,
                m: p("m")?,
                g: p("g")?,
                e_0: p("e_0")?,
                beta: component
                    .parameters
                    .get("beta")
                    .copied()
                    .flatten()
                    .unwrap_or(0.0),
                qbo: component
                    .parameters
                    .get("QBo")
                    .copied()
                    .flatten()
                    .unwrap_or(1.0),
                d: component
                    .parameters
                    .get("d")
                    .copied()
                    .flatten()
                    .unwrap_or(1.0),
                mb: component
                    .parameters
                    .get("mb")
                    .copied()
                    .flatten()
                    .unwrap_or(0.0),
                qb1o: component
                    .parameters
                    .get("QB1o")
                    .copied()
                    .flatten()
                    .unwrap_or(1.0),
                d1: component
                    .parameters
                    .get("d1")
                    .copied()
                    .flatten()
                    .unwrap_or(1.0),
                mb1: component
                    .parameters
                    .get("mb1")
                    .copied()
                    .flatten()
                    .unwrap_or(0.0),
            };
            MultiOscillatorGruneisen::new_with_atom_count(reference, parameters, p("n")?)
                .map(ThermalModel::MultiOscillatorGruneisen)
        }
        "thermal_modified_tait" => {
            check_type(component, "ThermalModifiedTait")?;
            let IsothermalModel::ModifiedTait(reference) = reference else {
                return Err(
                    "thermal_modified_tait requires a modified_tait reference EOS".to_owned(),
                );
            };
            ThermalModifiedTait::new(reference, p("Tr")?, p("theta")?, p("alpha0")?, p("n")?)
                .map(ThermalModel::ThermalModifiedTait)
        }
        _ => return Err(format!("unknown thermal model {model:?}")),
    };
    built.map_err(|error: EosError| error.to_string())
}
