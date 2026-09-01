//! Loading and executable evaluation of Peritheos `.eosmat` material files.
//!
//! The loader accepts canonical Peritheos format 3 and legacy Dioptas format
//! 2 documents. JSON extensions are retained in [`Material::document`] and
//! [`EosRecord::document`], while equation construction is restricted to the
//! built-in model registry.

use std::collections::HashMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use serde::Deserialize;
use serde_json::Value;

use crate::isothermal::{
    Holzapfel, ModifiedTait, Murnaghan, NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2,
    BM3, BM4,
};
use crate::thermal::{
    AsymptoticPowerLawMieGruneisenDebye, DebyeTemperatureLaw, LinearThermalPressure,
    LogVolumeThermalPressure, MieGruneisenDebye, MieGruneisenEinstein, MultiOscillatorGruneisen,
    ReferenceStateEos, ReferenceVolumeLaw, SokolovaParameters, ThermalExpansionLaw,
    ThermalModifiedTait, ThermalReferenceState,
};
use crate::{EosError, EosResult, IsothermalEos, ThermalEos};

const FORMAT: &str = "peritheos.material";
const FORMAT_VERSION: u64 = 3;
const LEGACY_FORMAT_VERSION: u64 = 2;
const CELL_ANGSTROM3_TO_FORMULA_MOLAR_J_PER_BAR: f64 = 0.060_221_407_6;

/// Errors encountered while decoding or constructing an `.eosmat` material.
#[derive(Debug)]
#[non_exhaustive]
pub enum EosmatError {
    /// The file could not be opened or read.
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

impl Display for EosmatError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "could not read eosmat file: {error}"),
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

/// Runtime-dispatched built-in thermal model.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum ThermalModel {
    /// Tange-type asymptotic-power-law Mie--Gruneisen--Debye EOS.
    AsymptoticPowerLawMieGruneisenDebye(AsymptoticPowerLawMieGruneisenDebye<IsothermalModel>),
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

    fn bulk_modulus(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.bulk_modulus(volume, temperature, 1.0e-6))
    }

    fn volume(&self, pressure: f64, temperature: f64) -> EosResult<f64> {
        dispatch_thermal!(self, model => model.volume(pressure, temperature))
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
}

impl LoadedEos {
    /// Whether the record has a thermal correction.
    #[must_use]
    pub const fn is_thermal(&self) -> bool {
        matches!(self, Self::Thermal(_))
    }

    /// Reference-isotherm model identifier.
    #[must_use]
    pub const fn isothermal_model_identifier(&self) -> &'static str {
        match self {
            Self::Isothermal(model) => model.model_identifier(),
            Self::Thermal(model) => match model {
                ThermalModel::AsymptoticPowerLawMieGruneisenDebye(value) => {
                    value.rt_eos.model_identifier()
                }
                ThermalModel::LinearThermalPressure(value) => value.rt_eos.model_identifier(),
                ThermalModel::LogVolumeThermalPressure(value) => value.rt_eos.model_identifier(),
                ThermalModel::MieGruneisenDebye(value) => value.rt_eos.model_identifier(),
                ThermalModel::MieGruneisenEinstein(value) => value.rt_eos.model_identifier(),
                ThermalModel::MultiOscillatorGruneisen(value) => value.rt_eos.model_identifier(),
                ThermalModel::ThermalModifiedTait(_) => "modified_tait",
                ThermalModel::ThermalReferenceState(value) => value.rt_eos.model_identifier(),
            },
        }
    }

    /// Thermal model identifier, if present.
    #[must_use]
    pub const fn thermal_model_identifier(&self) -> Option<&'static str> {
        match self {
            Self::Isothermal(_) => None,
            Self::Thermal(model) => Some(model.model_identifier()),
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
    /// Reference temperature in kelvin.
    pub reference_temperature: f64,
    /// Original record including all extension fields.
    pub document: Value,
    volume_scale: f64,
}

impl EosRecord {
    /// Reference volume in the file's conventional-cell volume unit.
    #[must_use]
    pub fn reference_volume(&self) -> f64 {
        let model_volume = match self.eos {
            LoadedEos::Isothermal(model) => model.reference_volume(),
            LoadedEos::Thermal(model) => model.reference_volume(),
        };
        model_volume / self.volume_scale
    }

    /// Pressure in `GPa` at conventional-cell `volume` and `temperature`.
    ///
    /// Isothermal records ignore `temperature`.
    ///
    /// # Errors
    ///
    /// Returns an error when the state is invalid or model evaluation fails.
    pub fn pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = volume * self.volume_scale;
        match self.eos {
            LoadedEos::Isothermal(model) => model.pressure(volume),
            LoadedEos::Thermal(model) => model.pressure(volume, temperature),
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

    /// Return the preferred record, or the first record when no default is marked.
    #[must_use]
    pub fn default_record(&self) -> Option<&EosRecord> {
        self.eos_records
            .iter()
            .find(|record| record.is_default)
            .or_else(|| self.eos_records.first())
    }
}

#[derive(Debug, Deserialize)]
struct RawMaterial {
    format: Option<String>,
    format_version: Option<u64>,
    identifier: Option<String>,
    name: Option<String>,
    formula: Option<String>,
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
    eos: Option<RawComponent>,
    thermal: Option<RawComponent>,
    temperature_ref: Option<f64>,
    volume: Option<RawVolume>,
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
    parameters: HashMap<String, f64>,
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

/// Construct an executable material from a decoded JSON value.
///
/// # Errors
///
/// Returns an error for an unsupported format or an invalid EOS record.
pub fn material_from_value(document: Value) -> Result<Material, EosmatError> {
    let raw: RawMaterial = serde_json::from_value(document.clone())?;
    let version = raw
        .format_version
        .ok_or_else(|| EosmatError::InvalidDocument("missing format_version".to_owned()))?;
    let canonical = raw.format.as_deref() == Some(FORMAT) && version == FORMAT_VERSION;
    let legacy = raw.format.is_none() && version == LEGACY_FORMAT_VERSION;
    if !canonical && !legacy {
        return Err(EosmatError::InvalidDocument(format!(
            "supported formats are {FORMAT} version {FORMAT_VERSION} and legacy Dioptas version {LEGACY_FORMAT_VERSION}"
        )));
    }

    let name = required_string(raw.name, "name")?;
    let formula = required_string(raw.formula, "formula")?;
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
        )?);
    }

    Ok(Material {
        identifier,
        name,
        formula,
        eos_records,
        document,
    })
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
    );
    result.map_err(|reason| EosmatError::InvalidRecord { identifier, reason })
}

fn build_record(
    raw: &RawRecord,
    identifier: String,
    label: String,
    document: Value,
    formula_units_per_cell: Option<f64>,
) -> Result<EosRecord, String> {
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

    let reference_component = raw.eos.as_ref().ok_or("missing eos component")?;
    let reference = build_isothermal(reference_component, volume_scale)?;
    let eos = match raw.thermal.as_ref() {
        None => LoadedEos::Isothermal(reference),
        Some(thermal) => LoadedEos::Thermal(build_thermal(thermal, reference)?),
    };
    let reference_temperature = raw.temperature_ref.unwrap_or_else(|| match eos {
        LoadedEos::Isothermal(_) => 300.0,
        LoadedEos::Thermal(model) => model.reference_temperature(),
    });
    if !reference_temperature.is_finite() || reference_temperature <= 0.0 {
        return Err("temperature_ref must be positive and finite".to_owned());
    }

    Ok(EosRecord {
        identifier,
        label,
        is_default: raw.default,
        eos,
        reference_temperature,
        document,
        volume_scale,
    })
}

fn is_molar_volume_model(model: &str) -> bool {
    matches!(
        model,
        "mie_gruneisen_debye"
            | "mie_gruneisen_einstein"
            | "asymptotic_power_law_mie_gruneisen_debye"
            | "multi_oscillator_gruneisen_thermal_pressure"
            | "thermal_modified_tait"
    )
}

fn parameter(component: &RawComponent, name: &str) -> Result<f64, String> {
    component
        .parameters
        .get(name)
        .copied()
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
            _ => return Err(format!("unknown isothermal type {model_type:?}")),
        }
    };
    Ok(identifier)
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
    let model = component_model_identifier(component, true)?;
    let p = |name| parameter(component, name);
    let built = match model {
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
                value => return Err(format!("unknown reference_volume_law {value:?}")),
            };
            ThermalReferenceState::new(
                reference,
                p("Tr")?,
                p("alpha0")?,
                p("dK_dT")?,
                component.parameters.get("alpha1").copied().unwrap_or(0.0),
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
                beta: component.parameters.get("beta").copied().unwrap_or(0.0),
                qbo: component.parameters.get("QBo").copied().unwrap_or(1.0),
                d: component.parameters.get("d").copied().unwrap_or(1.0),
                mb: component.parameters.get("mb").copied().unwrap_or(0.0),
                qb1o: component.parameters.get("QB1o").copied().unwrap_or(1.0),
                d1: component.parameters.get("d1").copied().unwrap_or(1.0),
                mb1: component.parameters.get("mb1").copied().unwrap_or(0.0),
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
