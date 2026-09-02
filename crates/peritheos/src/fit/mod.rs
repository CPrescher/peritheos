//! Fit Peritheos models and propagate parameter uncertainty in pure Rust.
//!
//! The highest-level entry points accept physical observations plus a model
//! factory. They assemble weighted residuals, bounds, and optional latent
//! state variables internally. Lower-level least-squares and covariance
//! routines remain available for custom models.
//!
//! # Fit an isothermal EOS
//!
//! This complete example estimates `K0` and `K0_prime` while keeping `V0`
//! fixed. The factory's parameter order is the order used by the initial
//! values, bounds, result, and covariance matrix.
//!
//! ```
//! use peritheos::isothermal::BM3;
//! use peritheos::fit::{
//!     fit_isothermal_eos, FitError, IsothermalObservations, SolverOptions,
//! };
//!
//! let pressure = [39.1, 25.3, 15.2, 7.7, 2.2];
//! let volume = [8.2, 8.6, 9.0, 9.4, 9.8];
//! let pressure_sigma = [0.2; 5];
//! let result = fit_isothermal_eos(
//!     IsothermalObservations {
//!         pressure: &pressure,
//!         volume: &volume,
//!         pressure_sigma: &pressure_sigma,
//!         volume_sigma: None,
//!         observation_cholesky: None,
//!     },
//!     &[140.0, 4.0],
//!     &[50.0, 1.0],
//!     &[300.0, 10.0],
//!     SolverOptions::default(),
//!     |parameters| {
//!         BM3::new(10.0, parameters[0], parameters[1])
//!             .map_err(FitError::from)
//!     },
//! )?;
//!
//! assert_eq!(result.parameters.len(), 2);
//! assert_eq!(result.covariance.len(), 4);
//! assert_eq!(result.standard_errors.len(), 2);
//! assert_eq!(result.predicted_pressure.len(), pressure.len());
//! # Ok::<(), FitError>(())
//! ```
//!
//! Set [`IsothermalObservations::volume_sigma`] to fit adjusted (latent)
//! volumes. For P-V-T data use [`fit_thermal_eos`] and
//! [`ThermalObservations`]; [`fit_joint_eos`] uses one ordered parameter slice
//! to reconstruct both the reference and thermal components.
//!
//! # Propagate fitted uncertainty
//!
//! [`parameter_covariance`] converts a final residual Jacobian to parameter
//! covariance. [`propagate_model_uncertainty`] then evaluates a model
//! Jacobian and applies the delta method. Use
//! [`monte_carlo_model_uncertainty`] when nonlinearity or invalid sampled
//! states matter.
//!
//! ```
//! use peritheos::{isothermal::BM3, IsothermalEos};
//! use peritheos::fit::{propagate_model_uncertainty, FitError};
//!
//! let parameters = [160.0, 4.0];
//! let covariance = [4.0, -0.05, -0.05, 0.04];
//! let volumes = [9.5, 9.0];
//! let propagated = propagate_model_uncertainty(
//!     &parameters,
//!     &covariance,
//!     &[0.0; 2],
//!     1.0e-5,
//!     true,
//!     |parameters| {
//!         let eos = BM3::new(10.0, parameters[0], parameters[1])
//!             .map_err(FitError::from)?;
//!         volumes
//!             .iter()
//!             .map(|&volume| {
//!                 eos.pressure(volume)
//!                     .map_err(FitError::from)
//!             })
//!             .collect()
//!     },
//! )?;
//!
//! assert_eq!(propagated.model.nominal.len(), 2);
//! assert!(propagated.propagation.variance.iter().all(|value| *value >= 0.0));
//! # Ok::<(), FitError>(())
//! ```
//!
//! Robust losses are selected with [`Loss`] through [`SolverOptions`]. All
//! errors are reported as [`FitError`], so factories can use
//! [`FitError::from`] to retain the original EOS error as a source.

mod eos;

pub use eos::{
    fit_isothermal_eos, fit_joint_eos, fit_thermal_eos, fit_thermal_eos_by, EosFitResult,
    IsothermalObservations, ThermalObservations,
};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::error::Error;
use std::fmt::{Display, Formatter};

use crate::EosError;

const EPSILON_SQRT: f64 = 1.490_116_119_384_765_6e-8;
const DEFAULT_TOLERANCE: f64 = 1.0e-8;

/// Machine-readable category for a [`FitError`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FitErrorKind {
    /// Observations, dimensions, bounds, or solver options are invalid.
    InvalidInput,
    /// The model or user callback failed during evaluation.
    Evaluation,
    /// An EOS failed during model construction or evaluation.
    EosEvaluation,
    /// A required linear system is singular.
    SingularSystem,
}

/// Error returned by fitting and covariance kernels.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum FitError {
    /// Observations, dimensions, bounds, or solver options are invalid.
    InvalidInput(String),
    /// A user-supplied model or callback failed and supplied diagnostic text.
    Evaluation(String),
    /// An EOS failed during model construction or evaluation.
    EosEvaluation(EosError),
    /// A required linear system is singular.
    SingularSystem,
}

impl FitError {
    /// Return the machine-readable error category.
    #[must_use]
    pub const fn kind(&self) -> FitErrorKind {
        match self {
            Self::InvalidInput(_) => FitErrorKind::InvalidInput,
            Self::Evaluation(_) => FitErrorKind::Evaluation,
            Self::EosEvaluation(_) => FitErrorKind::EosEvaluation,
            Self::SingularSystem => FitErrorKind::SingularSystem,
        }
    }

    /// Return a stable, language-independent error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidInput(_) => "fit.invalid_input",
            Self::Evaluation(_) => "fit.evaluation_failed",
            Self::EosEvaluation(error) => error.code(),
            Self::SingularSystem => "fit.singular_system",
        }
    }
}

impl Display for FitError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(message) | Self::Evaluation(message) => formatter.write_str(message),
            Self::EosEvaluation(error) => write!(formatter, "EOS evaluation failed: {error}"),
            Self::SingularSystem => formatter.write_str("linear system is singular"),
        }
    }
}

impl Error for FitError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::EosEvaluation(error) => Some(error),
            Self::InvalidInput(_) | Self::Evaluation(_) | Self::SingularSystem => None,
        }
    }
}

impl From<EosError> for FitError {
    fn from(error: EosError) -> Self {
        Self::EosEvaluation(error)
    }
}

/// Robust loss applied to squared, scaled residuals.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Loss {
    #[default]
    Linear,
    SoftL1,
    Huber,
    Cauchy,
    Arctan,
}

impl Loss {
    /// Parse a public Python loss name.
    ///
    /// # Errors
    ///
    /// Returns [`FitError::InvalidInput`] for an unsupported name.
    pub fn from_name(name: &str) -> Result<Self, FitError> {
        match name {
            "linear" => Ok(Self::Linear),
            "soft_l1" => Ok(Self::SoftL1),
            "huber" => Ok(Self::Huber),
            "cauchy" => Ok(Self::Cauchy),
            "arctan" => Ok(Self::Arctan),
            _ => Err(FitError::InvalidInput(format!(
                "unsupported loss function: {name}"
            ))),
        }
    }

    fn rho(self, squared: f64) -> f64 {
        match self {
            Self::Linear => squared,
            Self::SoftL1 => 2.0 * ((1.0 + squared).sqrt() - 1.0),
            Self::Huber if squared <= 1.0 => squared,
            Self::Huber => 2.0 * squared.sqrt() - 1.0,
            Self::Cauchy => squared.ln_1p(),
            Self::Arctan => squared.atan(),
        }
    }

    fn derivative(self, squared: f64) -> f64 {
        match self {
            Self::Linear => 1.0,
            Self::SoftL1 => 1.0 / (1.0 + squared).sqrt(),
            Self::Huber if squared <= 1.0 => 1.0,
            Self::Huber => 1.0 / squared.sqrt(),
            Self::Cauchy => 1.0 / (1.0 + squared),
            Self::Arctan => 1.0 / (1.0 + squared * squared),
        }
    }

    fn second_derivative(self, squared: f64) -> f64 {
        match self {
            Self::Linear => 0.0,
            Self::SoftL1 => -0.5 / (1.0 + squared).powf(1.5),
            Self::Huber if squared <= 1.0 => 0.0,
            Self::Huber => -0.5 / squared.powf(1.5),
            Self::Cauchy => -1.0 / (1.0 + squared).powi(2),
            Self::Arctan => -2.0 * squared / (1.0 + squared * squared).powi(2),
        }
    }
}

/// Controls the bounded nonlinear least-squares solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolverOptions {
    pub loss: Loss,
    pub f_scale: f64,
    pub max_evaluations: Option<usize>,
    pub function_tolerance: f64,
    pub gradient_tolerance: f64,
    pub step_tolerance: f64,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            loss: Loss::Linear,
            f_scale: 1.0,
            max_evaluations: None,
            function_tolerance: DEFAULT_TOLERANCE,
            gradient_tolerance: DEFAULT_TOLERANCE,
            step_tolerance: DEFAULT_TOLERANCE,
        }
    }
}

/// Result and diagnostics of a least-squares solve.
#[derive(Clone, Debug, PartialEq)]
pub struct SolverResult {
    pub parameters: Vec<f64>,
    pub residuals: Vec<f64>,
    /// Row-major residual Jacobian at the returned parameters.
    pub jacobian: Vec<f64>,
    pub residual_count: usize,
    pub cost: f64,
    pub optimality: f64,
    pub success: bool,
    pub status: i32,
    pub message: String,
    pub function_evaluations: usize,
    pub jacobian_evaluations: usize,
}

/// Block structure of an errors-in-variables least-squares problem.
///
/// Variables are ordered as global parameters followed by one point-sized
/// block per adjusted coordinate. Residuals use the same point-sized block
/// ordering, beginning with pressure residuals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StructuredLayout {
    pub global_parameter_count: usize,
    pub point_count: usize,
    pub latent_coordinate_count: usize,
}

/// Variance and optional full output covariance from delta-method propagation.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearPropagation {
    pub variance: Vec<f64>,
    pub covariance: Option<Vec<f64>>,
}

/// Summary statistics for accepted Monte Carlo output samples.
#[derive(Clone, Debug, PartialEq)]
pub struct MonteCarloSummary {
    pub standard_error: Vec<f64>,
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
    pub covariance: Option<Vec<f64>>,
}

/// Nominal model output and its finite-difference parameter Jacobian.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelJacobian {
    /// Model output at the unperturbed parameter vector.
    pub nominal: Vec<f64>,
    /// Row-major output-by-parameter Jacobian.
    pub jacobian: Vec<f64>,
}

/// Complete local linear model-uncertainty result.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelLinearPropagation {
    /// Model output and finite-difference parameter Jacobian.
    pub model: ModelJacobian,
    /// Propagated output variance and optional covariance.
    pub propagation: LinearPropagation,
}

/// Controls deterministic native Monte Carlo model evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MonteCarloOptions {
    /// Number of accepted output samples.
    pub sample_count: usize,
    /// Maximum attempted samples before returning an evaluation error.
    pub max_attempts: usize,
    /// Central confidence interval probability.
    pub confidence: f64,
    /// Whether to retain the full output covariance.
    pub full_covariance: bool,
    /// Deterministic native RNG seed.
    pub seed: u64,
}

impl Default for MonteCarloOptions {
    fn default() -> Self {
        Self {
            sample_count: 5_000,
            max_attempts: 100_000,
            confidence: 0.95,
            full_covariance: false,
            seed: 0,
        }
    }
}

/// Model-aware native Monte Carlo uncertainty result.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelMonteCarlo {
    /// Model output at the mean parameters and states.
    pub nominal: Vec<f64>,
    /// Statistics of accepted sampled model outputs.
    pub summary: MonteCarloSummary,
    /// Number of sampled parameter/state vectors attempted.
    pub attempted_samples: usize,
    /// Fraction of attempts rejected because model evaluation failed.
    pub rejected_fraction: f64,
}

/// Evaluate a model and its parameter Jacobian using adaptive centered differences.
///
/// Perturbation scales follow the Python compatibility convention:
/// `relative_step * max(abs(parameter), standard_error, 1)`. If only one side
/// of a perturbation is valid, a one-sided derivative is used.
///
/// # Errors
///
/// Returns an error for inconsistent or non-finite inputs, changed output
/// dimensions, or a parameter that cannot be perturbed on either side.
pub fn finite_difference_model_jacobian<F>(
    parameters: &[f64],
    parameter_covariance: &[f64],
    relative_step: f64,
    evaluate: F,
) -> Result<ModelJacobian, FitError>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, FitError>,
{
    let parameter_count = parameters.len();
    if parameter_count == 0
        || parameter_covariance.len() != parameter_count * parameter_count
        || !relative_step.is_finite()
        || relative_step <= 0.0
        || parameters
            .iter()
            .chain(parameter_covariance)
            .any(|value| !value.is_finite())
    {
        return Err(FitError::InvalidInput(
            "invalid model Jacobian parameters, covariance, or step".to_owned(),
        ));
    }
    positive_semidefinite_factor(parameter_covariance, parameter_count)?;
    let nominal = evaluate(parameters)?;
    if nominal.is_empty() || nominal.iter().any(|value| !value.is_finite()) {
        return Err(FitError::Evaluation(
            "nominal model output must be non-empty and finite".to_owned(),
        ));
    }
    let output_count = nominal.len();
    let mut jacobian = vec![0.0; output_count * parameter_count];
    for parameter in 0..parameter_count {
        let standard_error = parameter_covariance[parameter * parameter_count + parameter]
            .max(0.0)
            .sqrt();
        let step = relative_step * parameters[parameter].abs().max(standard_error).max(1.0);
        let mut plus_parameters = parameters.to_vec();
        let mut minus_parameters = parameters.to_vec();
        plus_parameters[parameter] += step;
        minus_parameters[parameter] -= step;
        let plus = evaluate(&plus_parameters).ok().filter(|values| {
            values.len() == output_count && values.iter().all(|value| value.is_finite())
        });
        let minus = evaluate(&minus_parameters).ok().filter(|values| {
            values.len() == output_count && values.iter().all(|value| value.is_finite())
        });
        for output in 0..output_count {
            jacobian[output * parameter_count + parameter] = match (&plus, &minus) {
                (Some(plus), Some(minus)) => (plus[output] - minus[output]) / (2.0 * step),
                (Some(plus), None) => (plus[output] - nominal[output]) / step,
                (None, Some(minus)) => (nominal[output] - minus[output]) / step,
                (None, None) => {
                    return Err(FitError::Evaluation(format!(
                        "could not perturb model parameter at index {parameter}"
                    )));
                }
            };
        }
    }
    Ok(ModelJacobian { nominal, jacobian })
}

/// Evaluate and propagate local linear parameter uncertainty for a Rust model.
///
/// `state_variance` contains any independently computed state-variable
/// contribution for each output. Pass zeros when only parameter uncertainty
/// is required.
///
/// # Errors
///
/// Returns an error from finite differences or covariance propagation.
pub fn propagate_model_uncertainty<F>(
    parameters: &[f64],
    parameter_covariance: &[f64],
    state_variance: &[f64],
    relative_step: f64,
    full_covariance: bool,
    evaluate: F,
) -> Result<ModelLinearPropagation, FitError>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, FitError>,
{
    let model = finite_difference_model_jacobian(
        parameters,
        parameter_covariance,
        relative_step,
        evaluate,
    )?;
    let propagation = propagate_linear_uncertainty(
        &model.jacobian,
        model.nominal.len(),
        parameters.len(),
        parameter_covariance,
        state_variance,
        full_covariance,
    )?;
    Ok(ModelLinearPropagation { model, propagation })
}

/// Sample parameters and independent states, evaluate a model, and summarize accepted outputs.
///
/// Parameter samples use a multivariate normal distribution with the supplied
/// positive-semidefinite covariance. State samples are independent normals;
/// pass empty state slices when the evaluator has no uncertain states. Native
/// seeded results are deterministic within this backend but intentionally do
/// not reproduce `NumPy`'s random stream.
///
/// # Errors
///
/// Returns an error for invalid dimensions or covariance, a changed model
/// output shape, or too few valid model evaluations before `max_attempts`.
#[allow(clippy::too_many_arguments)]
pub fn monte_carlo_model_uncertainty<F>(
    parameter_means: &[f64],
    parameter_covariance: &[f64],
    state_means: &[f64],
    state_sigmas: &[f64],
    options: MonteCarloOptions,
    evaluate: F,
) -> Result<ModelMonteCarlo, FitError>
where
    F: Fn(&[f64], &[f64]) -> Result<Vec<f64>, FitError>,
{
    let parameter_count = parameter_means.len();
    if parameter_count == 0
        || parameter_covariance.len() != parameter_count * parameter_count
        || state_means.len() != state_sigmas.len()
        || options.sample_count < 2
        || options.max_attempts < options.sample_count
        || !options.confidence.is_finite()
        || !(0.0..1.0).contains(&options.confidence)
        || parameter_means
            .iter()
            .chain(parameter_covariance)
            .chain(state_means)
            .chain(state_sigmas)
            .any(|value| !value.is_finite())
        || state_sigmas.iter().any(|value| *value <= 0.0)
    {
        return Err(FitError::InvalidInput(
            "invalid Monte Carlo parameters, states, or options".to_owned(),
        ));
    }
    let factor = positive_semidefinite_factor(parameter_covariance, parameter_count)?;
    let nominal = evaluate(parameter_means, state_means)?;
    if nominal.is_empty() || nominal.iter().any(|value| !value.is_finite()) {
        return Err(FitError::Evaluation(
            "nominal Monte Carlo model output must be non-empty and finite".to_owned(),
        ));
    }
    let output_count = nominal.len();
    let mut rng = StdRng::seed_from_u64(options.seed);
    let mut samples = Vec::with_capacity(options.sample_count * output_count);
    let mut attempts = 0;
    while samples.len() / output_count < options.sample_count && attempts < options.max_attempts {
        attempts += 1;
        let parameter_normals = standard_normal_values(&mut rng, parameter_count);
        let mut parameters = parameter_means.to_vec();
        for row in 0..parameter_count {
            parameters[row] += (0..=row)
                .map(|column| factor[row * parameter_count + column] * parameter_normals[column])
                .sum::<f64>();
        }
        let state_normals = standard_normal_values(&mut rng, state_means.len());
        let states = state_means
            .iter()
            .zip(state_sigmas)
            .zip(state_normals)
            .map(|((&mean, &sigma), normal)| mean + sigma * normal)
            .collect::<Vec<_>>();
        if let Ok(values) = evaluate(&parameters, &states) {
            if values.len() == output_count && values.iter().all(|value| value.is_finite()) {
                samples.extend(values);
            }
        }
    }
    if samples.len() / output_count < options.sample_count {
        return Err(FitError::Evaluation(format!(
            "accepted fewer than {} valid Monte Carlo samples in {} attempts",
            options.sample_count, attempts
        )));
    }
    let summary = summarize_monte_carlo(
        &samples,
        options.sample_count,
        output_count,
        options.confidence,
        options.full_covariance,
    )?;
    let accepted = options.sample_count;
    let rejected_fraction = f64::from(u32::try_from(attempts - accepted).unwrap_or(u32::MAX))
        / f64::from(u32::try_from(attempts).unwrap_or(u32::MAX));
    Ok(ModelMonteCarlo {
        nominal,
        summary,
        attempted_samples: attempts,
        rejected_fraction,
    })
}

fn positive_semidefinite_factor(covariance: &[f64], size: usize) -> Result<Vec<f64>, FitError> {
    let scale = covariance
        .iter()
        .map(|value| value.abs())
        .fold(1.0, f64::max);
    let tolerance = 1.0e-12 * scale;
    let mut factor = vec![0.0; size * size];
    for row in 0..size {
        for column in 0..=row {
            let mirrored = covariance[column * size + row];
            let value = covariance[row * size + column];
            if (value - mirrored).abs() > tolerance {
                return Err(FitError::InvalidInput(
                    "parameter covariance must be symmetric".to_owned(),
                ));
            }
            let correction = (0..column)
                .map(|index| factor[row * size + index] * factor[column * size + index])
                .sum::<f64>();
            let remainder = value - correction;
            if row == column {
                if remainder < -tolerance {
                    return Err(FitError::InvalidInput(
                        "parameter covariance must be positive semidefinite".to_owned(),
                    ));
                }
                factor[row * size + column] = remainder.max(0.0).sqrt();
            } else if factor[column * size + column] > tolerance {
                factor[row * size + column] = remainder / factor[column * size + column];
            } else if remainder.abs() > tolerance {
                return Err(FitError::InvalidInput(
                    "parameter covariance must be positive semidefinite".to_owned(),
                ));
            }
        }
    }
    Ok(factor)
}

fn standard_normal_values(rng: &mut StdRng, count: usize) -> Vec<f64> {
    let mut values = Vec::with_capacity(count);
    while values.len() < count {
        let first = rng.gen_range(f64::MIN_POSITIVE..1.0);
        let second = rng.gen_range(0.0..1.0);
        let radius = (-2.0 * first.ln()).sqrt();
        let angle = std::f64::consts::TAU * second;
        values.push(radius * angle.cos());
        if values.len() < count {
            values.push(radius * angle.sin());
        }
    }
    values
}

/// Propagate parameter covariance through an output-by-parameter Jacobian.
///
/// Independent state-variable contributions are supplied as a variance per
/// output and are added to the diagonal of the output covariance.
///
/// # Errors
///
/// Returns an error for inconsistent dimensions or non-finite inputs.
pub fn propagate_linear_uncertainty(
    jacobian: &[f64],
    output_count: usize,
    parameter_count: usize,
    parameter_covariance: &[f64],
    state_variance: &[f64],
    full_covariance: bool,
) -> Result<LinearPropagation, FitError> {
    if output_count == 0
        || jacobian.len() != output_count * parameter_count
        || parameter_covariance.len() != parameter_count * parameter_count
        || state_variance.len() != output_count
    {
        return Err(FitError::InvalidInput(
            "invalid dimensions for linear uncertainty propagation".to_owned(),
        ));
    }
    if jacobian
        .iter()
        .chain(parameter_covariance)
        .chain(state_variance)
        .any(|value| !value.is_finite())
    {
        return Err(FitError::InvalidInput(
            "linear uncertainty inputs must be finite".to_owned(),
        ));
    }
    if state_variance.iter().any(|value| *value < 0.0) {
        return Err(FitError::InvalidInput(
            "state variances must be nonnegative".to_owned(),
        ));
    }
    positive_semidefinite_factor(parameter_covariance, parameter_count)?;
    let mut transformed = vec![0.0; output_count * parameter_count];
    for output in 0..output_count {
        for right in 0..parameter_count {
            transformed[output * parameter_count + right] = (0..parameter_count)
                .map(|left| {
                    jacobian[output * parameter_count + left]
                        * parameter_covariance[left * parameter_count + right]
                })
                .sum();
        }
    }
    let mut variance = vec![0.0; output_count];
    for output in 0..output_count {
        let parameter_variance = (0..parameter_count)
            .map(|parameter| {
                transformed[output * parameter_count + parameter]
                    * jacobian[output * parameter_count + parameter]
            })
            .sum::<f64>();
        variance[output] = (parameter_variance + state_variance[output]).max(0.0);
    }
    let covariance = full_covariance.then(|| {
        let mut output_covariance = vec![0.0; output_count * output_count];
        for left in 0..output_count {
            for right in 0..output_count {
                output_covariance[left * output_count + right] = (0..parameter_count)
                    .map(|parameter| {
                        transformed[left * parameter_count + parameter]
                            * jacobian[right * parameter_count + parameter]
                    })
                    .sum();
            }
            output_covariance[left * output_count + left] += state_variance[left];
        }
        output_covariance
    });
    Ok(LinearPropagation {
        variance,
        covariance,
    })
}

/// Summarize row-major Monte Carlo samples using sample variance (`ddof=1`).
///
/// Quantiles use the same linear interpolation convention as `NumPy`'s default.
///
/// # Errors
///
/// Returns an error for fewer than two samples, inconsistent dimensions,
/// non-finite samples, or a confidence level outside `(0, 1)`.
pub fn summarize_monte_carlo(
    samples: &[f64],
    sample_count: usize,
    output_count: usize,
    confidence: f64,
    full_covariance: bool,
) -> Result<MonteCarloSummary, FitError> {
    if sample_count < 2
        || output_count == 0
        || samples.len() != sample_count * output_count
        || !confidence.is_finite()
        || !(0.0..1.0).contains(&confidence)
        || samples.iter().any(|value| !value.is_finite())
    {
        return Err(FitError::InvalidInput(
            "invalid Monte Carlo samples or confidence".to_owned(),
        ));
    }
    let sample_count_f64 = f64::from(u32::try_from(sample_count).unwrap_or(u32::MAX));
    let denominator = f64::from(u32::try_from(sample_count - 1).unwrap_or(u32::MAX));
    let means: Vec<_> = (0..output_count)
        .map(|output| {
            (0..sample_count)
                .map(|sample| samples[sample * output_count + output])
                .sum::<f64>()
                / sample_count_f64
        })
        .collect();
    let mut standard_error = vec![0.0; output_count];
    for output in 0..output_count {
        standard_error[output] = ((0..sample_count)
            .map(|sample| (samples[sample * output_count + output] - means[output]).powi(2))
            .sum::<f64>()
            / denominator)
            .sqrt();
    }
    let tail = (1.0 - confidence) / 2.0;
    let mut lower = vec![0.0; output_count];
    let mut upper = vec![0.0; output_count];
    for output in 0..output_count {
        let mut column: Vec<_> = (0..sample_count)
            .map(|sample| samples[sample * output_count + output])
            .collect();
        column.sort_by(f64::total_cmp);
        lower[output] = linear_quantile(&column, tail);
        upper[output] = linear_quantile(&column, 1.0 - tail);
    }
    let covariance = full_covariance.then(|| {
        let mut matrix = vec![0.0; output_count * output_count];
        for left in 0..output_count {
            for right in 0..=left {
                let value = (0..sample_count)
                    .map(|sample| {
                        (samples[sample * output_count + left] - means[left])
                            * (samples[sample * output_count + right] - means[right])
                    })
                    .sum::<f64>()
                    / denominator;
                matrix[left * output_count + right] = value;
                matrix[right * output_count + left] = value;
            }
        }
        matrix
    });
    Ok(MonteCarloSummary {
        standard_error,
        lower,
        upper,
        covariance,
    })
}

fn linear_quantile(sorted: &[f64], probability: f64) -> f64 {
    let last_index = sorted.len() - 1;
    let location = probability * f64::from(u32::try_from(last_index).unwrap_or(u32::MAX));
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let lower_index = location.floor() as usize;
    let upper_index = lower_index.saturating_add(1).min(last_index);
    let fraction = location - f64::from(u32::try_from(lower_index).unwrap_or(u32::MAX));
    sorted[lower_index] + fraction * (sorted[upper_index] - sorted[lower_index])
}

/// Solve a bounded nonlinear least-squares problem.
///
/// # Errors
///
/// Returns an error for inconsistent inputs, invalid residual evaluations, or
/// a linear algebra failure that cannot be regularized.
pub fn least_squares<F>(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    evaluate: F,
) -> Result<SolverResult, FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    least_squares_internal(initial, lower, upper, options, None, evaluate)
}

/// Solve a least-squares problem with observation-local latent coordinates.
///
/// This uses simultaneous finite-difference coloring for each latent coordinate
/// and solves damped normal equations through a block Schur complement.
///
/// # Errors
///
/// Returns an error for an inconsistent layout, invalid inputs or residuals,
/// or a linear algebra failure that cannot be regularized.
pub fn least_squares_structured<F>(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    layout: StructuredLayout,
    evaluate: F,
) -> Result<SolverResult, FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    if layout.global_parameter_count == 0
        || layout.point_count == 0
        || layout.latent_coordinate_count == 0
        || initial.len()
            != layout.global_parameter_count + layout.point_count * layout.latent_coordinate_count
    {
        return Err(FitError::InvalidInput(
            "structured layout does not match the parameter vector".to_owned(),
        ));
    }
    least_squares_internal(initial, lower, upper, options, Some(layout), evaluate)
}

#[allow(clippy::too_many_lines)]
fn least_squares_internal<F>(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    layout: Option<StructuredLayout>,
    mut evaluate: F,
) -> Result<SolverResult, FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    validate_problem(initial, lower, upper, options)?;
    let parameter_count = initial.len();
    let maximum_evaluations = options
        .max_evaluations
        .unwrap_or(100 * parameter_count.max(1));
    let mut evaluations = 1;
    let mut jacobian_evaluations = 0;
    let mut parameters = initial.to_vec();
    let mut residuals = evaluate(&parameters)?;
    validate_residuals(&residuals, None)?;
    let residual_count = residuals.len();
    if layout.is_some_and(|structure| {
        residual_count != structure.point_count * (1 + structure.latent_coordinate_count)
    }) {
        return Err(FitError::InvalidInput(
            "structured layout does not match the residual vector".to_owned(),
        ));
    }
    let mut cost = robust_cost(&residuals, options.loss, options.f_scale);
    let mut damping = 1.0e-3;
    let mut status = 0;
    let mut message = "The maximum number of function evaluations is exceeded.".to_owned();
    let mut success = false;
    let mut final_jacobian = Vec::new();
    let mut optimality = f64::INFINITY;

    while evaluations < maximum_evaluations {
        let (jacobian, _used) = finite_difference_for_layout(
            &parameters,
            lower,
            upper,
            &residuals,
            layout,
            &mut evaluate,
        )?;
        jacobian_evaluations += 1;
        final_jacobian.clone_from(&jacobian);
        let weights = robust_weights(&residuals, options.loss, options.f_scale);
        let (normal, gradient) = normal_equations_for_layout(
            &jacobian,
            &residuals,
            &weights,
            residual_count,
            parameter_count,
            layout,
        );
        optimality = projected_gradient_norm(&gradient, &parameters, lower, upper);
        if optimality <= options.gradient_tolerance {
            status = 1;
            "`gtol` termination condition is satisfied.".clone_into(&mut message);
            success = true;
            break;
        }

        let mut accepted = false;
        let mut accepted_step_norm = 0.0;
        let mut smallest_step_norm = f64::INFINITY;
        let mut computed_step = false;
        let previous_cost = cost;
        for _ in 0..24 {
            let Ok(step) = solve_normal_step(&normal, &gradient, damping) else {
                damping *= 10.0;
                continue;
            };
            computed_step = true;
            let mut candidate = parameters.clone();
            for index in 0..parameter_count {
                candidate[index] =
                    (parameters[index] + step[index]).clamp(lower[index], upper[index]);
            }
            accepted_step_norm = scaled_norm_difference(&candidate, &parameters);
            smallest_step_norm = smallest_step_norm.min(accepted_step_norm);
            if accepted_step_norm == 0.0 {
                damping *= 10.0;
                continue;
            }
            if evaluations >= maximum_evaluations {
                break;
            }
            let candidate_residuals = evaluate(&candidate)?;
            evaluations += 1;
            validate_residuals(&candidate_residuals, Some(residual_count))?;
            let candidate_cost = robust_cost(&candidate_residuals, options.loss, options.f_scale);
            if candidate_cost < cost {
                parameters = candidate;
                residuals = candidate_residuals;
                cost = candidate_cost;
                damping = (damping * 0.3).max(1.0e-15);
                accepted = true;
                break;
            }
            damping = (damping * 10.0).min(1.0e30);
        }

        if !accepted {
            if evaluations >= maximum_evaluations {
                break;
            }
            let parameter_norm = parameters
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();
            let step_tolerance = options.step_tolerance * (options.step_tolerance + parameter_norm);
            if computed_step && smallest_step_norm <= step_tolerance {
                status = 3;
                "`xtol` termination condition is satisfied.".clone_into(&mut message);
                success = true;
            } else {
                status = 0;
                "No acceptable optimization step could be computed.".clone_into(&mut message);
            }
            break;
        }
        let parameter_norm = parameters
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        let step_converged = accepted_step_norm
            <= options.step_tolerance * (options.step_tolerance + parameter_norm);
        let function_converged =
            (previous_cost - cost).abs() <= options.function_tolerance * previous_cost.max(1.0);
        if step_converged || function_converged {
            status = if step_converged && function_converged {
                4
            } else if function_converged {
                2
            } else {
                3
            };
            match status {
                4 => "Both `ftol` and `xtol` termination conditions are satisfied.",
                2 => "`ftol` termination condition is satisfied.",
                _ => "`xtol` termination condition is satisfied.",
            }
            .clone_into(&mut message);
            success = true;
            break;
        }
    }

    if final_jacobian.is_empty() || success {
        let (jacobian, _used) = finite_difference_for_layout(
            &parameters,
            lower,
            upper,
            &residuals,
            layout,
            &mut evaluate,
        )?;
        jacobian_evaluations += 1;
        final_jacobian = jacobian;
        let weights = robust_weights(&residuals, options.loss, options.f_scale);
        let (_, gradient) = normal_equations_for_layout(
            &final_jacobian,
            &residuals,
            &weights,
            residual_count,
            parameter_count,
            layout,
        );
        optimality = projected_gradient_norm(&gradient, &parameters, lower, upper);
    }

    let final_jacobian = robust_result_jacobian(
        final_jacobian,
        &residuals,
        parameter_count,
        options.loss,
        options.f_scale,
    );
    Ok(SolverResult {
        parameters,
        residuals,
        jacobian: final_jacobian,
        residual_count,
        cost,
        optimality,
        success,
        status,
        message,
        function_evaluations: evaluations,
        jacobian_evaluations,
    })
}

/// Return `(J^T J)^+` for the leading parameter columns of a Jacobian.
///
/// Latent variables are profiled out through the Schur complement when
/// `parameter_count` is smaller than the total number of columns.
///
/// # Errors
///
/// Returns an error for inconsistent dimensions or a singular information
/// matrix that cannot be regularized.
pub fn parameter_covariance(
    jacobian: &[f64],
    rows: usize,
    columns: usize,
    parameter_count: usize,
) -> Result<Vec<f64>, FitError> {
    if rows == 0
        || columns == 0
        || parameter_count == 0
        || parameter_count > columns
        || jacobian.len() != rows * columns
    {
        return Err(FitError::InvalidInput(
            "invalid Jacobian dimensions for covariance".to_owned(),
        ));
    }
    let mut information = vec![0.0; columns * columns];
    for row in 0..rows {
        for left in 0..columns {
            let value = jacobian[row * columns + left];
            for right in 0..=left {
                information[left * columns + right] += value * jacobian[row * columns + right];
            }
        }
    }
    for left in 0..columns {
        for right in 0..left {
            information[right * columns + left] = information[left * columns + right];
        }
    }
    let profiled = if parameter_count == columns {
        information
    } else {
        schur_complement(&information, columns, parameter_count)?
    };
    symmetric_pseudoinverse(&profiled, parameter_count, rows.max(columns))
}

/// Return the profiled model-parameter covariance for a structured fit.
///
/// This preserves the observation-local latent-coordinate blocks instead of
/// constructing a dense latent information matrix.
///
/// # Errors
///
/// Returns an error for inconsistent dimensions or a non-positive-semidefinite
/// information block.
pub fn parameter_covariance_structured(
    jacobian: &[f64],
    layout: StructuredLayout,
) -> Result<Vec<f64>, FitError> {
    let global_count = layout.global_parameter_count;
    let point_count = layout.point_count;
    let latent_count = layout.latent_coordinate_count;
    let rows = point_count * (1 + latent_count);
    let columns = global_count + point_count * latent_count;
    if global_count == 0
        || point_count == 0
        || latent_count == 0
        || jacobian.len() != rows * columns
    {
        return Err(FitError::InvalidInput(
            "invalid structured Jacobian dimensions for covariance".to_owned(),
        ));
    }

    let mut profiled = vec![0.0; global_count * global_count];
    for row in 0..rows {
        for left in 0..global_count {
            let left_value = jacobian[row * columns + left];
            for right in 0..=left {
                profiled[left * global_count + right] +=
                    left_value * jacobian[row * columns + right];
            }
        }
    }
    for left in 0..global_count {
        for right in 0..left {
            profiled[right * global_count + left] = profiled[left * global_count + right];
        }
    }

    for point in 0..point_count {
        let mut local = vec![0.0; latent_count * latent_count];
        let mut cross = vec![0.0; latent_count * global_count];
        for component in 0..=latent_count {
            let row = component * point_count + point;
            for latent_left in 0..latent_count {
                let left_column = global_count + latent_left * point_count + point;
                let left_value = jacobian[row * columns + left_column];
                for global in 0..global_count {
                    cross[latent_left * global_count + global] +=
                        left_value * jacobian[row * columns + global];
                }
                for latent_right in 0..=latent_left {
                    let right_column = global_count + latent_right * point_count + point;
                    local[latent_left * latent_count + latent_right] +=
                        left_value * jacobian[row * columns + right_column];
                }
            }
        }
        for left in 0..latent_count {
            for right in 0..left {
                local[right * latent_count + left] = local[left * latent_count + right];
            }
        }
        let local_inverse = symmetric_pseudoinverse(&local, latent_count, 1 + latent_count)?;
        for left in 0..global_count {
            for right in 0..global_count {
                let correction = (0..latent_count)
                    .map(|local_left| {
                        (0..latent_count)
                            .map(|local_right| {
                                cross[local_left * global_count + left]
                                    * local_inverse[local_left * latent_count + local_right]
                                    * cross[local_right * global_count + right]
                            })
                            .sum::<f64>()
                    })
                    .sum::<f64>();
                profiled[left * global_count + right] -= correction;
            }
        }
    }
    symmetric_pseudoinverse(&profiled, global_count, rows.max(columns))
}

fn validate_problem(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
) -> Result<(), FitError> {
    if initial.is_empty() || lower.len() != initial.len() || upper.len() != initial.len() {
        return Err(FitError::InvalidInput(
            "initial values and bounds must have the same non-zero length".to_owned(),
        ));
    }
    for index in 0..initial.len() {
        if !initial[index].is_finite()
            || lower[index].is_nan()
            || upper[index].is_nan()
            || lower[index] >= upper[index]
            || initial[index] < lower[index]
            || initial[index] > upper[index]
        {
            return Err(FitError::InvalidInput(format!(
                "invalid initial value or bounds at parameter {index}"
            )));
        }
    }
    if !options.f_scale.is_finite() || options.f_scale <= 0.0 {
        return Err(FitError::InvalidInput(
            "f_scale must be positive and finite".to_owned(),
        ));
    }
    if matches!(options.max_evaluations, Some(0)) {
        return Err(FitError::InvalidInput(
            "max_evaluations must be positive".to_owned(),
        ));
    }
    for tolerance in [
        options.function_tolerance,
        options.gradient_tolerance,
        options.step_tolerance,
    ] {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(FitError::InvalidInput(
                "solver tolerances must be positive and finite".to_owned(),
            ));
        }
    }
    Ok(())
}

fn validate_residuals(residuals: &[f64], expected: Option<usize>) -> Result<(), FitError> {
    if residuals.is_empty()
        || expected.is_some_and(|count| count != residuals.len())
        || residuals.iter().any(|value| !value.is_finite())
    {
        return Err(FitError::Evaluation(
            "residual evaluation returned an invalid vector".to_owned(),
        ));
    }
    Ok(())
}

fn robust_cost(residuals: &[f64], loss: Loss, scale: f64) -> f64 {
    let scale_squared = scale * scale;
    0.5 * scale_squared
        * residuals
            .iter()
            .map(|value| loss.rho(value * value / scale_squared))
            .sum::<f64>()
}

fn robust_weights(residuals: &[f64], loss: Loss, scale: f64) -> Vec<f64> {
    let scale_squared = scale * scale;
    residuals
        .iter()
        .map(|value| {
            loss.derivative(value * value / scale_squared)
                .max(f64::EPSILON)
        })
        .collect()
}

fn robust_result_jacobian(
    mut jacobian: Vec<f64>,
    residuals: &[f64],
    columns: usize,
    loss: Loss,
    scale: f64,
) -> Vec<f64> {
    let scale_squared = scale * scale;
    for (row, residual) in residuals.iter().enumerate() {
        let squared = residual * residual / scale_squared;
        let jacobian_scale = (loss.derivative(squared)
            + 2.0 * loss.second_derivative(squared) * squared)
            .max(f64::EPSILON)
            .sqrt();
        for column in 0..columns {
            jacobian[row * columns + column] *= jacobian_scale;
        }
    }
    jacobian
}

fn finite_difference_for_layout<F>(
    parameters: &[f64],
    lower: &[f64],
    upper: &[f64],
    residuals: &[f64],
    layout: Option<StructuredLayout>,
    evaluate: &mut F,
) -> Result<(Vec<f64>, usize), FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    match layout {
        Some(structure) => {
            finite_difference_structured(parameters, lower, upper, residuals, structure, evaluate)
        }
        None => {
            finite_difference_jacobian(parameters, lower, upper, residuals, usize::MAX, evaluate)
        }
    }
}

fn finite_difference_column<F>(
    parameters: &[f64],
    lower: &[f64],
    upper: &[f64],
    residuals: &[f64],
    column: usize,
    evaluate: &mut F,
) -> Result<(Vec<f64>, f64, usize), FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    let rows = residuals.len();
    let base_step = EPSILON_SQRT * parameters[column].abs().max(1.0);
    let forward_step = base_step.min(upper[column] - parameters[column]);
    let backward_step = base_step.min(parameters[column] - lower[column]);
    if forward_step > 0.0 && backward_step > 0.0 {
        let mut plus = parameters.to_vec();
        plus[column] += forward_step;
        let plus_residuals = evaluate(&plus)?;
        validate_residuals(&plus_residuals, Some(rows))?;
        let mut minus = parameters.to_vec();
        minus[column] -= backward_step;
        let minus_residuals = evaluate(&minus)?;
        validate_residuals(&minus_residuals, Some(rows))?;
        Ok((
            plus_residuals
                .iter()
                .zip(minus_residuals)
                .map(|(plus, minus)| plus - minus)
                .collect(),
            forward_step + backward_step,
            2,
        ))
    } else if forward_step > 0.0 {
        let mut plus = parameters.to_vec();
        plus[column] += forward_step;
        let plus_residuals = evaluate(&plus)?;
        validate_residuals(&plus_residuals, Some(rows))?;
        Ok((
            plus_residuals
                .iter()
                .zip(residuals)
                .map(|(plus, base)| plus - base)
                .collect(),
            forward_step,
            1,
        ))
    } else if backward_step > 0.0 {
        let mut minus = parameters.to_vec();
        minus[column] -= backward_step;
        let minus_residuals = evaluate(&minus)?;
        validate_residuals(&minus_residuals, Some(rows))?;
        Ok((
            residuals
                .iter()
                .zip(minus_residuals)
                .map(|(base, minus)| base - minus)
                .collect(),
            backward_step,
            1,
        ))
    } else {
        Err(FitError::InvalidInput(format!(
            "parameter {column} is fixed by equal bounds"
        )))
    }
}

fn finite_difference_jacobian<F>(
    parameters: &[f64],
    lower: &[f64],
    upper: &[f64],
    residuals: &[f64],
    evaluation_budget: usize,
    evaluate: &mut F,
) -> Result<(Vec<f64>, usize), FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    let rows = residuals.len();
    let columns = parameters.len();
    let mut jacobian = vec![0.0; rows * columns];
    let mut used = 0;
    for column in 0..columns {
        if used >= evaluation_budget {
            return Err(FitError::Evaluation(
                "function evaluation budget exhausted while calculating Jacobian".to_owned(),
            ));
        }
        let (difference, denominator, column_used) =
            finite_difference_column(parameters, lower, upper, residuals, column, evaluate)?;
        used += column_used;
        if used > evaluation_budget {
            return Err(FitError::Evaluation(
                "function evaluation budget exhausted while calculating Jacobian".to_owned(),
            ));
        }
        for row in 0..rows {
            jacobian[row * columns + column] = difference[row] / denominator;
        }
    }
    Ok((jacobian, used))
}

fn finite_difference_structured<F>(
    parameters: &[f64],
    lower: &[f64],
    upper: &[f64],
    residuals: &[f64],
    layout: StructuredLayout,
    evaluate: &mut F,
) -> Result<(Vec<f64>, usize), FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    let rows = residuals.len();
    let columns = parameters.len();
    let mut jacobian = vec![0.0; rows * columns];
    let mut used = 0;
    for column in 0..layout.global_parameter_count {
        let (difference, denominator, column_used) =
            finite_difference_column(parameters, lower, upper, residuals, column, evaluate)?;
        used += column_used;
        for row in 0..rows {
            jacobian[row * columns + column] = difference[row] / denominator;
        }
    }

    for coordinate in 0..layout.latent_coordinate_count {
        let mut plus = parameters.to_vec();
        let mut minus = parameters.to_vec();
        let mut denominators = vec![0.0; layout.point_count];
        for (point, denominator) in denominators.iter_mut().enumerate() {
            let column = layout.global_parameter_count + coordinate * layout.point_count + point;
            let base_step = EPSILON_SQRT * parameters[column].abs().max(1.0);
            let forward_step = base_step.min(upper[column] - parameters[column]);
            let backward_step = base_step.min(parameters[column] - lower[column]);
            *denominator = forward_step + backward_step;
            if *denominator <= 0.0 {
                return Err(FitError::InvalidInput(format!(
                    "parameter {column} is fixed by equal bounds"
                )));
            }
            plus[column] += forward_step;
            minus[column] -= backward_step;
        }
        let plus_residuals = evaluate(&plus)?;
        let minus_residuals = evaluate(&minus)?;
        used += 2;
        validate_residuals(&plus_residuals, Some(rows))?;
        validate_residuals(&minus_residuals, Some(rows))?;
        for (point, denominator) in denominators.iter().enumerate() {
            let column = layout.global_parameter_count + coordinate * layout.point_count + point;
            for component in 0..=layout.latent_coordinate_count {
                let row = component * layout.point_count + point;
                jacobian[row * columns + column] =
                    (plus_residuals[row] - minus_residuals[row]) / denominator;
            }
        }
    }
    Ok((jacobian, used))
}

#[derive(Clone, Debug)]
enum NormalSystem {
    Dense {
        matrix: Vec<f64>,
        size: usize,
    },
    Structured {
        global: Vec<f64>,
        cross: Vec<f64>,
        local: Vec<f64>,
        layout: StructuredLayout,
    },
}

fn normal_equations_for_layout(
    jacobian: &[f64],
    residuals: &[f64],
    weights: &[f64],
    rows: usize,
    columns: usize,
    layout: Option<StructuredLayout>,
) -> (NormalSystem, Vec<f64>) {
    if let Some(structure) = layout {
        structured_normal_equations(jacobian, residuals, weights, columns, structure)
    } else {
        let (matrix, gradient) = normal_equations(jacobian, residuals, weights, rows, columns);
        (
            NormalSystem::Dense {
                matrix,
                size: columns,
            },
            gradient,
        )
    }
}

fn structured_normal_equations(
    jacobian: &[f64],
    residuals: &[f64],
    weights: &[f64],
    columns: usize,
    layout: StructuredLayout,
) -> (NormalSystem, Vec<f64>) {
    let global_count = layout.global_parameter_count;
    let point_count = layout.point_count;
    let latent_count = layout.latent_coordinate_count;
    let mut global = vec![0.0; global_count * global_count];
    let mut cross = vec![0.0; point_count * global_count * latent_count];
    let mut local = vec![0.0; point_count * latent_count * latent_count];
    let mut gradient = vec![0.0; columns];
    for row in 0..residuals.len() {
        let point = row % point_count;
        let weight = weights[row];
        for left in 0..global_count {
            let left_value = jacobian[row * columns + left];
            gradient[left] += weight * left_value * residuals[row];
            for right in 0..=left {
                global[left * global_count + right] +=
                    weight * left_value * jacobian[row * columns + right];
            }
            for latent in 0..latent_count {
                let latent_column = global_count + latent * point_count + point;
                cross[(point * global_count + left) * latent_count + latent] +=
                    weight * left_value * jacobian[row * columns + latent_column];
            }
        }
        for left in 0..latent_count {
            let left_column = global_count + left * point_count + point;
            let left_value = jacobian[row * columns + left_column];
            gradient[left_column] += weight * left_value * residuals[row];
            for right in 0..=left {
                let right_column = global_count + right * point_count + point;
                local[(point * latent_count + left) * latent_count + right] +=
                    weight * left_value * jacobian[row * columns + right_column];
            }
        }
    }
    for left in 0..global_count {
        for right in 0..left {
            global[right * global_count + left] = global[left * global_count + right];
        }
    }
    for point in 0..point_count {
        for left in 0..latent_count {
            for right in 0..left {
                local[(point * latent_count + right) * latent_count + left] =
                    local[(point * latent_count + left) * latent_count + right];
            }
        }
    }
    (
        NormalSystem::Structured {
            global,
            cross,
            local,
            layout,
        },
        gradient,
    )
}

fn solve_normal_step(
    system: &NormalSystem,
    gradient: &[f64],
    damping: f64,
) -> Result<Vec<f64>, FitError> {
    match system {
        NormalSystem::Dense { matrix, size } => {
            let mut damped = matrix.clone();
            for index in 0..*size {
                let diagonal = matrix[index * size + index].abs();
                damped[index * size + index] +=
                    damping * if diagonal > 0.0 { diagonal } else { 1.0 };
            }
            let right_hand_side: Vec<_> = gradient.iter().map(|value| -value).collect();
            solve_scaled_linear_system(&damped, &right_hand_side, *size)
        }
        NormalSystem::Structured {
            global,
            cross,
            local,
            layout,
        } => solve_structured_normal_step(global, cross, local, gradient, damping, *layout),
    }
}

fn solve_structured_normal_step(
    global: &[f64],
    cross: &[f64],
    local: &[f64],
    gradient: &[f64],
    damping: f64,
    layout: StructuredLayout,
) -> Result<Vec<f64>, FitError> {
    let global_count = layout.global_parameter_count;
    let point_count = layout.point_count;
    let latent_count = layout.latent_coordinate_count;
    let mut schur = global.to_vec();
    let mut global_rhs: Vec<_> = gradient[..global_count]
        .iter()
        .map(|value| -value)
        .collect();
    for global_index in 0..global_count {
        let diagonal = global[global_index * global_count + global_index].abs();
        schur[global_index * global_count + global_index] +=
            damping * if diagonal > 0.0 { diagonal } else { 1.0 };
    }
    let mut damped_locals = Vec::with_capacity(point_count);
    for point in 0..point_count {
        let start = point * latent_count * latent_count;
        let mut local_matrix = local[start..start + latent_count * latent_count].to_vec();
        for latent in 0..latent_count {
            let diagonal = local_matrix[latent * latent_count + latent].abs();
            local_matrix[latent * latent_count + latent] +=
                damping * if diagonal > 0.0 { diagonal } else { 1.0 };
        }
        let local_gradient: Vec<_> = (0..latent_count)
            .map(|latent| gradient[global_count + latent * point_count + point])
            .collect();
        let solved_gradient =
            solve_scaled_linear_system(&local_matrix, &local_gradient, latent_count)?;
        for global_index in 0..global_count {
            global_rhs[global_index] += (0..latent_count)
                .map(|latent| {
                    cross[(point * global_count + global_index) * latent_count + latent]
                        * solved_gradient[latent]
                })
                .sum::<f64>();
        }
        for right_global in 0..global_count {
            let right_hand_side: Vec<_> = (0..latent_count)
                .map(|latent| cross[(point * global_count + right_global) * latent_count + latent])
                .collect();
            let solved_cross =
                solve_scaled_linear_system(&local_matrix, &right_hand_side, latent_count)?;
            for left_global in 0..global_count {
                schur[left_global * global_count + right_global] -= (0..latent_count)
                    .map(|latent| {
                        cross[(point * global_count + left_global) * latent_count + latent]
                            * solved_cross[latent]
                    })
                    .sum::<f64>();
            }
        }
        damped_locals.push(local_matrix);
    }
    let global_step = solve_scaled_linear_system(&schur, &global_rhs, global_count)?;
    let mut step = vec![0.0; global_count + point_count * latent_count];
    step[..global_count].copy_from_slice(&global_step);
    for point in 0..point_count {
        let local_rhs: Vec<_> = (0..latent_count)
            .map(|latent| {
                -gradient[global_count + latent * point_count + point]
                    - (0..global_count)
                        .map(|global_index| {
                            cross[(point * global_count + global_index) * latent_count + latent]
                                * global_step[global_index]
                        })
                        .sum::<f64>()
            })
            .collect();
        let local_step =
            solve_scaled_linear_system(&damped_locals[point], &local_rhs, latent_count)?;
        for latent in 0..latent_count {
            step[global_count + latent * point_count + point] = local_step[latent];
        }
    }
    Ok(step)
}

fn normal_equations(
    jacobian: &[f64],
    residuals: &[f64],
    weights: &[f64],
    rows: usize,
    columns: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut normal = vec![0.0; columns * columns];
    let mut gradient = vec![0.0; columns];
    for row in 0..rows {
        let weight = weights[row];
        for left in 0..columns {
            let left_value = jacobian[row * columns + left];
            gradient[left] += weight * left_value * residuals[row];
            for right in 0..=left {
                normal[left * columns + right] +=
                    weight * left_value * jacobian[row * columns + right];
            }
        }
    }
    for left in 0..columns {
        for right in 0..left {
            normal[right * columns + left] = normal[left * columns + right];
        }
    }
    (normal, gradient)
}

fn projected_gradient_norm(
    gradient: &[f64],
    parameters: &[f64],
    lower: &[f64],
    upper: &[f64],
) -> f64 {
    gradient
        .iter()
        .enumerate()
        .map(|(index, &value)| {
            if (parameters[index] <= lower[index] && value > 0.0)
                || (parameters[index] >= upper[index] && value < 0.0)
            {
                0.0
            } else {
                value.abs()
            }
        })
        .fold(0.0, f64::max)
}

fn scaled_norm_difference(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| {
            let scale = right.abs().max(1.0);
            ((left - right) / scale).powi(2)
        })
        .sum::<f64>()
        .sqrt()
}

fn solve_linear_system(
    matrix: &[f64],
    right_hand_side: &[f64],
    size: usize,
) -> Result<Vec<f64>, FitError> {
    if matrix.len() != size * size || right_hand_side.len() != size {
        return Err(FitError::InvalidInput(
            "invalid linear system dimensions".to_owned(),
        ));
    }
    let mut matrix = matrix.to_vec();
    let mut result = right_hand_side.to_vec();
    for pivot in 0..size {
        let mut selected = pivot;
        let mut selected_value = matrix[pivot * size + pivot].abs();
        for row in pivot + 1..size {
            let value = matrix[row * size + pivot].abs();
            if value > selected_value {
                selected = row;
                selected_value = value;
            }
        }
        let scale = matrix
            .iter()
            .map(|value| value.abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let dimension_scale = f64::from(u32::try_from(size).unwrap_or(u32::MAX));
        if selected_value <= f64::EPSILON * dimension_scale * scale {
            return Err(FitError::SingularSystem);
        }
        if selected != pivot {
            for column in 0..size {
                matrix.swap(pivot * size + column, selected * size + column);
            }
            result.swap(pivot, selected);
        }
        for row in pivot + 1..size {
            let factor = matrix[row * size + pivot] / matrix[pivot * size + pivot];
            matrix[row * size + pivot] = 0.0;
            for column in pivot + 1..size {
                matrix[row * size + column] -= factor * matrix[pivot * size + column];
            }
            result[row] -= factor * result[pivot];
        }
    }
    for row in (0..size).rev() {
        for column in row + 1..size {
            result[row] -= matrix[row * size + column] * result[column];
        }
        result[row] /= matrix[row * size + row];
    }
    Ok(result)
}

fn solve_scaled_linear_system(
    matrix: &[f64],
    right_hand_side: &[f64],
    size: usize,
) -> Result<Vec<f64>, FitError> {
    if matrix.len() != size * size || right_hand_side.len() != size {
        return Err(FitError::InvalidInput(
            "invalid scaled linear system dimensions".to_owned(),
        ));
    }
    let scales = (0..size)
        .map(|index| {
            let diagonal = matrix[index * size + index].abs();
            if diagonal > 0.0 {
                diagonal.sqrt()
            } else {
                1.0
            }
        })
        .collect::<Vec<_>>();
    let mut scaled_matrix = vec![0.0; size * size];
    for row in 0..size {
        for column in 0..size {
            scaled_matrix[row * size + column] =
                matrix[row * size + column] / (scales[row] * scales[column]);
        }
    }
    let scaled_rhs = right_hand_side
        .iter()
        .zip(&scales)
        .map(|(value, scale)| value / scale)
        .collect::<Vec<_>>();
    let scaled_solution = solve_linear_system(&scaled_matrix, &scaled_rhs, size)?;
    Ok(scaled_solution
        .into_iter()
        .zip(scales)
        .map(|(value, scale)| value / scale)
        .collect())
}

fn schur_complement(
    information: &[f64],
    size: usize,
    parameter_count: usize,
) -> Result<Vec<f64>, FitError> {
    let latent_count = size - parameter_count;
    let mut latent = vec![0.0; latent_count * latent_count];
    let mut cross = vec![0.0; latent_count * parameter_count];
    let mut parameter = vec![0.0; parameter_count * parameter_count];
    for row in 0..parameter_count {
        for column in 0..parameter_count {
            parameter[row * parameter_count + column] = information[row * size + column];
        }
    }
    for row in 0..latent_count {
        for column in 0..latent_count {
            latent[row * latent_count + column] =
                information[(parameter_count + row) * size + parameter_count + column];
        }
        for column in 0..parameter_count {
            cross[row * parameter_count + column] =
                information[(parameter_count + row) * size + column];
        }
    }
    for column in 0..parameter_count {
        let right_hand_side: Vec<_> = (0..latent_count)
            .map(|row| cross[row * parameter_count + column])
            .collect();
        let solved = solve_scaled_linear_system(&latent, &right_hand_side, latent_count)?;
        for left in 0..parameter_count {
            let correction = (0..latent_count)
                .map(|row| cross[row * parameter_count + left] * solved[row])
                .sum::<f64>();
            parameter[left * parameter_count + column] -= correction;
        }
    }
    Ok(parameter)
}

#[allow(clippy::too_many_lines)]
fn symmetric_pseudoinverse(
    matrix: &[f64],
    size: usize,
    effective_dimension: usize,
) -> Result<Vec<f64>, FitError> {
    if size == 0 || matrix.len() != size * size || matrix.iter().any(|value| !value.is_finite()) {
        return Err(FitError::InvalidInput(
            "invalid symmetric matrix for pseudoinverse".to_owned(),
        ));
    }
    let mut diagonalized = matrix.to_vec();
    for row in 0..size {
        for column in 0..row {
            let value =
                0.5 * (diagonalized[row * size + column] + diagonalized[column * size + row]);
            diagonalized[row * size + column] = value;
            diagonalized[column * size + row] = value;
        }
    }
    let mut eigenvectors = vec![0.0; size * size];
    for index in 0..size {
        eigenvectors[index * size + index] = 1.0;
    }
    let scale = diagonalized
        .iter()
        .map(|value| value.abs())
        .fold(0.0, f64::max);
    let convergence_tolerance =
        f64::EPSILON * f64::from(u32::try_from(size).unwrap_or(u32::MAX)) * scale;
    for _ in 0..(64 * size * size).max(1) {
        let mut pivot_row = 0;
        let mut pivot_column = 0;
        let mut largest = 0.0;
        for row in 0..size {
            for column in row + 1..size {
                let value = diagonalized[row * size + column].abs();
                if value > largest {
                    largest = value;
                    pivot_row = row;
                    pivot_column = column;
                }
            }
        }
        if largest <= convergence_tolerance {
            break;
        }
        let p = pivot_row;
        let q = pivot_column;
        let off_diagonal = diagonalized[p * size + q];
        let tau = (diagonalized[q * size + q] - diagonalized[p * size + p]) / (2.0 * off_diagonal);
        let tangent = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let cosine = 1.0 / (1.0 + tangent * tangent).sqrt();
        let rotation_sine = tangent * cosine;
        for index in 0..size {
            if index != p && index != q {
                let value_p = diagonalized[index * size + p];
                let value_q = diagonalized[index * size + q];
                let rotated_p = cosine * value_p - rotation_sine * value_q;
                let rotated_q = rotation_sine * value_p + cosine * value_q;
                diagonalized[index * size + p] = rotated_p;
                diagonalized[p * size + index] = rotated_p;
                diagonalized[index * size + q] = rotated_q;
                diagonalized[q * size + index] = rotated_q;
            }
        }
        let diagonal_p = diagonalized[p * size + p];
        let diagonal_q = diagonalized[q * size + q];
        diagonalized[p * size + p] = diagonal_p - tangent * off_diagonal;
        diagonalized[q * size + q] = diagonal_q + tangent * off_diagonal;
        diagonalized[p * size + q] = 0.0;
        diagonalized[q * size + p] = 0.0;
        for row in 0..size {
            let value_p = eigenvectors[row * size + p];
            let value_q = eigenvectors[row * size + q];
            eigenvectors[row * size + p] = cosine * value_p - rotation_sine * value_q;
            eigenvectors[row * size + q] = rotation_sine * value_p + cosine * value_q;
        }
    }

    let largest_eigenvalue = (0..size)
        .map(|index| diagonalized[index * size + index].abs())
        .fold(0.0, f64::max);
    let rank_tolerance = f64::EPSILON
        * f64::from(u32::try_from(effective_dimension).unwrap_or(u32::MAX))
        * largest_eigenvalue;
    let mut inverse = vec![0.0; size * size];
    for eigenvalue_index in 0..size {
        let eigenvalue = diagonalized[eigenvalue_index * size + eigenvalue_index];
        if eigenvalue < -rank_tolerance {
            return Err(FitError::InvalidInput(
                "information matrix must be positive semidefinite".to_owned(),
            ));
        }
        if eigenvalue <= rank_tolerance {
            continue;
        }
        for row in 0..size {
            for column in 0..size {
                inverse[row * size + column] += eigenvectors[row * size + eigenvalue_index]
                    * eigenvectors[column * size + eigenvalue_index]
                    / eigenvalue;
            }
        }
    }
    Ok(inverse)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::error::Error;

    use crate::EosError;

    use super::*;

    #[test]
    fn eos_failures_retain_their_code_and_source() {
        let source = EosError::InvalidState {
            name: "volume",
            reason: "must be positive and finite",
        };
        let error = FitError::from(source.clone());

        assert_eq!(error.kind(), FitErrorKind::EosEvaluation);
        assert_eq!(error.code(), "eos.invalid_state");
        assert_eq!(
            error.source().and_then(|value| value.downcast_ref()),
            Some(&source)
        );
    }
    #[test]
    fn bounded_solver_recovers_rosenbrock_minimum() {
        let result = least_squares(
            &[-1.2, 1.0],
            &[-2.0, -1.0],
            &[2.0, 3.0],
            SolverOptions {
                max_evaluations: Some(2_000),
                ..SolverOptions::default()
            },
            |parameters| {
                Ok(vec![
                    10.0 * (parameters[1] - parameters[0] * parameters[0]),
                    1.0 - parameters[0],
                ])
            },
        )
        .unwrap();

        assert!(result.success, "{}", result.message);
        assert!((result.parameters[0] - 1.0).abs() < 1.0e-5);
        assert!((result.parameters[1] - 1.0).abs() < 1.0e-5);
    }

    #[test]
    fn scaled_solver_handles_jacobian_columns_with_different_units() {
        let result = least_squares(
            &[0.0, 0.0],
            &[-10.0, -10.0],
            &[10.0, 10.0],
            SolverOptions {
                max_evaluations: Some(1_000),
                ..SolverOptions::default()
            },
            |parameters| {
                Ok(vec![
                    1.0e-9 * (parameters[0] - 2.0),
                    1.0e9 * (parameters[1] - 3.0),
                ])
            },
        )
        .unwrap();

        assert!(result.success, "{}", result.message);
        assert!((result.parameters[0] - 2.0).abs() < 1.0e-8, "{result:?}");
        assert!((result.parameters[1] - 3.0).abs() < 1.0e-8, "{result:?}");
        assert!(result.cost < 1.0e-10);
    }

    #[test]
    fn upper_bound_is_active_at_constrained_solution() {
        let result = least_squares(
            &[0.0],
            &[-1.0],
            &[2.0],
            SolverOptions::default(),
            |parameters| Ok(vec![parameters[0] - 5.0]),
        )
        .unwrap();

        assert!(result.success);
        assert!((result.parameters[0] - 2.0).abs() < 1.0e-10);
    }

    #[test]
    fn robust_loss_reduces_outlier_influence() {
        let observations = [1.0, 1.1, 0.9, 20.0];
        let solve = |loss| {
            least_squares(
                &[0.0],
                &[f64::NEG_INFINITY],
                &[f64::INFINITY],
                SolverOptions {
                    loss,
                    max_evaluations: Some(1_000),
                    ..SolverOptions::default()
                },
                |parameters| {
                    Ok(observations
                        .iter()
                        .map(|observation| parameters[0] - observation)
                        .collect())
                },
            )
            .unwrap()
        };
        let linear = solve(Loss::Linear);
        let robust = solve(Loss::SoftL1);

        assert!((robust.parameters[0] - 1.0).abs() < (linear.parameters[0] - 1.0).abs());
    }

    #[test]
    fn covariance_matches_linear_closed_form() {
        let jacobian = [1.0, 2.0, 3.0];
        let covariance = parameter_covariance(&jacobian, 3, 1, 1).unwrap();
        assert!((covariance[0] - 1.0 / 14.0).abs() < 1.0e-14);

        let scaled_jacobian = [1.0e-9, 2.0e-9, 3.0e-9];
        let scaled_covariance = parameter_covariance(&scaled_jacobian, 3, 1, 1).unwrap();
        assert!((scaled_covariance[0] - 1.0e18 / 14.0).abs() < 1.0e4);
    }

    #[test]
    fn covariance_uses_moore_penrose_inverse_for_rank_loss() {
        let jacobian = [1.0, 2.0, 1.0, 2.0, 1.0, 2.0];
        let covariance = parameter_covariance(&jacobian, 3, 2, 2).unwrap();
        let expected = [1.0 / 75.0, 2.0 / 75.0, 2.0 / 75.0, 4.0 / 75.0];
        for (actual, expected) in covariance.iter().zip(expected) {
            assert!((actual - expected).abs() < 1.0e-12, "{covariance:?}");
        }
    }

    #[test]
    fn linear_uncertainty_matches_dense_matrix_product() {
        let result = propagate_linear_uncertainty(
            &[1.0, 2.0, -1.0, 0.5],
            2,
            2,
            &[4.0, 0.5, 0.5, 1.0],
            &[0.25, 0.0],
            true,
        )
        .unwrap();

        assert_eq!(result.variance, vec![10.25, 3.75]);
        assert_eq!(result.covariance.unwrap(), vec![10.25, -3.75, -3.75, 3.75]);
    }

    #[test]
    fn linear_uncertainty_supports_state_only_propagation() {
        let result = propagate_linear_uncertainty(&[], 2, 0, &[], &[0.25, 1.0], true).unwrap();

        assert_eq!(result.variance, vec![0.25, 1.0]);
        assert_eq!(result.covariance.unwrap(), vec![0.25, 0.0, 0.0, 1.0]);
    }

    #[test]
    fn linear_uncertainty_rejects_invalid_variances() {
        assert!(propagate_linear_uncertainty(&[1.0], 1, 1, &[-1.0], &[0.0], true).is_err());
        assert!(propagate_linear_uncertainty(&[1.0], 1, 1, &[1.0], &[-1.0], true).is_err());
    }

    #[test]
    fn monte_carlo_summary_uses_sample_statistics_and_linear_quantiles() {
        let summary = summarize_monte_carlo(
            &[1.0, 10.0, 2.0, 20.0, 3.0, 30.0, 4.0, 40.0],
            4,
            2,
            0.5,
            true,
        )
        .unwrap();

        let expected_error = (5.0_f64 / 3.0).sqrt();
        assert!((summary.standard_error[0] - expected_error).abs() < 1.0e-14);
        assert!((summary.standard_error[1] - 10.0 * expected_error).abs() < 1.0e-13);
        assert_eq!(summary.lower, vec![1.75, 17.5]);
        assert_eq!(summary.upper, vec![3.25, 32.5]);
        assert_eq!(
            summary.covariance.unwrap(),
            vec![5.0 / 3.0, 50.0 / 3.0, 50.0 / 3.0, 500.0 / 3.0]
        );
    }

    #[test]
    fn model_jacobian_and_linear_propagation_match_closed_form() {
        let parameters = [2.0, 3.0];
        let covariance = [0.04, 0.01, 0.01, 0.09];
        let evaluate =
            |values: &[f64]| Ok(vec![values[0] + 2.0 * values[1], values[0] * values[1]]);
        let result = propagate_model_uncertainty(
            &parameters,
            &covariance,
            &[0.0, 0.0],
            1.0e-6,
            true,
            evaluate,
        )
        .unwrap();

        assert_eq!(result.model.nominal, vec![8.0, 6.0]);
        let expected_jacobian = [1.0, 2.0, 3.0, 2.0];
        for (actual, expected) in result.model.jacobian.iter().zip(expected_jacobian) {
            assert!((actual - expected).abs() < 1.0e-8);
        }
        assert!((result.propagation.variance[0] - 0.44).abs() < 1.0e-8);
        assert!((result.propagation.variance[1] - 0.84).abs() < 1.0e-8);
        assert!(result.propagation.covariance.is_some());
    }

    #[test]
    fn native_monte_carlo_is_seeded_and_combines_parameter_and_state_variance() {
        let options = MonteCarloOptions {
            sample_count: 20_000,
            max_attempts: 20_000,
            confidence: 0.95,
            full_covariance: true,
            seed: 42,
        };
        let evaluate = |parameters: &[f64], states: &[f64]| Ok(vec![parameters[0] + states[0]]);
        let first =
            monte_carlo_model_uncertainty(&[1.0], &[4.0], &[2.0], &[3.0], options, evaluate)
                .unwrap();
        let second =
            monte_carlo_model_uncertainty(&[1.0], &[4.0], &[2.0], &[3.0], options, evaluate)
                .unwrap();

        assert_eq!(first, second);
        assert_eq!(first.nominal, vec![3.0]);
        assert_eq!(first.attempted_samples, options.sample_count);
        assert_eq!(first.rejected_fraction.to_bits(), 0.0_f64.to_bits());
        assert!((first.summary.standard_error[0] - 13.0_f64.sqrt()).abs() < 0.08);
    }

    #[test]
    fn native_monte_carlo_tracks_rejected_model_evaluations() {
        let result = monte_carlo_model_uncertainty(
            &[0.0],
            &[1.0],
            &[],
            &[],
            MonteCarloOptions {
                sample_count: 1_000,
                max_attempts: 4_000,
                seed: 7,
                ..MonteCarloOptions::default()
            },
            |parameters, _| {
                if parameters[0] < 0.0 {
                    Err(FitError::Evaluation("negative sample".to_owned()))
                } else {
                    Ok(vec![parameters[0]])
                }
            },
        )
        .unwrap();

        assert!(result.rejected_fraction > 0.4);
        assert!(result.rejected_fraction < 0.6);
        assert!(result.summary.lower[0] >= 0.0);
    }

    #[test]
    fn native_monte_carlo_accepts_singular_positive_semidefinite_covariance() {
        let result = monte_carlo_model_uncertainty(
            &[1.0, 1.0],
            &[1.0, 1.0, 1.0, 1.0],
            &[],
            &[],
            MonteCarloOptions {
                sample_count: 100,
                max_attempts: 100,
                seed: 9,
                ..MonteCarloOptions::default()
            },
            |parameters, _| Ok(vec![parameters[0] - parameters[1]]),
        )
        .unwrap();

        assert!(result.summary.standard_error[0] < 1.0e-14);
    }

    #[test]
    fn structured_latent_solver_matches_dense_with_fewer_evaluations() {
        let point_count = 80;
        let observed_coordinates: Vec<_> = (0..point_count)
            .map(|index| {
                let value = f64::from(u32::try_from(index).unwrap());
                0.5 + value / 100.0 + 0.002 * value.sin()
            })
            .collect();
        let observations: Vec<_> = (0..point_count)
            .map(|index| {
                let value = f64::from(u32::try_from(index).unwrap());
                2.0 * (0.5 + value / 100.0)
            })
            .collect();
        let mut initial = vec![1.8];
        initial.extend(&observed_coordinates);
        let mut lower = vec![0.1];
        lower.extend(vec![f64::MIN_POSITIVE; point_count]);
        let upper = vec![f64::INFINITY; 1 + point_count];

        let run = |structured: bool| {
            let calls = Cell::new(0_usize);
            let residual = |parameters: &[f64]| {
                calls.set(calls.get() + 1);
                let model = parameters[0];
                let latent = &parameters[1..];
                let mut values = Vec::with_capacity(2 * point_count);
                values.extend(
                    latent
                        .iter()
                        .zip(&observations)
                        .map(|(coordinate, observed)| (model * coordinate - observed) / 0.01),
                );
                values.extend(
                    latent
                        .iter()
                        .zip(&observed_coordinates)
                        .map(|(coordinate, observed)| (coordinate - observed) / 0.005),
                );
                Ok(values)
            };
            let options = SolverOptions {
                max_evaluations: Some(500),
                ..SolverOptions::default()
            };
            let result = if structured {
                least_squares_structured(
                    &initial,
                    &lower,
                    &upper,
                    options,
                    StructuredLayout {
                        global_parameter_count: 1,
                        point_count,
                        latent_coordinate_count: 1,
                    },
                    residual,
                )
            } else {
                least_squares(&initial, &lower, &upper, options, residual)
            }
            .unwrap();
            (result, calls.get())
        };

        let (dense, dense_calls) = run(false);
        let (structured, structured_calls) = run(true);
        assert!(dense.success && structured.success);
        assert!((dense.parameters[0] - structured.parameters[0]).abs() < 1.0e-8);
        assert!(structured_calls * 5 < dense_calls);

        let layout = StructuredLayout {
            global_parameter_count: 1,
            point_count,
            latent_coordinate_count: 1,
        };
        let structured_covariance =
            parameter_covariance_structured(&structured.jacobian, layout).unwrap();
        let dense_covariance = parameter_covariance(
            &structured.jacobian,
            structured.residual_count,
            structured.parameters.len(),
            1,
        )
        .unwrap();
        assert!((structured_covariance[0] - dense_covariance[0]).abs() < 1.0e-10);
    }
}
