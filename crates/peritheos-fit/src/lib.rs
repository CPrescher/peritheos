//! Numerical fitting primitives used by Peritheos.
//!
//! The solver deliberately owns its numerical conventions instead of adapting
//! a partial third-party API. It supports box constraints, the robust losses
//! exposed by the Python API, deterministic finite-difference Jacobians, and
//! the diagnostics needed to preserve `scipy.optimize.least_squares`-style
//! results.

mod eos;

pub use eos::{
    fit_isothermal_eos, fit_thermal_eos, fit_thermal_eos_by, EosFitResult, IsothermalObservations,
    ThermalObservations,
};

use std::error::Error;
use std::fmt::{Display, Formatter};

const EPSILON_SQRT: f64 = 1.490_116_119_384_765_6e-8;
const DEFAULT_TOLERANCE: f64 = 1.0e-8;

/// Error returned by fitting and covariance kernels.
#[derive(Clone, Debug, PartialEq)]
pub enum FitError {
    InvalidInput(String),
    Evaluation(String),
    SingularSystem,
}

impl Display for FitError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(message) | Self::Evaluation(message) => formatter.write_str(message),
            Self::SingularSystem => formatter.write_str("linear system is singular"),
        }
    }
}

impl Error for FitError {}

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
        || parameter_count == 0
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
        let previous_cost = cost;
        for _ in 0..24 {
            let Ok(step) = solve_normal_step(&normal, &gradient, damping) else {
                damping *= 10.0;
                continue;
            };
            let mut candidate = parameters.clone();
            for index in 0..parameter_count {
                candidate[index] =
                    (parameters[index] + step[index]).clamp(lower[index], upper[index]);
            }
            accepted_step_norm = scaled_norm_difference(&candidate, &parameters);
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
            status = 3;
            "`xtol` termination condition is satisfied.".clone_into(&mut message);
            success = true;
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
    regularized_inverse(&profiled, parameter_count)
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
                damped[index * size + index] += damping * matrix[index * size + index].max(1.0);
            }
            let right_hand_side: Vec<_> = gradient.iter().map(|value| -value).collect();
            solve_linear_system(&damped, &right_hand_side, *size)
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
        schur[global_index * global_count + global_index] +=
            damping * global[global_index * global_count + global_index].max(1.0);
    }
    let mut damped_locals = Vec::with_capacity(point_count);
    for point in 0..point_count {
        let start = point * latent_count * latent_count;
        let mut local_matrix = local[start..start + latent_count * latent_count].to_vec();
        for latent in 0..latent_count {
            local_matrix[latent * latent_count + latent] +=
                damping * local_matrix[latent * latent_count + latent].max(1.0);
        }
        let local_gradient: Vec<_> = (0..latent_count)
            .map(|latent| gradient[global_count + latent * point_count + point])
            .collect();
        let solved_gradient = solve_linear_system(&local_matrix, &local_gradient, latent_count)?;
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
            let solved_cross = solve_linear_system(&local_matrix, &right_hand_side, latent_count)?;
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
    let global_step = solve_linear_system(&schur, &global_rhs, global_count)?;
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
        let local_step = solve_linear_system(&damped_locals[point], &local_rhs, latent_count)?;
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
        let solved = solve_linear_system(&latent, &right_hand_side, latent_count)?;
        for left in 0..parameter_count {
            let correction = (0..latent_count)
                .map(|row| cross[row * parameter_count + left] * solved[row])
                .sum::<f64>();
            parameter[left * parameter_count + column] -= correction;
        }
    }
    Ok(parameter)
}

fn regularized_inverse(matrix: &[f64], size: usize) -> Result<Vec<f64>, FitError> {
    let scale = (0..size)
        .map(|index| matrix[index * size + index].abs())
        .fold(0.0, f64::max)
        .max(1.0);
    for regularization in [0.0, 1.0e-14, 1.0e-12, 1.0e-10, 1.0e-8] {
        let mut candidate = matrix.to_vec();
        for index in 0..size {
            candidate[index * size + index] += regularization * scale;
        }
        let mut inverse = vec![0.0; size * size];
        let mut valid = true;
        for column in 0..size {
            let mut unit = vec![0.0; size];
            unit[column] = 1.0;
            if let Ok(solution) = solve_linear_system(&candidate, &unit, size) {
                for row in 0..size {
                    inverse[row * size + column] = solution[row];
                }
            } else {
                valid = false;
                break;
            }
        }
        if valid {
            return Ok(inverse);
        }
    }
    Err(FitError::SingularSystem)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

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
    }
}
