//! Numerical fitting primitives used by Peritheos.
//!
//! The solver deliberately owns its numerical conventions instead of adapting
//! a partial third-party API. It supports box constraints, the robust losses
//! exposed by the Python API, deterministic finite-difference Jacobians, and
//! the diagnostics needed to preserve [`scipy.optimize.least_squares`]-style
//! results.

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

/// Solve a bounded nonlinear least-squares problem.
///
/// # Errors
///
/// Returns an error for inconsistent inputs, invalid residual evaluations, or
/// a linear algebra failure that cannot be regularized.
#[allow(clippy::too_many_lines)]
pub fn least_squares<F>(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
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
    let mut cost = robust_cost(&residuals, options.loss, options.f_scale);
    let mut damping = 1.0e-3;
    let mut status = 0;
    let mut message = "The maximum number of function evaluations is exceeded.".to_owned();
    let mut success = false;
    let mut final_jacobian = Vec::new();
    let mut optimality = f64::INFINITY;

    while evaluations < maximum_evaluations {
        let (jacobian, _used) = finite_difference_jacobian(
            &parameters,
            lower,
            upper,
            &residuals,
            usize::MAX,
            &mut evaluate,
        )?;
        jacobian_evaluations += 1;
        final_jacobian.clone_from(&jacobian);
        let weights = robust_weights(&residuals, options.loss, options.f_scale);
        let (normal, gradient) = normal_equations(
            &jacobian,
            &residuals,
            &weights,
            residual_count,
            parameter_count,
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
            let mut damped = normal.clone();
            for index in 0..parameter_count {
                let diagonal = normal[index * parameter_count + index].max(1.0);
                damped[index * parameter_count + index] += damping * diagonal;
            }
            let right_hand_side: Vec<_> = gradient.iter().map(|value| -value).collect();
            let Ok(step) = solve_linear_system(&damped, &right_hand_side, parameter_count) else {
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
        let (jacobian, _used) = finite_difference_jacobian(
            &parameters,
            lower,
            upper,
            &residuals,
            usize::MAX,
            &mut evaluate,
        )?;
        jacobian_evaluations += 1;
        final_jacobian = jacobian;
        let weights = robust_weights(&residuals, options.loss, options.f_scale);
        let (_, gradient) = normal_equations(
            &final_jacobian,
            &residuals,
            &weights,
            residual_count,
            parameter_count,
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
        let base_step = EPSILON_SQRT * parameters[column].abs().max(1.0);
        let forward_room = upper[column] - parameters[column];
        let backward_room = parameters[column] - lower[column];
        let forward_step = base_step.min(forward_room);
        let backward_step = base_step.min(backward_room);
        let (difference, denominator) =
            if forward_step > 0.0 && backward_step > 0.0 && used + 2 <= evaluation_budget {
                let mut plus = parameters.to_vec();
                plus[column] += forward_step;
                let plus_residuals = evaluate(&plus)?;
                used += 1;
                validate_residuals(&plus_residuals, Some(rows))?;
                let mut minus = parameters.to_vec();
                minus[column] -= backward_step;
                let minus_residuals = evaluate(&minus)?;
                used += 1;
                validate_residuals(&minus_residuals, Some(rows))?;
                (
                    plus_residuals
                        .iter()
                        .zip(minus_residuals)
                        .map(|(plus, minus)| plus - minus)
                        .collect::<Vec<_>>(),
                    forward_step + backward_step,
                )
            } else if forward_step > 0.0 {
                let mut plus = parameters.to_vec();
                plus[column] += forward_step;
                let plus_residuals = evaluate(&plus)?;
                used += 1;
                validate_residuals(&plus_residuals, Some(rows))?;
                (
                    plus_residuals
                        .iter()
                        .zip(residuals)
                        .map(|(plus, base)| plus - base)
                        .collect(),
                    forward_step,
                )
            } else if backward_step > 0.0 {
                let mut minus = parameters.to_vec();
                minus[column] -= backward_step;
                let minus_residuals = evaluate(&minus)?;
                used += 1;
                validate_residuals(&minus_residuals, Some(rows))?;
                (
                    residuals
                        .iter()
                        .zip(minus_residuals)
                        .map(|(base, minus)| base - minus)
                        .collect(),
                    backward_step,
                )
            } else {
                return Err(FitError::InvalidInput(format!(
                    "parameter {column} is fixed by equal bounds"
                )));
            };
        for row in 0..rows {
            jacobian[row * columns + column] = difference[row] / denominator;
        }
    }
    Ok((jacobian, used))
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
}
