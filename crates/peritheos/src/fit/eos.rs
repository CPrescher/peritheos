//! End-to-end native EOS fitting built on the least-squares kernel.

use crate::{IsothermalEos, ThermalEos};

use super::{
    least_squares, least_squares_structured, parameter_covariance, parameter_covariance_structured,
    FitError, SolverOptions, SolverResult, StructuredLayout,
};

/// Isothermal pressure-volume observations and their error model.
#[derive(Clone, Copy, Debug)]
pub struct IsothermalObservations<'a> {
    /// Observed pressures, flattened in row-major order.
    pub pressure: &'a [f64],
    /// Measured volumes in the same order as `pressure`.
    pub volume: &'a [f64],
    /// Pressure standard deviations. One value is required per observation.
    pub pressure_sigma: &'a [f64],
    /// Optional volume standard deviations. Their presence makes volume latent.
    pub volume_sigma: Option<&'a [f64]>,
    /// Optional per-observation lower Cholesky factors in row-major order.
    ///
    /// When present, the factors have shape `(point_count, 2, 2)` and replace
    /// the independent standard deviations when whitening residuals.
    pub observation_cholesky: Option<&'a [f64]>,
}

/// Thermal pressure-volume-temperature observations and their error model.
#[derive(Clone, Copy, Debug)]
pub struct ThermalObservations<'a> {
    /// Observed pressures, flattened in row-major order.
    pub pressure: &'a [f64],
    /// Measured volumes in the same order as `pressure`.
    pub volume: &'a [f64],
    /// Measured temperatures in the same order as `pressure`.
    pub temperature: &'a [f64],
    /// Pressure standard deviations. One value is required per observation.
    pub pressure_sigma: &'a [f64],
    /// Optional volume standard deviations. Their presence makes volume latent.
    pub volume_sigma: Option<&'a [f64]>,
    /// Optional temperature standard deviations. Their presence makes
    /// temperature latent.
    pub temperature_sigma: Option<&'a [f64]>,
    /// Optional per-observation lower Cholesky factors in row-major order.
    ///
    /// When present, the factors have shape `(point_count, 3, 3)` and replace
    /// the independent standard deviations when whitening residuals.
    pub observation_cholesky: Option<&'a [f64]>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct FitStatistics {
    chi_square: f64,
    reduced_chi_square: Option<f64>,
    degrees_of_freedom: i64,
    aic: f64,
    bic: f64,
}

/// Interpreted EOS fit result with parameter uncertainty and adjusted states.
///
/// Matrix fields use row-major model-parameter order. [`Self::solver`] retains
/// the complete optimizer vector and Jacobian for advanced diagnostics, while
/// [`Self::parameters`] contains only the parameters passed to the model
/// factory.
#[derive(Clone, Debug, PartialEq)]
pub struct EosFitResult {
    /// Final model parameters in factory order, excluding latent states.
    pub parameters: Vec<f64>,
    /// Profiled, unscaled covariance of the model parameters.
    pub covariance: Vec<f64>,
    /// Parameter standard errors derived from the covariance diagonal.
    pub standard_errors: Vec<f64>,
    /// Parameter correlation matrix.
    pub correlation: Vec<f64>,
    /// Final volumes used for prediction; measured values when not latent.
    pub adjusted_volume: Vec<f64>,
    /// Final temperatures used for prediction; absent for isothermal fits.
    pub adjusted_temperature: Option<Vec<f64>>,
    /// Pressure predicted at the final adjusted coordinates.
    pub predicted_pressure: Vec<f64>,
    /// Sum of squared weighted residuals.
    pub chi_square: f64,
    /// Chi-square divided by the degrees of freedom, when positive.
    pub reduced_chi_square: Option<f64>,
    /// Residual count minus all fitted variables, including latent states.
    pub degrees_of_freedom: i64,
    /// Akaike information criterion using all fitted variables.
    pub aic: f64,
    /// Bayesian information criterion using all fitted variables.
    pub bic: f64,
    /// Complete least-squares variables and solver diagnostics.
    pub solver: SolverResult,
}

fn fit_statistics(solver: &SolverResult) -> Result<FitStatistics, FitError> {
    let residual_count = u32::try_from(solver.residual_count).map_err(|_| {
        FitError::InvalidInput("residual count exceeds the statistics limit".to_owned())
    })?;
    let variable_count = u32::try_from(solver.parameters.len()).map_err(|_| {
        FitError::InvalidInput("variable count exceeds the statistics limit".to_owned())
    })?;
    let positive_degrees_of_freedom = residual_count
        .checked_sub(variable_count)
        .filter(|value| *value > 0);
    let degrees_of_freedom = i64::from(residual_count) - i64::from(variable_count);
    let chi_square: f64 = solver.residuals.iter().map(|value| value * value).sum();
    let reduced_chi_square = positive_degrees_of_freedom.map(|value| chi_square / f64::from(value));
    let residual_count = f64::from(residual_count);
    let variable_count = f64::from(variable_count);
    let log_variance = (chi_square / residual_count).max(f64::MIN_POSITIVE).ln();
    let aic = residual_count * log_variance + 2.0 * variable_count;
    let bic = residual_count * log_variance + variable_count * residual_count.ln();
    Ok(FitStatistics {
        chi_square,
        reduced_chi_square,
        degrees_of_freedom,
        aic,
        bic,
    })
}

fn interpret_result(
    solver: SolverResult,
    predicted_pressure: Vec<f64>,
    adjusted_volume: Vec<f64>,
    adjusted_temperature: Option<Vec<f64>>,
    parameter_count: usize,
    point_count: usize,
    latent_coordinate_count: usize,
) -> Result<EosFitResult, FitError> {
    let column_count = solver.parameters.len();
    let covariance = if latent_coordinate_count == 0 {
        parameter_covariance(
            &solver.jacobian,
            solver.residual_count,
            column_count,
            parameter_count,
        )
    } else {
        parameter_covariance_structured(
            &solver.jacobian,
            StructuredLayout {
                global_parameter_count: parameter_count,
                point_count,
                latent_coordinate_count,
            },
        )
    }?;
    let standard_errors = (0..parameter_count)
        .map(|index| covariance[index * parameter_count + index].max(0.0).sqrt())
        .collect::<Vec<_>>();
    let mut correlation = vec![0.0; parameter_count * parameter_count];
    for row in 0..parameter_count {
        for column in 0..parameter_count {
            let denominator = standard_errors[row] * standard_errors[column];
            if denominator > 0.0 {
                correlation[row * parameter_count + column] =
                    covariance[row * parameter_count + column] / denominator;
            }
        }
    }
    let statistics = fit_statistics(&solver)?;
    Ok(EosFitResult {
        parameters: solver.parameters[..parameter_count].to_vec(),
        covariance,
        standard_errors,
        correlation,
        adjusted_volume,
        adjusted_temperature,
        predicted_pressure,
        chi_square: statistics.chi_square,
        reduced_chi_square: statistics.reduced_chi_square,
        degrees_of_freedom: statistics.degrees_of_freedom,
        aic: statistics.aic,
        bic: statistics.bic,
        solver,
    })
}

fn validate_common(
    pressure: &[f64],
    volume: &[f64],
    pressure_sigma: &[f64],
    point_count: usize,
) -> Result<(), FitError> {
    if point_count == 0 || volume.len() != point_count || pressure_sigma.len() != point_count {
        return Err(FitError::InvalidInput(
            "EOS observation arrays must have the same non-zero length".to_owned(),
        ));
    }
    if pressure.iter().any(|value| !value.is_finite())
        || volume
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        || pressure_sigma
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(FitError::InvalidInput(
            "EOS observations and uncertainties must be finite and valid".to_owned(),
        ));
    }
    Ok(())
}

fn validate_optional_sigma(
    sigma: Option<&[f64]>,
    point_count: usize,
    name: &str,
) -> Result<(), FitError> {
    if let Some(values) = sigma {
        if values.len() != point_count
            || values
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(FitError::InvalidInput(format!(
                "{name} must contain one positive finite value per observation"
            )));
        }
    }
    Ok(())
}

fn validate_cholesky(
    factors: Option<&[f64]>,
    point_count: usize,
    component_count: usize,
) -> Result<(), FitError> {
    let Some(factors) = factors else {
        return Ok(());
    };
    if factors.len() != point_count * component_count * component_count
        || factors.iter().any(|value| !value.is_finite())
    {
        return Err(FitError::InvalidInput(
            "observation Cholesky factors have invalid dimensions or values".to_owned(),
        ));
    }
    for point in 0..point_count {
        let offset = point * component_count * component_count;
        for row in 0..component_count {
            if factors[offset + row * component_count + row] <= 0.0 {
                return Err(FitError::InvalidInput(
                    "observation Cholesky factors must have positive diagonals".to_owned(),
                ));
            }
            for column in row + 1..component_count {
                if factors[offset + row * component_count + column] != 0.0 {
                    return Err(FitError::InvalidInput(
                        "observation Cholesky factors must be lower triangular".to_owned(),
                    ));
                }
            }
        }
    }
    Ok(())
}

fn latent_slice(
    parameters: &[f64],
    measured: &[f64],
    adjusted: bool,
    offset: &mut usize,
) -> Vec<f64> {
    if adjusted {
        let values = parameters[*offset..*offset + measured.len()].to_vec();
        *offset += measured.len();
        values
    } else {
        measured.to_vec()
    }
}

fn whiten_correlated(raw_components: &[&[f64]], cholesky: &[f64], point_count: usize) -> Vec<f64> {
    let component_count = raw_components.len();
    let mut residuals = vec![0.0; component_count * point_count];
    let mut solved = vec![0.0; component_count];
    for point in 0..point_count {
        let factor_offset = point * component_count * component_count;
        for row in 0..component_count {
            let mut value = raw_components[row][point];
            for column in 0..row {
                value -= cholesky[factor_offset + row * component_count + column] * solved[column];
            }
            solved[row] = value / cholesky[factor_offset + row * component_count + row];
            residuals[row * point_count + point] = solved[row];
        }
    }
    residuals
}

#[allow(clippy::too_many_arguments)]
fn solve_native<F>(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    parameter_count: usize,
    point_count: usize,
    latent_coordinate_count: usize,
    evaluate: F,
) -> Result<SolverResult, FitError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, FitError>,
{
    if parameter_count == 0
        || initial.len() != parameter_count + point_count * latent_coordinate_count
    {
        return Err(FitError::InvalidInput(
            "EOS parameter and latent-coordinate dimensions are inconsistent".to_owned(),
        ));
    }
    if latent_coordinate_count == 0 {
        least_squares(initial, lower, upper, options, evaluate)
    } else {
        least_squares_structured(
            initial,
            lower,
            upper,
            options,
            StructuredLayout {
                global_parameter_count: parameter_count,
                point_count,
                latent_coordinate_count,
            },
            evaluate,
        )
    }
}

type OptimizerVectors = (Vec<f64>, Vec<f64>, Vec<f64>);

fn optimizer_vectors(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    latent_coordinates: &[&[f64]],
) -> Result<OptimizerVectors, FitError> {
    if initial.is_empty() || lower.len() != initial.len() || upper.len() != initial.len() {
        return Err(FitError::InvalidInput(
            "EOS initial parameters and bounds must have the same non-zero length".to_owned(),
        ));
    }
    let latent_count: usize = latent_coordinates.iter().map(|values| values.len()).sum();
    let mut optimizer_initial = Vec::with_capacity(initial.len() + latent_count);
    let mut optimizer_lower = Vec::with_capacity(lower.len() + latent_count);
    let mut optimizer_upper = Vec::with_capacity(upper.len() + latent_count);
    optimizer_initial.extend_from_slice(initial);
    optimizer_lower.extend_from_slice(lower);
    optimizer_upper.extend_from_slice(upper);
    for values in latent_coordinates {
        optimizer_initial.extend_from_slice(values);
        optimizer_lower.extend(std::iter::repeat_n(f64::MIN_POSITIVE, values.len()));
        optimizer_upper.extend(std::iter::repeat_n(f64::INFINITY, values.len()));
    }
    Ok((optimizer_initial, optimizer_lower, optimizer_upper))
}

/// Fit an isothermal EOS without leaving Rust during residual evaluation.
///
/// `factory` constructs the EOS from the fitted model parameters. Latent
/// volumes and their positivity bounds are assembled internally when
/// `volume_sigma` is present.
///
/// # Errors
///
/// Returns [`FitError`] for inconsistent inputs, model construction or
/// evaluation failures, and solver failures.
#[allow(clippy::too_many_lines)]
pub fn fit_isothermal_eos<M, F>(
    observations: IsothermalObservations<'_>,
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    factory: F,
) -> Result<EosFitResult, FitError>
where
    M: IsothermalEos,
    F: Fn(&[f64]) -> Result<M, FitError>,
{
    let parameter_count = initial.len();
    let point_count = observations.pressure.len();
    validate_common(
        observations.pressure,
        observations.volume,
        observations.pressure_sigma,
        point_count,
    )?;
    validate_optional_sigma(observations.volume_sigma, point_count, "volume_sigma")?;
    validate_cholesky(observations.observation_cholesky, point_count, 2)?;
    if observations.observation_cholesky.is_some() && observations.volume_sigma.is_none() {
        return Err(FitError::InvalidInput(
            "correlated pressure-volume errors require latent volumes".to_owned(),
        ));
    }
    let latent_count = usize::from(observations.volume_sigma.is_some());
    let latent_coordinates = observations
        .volume_sigma
        .map_or_else(Vec::new, |_| vec![observations.volume]);
    let (optimizer_initial, optimizer_lower, optimizer_upper) =
        optimizer_vectors(initial, lower, upper, &latent_coordinates)?;
    let evaluate = |parameters: &[f64]| {
        let model = factory(&parameters[..parameter_count])?;
        let mut offset = parameter_count;
        let volume = latent_slice(
            parameters,
            observations.volume,
            observations.volume_sigma.is_some(),
            &mut offset,
        );
        let predicted = volume
            .iter()
            .map(|value| model.pressure(*value).map_err(FitError::from))
            .collect::<Result<Vec<_>, _>>()?;
        let pressure_residual = predicted
            .iter()
            .zip(observations.pressure)
            .map(|(predicted, observed)| predicted - observed)
            .collect::<Vec<_>>();
        if let Some(cholesky) = observations.observation_cholesky {
            let volume_residual = volume
                .iter()
                .zip(observations.volume)
                .map(|(adjusted, measured)| adjusted - measured)
                .collect::<Vec<_>>();
            return Ok(whiten_correlated(
                &[&pressure_residual, &volume_residual],
                cholesky,
                point_count,
            ));
        }
        let mut residuals = pressure_residual
            .iter()
            .zip(observations.pressure_sigma)
            .map(|(residual, sigma)| residual / sigma)
            .collect::<Vec<_>>();
        if let Some(sigma) = observations.volume_sigma {
            residuals.extend(
                volume
                    .iter()
                    .zip(observations.volume)
                    .zip(sigma)
                    .map(|((adjusted, measured), sigma)| (adjusted - measured) / sigma),
            );
        }
        Ok(residuals)
    };
    let solver = solve_native(
        &optimizer_initial,
        &optimizer_lower,
        &optimizer_upper,
        options,
        parameter_count,
        point_count,
        latent_count,
        evaluate,
    )?;
    let model = factory(&solver.parameters[..parameter_count])?;
    let volume = if observations.volume_sigma.is_some() {
        &solver.parameters[parameter_count..parameter_count + point_count]
    } else {
        observations.volume
    };
    let predicted_pressure = volume
        .iter()
        .map(|value| model.pressure(*value).map_err(FitError::from))
        .collect::<Result<Vec<_>, _>>()?;
    let adjusted_volume = volume.to_vec();
    interpret_result(
        solver,
        predicted_pressure,
        adjusted_volume,
        None,
        parameter_count,
        point_count,
        latent_count,
    )
}

/// Fit a thermal EOS without leaving Rust during residual evaluation.
///
/// `factory` constructs the EOS from the fitted model parameters. Latent
/// volume and temperature blocks and their positivity bounds are assembled
/// internally when the corresponding standard deviations are present.
///
/// # Errors
///
/// Returns [`FitError`] for inconsistent inputs, model construction or
/// evaluation failures, and solver failures.
pub fn fit_thermal_eos<M, F>(
    observations: ThermalObservations<'_>,
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    factory: F,
) -> Result<EosFitResult, FitError>
where
    M: ThermalEos,
    F: Fn(&[f64]) -> Result<M, FitError>,
{
    fit_thermal_eos_by(
        observations,
        initial,
        lower,
        upper,
        options,
        factory,
        |model, volume, temperature| model.pressure(volume, temperature).map_err(FitError::from),
    )
}

/// Jointly fit reference-isotherm and thermal parameters in one Rust model factory.
///
/// Rust represents the combined parameter set as one ordered slice rather
/// than using Python's dotted parameter names. The factory receives that full
/// slice and may reconstruct both the reference EOS and the thermal model.
/// This is the Rust counterpart of Python's `fit_joint_eos` convenience API.
///
/// # Errors
///
/// Returns [`FitError`] for inconsistent inputs, model construction or
/// evaluation failures, and solver failures.
pub fn fit_joint_eos<M, F>(
    observations: ThermalObservations<'_>,
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    factory: F,
) -> Result<EosFitResult, FitError>
where
    M: ThermalEos,
    F: Fn(&[f64]) -> Result<M, FitError>,
{
    fit_thermal_eos(observations, initial, lower, upper, options, factory)
}

/// Fit a thermal pressure model through an explicit Rust pressure evaluator.
///
/// This variant supports type-erased Rust model enums which cannot implement
/// [`ThermalEos`] because its reference-EOS associated type differs between
/// enum variants. Ordinary concrete thermal models should use
/// [`fit_thermal_eos`] instead.
///
/// # Errors
///
/// Returns [`FitError`] for inconsistent inputs, model construction or
/// evaluation failures, and solver failures.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn fit_thermal_eos_by<M, F, P>(
    observations: ThermalObservations<'_>,
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    options: SolverOptions,
    factory: F,
    pressure: P,
) -> Result<EosFitResult, FitError>
where
    F: Fn(&[f64]) -> Result<M, FitError>,
    P: Fn(&M, f64, f64) -> Result<f64, FitError>,
{
    let parameter_count = initial.len();
    let point_count = observations.pressure.len();
    validate_common(
        observations.pressure,
        observations.volume,
        observations.pressure_sigma,
        point_count,
    )?;
    if observations.temperature.len() != point_count
        || observations
            .temperature
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(FitError::InvalidInput(
            "temperature must contain one positive finite value per observation".to_owned(),
        ));
    }
    validate_optional_sigma(observations.volume_sigma, point_count, "volume_sigma")?;
    validate_optional_sigma(
        observations.temperature_sigma,
        point_count,
        "temperature_sigma",
    )?;
    validate_cholesky(observations.observation_cholesky, point_count, 3)?;
    if observations.observation_cholesky.is_some()
        && (observations.volume_sigma.is_none() || observations.temperature_sigma.is_none())
    {
        return Err(FitError::InvalidInput(
            "correlated pressure-volume-temperature errors require latent coordinates".to_owned(),
        ));
    }
    let latent_count = usize::from(observations.volume_sigma.is_some())
        + usize::from(observations.temperature_sigma.is_some());
    let mut latent_coordinates = Vec::with_capacity(latent_count);
    if observations.volume_sigma.is_some() {
        latent_coordinates.push(observations.volume);
    }
    if observations.temperature_sigma.is_some() {
        latent_coordinates.push(observations.temperature);
    }
    let (optimizer_initial, optimizer_lower, optimizer_upper) =
        optimizer_vectors(initial, lower, upper, &latent_coordinates)?;
    let evaluate = |parameters: &[f64]| {
        let model = factory(&parameters[..parameter_count])?;
        let mut offset = parameter_count;
        let volume = latent_slice(
            parameters,
            observations.volume,
            observations.volume_sigma.is_some(),
            &mut offset,
        );
        let temperature = latent_slice(
            parameters,
            observations.temperature,
            observations.temperature_sigma.is_some(),
            &mut offset,
        );
        let predicted = volume
            .iter()
            .zip(&temperature)
            .map(|(volume, temperature)| pressure(&model, *volume, *temperature))
            .collect::<Result<Vec<_>, _>>()?;
        let pressure_residual = predicted
            .iter()
            .zip(observations.pressure)
            .map(|(predicted, observed)| predicted - observed)
            .collect::<Vec<_>>();
        if let Some(cholesky) = observations.observation_cholesky {
            let volume_residual = volume
                .iter()
                .zip(observations.volume)
                .map(|(adjusted, measured)| adjusted - measured)
                .collect::<Vec<_>>();
            let temperature_residual = temperature
                .iter()
                .zip(observations.temperature)
                .map(|(adjusted, measured)| adjusted - measured)
                .collect::<Vec<_>>();
            return Ok(whiten_correlated(
                &[&pressure_residual, &volume_residual, &temperature_residual],
                cholesky,
                point_count,
            ));
        }
        let mut residuals = pressure_residual
            .iter()
            .zip(observations.pressure_sigma)
            .map(|(residual, sigma)| residual / sigma)
            .collect::<Vec<_>>();
        if let Some(sigma) = observations.volume_sigma {
            residuals.extend(
                volume
                    .iter()
                    .zip(observations.volume)
                    .zip(sigma)
                    .map(|((adjusted, measured), sigma)| (adjusted - measured) / sigma),
            );
        }
        if let Some(sigma) = observations.temperature_sigma {
            residuals.extend(
                temperature
                    .iter()
                    .zip(observations.temperature)
                    .zip(sigma)
                    .map(|((adjusted, measured), sigma)| (adjusted - measured) / sigma),
            );
        }
        Ok(residuals)
    };
    let solver = solve_native(
        &optimizer_initial,
        &optimizer_lower,
        &optimizer_upper,
        options,
        parameter_count,
        point_count,
        latent_count,
        evaluate,
    )?;
    let model = factory(&solver.parameters[..parameter_count])?;
    let mut offset = parameter_count;
    let volume = latent_slice(
        &solver.parameters,
        observations.volume,
        observations.volume_sigma.is_some(),
        &mut offset,
    );
    let temperature = latent_slice(
        &solver.parameters,
        observations.temperature,
        observations.temperature_sigma.is_some(),
        &mut offset,
    );
    let predicted_pressure = volume
        .iter()
        .zip(&temperature)
        .map(|(volume, temperature)| pressure(&model, *volume, *temperature))
        .collect::<Result<Vec<_>, _>>()?;
    interpret_result(
        solver,
        predicted_pressure,
        volume,
        Some(temperature),
        parameter_count,
        point_count,
        latent_count,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isothermal::BM3;
    use crate::thermal::MieGruneisenDebye;

    #[test]
    fn isothermal_fit_recovers_model_without_callback_runtime() {
        let expected = BM3::new(10.0, 120.0, 4.3).unwrap();
        let volume = [8.0, 8.5, 9.0, 9.5, 10.0];
        let pressure = volume
            .iter()
            .map(|value| expected.pressure(*value).unwrap())
            .collect::<Vec<_>>();
        let sigma = [0.1; 5];
        let result = fit_isothermal_eos(
            IsothermalObservations {
                pressure: &pressure,
                volume: &volume,
                pressure_sigma: &sigma,
                volume_sigma: None,
                observation_cholesky: None,
            },
            &[110.0, 4.0],
            &[50.0, f64::NEG_INFINITY],
            &[200.0, f64::INFINITY],
            SolverOptions::default(),
            |parameters| BM3::new(10.0, parameters[0], parameters[1]).map_err(FitError::from),
        )
        .unwrap();

        assert!(result.solver.success, "{}", result.solver.message);
        assert!((result.parameters[0] - 120.0).abs() < 1.0e-7);
        assert!((result.parameters[1] - 4.3).abs() < 1.0e-7);
        assert_eq!(result.parameters.len(), 2);
        assert_eq!(result.adjusted_volume, volume);
        assert_eq!(result.adjusted_temperature, None);
        assert_eq!(result.covariance.len(), 4);
        assert_eq!(result.standard_errors.len(), 2);
        assert_eq!(result.correlation.len(), 4);
        assert!((result.standard_errors[0].powi(2) - result.covariance[0]).abs() < 1.0e-12);
        assert!((result.standard_errors[1].powi(2) - result.covariance[3]).abs() < 1.0e-12);
        assert!((result.correlation[0] - 1.0).abs() < 1.0e-12);
        assert!((result.correlation[3] - 1.0).abs() < 1.0e-12);
        assert_eq!(result.degrees_of_freedom, 3);
        assert!(result.chi_square < 1.0e-10);
        assert!(result.reduced_chi_square.unwrap() < 1.0e-10);
        assert!(result.aic.is_finite());
        assert!(result.bic.is_finite());
        assert_eq!(result.predicted_pressure.len(), volume.len());
    }

    #[test]
    fn isothermal_fit_separates_latent_volumes_from_model_parameters() {
        let expected = BM3::new(10.0, 120.0, 4.3).unwrap();
        let volume = [8.0, 8.5, 9.0, 9.5, 10.0];
        let pressure = volume
            .iter()
            .map(|value| expected.pressure(*value).unwrap())
            .collect::<Vec<_>>();
        let pressure_sigma = [0.1; 5];
        let volume_sigma = [0.01; 5];
        let result = fit_isothermal_eos(
            IsothermalObservations {
                pressure: &pressure,
                volume: &volume,
                pressure_sigma: &pressure_sigma,
                volume_sigma: Some(&volume_sigma),
                observation_cholesky: None,
            },
            &[120.0, 4.3],
            &[50.0, 1.0],
            &[200.0, 10.0],
            SolverOptions::default(),
            |parameters| BM3::new(10.0, parameters[0], parameters[1]).map_err(FitError::from),
        )
        .unwrap();

        assert_eq!(result.parameters.len(), 2);
        assert_eq!(result.solver.parameters.len(), 7);
        assert_eq!(result.adjusted_volume.len(), volume.len());
        assert_eq!(result.degrees_of_freedom, 3);
        assert_eq!(result.covariance.len(), 4);
    }

    #[test]
    fn correlated_whitening_uses_component_major_residual_order() {
        let raw_pressure = [2.0, 4.0];
        let raw_volume = [3.0, 6.0];
        let factors = [2.0, 0.0, 1.0, 1.0, 4.0, 0.0, 2.0, 2.0];
        let whitened = whiten_correlated(&[&raw_pressure, &raw_volume], &factors, 2);
        assert_eq!(whitened, vec![1.0, 1.0, 2.0, 2.0]);
    }

    #[test]
    fn joint_fit_factory_updates_reference_and_thermal_parameters() {
        let expected = MieGruneisenDebye::new(
            BM3::new(1.0, 160.0, 4.0).unwrap(),
            300.0,
            800.0,
            1.5,
            1.0,
            2.0,
        )
        .unwrap();
        let volume = [0.75, 0.8, 0.85, 0.9, 0.95, 1.0];
        let temperature = [400.0, 700.0, 1_000.0, 1_400.0, 1_900.0, 2_500.0];
        let pressure = volume
            .iter()
            .zip(temperature)
            .map(|(&volume, temperature)| expected.pressure(volume, temperature).unwrap())
            .collect::<Vec<_>>();
        let sigma = [0.01; 6];
        let result = fit_joint_eos(
            ThermalObservations {
                pressure: &pressure,
                volume: &volume,
                temperature: &temperature,
                pressure_sigma: &sigma,
                volume_sigma: None,
                temperature_sigma: None,
                observation_cholesky: None,
            },
            &[150.0, 1.3],
            &[100.0, 0.5],
            &[220.0, 3.0],
            SolverOptions::default(),
            |parameters| {
                let reference = BM3::new(1.0, parameters[0], 4.0).map_err(FitError::from)?;
                MieGruneisenDebye::new(reference, 300.0, 800.0, parameters[1], 1.0, 2.0)
                    .map_err(FitError::from)
            },
        )
        .unwrap();

        assert!(result.solver.success, "{}", result.solver.message);
        assert!((result.parameters[0] - 160.0).abs() < 1.0e-6);
        assert!((result.parameters[1] - 1.5).abs() < 1.0e-6);
        assert_eq!(result.adjusted_volume, volume);
        assert_eq!(result.adjusted_temperature, Some(temperature.to_vec()));
        assert_eq!(result.degrees_of_freedom, 4);
    }
}
