use crate::validation::{finite_state, positive_state};
use crate::{EosError, EosResult, IsothermalEos};

const MAX_BRACKET_ITERATIONS: usize = 160;
const MAX_BISECTION_ITERATIONS: usize = 256;
const RELATIVE_ROOT_TOLERANCE: f64 = 1.0e-12;
const RESIDUAL_TOLERANCE: f64 = 1.0e-8;

pub(crate) fn solve_volume<E>(eos: &E, pressure: f64) -> EosResult<f64>
where
    E: IsothermalEos + ?Sized,
{
    solve_volume_function(
        |volume| eos.pressure(volume),
        pressure,
        eos.reference_volume(),
    )
}

#[allow(clippy::similar_names)]
pub(crate) fn solve_volume_function<F>(
    mut pressure_function: F,
    pressure: f64,
    reference_volume: f64,
) -> EosResult<f64>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    let target = finite_state(pressure, "pressure")?;
    let reference_volume = positive_state(reference_volume, "reference volume")?;
    let reference_pressure = pressure_function(reference_volume)?;
    let reference_residual = reference_pressure - target;
    let pressure_tolerance = 1.0e-10 * target.abs().max(reference_pressure.abs()).max(1.0);
    if reference_residual.abs() <= pressure_tolerance {
        return Ok(reference_volume);
    }

    let (lower, lower_residual, upper, upper_residual) = if reference_residual < 0.0 {
        let upper = reference_volume;
        let upper_residual = reference_residual;
        let mut lower = reference_volume;
        let mut lower_residual = reference_residual;
        let minimum_volume = reference_volume * 1.0e-14;

        for _ in 0..MAX_BRACKET_ITERATIONS {
            lower *= 0.8;
            if lower <= minimum_volume {
                break;
            }
            lower_residual = pressure_function(lower)? - target;
            if lower_residual >= 0.0 {
                break;
            }
        }
        if lower_residual < 0.0 {
            return Err(EosError::BracketingFailed);
        }
        (lower, lower_residual, upper, upper_residual)
    } else {
        let lower = reference_volume;
        let lower_residual = reference_residual;
        let mut upper = reference_volume;
        let mut upper_residual = reference_residual;
        let maximum_volume = reference_volume * 1.0e4;

        for _ in 0..MAX_BRACKET_ITERATIONS {
            upper *= 1.05;
            if upper >= maximum_volume {
                break;
            }
            upper_residual = pressure_function(upper)? - target;
            if upper_residual <= 0.0 {
                break;
            }
        }
        if upper_residual > 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        (lower, lower_residual, upper, upper_residual)
    };

    let absolute_tolerance = (f64::EPSILON * reference_volume).max(1.0e-14);
    let result = bisect_bracket(
        &mut pressure_function,
        target,
        lower,
        lower_residual,
        upper,
        upper_residual,
        absolute_tolerance,
        pressure_tolerance,
    )?;
    check_residual(&mut pressure_function, target, result)
}

#[allow(clippy::similar_names)]
pub(crate) fn solve_temperature_function<F>(
    mut pressure_function: F,
    pressure: f64,
    reference_temperature: f64,
) -> EosResult<f64>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    let target = finite_state(pressure, "pressure")?;
    let reference_temperature = positive_state(reference_temperature, "temperature")?;
    let reference_pressure = pressure_function(reference_temperature)?;
    let reference_residual = reference_pressure - target;
    let pressure_tolerance = 1.0e-10 * target.abs().max(reference_pressure.abs()).max(1.0);
    if reference_residual.abs() <= pressure_tolerance {
        return Ok(reference_temperature);
    }

    let minimum_temperature = reference_temperature * 1.0e-14;
    let maximum_temperature = reference_temperature * 1.0e8;
    let mut lower = reference_temperature;
    let mut upper = reference_temperature;
    let mut lower_residual = reference_residual;
    let mut upper_residual = reference_residual;
    let mut lower_active = true;
    let mut upper_active = true;
    let mut brackets = Vec::with_capacity(2);

    for _ in 0..MAX_BRACKET_ITERATIONS {
        if lower_active {
            let next_lower = (lower * 0.8).max(minimum_temperature);
            let next_residual = pressure_function(next_lower)? - target;
            if next_residual * lower_residual <= 0.0 {
                brackets.push((next_lower, next_residual, lower, lower_residual));
            }
            lower = next_lower;
            lower_residual = next_residual;
            lower_active = lower > minimum_temperature;
        }

        if upper_active {
            let next_upper = (upper * 1.25).min(maximum_temperature);
            let next_residual = pressure_function(next_upper)? - target;
            if next_residual * upper_residual <= 0.0 {
                brackets.push((upper, upper_residual, next_upper, next_residual));
            }
            upper = next_upper;
            upper_residual = next_residual;
            upper_active = upper < maximum_temperature;
        }

        if !brackets.is_empty() || (!lower_active && !upper_active) {
            break;
        }
    }

    if brackets.is_empty() {
        return Err(EosError::OutsideInvertibleRange);
    }

    let absolute_tolerance = f64::EPSILON * reference_temperature.max(1.0);
    let mut roots = Vec::with_capacity(brackets.len());
    for (bracket_lower, lower_value, bracket_upper, upper_value) in brackets {
        roots.push(bisect_bracket(
            &mut pressure_function,
            target,
            bracket_lower,
            lower_value,
            bracket_upper,
            upper_value,
            absolute_tolerance,
            pressure_tolerance,
        )?);
    }
    let result = roots
        .into_iter()
        .min_by(|left, right| {
            (left / reference_temperature)
                .ln()
                .abs()
                .total_cmp(&(right / reference_temperature).ln().abs())
        })
        .ok_or(EosError::ConvergenceFailed)?;
    check_residual(&mut pressure_function, target, result)
}

#[allow(clippy::too_many_arguments, clippy::similar_names)]
fn bisect_bracket<F>(
    pressure_function: &mut F,
    target: f64,
    mut lower: f64,
    mut lower_residual: f64,
    mut upper: f64,
    upper_residual: f64,
    absolute_tolerance: f64,
    pressure_tolerance: f64,
) -> EosResult<f64>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    if lower_residual * upper_residual > 0.0 {
        return Err(EosError::BracketingFailed);
    }
    let lower_is_positive = lower_residual >= 0.0;
    for _ in 0..MAX_BISECTION_ITERATIONS {
        let result = lower + 0.5 * (upper - lower);
        let residual = pressure_function(result)? - target;
        let value_tolerance = absolute_tolerance + RELATIVE_ROOT_TOLERANCE * result.abs();
        if residual.abs() <= pressure_tolerance || 0.5 * (upper - lower) <= value_tolerance {
            return Ok(result);
        }
        if (residual >= 0.0) == lower_is_positive {
            lower = result;
            lower_residual = residual;
        } else {
            upper = result;
        }
    }
    debug_assert!(lower_residual.is_finite());
    Err(EosError::ConvergenceFailed)
}

fn check_residual<F>(pressure_function: &mut F, target: f64, result: f64) -> EosResult<f64>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    let residual = (pressure_function(result)? - target).abs();
    if residual > RESIDUAL_TOLERANCE * target.abs().max(1.0) {
        Err(EosError::ConvergenceFailed)
    } else {
        Ok(result)
    }
}
