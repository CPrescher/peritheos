use crate::validation::finite_state;
use crate::{EosError, EosResult, IsothermalEos};

const MAX_BRACKET_ITERATIONS: usize = 160;
const MAX_BISECTION_ITERATIONS: usize = 256;
const RELATIVE_ROOT_TOLERANCE: f64 = 1.0e-12;
const RESIDUAL_TOLERANCE: f64 = 1.0e-8;

#[allow(clippy::similar_names)]
pub(crate) fn solve_volume<E>(eos: &E, pressure: f64) -> EosResult<f64>
where
    E: IsothermalEos + ?Sized,
{
    let target = finite_state(pressure, "pressure")?;
    let reference_volume = eos.reference_volume();
    let reference_pressure = eos.pressure(reference_volume)?;
    let reference_residual = reference_pressure - target;
    let pressure_tolerance = 1.0e-10 * target.abs().max(reference_pressure.abs()).max(1.0);
    if reference_residual.abs() <= pressure_tolerance {
        return Ok(reference_volume);
    }

    let (mut lower, mut lower_residual, mut upper, mut upper_residual) = if reference_residual < 0.0
    {
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
            lower_residual = eos.pressure(lower)? - target;
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
            upper_residual = eos.pressure(upper)? - target;
            if upper_residual <= 0.0 {
                break;
            }
        }
        if upper_residual > 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        (lower, lower_residual, upper, upper_residual)
    };

    if lower_residual < 0.0 || upper_residual > 0.0 {
        return Err(EosError::BracketingFailed);
    }

    let absolute_tolerance = (f64::EPSILON * reference_volume).max(1.0e-14);
    let mut result = lower + 0.5 * (upper - lower);
    for _ in 0..MAX_BISECTION_ITERATIONS {
        result = lower + 0.5 * (upper - lower);
        let residual = eos.pressure(result)? - target;
        let volume_tolerance = absolute_tolerance + RELATIVE_ROOT_TOLERANCE * result.abs();

        if residual.abs() <= pressure_tolerance || 0.5 * (upper - lower) <= volume_tolerance {
            break;
        }
        if residual > 0.0 {
            lower = result;
            lower_residual = residual;
        } else {
            upper = result;
            upper_residual = residual;
        }
    }

    let residual = (eos.pressure(result)? - target).abs();
    if residual > RESIDUAL_TOLERANCE * target.abs().max(1.0) {
        return Err(EosError::ConvergenceFailed);
    }

    // Keep the sign variables live through the final iteration so accidental
    // loss of the physical bracket is visible to debug assertions.
    debug_assert!(lower_residual >= 0.0 && upper_residual <= 0.0);
    Ok(result)
}
