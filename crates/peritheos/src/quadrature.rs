use crate::{EosError, EosResult};

const INITIAL_PANELS: u32 = 32;
const MAX_DEPTH: usize = 24;
const ABSOLUTE_TOLERANCE: f64 = 1.0e-10;
const RELATIVE_TOLERANCE: f64 = 1.0e-10;

#[allow(clippy::float_cmp)]
pub(crate) fn integrate<F>(mut function: F, lower: f64, upper: f64) -> EosResult<f64>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    if !lower.is_finite() || !upper.is_finite() {
        return Err(EosError::InvalidState {
            name: "integration bounds",
            reason: "must be finite",
        });
    }
    if lower == upper {
        return Ok(0.0);
    }
    if lower > upper {
        return integrate(function, upper, lower).map(|value| -value);
    }

    let panel_width = (upper - lower) / f64::from(INITIAL_PANELS);
    let mut total = 0.0;
    for panel in 0..INITIAL_PANELS {
        let a = lower + f64::from(panel) * panel_width;
        let b = if panel + 1 == INITIAL_PANELS {
            upper
        } else {
            a + panel_width
        };
        let midpoint = a + 0.5 * (b - a);
        let fa = finite_integrand(function(a)?)?;
        let fm = finite_integrand(function(midpoint)?)?;
        let fb = finite_integrand(function(b)?)?;
        let whole = simpson(a, b, fa, fm, fb);
        let tolerance =
            (ABSOLUTE_TOLERANCE + RELATIVE_TOLERANCE * whole.abs()) / f64::from(INITIAL_PANELS);
        total += adaptive_simpson(&mut function, a, b, fa, fm, fb, whole, tolerance, MAX_DEPTH)?;
    }
    finite_integrand(total)
}

#[allow(clippy::too_many_arguments)]
fn adaptive_simpson<F>(
    function: &mut F,
    a: f64,
    b: f64,
    fa: f64,
    fm: f64,
    fb: f64,
    whole: f64,
    tolerance: f64,
    depth: usize,
) -> EosResult<f64>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    let midpoint = a + 0.5 * (b - a);
    let left_midpoint = a + 0.5 * (midpoint - a);
    let right_midpoint = midpoint + 0.5 * (b - midpoint);
    let f_left_midpoint = finite_integrand(function(left_midpoint)?)?;
    let f_right_midpoint = finite_integrand(function(right_midpoint)?)?;
    let left = simpson(a, midpoint, fa, f_left_midpoint, fm);
    let right = simpson(midpoint, b, fm, f_right_midpoint, fb);
    let correction = left + right - whole;

    if correction.abs() <= 15.0 * tolerance {
        return finite_integrand(left + right + correction / 15.0);
    }
    if depth == 0 {
        // At machine precision, finite-difference integrands can stop reducing
        // their local error estimate. Preserve the best corrected estimate;
        // callers still apply independent state and residual validation.
        return finite_integrand(left + right + correction / 15.0);
    }

    let left_integral = adaptive_simpson(
        function,
        a,
        midpoint,
        fa,
        f_left_midpoint,
        fm,
        left,
        tolerance / 2.0,
        depth - 1,
    )?;
    let right_integral = adaptive_simpson(
        function,
        midpoint,
        b,
        fm,
        f_right_midpoint,
        fb,
        right,
        tolerance / 2.0,
        depth - 1,
    )?;
    finite_integrand(left_integral + right_integral)
}

fn simpson(a: f64, b: f64, fa: f64, fm: f64, fb: f64) -> f64 {
    (b - a) * (fa + 4.0 * fm + fb) / 6.0
}

fn finite_integrand(value: f64) -> EosResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EosError::NonFiniteResult)
    }
}

#[cfg(test)]
mod tests {
    use super::integrate;

    #[test]
    fn integrates_polynomial_and_reversed_bounds() {
        let forward = integrate(|x| Ok(x * x), 0.0, 1.0).unwrap();
        let reverse = integrate(|x| Ok(x * x), 1.0, 0.0).unwrap();
        assert!((forward - 1.0 / 3.0).abs() < 1.0e-13);
        assert!((reverse + 1.0 / 3.0).abs() < 1.0e-13);
    }

    #[test]
    fn initial_panels_find_a_boundary_localized_integrand() {
        let value = integrate(|x| Ok(x * x * (-x).exp()), 0.0, 150.0).unwrap();
        assert!((value - 2.0).abs() < 1.0e-11);
    }
}
