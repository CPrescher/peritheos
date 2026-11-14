//! Analytic pressure Jacobians for the canonical BM3 and q-compromise laws.

use crate::isothermal::BM3;
use crate::thermal::{debye_function_3, QCompromiseDebye, GAS_CONSTANT};
use crate::validation::positive_state;

use super::{FitError, PressureDerivatives};

/// Differentiate BM3 pressure in the requested canonical parameter order.
///
/// Names are `V0`, `K0` and `K0_prime`; subsets and permutations are allowed.
/// The factory must use these parameters directly, without transformations.
/// Pressure and volume units are those of the model.
///
/// # Errors
/// Returns an error for an invalid volume, unknown name or non-finite result.
pub fn bm3_pressure_derivatives<S: AsRef<str>>(
    model: &BM3,
    volume: f64,
    parameter_names: &[S],
) -> Result<PressureDerivatives, FitError> {
    let volume = positive_state(volume, "volume")?;
    let eta = (model.v0 / volume).cbrt();
    let h = eta * eta;
    let eta5 = h * h * eta;
    let strain = (h - 1.0) / 2.0;
    let modulus = eta5
        * (model.k0
            + (3.0 * model.k0 * model.k0_prime - 5.0 * model.k0) * strain
            + 13.5 * (model.k0 * model.k0_prime - 4.0 * model.k0) * strain * strain);
    let pressure_per_k0 = 3.0 * strain * eta5 * (1.0 + 1.5 * (model.k0_prime - 4.0) * strain);
    let parameters = parameter_names
        .iter()
        .map(|name| match name.as_ref() {
            "V0" => Ok(modulus / model.v0),
            "K0" => Ok(pressure_per_k0),
            "K0_prime" => Ok(4.5 * model.k0 * strain * strain * eta5),
            name => Err(FitError::InvalidInput(format!(
                "unknown BM3 derivative parameter {name:?}"
            ))),
        })
        .collect::<Result<Vec<_>, _>>()?;
    checked(PressureDerivatives {
        parameters,
        volume: -modulus / volume,
        temperature: 0.0,
    })
}

// Match the pressure kernel's stable small/large argument approximations.
fn debye_energy_derivatives(
    theta: f64,
    temperature: f64,
    n: f64,
) -> Result<(f64, f64, f64), FitError> {
    let y = theta / temperature;
    let debye = debye_function_3(y)?;
    let derivative = if y < 1.0e-3 {
        -3.0 / 8.0 + y / 10.0 - y.powi(3) / 420.0
    } else if y > 150.0 {
        -3.0 * debye / y
    } else {
        3.0 / y.exp_m1() - 3.0 * debye / y
    };
    let factor = 3.0 * n * GAS_CONSTANT;
    Ok((
        factor * temperature * debye,
        factor * derivative,
        factor * (debye - y * derivative),
    ))
}

/// Differentiate q-compromise Debye pressure with a BM3 reference isotherm.
///
/// Supports `rt_eos.V0`, `rt_eos.K0`, `rt_eos.K0_prime`, `Tr`, `theta0`,
/// `gamma0` and `n`, in any subset/order. Volumes are J/bar/mol, pressures `GPa`,
/// temperatures kelvin and `n` counts atoms per formula unit. The derivative
/// of reference volume includes the thermal-pressure denominator.
///
/// # Errors
/// Returns an error for invalid states, unknown names or non-finite results.
pub fn q_compromise_bm3_pressure_derivatives<S: AsRef<str>>(
    model: &QCompromiseDebye<BM3>,
    volume: f64,
    temperature: f64,
    parameter_names: &[S],
) -> Result<PressureDerivatives, FitError> {
    let temperature = positive_state(temperature, "temperature")?;
    let cold = bm3_pressure_derivatives(&model.rt_eos, volume, &["V0", "K0", "K0_prime"])?;
    let (energy, theta_derivative, temperature_derivative) =
        debye_energy_derivatives(model.theta0, temperature, model.n)?;
    let (reference, reference_theta_derivative, reference_temperature_derivative) =
        debye_energy_derivatives(model.theta0, model.tr, model.n)?;
    let denominator = model.rt_eos.v0 * 1.0e4;
    let factor = model.gamma0 / denominator;
    let thermal_pressure = model.gamma0 * (energy - reference) / model.rt_eos.v0 / 1.0e4;
    let parameters = parameter_names
        .iter()
        .map(|name| match name.as_ref() {
            "rt_eos.V0" => Ok(cold.parameters[0] - thermal_pressure / model.rt_eos.v0),
            "rt_eos.K0" => Ok(cold.parameters[1]),
            "rt_eos.K0_prime" => Ok(cold.parameters[2]),
            "Tr" => Ok(-factor * reference_temperature_derivative),
            "theta0" => Ok(factor * (theta_derivative - reference_theta_derivative)),
            "gamma0" => Ok((energy - reference) / denominator),
            "n" => Ok(thermal_pressure / model.n),
            name => Err(FitError::InvalidInput(format!(
                "unknown q-compromise BM3 derivative parameter {name:?}"
            ))),
        })
        .collect::<Result<Vec<_>, _>>()?;
    checked(PressureDerivatives {
        parameters,
        volume: cold.volume,
        temperature: factor * temperature_derivative,
    })
}

fn checked(derivatives: PressureDerivatives) -> Result<PressureDerivatives, FitError> {
    if derivatives
        .parameters
        .iter()
        .chain([&derivatives.volume, &derivatives.temperature])
        .any(|value| !value.is_finite())
    {
        return Err(FitError::InvalidInput(
            "pressure derivatives must be finite".to_owned(),
        ));
    }
    Ok(derivatives)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{IsothermalEos, ThermalEos};

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2.0e-8 * expected.abs().max(1.0e-7),
            "{actual:e} vs {expected:e}"
        );
    }

    fn difference(evaluate: impl Fn(f64) -> f64, at: f64) -> f64 {
        let step = 1.0e-4 * at.abs().max(1.0);
        (evaluate(at - 2.0 * step) - 8.0 * evaluate(at - step) + 8.0 * evaluate(at + step)
            - evaluate(at + 2.0 * step))
            / (12.0 * step)
    }

    #[test]
    fn bm3_partials_match_independent_pressure_differences() {
        for ratio in [0.65, 0.95, 1.0, 1.1] {
            let values = [11.3, 171.0, 4.7];
            let model = BM3::new(values[0], values[1], values[2]).unwrap();
            let volume = ratio * model.v0;
            let gradients =
                bm3_pressure_derivatives(&model, volume, &["K0_prime", "V0", "K0"]).unwrap();
            for (column, index) in [2, 0, 1].into_iter().enumerate() {
                let numerical = difference(
                    |value| {
                        let mut changed = values;
                        changed[index] = value;
                        BM3::new(changed[0], changed[1], changed[2])
                            .unwrap()
                            .pressure(volume)
                            .unwrap()
                    },
                    values[index],
                );
                // At V0 the K0 and Kprime derivatives vanish exactly.
                if numerical.abs() < 1.0e-10 {
                    assert!(gradients.parameters[column].abs() < 1.0e-10);
                } else {
                    close(gradients.parameters[column], numerical);
                }
            }
            close(
                gradients.volume,
                difference(|v| model.pressure(v).unwrap(), volume),
            );
        }
    }

    #[test]
    fn reference_state_derivatives_have_exact_limits() {
        let model = BM3::new(11.3, 171.0, 4.7).unwrap();
        let gradient =
            bm3_pressure_derivatives(&model, model.v0, &["V0", "K0", "K0_prime"]).unwrap();
        assert_eq!(gradient.parameters, [171.0 / 11.3, 0.0, 0.0]);
        assert_eq!(gradient.volume, -171.0 / 11.3);
    }

    #[test]
    fn q_compromise_partials_include_reference_denominator_and_zero_gamma() {
        let names = [
            "theta0",
            "rt_eos.V0",
            "Tr",
            "gamma0",
            "rt_eos.K0_prime",
            "n",
            "rt_eos.K0",
        ];
        for gamma in [0.0, 1.1, -0.3] {
            let values = [820.0, 11.3, 298.0, gamma, 4.7, 20.0, 171.0];
            let make = |p: &[f64]| {
                QCompromiseDebye::new(BM3::new(p[1], p[6], p[4]).unwrap(), p[2], p[0], p[3], p[5])
                    .unwrap()
            };
            let model = make(&values);
            for temperature in [298.0, 900.0, 2500.0] {
                let volume = 10.1;
                let gradient =
                    q_compromise_bm3_pressure_derivatives(&model, volume, temperature, &names)
                        .unwrap();
                for index in 0..names.len() {
                    let numerical = difference(
                        |value| {
                            let mut changed = values;
                            changed[index] = value;
                            make(&changed).pressure(volume, temperature).unwrap()
                        },
                        values[index],
                    );
                    // Subtraction at Tr and zero gamma produces exact zero partials.
                    if numerical.abs() < 1.0e-10 {
                        assert!(gradient.parameters[index].abs() < 1.0e-10);
                    } else {
                        close(gradient.parameters[index], numerical);
                    }
                }
                close(
                    gradient.volume,
                    difference(|v| model.pressure(v, temperature).unwrap(), volume),
                );
                if gamma == 0.0 {
                    assert_eq!(gradient.temperature, 0.0);
                } else {
                    close(
                        gradient.temperature,
                        difference(|t| model.pressure(volume, t).unwrap(), temperature),
                    );
                }
            }
        }
    }

    #[test]
    fn debye_energy_partials_follow_all_kernel_branches() {
        for theta in [0.1, 800.0, 100_000.0] {
            let temperature = 300.0;
            let n = 20.0;
            let (_, dtheta, dt) = debye_energy_derivatives(theta, temperature, n).unwrap();
            let energy =
                |th: f64, t: f64| 3.0 * n * GAS_CONSTANT * t * debye_function_3(th / t).unwrap();
            close(dtheta, difference(|th| energy(th, temperature), theta));
            close(dt, difference(|t| energy(theta, t), temperature));
        }
    }

    #[test]
    fn derivative_helpers_reject_invalid_states_and_unknown_names() {
        let cold = BM3::new(11.3, 171.0, 4.7).unwrap();
        assert!(bm3_pressure_derivatives(&cold, 0.0, &["K0"]).is_err());
        assert!(bm3_pressure_derivatives(&cold, 10.0, &["Kp"]).is_err());
        let thermal = QCompromiseDebye::new(cold, 300.0, 800.0, 1.0, 20.0).unwrap();
        assert!(q_compromise_bm3_pressure_derivatives(&thermal, 10.0, 0.0, &["theta0"]).is_err());
        assert!(q_compromise_bm3_pressure_derivatives(&thermal, 10.0, 800.0, &["q"]).is_err());
    }
}
