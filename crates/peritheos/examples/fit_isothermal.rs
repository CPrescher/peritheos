use peritheos::fit::{
    fit_isothermal_eos, parameter_covariance, FitError, IsothermalObservations, Loss, SolverOptions,
};
use peritheos::{isothermal::BM3, IsothermalEos};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let truth = BM3::new(10.0, 160.0, 4.2)?;
    let volume = [9.8, 9.5, 9.2, 8.9, 8.6, 8.3];
    let pressure = volume
        .iter()
        .map(|&value| truth.pressure(value))
        .collect::<Result<Vec<_>, _>>()?;
    let pressure_sigma = [0.2; 6];

    let result = fit_isothermal_eos(
        IsothermalObservations {
            pressure: &pressure,
            volume: &volume,
            pressure_sigma: &pressure_sigma,
            volume_sigma: None,
            observation_cholesky: None,
        },
        &[140.0, 4.0],
        &[50.0, 1.0],
        &[300.0, 10.0],
        SolverOptions {
            loss: Loss::SoftL1,
            ..SolverOptions::default()
        },
        |parameters| {
            BM3::new(10.0, parameters[0], parameters[1])
                .map_err(|error| FitError::Evaluation(error.to_string()))
        },
    )?;

    let parameter_count = 2;
    let covariance = parameter_covariance(
        &result.solver.jacobian,
        result.solver.residual_count,
        result.solver.parameters.len(),
        parameter_count,
    )?;

    println!("success: {}", result.solver.message);
    println!("K0 = {:.3} GPa", result.solver.parameters[0]);
    println!("K0' = {:.3}", result.solver.parameters[1]);
    println!(
        "standard errors: {:.3}, {:.3}",
        covariance[0].sqrt(),
        covariance[3].sqrt()
    );
    Ok(())
}
