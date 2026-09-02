use peritheos::fit::{fit_joint_eos, FitError, SolverOptions, ThermalObservations};
use peritheos::{isothermal::BM3, thermal::LinearThermalPressure, ThermalEos};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let truth = LinearThermalPressure::new(BM3::new(10.0, 160.0, 4.2)?, 300.0, 0.005)?;
    let volume = [9.8, 9.5, 9.2, 8.9, 9.7, 9.4, 9.1, 8.8];
    let temperature = [
        300.0, 300.0, 300.0, 300.0, 1_000.0, 1_200.0, 1_400.0, 1_600.0,
    ];
    let pressure = volume
        .iter()
        .zip(temperature)
        .map(|(&volume, temperature)| truth.pressure(volume, temperature))
        .collect::<Result<Vec<_>, _>>()?;
    let pressure_sigma = [0.1; 8];

    // Factory order is K0, K0_prime, alpha_KT. V0 and Tr stay fixed.
    let result = fit_joint_eos(
        ThermalObservations {
            pressure: &pressure,
            volume: &volume,
            temperature: &temperature,
            pressure_sigma: &pressure_sigma,
            volume_sigma: None,
            temperature_sigma: None,
            observation_cholesky: None,
        },
        &[140.0, 4.0, 0.003],
        &[50.0, 1.0, 0.0],
        &[300.0, 10.0, 0.02],
        SolverOptions::default(),
        |parameters| {
            let reference = BM3::new(10.0, parameters[0], parameters[1])
                .map_err(|error| FitError::Evaluation(error.to_string()))?;
            LinearThermalPressure::new(reference, 300.0, parameters[2])
                .map_err(|error| FitError::Evaluation(error.to_string()))
        },
    )?;

    println!("K0 = {:.3} GPa", result.parameters[0]);
    println!("K0' = {:.3}", result.parameters[1]);
    println!("alpha_KT = {:.6} GPa K^-1", result.parameters[2]);
    Ok(())
}
