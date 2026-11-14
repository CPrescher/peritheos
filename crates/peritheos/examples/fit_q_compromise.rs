//! Analytic joint fitting with latent volume and temperature observations.
use peritheos::fit::{
    fit_thermal_eos_with_derivatives, q_compromise_bm3_pressure_derivatives, FitError,
    SolverOptions, ThermalObservations,
};
use peritheos::{isothermal::BM3, thermal::QCompromiseDebye, ThermalEos};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Volumes: J/bar/mol; pressures: GPa; temperatures: K; n: atoms/formula unit.
    let truth = QCompromiseDebye::new(BM3::new(11.3, 171.0, 4.7)?, 298.0, 820.0, 1.1, 20.0)?;
    let volume = [11.0, 10.0, 9.0, 11.0, 10.0, 9.0, 11.0, 10.0, 9.0];
    let temperature = [
        300.0, 300.0, 300.0, 900.0, 900.0, 900.0, 1800.0, 1800.0, 1800.0,
    ];
    let pressure = volume
        .iter()
        .zip(temperature)
        .enumerate()
        .map(|(index, (&v, t))| {
            // Explicit deterministic synthetic measurement perturbation.
            Ok(truth.pressure(v, t)? + if index % 2 == 0 { 0.01 } else { -0.01 })
        })
        .collect::<Result<Vec<_>, peritheos::EosError>>()?;
    let result = fit_thermal_eos_with_derivatives(
        ThermalObservations {
            pressure: &pressure,
            volume: &volume,
            temperature: &temperature,
            pressure_sigma: &[0.04; 9],
            volume_sigma: Some(&[0.008; 9]),
            temperature_sigma: Some(&[4.0; 9]),
            observation_cholesky: None,
        },
        &[163.0, 750.0, 1.0],
        &[100.0, 100.0, 0.2],
        &[220.0, 1800.0, 2.0],
        SolverOptions::default(),
        |p| {
            let reference = BM3::new(11.3, p[0], 4.7)?;
            QCompromiseDebye::new(reference, 298.0, p[1], p[2], 20.0).map_err(FitError::from)
        },
        |model, v, t| {
            q_compromise_bm3_pressure_derivatives(model, v, t, &["rt_eos.K0", "theta0", "gamma0"])
        },
    )?;
    for (index, name) in ["K0", "theta0", "gamma0"].iter().enumerate() {
        println!(
            "{name} = {:.6} ± {:.6}",
            result.parameters[index], result.standard_errors[index]
        );
    }
    println!(
        "chi-square = {:.6}; {}",
        result.chi_square, result.solver.message
    );
    Ok(())
}
