use peritheos::fit::{propagate_model_uncertainty, FitError};
use peritheos::{isothermal::BM3, IsothermalEos};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let parameters = [160.0, 4.0];
    // Row-major covariance of K0 and K0_prime.
    let covariance = [4.0, -0.05, -0.05, 0.04];
    let volumes = [9.5, 9.0, 8.5];

    let result = propagate_model_uncertainty(
        &parameters,
        &covariance,
        &[0.0; 3],
        1.0e-5,
        true,
        |parameters| {
            let eos = BM3::new(10.0, parameters[0], parameters[1])
                .map_err(|error| FitError::Evaluation(error.to_string()))?;
            volumes
                .iter()
                .map(|&volume| {
                    eos.pressure(volume)
                        .map_err(|error| FitError::Evaluation(error.to_string()))
                })
                .collect()
        },
    )?;

    println!("volume\tpressure (GPa)\tstandard uncertainty (GPa)");
    for ((volume, pressure), variance) in volumes
        .iter()
        .zip(&result.model.nominal)
        .zip(&result.propagation.variance)
    {
        println!("{volume:.3}\t{pressure:.3}\t\t{:.3}", variance.sqrt());
    }
    Ok(())
}
