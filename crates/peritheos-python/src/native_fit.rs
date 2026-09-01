//! Thin `PyO3` adapters for end-to-end native EOS fitting.

use numpy::PyReadonlyArrayDyn;
use peritheos_core::isothermal::{
    Holzapfel, ModifiedTait, Murnaghan, NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2,
    BM3, BM4,
};
use peritheos_core::thermal::{
    AsymptoticPowerLawMieGruneisenDebye, LinearThermalPressure, LogVolumeThermalPressure,
    MieGruneisenDebye, MieGruneisenEinstein, MultiOscillatorGruneisen, SokolovaParameters,
    ThermalModifiedTait, ThermalReferenceState,
};
use peritheos_core::{EosResult, ThermalEos};
use peritheos_fit::{
    fit_isothermal_eos, fit_thermal_eos_by, FitError, IsothermalObservations, Loss, SolverOptions,
    StructuredLayout, ThermalObservations,
};
use pyo3::prelude::*;

use super::{
    to_python_fit_error, PyLeastSquaresResult, PyRtEos, PyThermalEos, RtModel, ThermalModel,
};

fn array_values(array: PyReadonlyArrayDyn<'_, f64>) -> Vec<f64> {
    array.as_array().iter().copied().collect()
}

fn optional_array_values(array: Option<PyReadonlyArrayDyn<'_, f64>>) -> Option<Vec<f64>> {
    array.map(array_values)
}

fn ensure_names(names: &[String], allowed: &[&str], allow_reference: bool) -> Result<(), FitError> {
    if let Some(name) = names.iter().find(|name| {
        !(allowed.contains(&name.as_str()) || allow_reference && name.starts_with("rt_eos."))
    }) {
        return Err(FitError::InvalidInput(format!(
            "parameter {name:?} is unavailable for the native EOS"
        )));
    }
    Ok(())
}

fn value(names: &[String], values: &[f64], name: &str, current: f64) -> f64 {
    names
        .iter()
        .position(|candidate| candidate == name)
        .map_or(current, |index| values[index])
}

impl RtModel {
    #[allow(clippy::too_many_lines)]
    fn with_parameters(self, names: &[String], values: &[f64]) -> Result<Self, FitError> {
        if names.len() != values.len() {
            return Err(FitError::InvalidInput(
                "native EOS parameter names and values have different lengths".to_owned(),
            ));
        }
        let model = match self {
            Self::BM2(model) => {
                ensure_names(names, &["V0", "K0"], false)?;
                Self::BM2(
                    BM2::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::BM3(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime"], false)?;
                Self::BM3(
                    BM3::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::BM4(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime", "K0_double_prime"], false)?;
                Self::BM4(
                    BM4::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                        value(names, values, "K0_double_prime", model.k0_double_prime),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::Murnaghan(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime"], false)?;
                Self::Murnaghan(
                    Murnaghan::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::ModifiedTait(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime", "K0_double_prime"], false)?;
                Self::ModifiedTait(
                    ModifiedTait::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                        value(names, values, "K0_double_prime", model.k0_double_prime),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::NaturalStrain2(model) => {
                ensure_names(names, &["V0", "K0"], false)?;
                Self::NaturalStrain2(
                    NaturalStrain2::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::NaturalStrain3(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime"], false)?;
                Self::NaturalStrain3(
                    NaturalStrain3::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::NaturalStrain4(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime", "K0_double_prime"], false)?;
                Self::NaturalStrain4(
                    NaturalStrain4::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                        value(names, values, "K0_double_prime", model.k0_double_prime),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::Vinet(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime"], false)?;
                Self::Vinet(
                    Vinet::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::Holzapfel(model) => {
                ensure_names(names, &["V0", "K0", "K0_prime", "n", "Z"], false)?;
                Self::Holzapfel(
                    Holzapfel::new(
                        value(names, values, "V0", model.v0),
                        value(names, values, "K0", model.k0),
                        value(names, values, "K0_prime", model.k0_prime),
                        value(names, values, "n", model.n),
                        value(names, values, "Z", model.z),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
        };
        Ok(model)
    }
}

fn reference_updates(names: &[String], values: &[f64]) -> (Vec<String>, Vec<f64>) {
    names
        .iter()
        .zip(values)
        .filter_map(|(name, value)| {
            name.strip_prefix("rt_eos.")
                .map(|name| (name.to_owned(), *value))
        })
        .unzip()
}

impl ThermalModel {
    fn pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        match self {
            Self::AsymptoticPowerLawMieGruneisenDebye(model) => model.pressure(volume, temperature),
            Self::LinearThermalPressure(model) => model.pressure(volume, temperature),
            Self::LogVolumeThermalPressure(model) => model.pressure(volume, temperature),
            Self::MieGruneisenDebye(model) => model.pressure(volume, temperature),
            Self::MieGruneisenEinstein(model) => model.pressure(volume, temperature),
            Self::ThermalModifiedTait(model) => model.pressure(volume, temperature),
            Self::Sokolova2016(model) => model.pressure(volume, temperature),
            Self::ThermalReferenceState(model) => model.pressure(volume, temperature),
        }
    }

    #[allow(clippy::too_many_lines)]
    fn with_parameters(self, names: &[String], values: &[f64]) -> Result<Self, FitError> {
        if names.len() != values.len() {
            return Err(FitError::InvalidInput(
                "native EOS parameter names and values have different lengths".to_owned(),
            ));
        }
        let (reference_names, reference_values) = reference_updates(names, values);
        let model = match self {
            Self::AsymptoticPowerLawMieGruneisenDebye(model) => {
                ensure_names(names, &["Tr", "theta0", "gamma0", "a", "b", "n"], true)?;
                let reference = model
                    .rt_eos
                    .with_parameters(&reference_names, &reference_values)?;
                Self::AsymptoticPowerLawMieGruneisenDebye(
                    AsymptoticPowerLawMieGruneisenDebye::new(
                        reference,
                        value(names, values, "Tr", model.tr),
                        value(names, values, "theta0", model.theta0),
                        value(names, values, "gamma0", model.gamma0),
                        value(names, values, "a", model.a),
                        value(names, values, "b", model.b),
                        value(names, values, "n", model.n),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::LinearThermalPressure(model) => {
                ensure_names(names, &["Tr", "alpha_KT"], true)?;
                let reference = model
                    .rt_eos
                    .with_parameters(&reference_names, &reference_values)?;
                Self::LinearThermalPressure(
                    LinearThermalPressure::new(
                        reference,
                        value(names, values, "Tr", model.tr),
                        value(names, values, "alpha_KT", model.alpha_kt),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::LogVolumeThermalPressure(model) => {
                ensure_names(names, &["Tr", "alpha_KT_ref", "dK_dT_V"], true)?;
                let reference = model
                    .rt_eos
                    .with_parameters(&reference_names, &reference_values)?;
                Self::LogVolumeThermalPressure(
                    LogVolumeThermalPressure::new(
                        reference,
                        value(names, values, "Tr", model.tr),
                        value(names, values, "alpha_KT_ref", model.alpha_kt_ref),
                        value(names, values, "dK_dT_V", model.dk_dt_v),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::MieGruneisenDebye(model) => {
                ensure_names(names, &["Tr", "theta0", "gamma0", "q", "n"], true)?;
                let reference = model
                    .rt_eos
                    .with_parameters(&reference_names, &reference_values)?;
                Self::MieGruneisenDebye(
                    MieGruneisenDebye::new_with_temperature_law(
                        reference,
                        value(names, values, "Tr", model.tr),
                        value(names, values, "theta0", model.theta0),
                        value(names, values, "gamma0", model.gamma0),
                        value(names, values, "q", model.q),
                        value(names, values, "n", model.n),
                        model.debye_temperature_law,
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::MieGruneisenEinstein(model) => {
                ensure_names(names, &["Tr", "theta0", "gamma0", "q", "n"], true)?;
                let reference = model
                    .rt_eos
                    .with_parameters(&reference_names, &reference_values)?;
                Self::MieGruneisenEinstein(
                    MieGruneisenEinstein::new(
                        reference,
                        value(names, values, "Tr", model.tr),
                        value(names, values, "theta0", model.theta0),
                        value(names, values, "gamma0", model.gamma0),
                        value(names, values, "q", model.q),
                        value(names, values, "n", model.n),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::ThermalModifiedTait(model) => {
                ensure_names(names, &["Tr", "theta", "alpha0", "n"], true)?;
                let reference = RtModel::ModifiedTait(model.rt_eos)
                    .with_parameters(&reference_names, &reference_values)?;
                let RtModel::ModifiedTait(reference) = reference else {
                    unreachable!("the reference EOS variant cannot change")
                };
                Self::ThermalModifiedTait(
                    ThermalModifiedTait::new(
                        reference,
                        value(names, values, "Tr", model.tr),
                        value(names, values, "theta", model.theta),
                        value(names, values, "alpha0", model.alpha0),
                        value(names, values, "n", model.n),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::Sokolova2016(model) => {
                ensure_names(
                    names,
                    &[
                        "Tr", "QE1o", "mE1", "QE2o", "mE2", "delta", "t", "a_0", "m", "g", "e_0",
                        "beta", "QBo", "d", "mb", "QB1o", "d1", "mb1", "n",
                    ],
                    true,
                )?;
                let reference = model
                    .rt_eos
                    .with_parameters(&reference_names, &reference_values)?;
                let current = model.parameters;
                let parameters = SokolovaParameters {
                    tr: value(names, values, "Tr", current.tr),
                    qe1o: value(names, values, "QE1o", current.qe1o),
                    me1: value(names, values, "mE1", current.me1),
                    qe2o: value(names, values, "QE2o", current.qe2o),
                    me2: value(names, values, "mE2", current.me2),
                    delta: value(names, values, "delta", current.delta),
                    t: value(names, values, "t", current.t),
                    a_0: value(names, values, "a_0", current.a_0),
                    m: value(names, values, "m", current.m),
                    g: value(names, values, "g", current.g),
                    e_0: value(names, values, "e_0", current.e_0),
                    beta: value(names, values, "beta", current.beta),
                    qbo: value(names, values, "QBo", current.qbo),
                    d: value(names, values, "d", current.d),
                    mb: value(names, values, "mb", current.mb),
                    qb1o: value(names, values, "QB1o", current.qb1o),
                    d1: value(names, values, "d1", current.d1),
                    mb1: value(names, values, "mb1", current.mb1),
                };
                Self::Sokolova2016(
                    MultiOscillatorGruneisen::new_with_atom_count(
                        reference,
                        parameters,
                        value(names, values, "n", model.n),
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
            Self::ThermalReferenceState(model) => {
                ensure_names(names, &["Tr", "alpha0", "dK_dT", "alpha1"], true)?;
                let reference = model
                    .rt_eos
                    .with_parameters(&reference_names, &reference_values)?;
                Self::ThermalReferenceState(
                    ThermalReferenceState::new(
                        reference,
                        value(names, values, "Tr", model.tr),
                        value(names, values, "alpha0", model.alpha0),
                        value(names, values, "dK_dT", model.dk_dt),
                        value(names, values, "alpha1", model.alpha1),
                        model.thermal_expansion_law,
                        model.reference_volume_law,
                    )
                    .map_err(|error| FitError::Evaluation(error.to_string()))?,
                )
            }
        };
        Ok(model)
    }
}

#[pyfunction]
#[pyo3(signature = (
    model, parameter_names, initial, lower, upper, pressure, volume,
    pressure_sigma, volume_sigma=None, observation_cholesky=None,
    loss="linear", f_scale=1.0, max_nfev=None
))]
#[allow(clippy::too_many_arguments)]
pub(super) fn fit_rt_eos_native(
    py: Python<'_>,
    model: PyRef<'_, PyRtEos>,
    parameter_names: Vec<String>,
    initial: PyReadonlyArrayDyn<'_, f64>,
    lower: PyReadonlyArrayDyn<'_, f64>,
    upper: PyReadonlyArrayDyn<'_, f64>,
    pressure: PyReadonlyArrayDyn<'_, f64>,
    volume: PyReadonlyArrayDyn<'_, f64>,
    pressure_sigma: PyReadonlyArrayDyn<'_, f64>,
    volume_sigma: Option<PyReadonlyArrayDyn<'_, f64>>,
    observation_cholesky: Option<PyReadonlyArrayDyn<'_, f64>>,
    loss: &str,
    f_scale: f64,
    max_nfev: Option<usize>,
) -> PyResult<PyLeastSquaresResult> {
    let model = model.model;
    let initial = array_values(initial);
    let lower = array_values(lower);
    let upper = array_values(upper);
    let pressure = array_values(pressure);
    let volume = array_values(volume);
    let pressure_sigma = array_values(pressure_sigma);
    let volume_sigma = optional_array_values(volume_sigma);
    let observation_cholesky = optional_array_values(observation_cholesky);
    let global_parameter_count = initial.len();
    let structured_layout = volume_sigma.as_ref().map(|_| StructuredLayout {
        global_parameter_count,
        point_count: pressure.len(),
        latent_coordinate_count: 1,
    });
    let options = SolverOptions {
        loss: Loss::from_name(loss).map_err(to_python_fit_error)?,
        f_scale,
        max_evaluations: max_nfev,
        ..SolverOptions::default()
    };
    let result = py
        .detach(move || {
            fit_isothermal_eos(
                IsothermalObservations {
                    pressure: &pressure,
                    volume: &volume,
                    pressure_sigma: &pressure_sigma,
                    volume_sigma: volume_sigma.as_deref(),
                    observation_cholesky: observation_cholesky.as_deref(),
                },
                &initial,
                &lower,
                &upper,
                options,
                |parameters| model.with_parameters(&parameter_names, parameters),
            )
        })
        .map_err(to_python_fit_error)?;
    Ok(PyLeastSquaresResult {
        global_parameter_count,
        structured_layout,
        predicted_pressure: Some(result.predicted_pressure),
        result: result.solver,
    })
}

#[pyfunction]
#[pyo3(signature = (
    model, parameter_names, initial, lower, upper, pressure, volume,
    temperature, pressure_sigma, volume_sigma=None, temperature_sigma=None,
    observation_cholesky=None, loss="linear", f_scale=1.0, max_nfev=None
))]
#[allow(clippy::too_many_arguments)]
pub(super) fn fit_thermal_eos_native(
    py: Python<'_>,
    model: PyRef<'_, PyThermalEos>,
    parameter_names: Vec<String>,
    initial: PyReadonlyArrayDyn<'_, f64>,
    lower: PyReadonlyArrayDyn<'_, f64>,
    upper: PyReadonlyArrayDyn<'_, f64>,
    pressure: PyReadonlyArrayDyn<'_, f64>,
    volume: PyReadonlyArrayDyn<'_, f64>,
    temperature: PyReadonlyArrayDyn<'_, f64>,
    pressure_sigma: PyReadonlyArrayDyn<'_, f64>,
    volume_sigma: Option<PyReadonlyArrayDyn<'_, f64>>,
    temperature_sigma: Option<PyReadonlyArrayDyn<'_, f64>>,
    observation_cholesky: Option<PyReadonlyArrayDyn<'_, f64>>,
    loss: &str,
    f_scale: f64,
    max_nfev: Option<usize>,
) -> PyResult<PyLeastSquaresResult> {
    let model = model.model;
    let initial = array_values(initial);
    let lower = array_values(lower);
    let upper = array_values(upper);
    let pressure = array_values(pressure);
    let volume = array_values(volume);
    let temperature = array_values(temperature);
    let pressure_sigma = array_values(pressure_sigma);
    let volume_sigma = optional_array_values(volume_sigma);
    let temperature_sigma = optional_array_values(temperature_sigma);
    let observation_cholesky = optional_array_values(observation_cholesky);
    let global_parameter_count = initial.len();
    let latent_coordinate_count =
        usize::from(volume_sigma.is_some()) + usize::from(temperature_sigma.is_some());
    let structured_layout = (latent_coordinate_count > 0).then_some(StructuredLayout {
        global_parameter_count,
        point_count: pressure.len(),
        latent_coordinate_count,
    });
    let options = SolverOptions {
        loss: Loss::from_name(loss).map_err(to_python_fit_error)?,
        f_scale,
        max_evaluations: max_nfev,
        ..SolverOptions::default()
    };
    let result = py
        .detach(move || {
            fit_thermal_eos_by(
                ThermalObservations {
                    pressure: &pressure,
                    volume: &volume,
                    temperature: &temperature,
                    pressure_sigma: &pressure_sigma,
                    volume_sigma: volume_sigma.as_deref(),
                    temperature_sigma: temperature_sigma.as_deref(),
                    observation_cholesky: observation_cholesky.as_deref(),
                },
                &initial,
                &lower,
                &upper,
                options,
                |parameters| model.with_parameters(&parameter_names, parameters),
                |model, volume, temperature| {
                    model
                        .pressure(volume, temperature)
                        .map_err(|error| FitError::Evaluation(error.to_string()))
                },
            )
        })
        .map_err(to_python_fit_error)?;
    Ok(PyLeastSquaresResult {
        global_parameter_count,
        structured_layout,
        predicted_pressure: Some(result.predicted_pressure),
        result: result.solver,
    })
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn names(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    }

    #[test]
    fn rt_reconstruction_maps_all_parameter_conventions() {
        let bm4 = RtModel::BM4(BM4::new(10.0, 100.0, 4.0, -0.01).unwrap())
            .with_parameters(
                &names(&["V0", "K0", "K0_prime", "K0_double_prime"]),
                &[11.0, 120.0, 4.5, -0.02],
            )
            .unwrap();
        let RtModel::BM4(bm4) = bm4 else {
            panic!("model variant changed")
        };
        assert_eq!((bm4.v0, bm4.k0, bm4.k0_prime), (11.0, 120.0, 4.5));
        assert_eq!(bm4.k0_double_prime, -0.02);

        let holzapfel = RtModel::Holzapfel(Holzapfel::new(0.3414, 441.5, 3.9, 1.0, 6.0).unwrap())
            .with_parameters(
                &names(&["V0", "K0", "K0_prime", "n", "Z"]),
                &[0.35, 450.0, 4.0, 2.0, 8.0],
            )
            .unwrap();
        let RtModel::Holzapfel(holzapfel) = holzapfel else {
            panic!("model variant changed")
        };
        assert_eq!(
            (
                holzapfel.v0,
                holzapfel.k0,
                holzapfel.k0_prime,
                holzapfel.n,
                holzapfel.z,
            ),
            (0.35, 450.0, 4.0, 2.0, 8.0)
        );
    }

    #[test]
    fn thermal_reconstruction_updates_own_and_reference_parameters() {
        let reference = RtModel::BM3(BM3::new(1.0, 160.0, 4.0).unwrap());
        let model = ThermalModel::MieGruneisenEinstein(
            MieGruneisenEinstein::new(reference, 300.0, 800.0, 1.5, 1.0, 2.0).unwrap(),
        )
        .with_parameters(
            &names(&[
                "rt_eos.V0",
                "rt_eos.K0",
                "rt_eos.K0_prime",
                "Tr",
                "theta0",
                "gamma0",
                "q",
                "n",
            ]),
            &[1.1, 170.0, 4.2, 298.15, 900.0, 1.7, 1.2, 3.0],
        )
        .unwrap();
        let ThermalModel::MieGruneisenEinstein(model) = model else {
            panic!("model variant changed")
        };
        let RtModel::BM3(reference) = model.rt_eos else {
            panic!("reference variant changed")
        };
        assert_eq!(
            (reference.v0, reference.k0, reference.k0_prime),
            (1.1, 170.0, 4.2)
        );
        assert_eq!(
            (model.tr, model.theta0, model.gamma0, model.q, model.n),
            (298.15, 900.0, 1.7, 1.2, 3.0)
        );
    }

    #[test]
    fn sokolova_reconstruction_preserves_case_sensitive_public_names() {
        let reference = Holzapfel::new(0.3414, 441.5, 3.9, 1.0, 6.0).unwrap();
        let current = SokolovaParameters::reduced(
            298.15, 684.0, 0.564, 1561.0, 2.436, -0.506, 1.085, 0.0, 0.0, 0.0, 0.0,
        );
        let model = ThermalModel::Sokolova2016(
            MultiOscillatorGruneisen::new_with_atom_count(
                RtModel::Holzapfel(reference),
                current,
                1.0,
            )
            .unwrap(),
        )
        .with_parameters(
            &names(&["QE1o", "mE1", "QBo", "QB1o", "mb1"]),
            &[700.0, 0.6, 500.0, 1200.0, 0.4],
        )
        .unwrap();
        let ThermalModel::Sokolova2016(model) = model else {
            panic!("model variant changed")
        };
        assert_eq!(model.parameters.qe1o, 700.0);
        assert_eq!(model.parameters.me1, 0.6);
        assert_eq!(model.parameters.qbo, 500.0);
        assert_eq!(model.parameters.qb1o, 1200.0);
        assert_eq!(model.parameters.mb1, 0.4);
    }
}
