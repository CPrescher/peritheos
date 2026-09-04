#![allow(clippy::needless_pass_by_value, clippy::similar_names)]

mod native_fit;

use numpy::ndarray::ArrayD;
use numpy::{IntoPyArray, PyArrayDyn, PyReadonlyArrayDyn};
use peritheos::fit::{
    least_squares, least_squares_structured, parameter_covariance_structured,
    propagate_linear_uncertainty, summarize_monte_carlo, FitError, LinearPropagation, Loss,
    MonteCarloSummary, SolverOptions, StructuredLayout,
};
use peritheos::hugoniot::{Hugoniot, LinearUsUpHugoniot};
use peritheos::isothermal::{
    holzapfel_bulk_modulus_derivative_analytical, Holzapfel, ModifiedTait, Murnaghan,
    NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2, BM3, BM4,
};
use peritheos::thermal::{
    AsymptoticPowerLawMieGruneisenDebye, DebyeTemperatureLaw, DorogokupetsOganov2007,
    DorogokupetsOganov2007Parameters, LinearThermalPressure, LogVolumeThermalPressure,
    MieGruneisenDebye, MieGruneisenEinstein, MultiOscillatorGruneisen, ReferenceStateEos,
    ReferenceVolumeLaw, SecondOrderTaylorThermalPressure, SokolovaParameters, ThermalExpansionLaw,
    ThermalModifiedTait, ThermalReferenceState,
};
use peritheos::{CaloricEos, EosError, EosResult, IsothermalEos, ThermalEos};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyModule};
use rayon::prelude::*;

const PARALLEL_ELEMENTWISE_THRESHOLD: usize = 65_536;
const PARALLEL_ROOT_THRESHOLD: usize = 512;
const PARALLEL_THERMAL_THRESHOLD: usize = 2_048;

#[derive(Clone, Copy, Debug, PartialEq)]
enum RtModel {
    BM2(BM2),
    BM3(BM3),
    BM4(BM4),
    Murnaghan(Murnaghan),
    ModifiedTait(ModifiedTait),
    NaturalStrain2(NaturalStrain2),
    NaturalStrain3(NaturalStrain3),
    NaturalStrain4(NaturalStrain4),
    Vinet(Vinet),
    Holzapfel(Holzapfel),
}

impl RtModel {
    fn name(self) -> &'static str {
        match self {
            Self::BM2(_) => "BM2",
            Self::BM3(_) => "BM3",
            Self::BM4(_) => "BM4",
            Self::Murnaghan(_) => "Murnaghan",
            Self::ModifiedTait(_) => "ModifiedTait",
            Self::NaturalStrain2(_) => "NaturalStrain2",
            Self::NaturalStrain3(_) => "NaturalStrain3",
            Self::NaturalStrain4(_) => "NaturalStrain4",
            Self::Vinet(_) => "Vinet",
            Self::Holzapfel(_) => "Holzapfel",
        }
    }
}

impl IsothermalEos for RtModel {
    fn reference_volume(&self) -> f64 {
        match self {
            Self::BM2(model) => model.reference_volume(),
            Self::BM3(model) => model.reference_volume(),
            Self::BM4(model) => model.reference_volume(),
            Self::Murnaghan(model) => model.reference_volume(),
            Self::ModifiedTait(model) => model.reference_volume(),
            Self::NaturalStrain2(model) => model.reference_volume(),
            Self::NaturalStrain3(model) => model.reference_volume(),
            Self::NaturalStrain4(model) => model.reference_volume(),
            Self::Vinet(model) => model.reference_volume(),
            Self::Holzapfel(model) => model.reference_volume(),
        }
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::BM2(model) => model.pressure(volume),
            Self::BM3(model) => model.pressure(volume),
            Self::BM4(model) => model.pressure(volume),
            Self::Murnaghan(model) => model.pressure(volume),
            Self::ModifiedTait(model) => model.pressure(volume),
            Self::NaturalStrain2(model) => model.pressure(volume),
            Self::NaturalStrain3(model) => model.pressure(volume),
            Self::NaturalStrain4(model) => model.pressure(volume),
            Self::Vinet(model) => model.pressure(volume),
            Self::Holzapfel(model) => model.pressure(volume),
        }
    }

    fn bulk_modulus(&self, volume: f64) -> EosResult<f64> {
        match self {
            Self::BM2(model) => model.bulk_modulus(volume),
            Self::BM3(model) => model.bulk_modulus(volume),
            Self::BM4(model) => model.bulk_modulus(volume),
            Self::Murnaghan(model) => model.bulk_modulus(volume),
            Self::ModifiedTait(model) => model.bulk_modulus(volume),
            Self::NaturalStrain2(model) => model.bulk_modulus(volume),
            Self::NaturalStrain3(model) => model.bulk_modulus(volume),
            Self::NaturalStrain4(model) => model.bulk_modulus(volume),
            Self::Vinet(model) => model.bulk_modulus(volume),
            Self::Holzapfel(model) => model.bulk_modulus(volume),
        }
    }
}

impl ReferenceStateEos for RtModel {
    fn reference_bulk_modulus(&self) -> f64 {
        match self {
            Self::BM2(model) => model.k0,
            Self::BM3(model) => model.k0,
            Self::BM4(model) => model.k0,
            Self::Murnaghan(model) => model.k0,
            Self::ModifiedTait(model) => model.k0,
            Self::NaturalStrain2(model) => model.k0,
            Self::NaturalStrain3(model) => model.k0,
            Self::NaturalStrain4(model) => model.k0,
            Self::Vinet(model) => model.k0,
            Self::Holzapfel(model) => model.k0,
        }
    }

    fn with_reference_state(&self, volume: f64, bulk_modulus: f64) -> EosResult<Self> {
        match self {
            Self::BM2(_) => BM2::new(volume, bulk_modulus).map(Self::BM2),
            Self::BM3(model) => BM3::new(volume, bulk_modulus, model.k0_prime).map(Self::BM3),
            Self::BM4(model) => {
                BM4::new(volume, bulk_modulus, model.k0_prime, model.k0_double_prime).map(Self::BM4)
            }
            Self::Murnaghan(model) => {
                Murnaghan::new(volume, bulk_modulus, model.k0_prime).map(Self::Murnaghan)
            }
            Self::ModifiedTait(model) => {
                ModifiedTait::new(volume, bulk_modulus, model.k0_prime, model.k0_double_prime)
                    .map(Self::ModifiedTait)
            }
            Self::NaturalStrain2(_) => {
                NaturalStrain2::new(volume, bulk_modulus).map(Self::NaturalStrain2)
            }
            Self::NaturalStrain3(model) => {
                NaturalStrain3::new(volume, bulk_modulus, model.k0_prime).map(Self::NaturalStrain3)
            }
            Self::NaturalStrain4(model) => {
                NaturalStrain4::new(volume, bulk_modulus, model.k0_prime, model.k0_double_prime)
                    .map(Self::NaturalStrain4)
            }
            Self::Vinet(model) => Vinet::new(volume, bulk_modulus, model.k0_prime).map(Self::Vinet),
            Self::Holzapfel(model) => {
                Holzapfel::new(volume, bulk_modulus, model.k0_prime, model.n, model.z)
                    .map(Self::Holzapfel)
            }
        }
    }
}

/// Private native representation of a built-in isothermal EOS.
#[pyclass(
    name = "RtEos",
    frozen,
    module = "peritheos._rust",
    skip_from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq)]
struct PyRtEos {
    model: RtModel,
}

#[pymethods]
impl PyRtEos {
    #[staticmethod]
    fn bm2(v0: f64, k0: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::BM2(BM2::new(v0, k0).map_err(to_python_error)?),
        })
    }

    #[staticmethod]
    fn bm3(v0: f64, k0: f64, k0_prime: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::BM3(BM3::new(v0, k0, k0_prime).map_err(to_python_error)?),
        })
    }

    #[staticmethod]
    fn bm4(v0: f64, k0: f64, k0_prime: f64, k0_double_prime: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::BM4(
                BM4::new(v0, k0, k0_prime, k0_double_prime).map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn murnaghan(v0: f64, k0: f64, k0_prime: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::Murnaghan(Murnaghan::new(v0, k0, k0_prime).map_err(to_python_error)?),
        })
    }

    #[staticmethod]
    fn modified_tait(v0: f64, k0: f64, k0_prime: f64, k0_double_prime: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::ModifiedTait(
                ModifiedTait::new(v0, k0, k0_prime, k0_double_prime).map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn natural_strain2(v0: f64, k0: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::NaturalStrain2(NaturalStrain2::new(v0, k0).map_err(to_python_error)?),
        })
    }

    #[staticmethod]
    fn natural_strain3(v0: f64, k0: f64, k0_prime: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::NaturalStrain3(
                NaturalStrain3::new(v0, k0, k0_prime).map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn natural_strain4(v0: f64, k0: f64, k0_prime: f64, k0_double_prime: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::NaturalStrain4(
                NaturalStrain4::new(v0, k0, k0_prime, k0_double_prime).map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn vinet(v0: f64, k0: f64, k0_prime: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::Vinet(Vinet::new(v0, k0, k0_prime).map_err(to_python_error)?),
        })
    }

    #[staticmethod]
    fn holzapfel(v0: f64, k0: f64, k0_prime: f64, n: f64, z: f64) -> PyResult<Self> {
        Ok(Self {
            model: RtModel::Holzapfel(
                Holzapfel::new(v0, k0, k0_prime, n, z).map_err(to_python_error)?,
            ),
        })
    }

    #[getter]
    fn model_name(&self) -> &'static str {
        self.model.name()
    }

    #[getter]
    fn reference_volume(&self) -> f64 {
        self.model.reference_volume()
    }

    fn pressure_scalar(&self, volume: f64) -> PyResult<f64> {
        self.model.pressure(volume).map_err(to_python_error)
    }

    fn bulk_modulus_scalar(&self, volume: f64) -> PyResult<f64> {
        self.model.bulk_modulus(volume).map_err(to_python_error)
    }

    fn bulk_modulus_derivative_scalar(&self, volume: f64, epsilon: f64) -> PyResult<f64> {
        match self.model {
            RtModel::Holzapfel(model) => model
                .bulk_modulus_derivative(volume, epsilon)
                .map_err(to_python_error),
            _ => Err(python_unsupported_error(
                "bulk-modulus derivative is only defined for Holzapfel",
            )),
        }
    }

    fn volume_scalar(&self, pressure: f64) -> PyResult<f64> {
        self.model.volume(pressure).map_err(to_python_error)
    }

    fn pressure_array<'py>(
        &self,
        py: Python<'py>,
        volumes: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        map_array(py, volumes, PARALLEL_ELEMENTWISE_THRESHOLD, |volume| {
            self.model.pressure(volume)
        })
    }

    fn bulk_modulus_array<'py>(
        &self,
        py: Python<'py>,
        volumes: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        map_array(py, volumes, PARALLEL_ELEMENTWISE_THRESHOLD, |volume| {
            self.model.bulk_modulus(volume)
        })
    }

    fn bulk_modulus_derivative_array<'py>(
        &self,
        py: Python<'py>,
        volumes: PyReadonlyArrayDyn<'py, f64>,
        epsilon: f64,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let RtModel::Holzapfel(model) = self.model else {
            return Err(python_unsupported_error(
                "bulk-modulus derivative is only defined for Holzapfel",
            ));
        };
        map_array(py, volumes, PARALLEL_ELEMENTWISE_THRESHOLD, |volume| {
            model.bulk_modulus_derivative(volume, epsilon)
        })
    }

    fn volume_array<'py>(
        &self,
        py: Python<'py>,
        pressures: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        map_array(py, pressures, PARALLEL_ROOT_THRESHOLD, |pressure| {
            self.model.volume(pressure)
        })
    }

    fn __repr__(&self) -> String {
        format!("RtEos(model_name='{}')", self.model.name())
    }
}

/// Private native representation of a built-in shock Hugoniot.
#[pyclass(
    name = "HugoniotEos",
    frozen,
    module = "peritheos._rust",
    skip_from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq)]
struct PyHugoniotEos {
    model: LinearUsUpHugoniot,
}

impl PyHugoniotEos {
    fn evaluate_model(&self, quantity: &str, value: f64) -> EosResult<f64> {
        match quantity {
            "pressure" => self.model.pressure(value),
            "volume" => self.model.volume(value),
            "particle_velocity" => self.model.particle_velocity(value),
            "shock_velocity" => self.model.shock_velocity(value),
            "density" => self.model.density(value),
            "specific_internal_energy_change" => self.model.specific_internal_energy_change(value),
            "tangent_modulus" => self.model.tangent_modulus(value),
            "shock_velocity_from_particle_velocity" => {
                self.model.shock_velocity_from_particle_velocity(value)
            }
            "pressure_from_particle_velocity" => self.model.pressure_from_particle_velocity(value),
            "volume_from_particle_velocity" => self.model.volume_from_particle_velocity(value),
            _ => Err(EosError::InvalidState {
                name: "quantity",
                reason: "unknown shock Hugoniot quantity",
            }),
        }
    }
}

#[pymethods]
impl PyHugoniotEos {
    #[staticmethod]
    fn linear_us_up(v0: f64, rho0: f64, c0: f64, s: f64, p0: f64) -> PyResult<Self> {
        Ok(Self {
            model: LinearUsUpHugoniot::new(v0, rho0, c0, s, p0).map_err(to_python_error)?,
        })
    }

    fn evaluate_scalar(&self, quantity: &str, value: f64) -> PyResult<f64> {
        self.evaluate_model(quantity, value)
            .map_err(to_python_error)
    }

    fn evaluate_array<'py>(
        &self,
        py: Python<'py>,
        quantity: &str,
        values: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        map_array(py, values, PARALLEL_ELEMENTWISE_THRESHOLD, |value| {
            self.evaluate_model(quantity, value)
        })
    }

    fn __repr__(&self) -> String {
        format!("HugoniotEos({:?})", self.model)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum ThermalModel {
    AsymptoticPowerLawMieGruneisenDebye(AsymptoticPowerLawMieGruneisenDebye<RtModel>),
    DorogokupetsOganov2007(DorogokupetsOganov2007<RtModel>),
    LinearThermalPressure(LinearThermalPressure<RtModel>),
    LogVolumeThermalPressure(LogVolumeThermalPressure<RtModel>),
    SecondOrderTaylorThermalPressure(SecondOrderTaylorThermalPressure<RtModel>),
    MieGruneisenDebye(MieGruneisenDebye<RtModel>),
    MieGruneisenEinstein(MieGruneisenEinstein<RtModel>),
    ThermalModifiedTait(ThermalModifiedTait),
    Sokolova2016(MultiOscillatorGruneisen<RtModel>),
    ThermalReferenceState(ThermalReferenceState<RtModel>),
}

impl ThermalModel {
    fn name(self) -> &'static str {
        match self {
            Self::AsymptoticPowerLawMieGruneisenDebye(_) => "AsymptoticPowerLawMieGruneisenDebye",
            Self::DorogokupetsOganov2007(_) => "DorogokupetsOganov2007",
            Self::LinearThermalPressure(_) => "LinearThermalPressure",
            Self::LogVolumeThermalPressure(_) => "LogVolumeThermalPressure",
            Self::SecondOrderTaylorThermalPressure(_) => "SecondOrderTaylorThermalPressure",
            Self::MieGruneisenDebye(_) => "MieGruneisenDebye",
            Self::MieGruneisenEinstein(_) => "MieGruneisenEinstein",
            Self::ThermalModifiedTait(_) => "ThermalModifiedTait",
            Self::Sokolova2016(_) => "MultiOscillatorGruneisenThermalEOS",
            Self::ThermalReferenceState(_) => "ThermalReferenceStateEOS",
        }
    }

    fn evaluate(self, quantity: &str, first: f64, second: f64) -> PyResult<f64> {
        match self {
            Self::AsymptoticPowerLawMieGruneisenDebye(model) => {
                evaluate_asymptotic_mie_quantity(&model, quantity, first, second)
            }
            Self::DorogokupetsOganov2007(model) => {
                evaluate_thermal_quantity(&model, quantity, first, second)
            }
            Self::LinearThermalPressure(model) => {
                evaluate_thermal_quantity(&model, quantity, first, second)
            }
            Self::LogVolumeThermalPressure(model) => {
                evaluate_thermal_quantity(&model, quantity, first, second)
            }
            Self::SecondOrderTaylorThermalPressure(model) => {
                evaluate_thermal_quantity(&model, quantity, first, second)
            }
            Self::MieGruneisenDebye(model) => {
                evaluate_mie_quantity(&model, quantity, first, second)
            }
            Self::MieGruneisenEinstein(model) => {
                evaluate_mie_quantity(&model, quantity, first, second)
            }
            Self::ThermalModifiedTait(model) => {
                evaluate_caloric_quantity(&model, quantity, first, second)
            }
            Self::Sokolova2016(model) => evaluate_thermal_quantity(&model, quantity, first, second),
            Self::ThermalReferenceState(model) => {
                evaluate_thermal_quantity(&model, quantity, first, second)
            }
        }
    }

    fn volume_with_dac_confinement(
        self,
        cold_pressure: f64,
        temperature: f64,
        f_dac: f64,
    ) -> PyResult<f64> {
        let result = match self {
            Self::AsymptoticPowerLawMieGruneisenDebye(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::DorogokupetsOganov2007(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::LinearThermalPressure(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::LogVolumeThermalPressure(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::SecondOrderTaylorThermalPressure(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::MieGruneisenDebye(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::MieGruneisenEinstein(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::ThermalModifiedTait(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::Sokolova2016(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
            Self::ThermalReferenceState(model) => {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            }
        };
        result.map_err(to_python_error)
    }
}

/// Private native representation of a built-in thermal EOS.
#[pyclass(
    name = "ThermalEos",
    frozen,
    module = "peritheos._rust",
    skip_from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq)]
struct PyThermalEos {
    model: ThermalModel,
}

#[pymethods]
impl PyThermalEos {
    #[staticmethod]
    #[pyo3(signature = (
        rt_eos, tr, theta0, gamma0, q, n,
        debye_temperature_law="integrated_gruneisen"
    ))]
    fn mie_gruneisen_debye(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        theta0: f64,
        gamma0: f64,
        q: f64,
        n: f64,
        debye_temperature_law: &str,
    ) -> PyResult<Self> {
        let law = match debye_temperature_law {
            "integrated_gruneisen" => DebyeTemperatureLaw::IntegratedGruneisen,
            "variable_exponent" => DebyeTemperatureLaw::VariableExponent,
            _ => {
                return Err(python_validation_error(
                    "debye_temperature_law must be 'integrated_gruneisen' or 'variable_exponent'",
                ));
            }
        };
        Ok(Self {
            model: ThermalModel::MieGruneisenDebye(
                MieGruneisenDebye::new_with_temperature_law(
                    rt_eos.model,
                    tr,
                    theta0,
                    gamma0,
                    q,
                    n,
                    law,
                )
                .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn mie_gruneisen_einstein(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        theta0: f64,
        gamma0: f64,
        q: f64,
        n: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            model: ThermalModel::MieGruneisenEinstein(
                MieGruneisenEinstein::new(rt_eos.model, tr, theta0, gamma0, q, n)
                    .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    fn asymptotic_power_law_mie_gruneisen_debye(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        theta0: f64,
        gamma0: f64,
        a: f64,
        b: f64,
        n: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            model: ThermalModel::AsymptoticPowerLawMieGruneisenDebye(
                AsymptoticPowerLawMieGruneisenDebye::new(rt_eos.model, tr, theta0, gamma0, a, b, n)
                    .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn linear_thermal_pressure(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        alpha_kt: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            model: ThermalModel::LinearThermalPressure(
                LinearThermalPressure::new(rt_eos.model, tr, alpha_kt).map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn log_volume_thermal_pressure(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        alpha_kt_ref: f64,
        dk_dt_v: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            model: ThermalModel::LogVolumeThermalPressure(
                LogVolumeThermalPressure::new(rt_eos.model, tr, alpha_kt_ref, dk_dt_v)
                    .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    fn second_order_taylor_thermal_pressure(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        eta0: f64,
        c0: f64,
        c1: f64,
        c2: f64,
        c3: f64,
        c4: f64,
        c5: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            model: ThermalModel::SecondOrderTaylorThermalPressure(
                SecondOrderTaylorThermalPressure::new(
                    rt_eos.model,
                    tr,
                    eta0,
                    c0,
                    c1,
                    c2,
                    c3,
                    c4,
                    c5,
                )
                .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    #[pyo3(signature = (
        rt_eos, tr, alpha0, dk_dt, alpha1=0.0,
        thermal_expansion_law="constant",
        reference_volume_law="integrated_expansivity"
    ))]
    #[allow(clippy::too_many_arguments)]
    fn thermal_reference_state(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        alpha0: f64,
        dk_dt: f64,
        alpha1: f64,
        thermal_expansion_law: &str,
        reference_volume_law: &str,
    ) -> PyResult<Self> {
        let expansion_law = match thermal_expansion_law {
            "constant" => ThermalExpansionLaw::Constant,
            "linear_temperature" => ThermalExpansionLaw::LinearTemperature,
            _ => {
                return Err(python_validation_error(
                    "thermal_expansion_law must be 'constant' or 'linear_temperature'",
                ));
            }
        };
        let volume_law = match reference_volume_law {
            "integrated_expansivity" => ReferenceVolumeLaw::IntegratedExpansivity,
            "linear_temperature" => ReferenceVolumeLaw::LinearTemperature,
            "berman" => ReferenceVolumeLaw::Berman,
            _ => {
                return Err(python_validation_error(
                    "reference_volume_law must be 'integrated_expansivity', 'linear_temperature', or 'berman'",
                ));
            }
        };
        Ok(Self {
            model: ThermalModel::ThermalReferenceState(
                ThermalReferenceState::new(
                    rt_eos.model,
                    tr,
                    alpha0,
                    dk_dt,
                    alpha1,
                    expansion_law,
                    volume_law,
                )
                .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    fn thermal_modified_tait(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        theta: f64,
        alpha0: f64,
        n: f64,
    ) -> PyResult<Self> {
        let RtModel::ModifiedTait(reference) = rt_eos.model else {
            return Err(python_configuration_error(
                "ThermalModifiedTait requires a ModifiedTait EOS",
            ));
        };
        Ok(Self {
            model: ThermalModel::ThermalModifiedTait(
                ThermalModifiedTait::new(reference, tr, theta, alpha0, n)
                    .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    fn dorogokupets_oganov_2007(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        theta_b1: f64,
        d_b1: f64,
        m_b1: f64,
        theta_b2: f64,
        d_b2: f64,
        m_b2: f64,
        theta_e1: f64,
        m_e1: f64,
        theta_e2: f64,
        m_e2: f64,
        gamma0: f64,
        gamma_inf: f64,
        beta: f64,
        anharmonic_a: f64,
        anharmonic_m: f64,
        electronic_e: f64,
        electronic_g: f64,
        defect_h: f64,
        defect_s: f64,
        n: f64,
    ) -> PyResult<Self> {
        let parameters = DorogokupetsOganov2007Parameters {
            tr,
            theta_b1,
            d_b1,
            m_b1,
            theta_b2,
            d_b2,
            m_b2,
            theta_e1,
            m_e1,
            theta_e2,
            m_e2,
            gamma0,
            gamma_inf,
            beta,
            anharmonic_a,
            anharmonic_m,
            electronic_e,
            electronic_g,
            defect_h,
            defect_s,
        };
        Ok(Self {
            model: ThermalModel::DorogokupetsOganov2007(
                DorogokupetsOganov2007::new(rt_eos.model, parameters, n)
                    .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::many_single_char_names)]
    #[pyo3(signature = (
        rt_eos, tr, qe1o, me1, qe2o, me2, delta, t, a_0, m, g, e_0,
        beta=0.0, qbo=1.0, d=1.0, mb=0.0, qb1o=1.0, d1=1.0, mb1=0.0,
        n=None
    ))]
    fn multi_oscillator_gruneisen(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        qe1o: f64,
        me1: f64,
        qe2o: f64,
        me2: f64,
        delta: f64,
        t: f64,
        a_0: f64,
        m: f64,
        g: f64,
        e_0: f64,
        beta: f64,
        qbo: f64,
        d: f64,
        mb: f64,
        qb1o: f64,
        d1: f64,
        mb1: f64,
        n: Option<f64>,
    ) -> PyResult<Self> {
        let parameters = SokolovaParameters {
            tr,
            qe1o,
            me1,
            qe2o,
            me2,
            delta,
            t,
            a_0,
            m,
            g,
            e_0,
            beta,
            qbo,
            d,
            mb,
            qb1o,
            d1,
            mb1,
        };
        let atom_count = n.or_else(|| match rt_eos.model {
            RtModel::Holzapfel(reference) => Some(reference.n),
            _ => None,
        });
        Ok(Self {
            model: ThermalModel::Sokolova2016(
                MultiOscillatorGruneisen::new_with_atom_count(
                    rt_eos.model,
                    parameters,
                    atom_count.ok_or_else(|| {
                        python_validation_error(
                            "n is required for the generic multi-oscillator model",
                        )
                    })?,
                )
                .map_err(to_python_error)?,
            ),
        })
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::many_single_char_names)]
    #[pyo3(signature = (
        rt_eos, tr, qe1o, me1, qe2o, me2, delta, t, a_0, m, g, e_0,
        beta=0.0, qbo=1.0, d=1.0, mb=0.0, qb1o=1.0, d1=1.0, mb1=0.0,
        n=None
    ))]
    fn sokolova2016(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        qe1o: f64,
        me1: f64,
        qe2o: f64,
        me2: f64,
        delta: f64,
        t: f64,
        a_0: f64,
        m: f64,
        g: f64,
        e_0: f64,
        beta: f64,
        qbo: f64,
        d: f64,
        mb: f64,
        qb1o: f64,
        d1: f64,
        mb1: f64,
        n: Option<f64>,
    ) -> PyResult<Self> {
        Self::multi_oscillator_gruneisen(
            rt_eos, tr, qe1o, me1, qe2o, me2, delta, t, a_0, m, g, e_0, beta, qbo, d, mb, qb1o, d1,
            mb1, n,
        )
    }

    #[getter]
    fn model_name(&self) -> &'static str {
        self.model.name()
    }

    fn evaluate_scalar(&self, quantity: &str, first: f64, second: f64) -> PyResult<f64> {
        self.model.evaluate(quantity, first, second)
    }

    fn evaluate_array<'py>(
        &self,
        py: Python<'py>,
        quantity: &str,
        first: PyReadonlyArrayDyn<'py, f64>,
        second: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let model = self.model;
        let quantity = quantity.to_owned();
        map_array2(
            py,
            first,
            second,
            PARALLEL_THERMAL_THRESHOLD,
            move |left, right| model.evaluate(&quantity, left, right),
        )
    }

    fn volume_with_dac_confinement_scalar(
        &self,
        cold_pressure: f64,
        temperature: f64,
        f_dac: f64,
    ) -> PyResult<f64> {
        self.model
            .volume_with_dac_confinement(cold_pressure, temperature, f_dac)
    }

    fn volume_with_dac_confinement_array<'py>(
        &self,
        py: Python<'py>,
        cold_pressures: PyReadonlyArrayDyn<'py, f64>,
        temperatures: PyReadonlyArrayDyn<'py, f64>,
        f_dac: f64,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let model = self.model;
        map_array2(
            py,
            cold_pressures,
            temperatures,
            PARALLEL_THERMAL_THRESHOLD,
            move |cold_pressure, temperature| {
                model.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            },
        )
    }

    fn __repr__(&self) -> String {
        format!("ThermalEos(model_name='{}')", self.model.name())
    }
}

fn evaluate_thermal_quantity<T: ThermalEos>(
    model: &T,
    quantity: &str,
    first: f64,
    second: f64,
) -> PyResult<f64> {
    let result = match quantity {
        "thermal_pressure" => model.thermal_pressure(first, second),
        "thermal_pressure_increment" => model.thermal_pressure_increment(first, second),
        "pressure" => model.pressure(first, second),
        "volume" => model.volume(first, second),
        "temperature" => model.temperature(first, second),
        "bulk_modulus" => model.bulk_modulus(first, second, 1.0e-6),
        "isothermal_compressibility" => model.isothermal_compressibility(first, second),
        "thermal_expansivity" => model.thermal_expansivity(first, second, 1.0e-5),
        _ => {
            return Err(python_unsupported_error(format!(
                "quantity '{quantity}' is unavailable for this thermal model"
            )));
        }
    };
    result.map_err(to_python_error)
}

fn evaluate_caloric_quantity<T: CaloricEos>(
    model: &T,
    quantity: &str,
    first: f64,
    second: f64,
) -> PyResult<f64> {
    match quantity {
        "molar_heat_capacity_v" | "heat_capacity_v" => model
            .molar_heat_capacity_v(first, second)
            .map_err(to_python_error),
        "molar_heat_capacity_p" | "heat_capacity_p" => model
            .molar_heat_capacity_p(first, second)
            .map_err(to_python_error),
        "gruneisen_parameter" => model
            .gruneisen_parameter(first, second)
            .map_err(to_python_error),
        "adiabatic_bulk_modulus" => model
            .adiabatic_bulk_modulus(first, second)
            .map_err(to_python_error),
        _ => evaluate_thermal_quantity(model, quantity, first, second),
    }
}

fn evaluate_mie_quantity<R, const DEBYE: bool>(
    model: &peritheos::thermal::MieGruneisen<R, DEBYE>,
    quantity: &str,
    first: f64,
    second: f64,
) -> PyResult<f64>
where
    R: IsothermalEos,
{
    match quantity {
        "characteristic_temperature" => model
            .characteristic_temperature(first)
            .map_err(to_python_error),
        "thermal_energy" | "thermal_internal_energy" => {
            model.thermal_energy(first, second).map_err(to_python_error)
        }
        "thermal_entropy" => model
            .thermal_entropy(first, second)
            .map_err(to_python_error),
        "vibrational_pressure" => model
            .vibrational_pressure(first, second)
            .map_err(to_python_error),
        "thermal_helmholtz_free_energy" => model
            .thermal_helmholtz_free_energy(first, second)
            .map_err(to_python_error),
        "thermal_enthalpy" => model
            .thermal_enthalpy(first, second)
            .map_err(to_python_error),
        "thermal_gibbs_free_energy" => model
            .thermal_gibbs_free_energy(first, second)
            .map_err(to_python_error),
        _ => evaluate_caloric_quantity(model, quantity, first, second),
    }
}

fn evaluate_asymptotic_mie_quantity<R: IsothermalEos>(
    model: &AsymptoticPowerLawMieGruneisenDebye<R>,
    quantity: &str,
    first: f64,
    second: f64,
) -> PyResult<f64> {
    match quantity {
        "characteristic_temperature" => model
            .characteristic_temperature(first)
            .map_err(to_python_error),
        "thermal_energy" | "thermal_internal_energy" => {
            model.thermal_energy(first, second).map_err(to_python_error)
        }
        "thermal_entropy" => model
            .thermal_entropy(first, second)
            .map_err(to_python_error),
        "vibrational_pressure" => model
            .vibrational_pressure(first, second)
            .map_err(to_python_error),
        "thermal_helmholtz_free_energy" => model
            .thermal_helmholtz_free_energy(first, second)
            .map_err(to_python_error),
        "thermal_enthalpy" => model
            .thermal_enthalpy(first, second)
            .map_err(to_python_error),
        "thermal_gibbs_free_energy" => model
            .thermal_gibbs_free_energy(first, second)
            .map_err(to_python_error),
        _ => evaluate_caloric_quantity(model, quantity, first, second),
    }
}

#[allow(clippy::needless_pass_by_value)]
fn map_array2<'py, F>(
    py: Python<'py>,
    first: PyReadonlyArrayDyn<'py, f64>,
    second: PyReadonlyArrayDyn<'py, f64>,
    parallel_threshold: usize,
    evaluator: F,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>>
where
    F: Fn(f64, f64) -> PyResult<f64> + Send + Sync,
{
    let first_view = first.as_array();
    let second_view = second.as_array();
    if first_view.shape() != second_view.shape() {
        return Err(python_validation_error(
            "native array inputs must have matching broadcasted shapes",
        ));
    }
    let shape = first_view.raw_dim();
    let values = if first_view.len() >= parallel_threshold {
        let inputs = first_view
            .iter()
            .copied()
            .zip(second_view.iter().copied())
            .collect::<Vec<_>>();
        py.detach(move || {
            inputs
                .into_par_iter()
                .map(|(left, right)| evaluator(left, right))
                .collect::<Vec<_>>()
        })
        .into_iter()
        .collect::<PyResult<Vec<_>>>()?
    } else {
        first_view
            .iter()
            .zip(second_view.iter())
            .map(|(&left, &right)| evaluator(left, right))
            .collect::<PyResult<Vec<_>>>()?
    };
    let output = ArrayD::from_shape_vec(shape, values)
        .map_err(|error| python_validation_error(error.to_string()))?;
    Ok(output.into_pyarray(py))
}

#[allow(clippy::needless_pass_by_value)]
fn map_array<'py, F>(
    py: Python<'py>,
    input: PyReadonlyArrayDyn<'py, f64>,
    parallel_threshold: usize,
    evaluator: F,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>>
where
    F: Fn(f64) -> EosResult<f64> + Send + Sync,
{
    let view = input.as_array();
    let shape = view.raw_dim();
    let values = if view.len() >= parallel_threshold {
        let input_values = view.iter().copied().collect::<Vec<_>>();
        py.detach(move || {
            input_values
                .into_par_iter()
                .map(evaluator)
                .collect::<Vec<_>>()
        })
        .into_iter()
        .collect::<EosResult<Vec<_>>>()
        .map_err(to_python_error)?
    } else {
        view.iter()
            .map(|&value| evaluator(value).map_err(to_python_error))
            .collect::<PyResult<Vec<_>>>()?
    };
    let output = ArrayD::from_shape_vec(shape, values)
        .map_err(|error| python_validation_error(error.to_string()))?;
    Ok(output.into_pyarray(py))
}

#[allow(clippy::needless_pass_by_value)]
fn to_python_error(error: EosError) -> PyErr {
    let class_name = if error.is_validation() {
        "EosValidationError"
    } else {
        "EosNumericalError"
    };
    let message = if matches!(error, EosError::OutsideInvertibleRange) {
        "pressure is outside the invertible expansion range".to_owned()
    } else {
        error.to_string()
    };
    peritheos_python_error(
        class_name,
        message,
        error.code(),
        Some("eos"),
        error.field(),
    )
}

fn peritheos_python_error(
    class_name: &str,
    message: String,
    code: &str,
    operation: Option<&str>,
    field: Option<&str>,
) -> PyErr {
    Python::attach(|py| {
        let converted = (|| -> PyResult<PyErr> {
            let exception_class = py.import("peritheos.errors")?.getattr(class_name)?;
            let keyword_arguments = PyDict::new(py);
            keyword_arguments.set_item("code", code)?;
            if let Some(operation) = operation {
                keyword_arguments.set_item("operation", operation)?;
            }
            if let Some(field) = field {
                keyword_arguments.set_item("field", field)?;
            }
            let instance = exception_class.call((message.as_str(),), Some(&keyword_arguments))?;
            Ok(PyErr::from_value(instance))
        })();
        converted.unwrap_or_else(|conversion_error| {
            PyRuntimeError::new_err(format!(
                "{message} (could not construct peritheos.errors.{class_name}: {conversion_error})"
            ))
        })
    })
}

fn python_validation_error(message: impl ToString) -> PyErr {
    peritheos_python_error(
        "ValidationError",
        message.to_string(),
        "validation.invalid_input",
        None,
        None,
    )
}

fn python_configuration_error(message: impl ToString) -> PyErr {
    peritheos_python_error(
        "ConfigurationError",
        message.to_string(),
        "configuration.invalid",
        None,
        None,
    )
}

fn python_fit_validation_error(message: impl ToString) -> PyErr {
    peritheos_python_error(
        "FitValidationError",
        message.to_string(),
        "fit.invalid_input",
        Some("fit"),
        None,
    )
}

fn python_unsupported_error(message: impl ToString) -> PyErr {
    peritheos_python_error(
        "UnsupportedOperationError",
        message.to_string(),
        "operation.unsupported",
        None,
        None,
    )
}

/// Private SciPy-compatible view of a native least-squares result.
#[pyclass(
    name = "LeastSquaresResult",
    frozen,
    module = "peritheos._rust",
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
struct PyLeastSquaresResult {
    result: peritheos::fit::SolverResult,
    global_parameter_count: usize,
    structured_layout: Option<StructuredLayout>,
    predicted_pressure: Option<Vec<f64>>,
    parameter_covariance: Option<Vec<f64>>,
}

#[pymethods]
impl PyLeastSquaresResult {
    #[getter]
    fn x<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        vector_array(py, self.result.parameters.clone())
    }

    #[getter]
    fn fun<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        vector_array(py, self.result.residuals.clone())
    }

    /// Pressure predictions are present for end-to-end native EOS fits.
    #[getter]
    fn predicted_pressure<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArrayDyn<f64>>> {
        self.predicted_pressure
            .clone()
            .map(|values| vector_array(py, values))
    }

    #[getter]
    fn jac<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let column_count = self.result.parameters.len();
        let output = ArrayD::from_shape_vec(
            numpy::ndarray::IxDyn(&[self.result.residual_count, column_count]),
            self.result.jacobian.clone(),
        )
        .map_err(|error| python_validation_error(error.to_string()))?;
        Ok(output.into_pyarray(py))
    }

    /// Leading model-parameter covariance after profiling latent coordinates.
    #[getter]
    fn parameter_covariance<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let covariance = if let Some(covariance) = &self.parameter_covariance {
            covariance.clone()
        } else {
            let column_count = self.result.parameters.len();
            if let Some(layout) = self.structured_layout {
                parameter_covariance_structured(&self.result.jacobian, layout)
            } else {
                peritheos::fit::parameter_covariance(
                    &self.result.jacobian,
                    self.result.residual_count,
                    column_count,
                    self.global_parameter_count,
                )
            }
            .map_err(to_python_fit_error)?
        };
        let output = ArrayD::from_shape_vec(
            numpy::ndarray::IxDyn(&[self.global_parameter_count, self.global_parameter_count]),
            covariance,
        )
        .map_err(|error| python_validation_error(error.to_string()))?;
        Ok(output.into_pyarray(py))
    }

    #[getter]
    fn cost(&self) -> f64 {
        self.result.cost
    }

    #[getter]
    fn optimality(&self) -> f64 {
        self.result.optimality
    }

    #[getter]
    fn success(&self) -> bool {
        self.result.success
    }

    #[getter]
    fn message(&self) -> &str {
        &self.result.message
    }

    #[getter]
    fn status(&self) -> i32 {
        self.result.status
    }

    #[getter]
    fn nfev(&self) -> usize {
        self.result.function_evaluations
    }

    #[getter]
    fn njev(&self) -> usize {
        self.result.jacobian_evaluations
    }
}

fn vector_array(py: Python<'_>, values: Vec<f64>) -> Bound<'_, PyArrayDyn<f64>> {
    ArrayD::from_shape_vec(numpy::ndarray::IxDyn(&[values.len()]), values)
        .expect("a vector always has a valid one-dimensional shape")
        .into_pyarray(py)
}

#[pyfunction]
#[pyo3(signature = (
    residual_function, initial, lower, upper, loss="linear", f_scale=1.0,
    max_nfev=None, global_parameter_count=None, point_count=None,
    latent_coordinate_count=None
))]
#[allow(clippy::too_many_arguments)]
fn fit_least_squares<'py>(
    py: Python<'py>,
    residual_function: Bound<'py, PyAny>,
    initial: PyReadonlyArrayDyn<'py, f64>,
    lower: PyReadonlyArrayDyn<'py, f64>,
    upper: PyReadonlyArrayDyn<'py, f64>,
    loss: &str,
    f_scale: f64,
    max_nfev: Option<usize>,
    global_parameter_count: Option<usize>,
    point_count: Option<usize>,
    latent_coordinate_count: Option<usize>,
) -> PyResult<PyLeastSquaresResult> {
    let initial = initial.as_array().iter().copied().collect::<Vec<_>>();
    let lower = lower.as_array().iter().copied().collect::<Vec<_>>();
    let upper = upper.as_array().iter().copied().collect::<Vec<_>>();
    let options = SolverOptions {
        loss: Loss::from_name(loss).map_err(to_python_fit_error)?,
        f_scale,
        max_evaluations: max_nfev,
        ..SolverOptions::default()
    };
    let evaluator = |parameters: &[f64]| {
        let argument = vector_array(py, parameters.to_vec());
        let result = residual_function
            .call1((argument,))
            .map_err(|error| FitError::Evaluation(error.to_string()))?;
        let array = result
            .extract::<PyReadonlyArrayDyn<'_, f64>>()
            .map_err(|error| FitError::Evaluation(error.to_string()))?;
        Ok(array.as_array().iter().copied().collect())
    };
    let layout = match (global_parameter_count, point_count, latent_coordinate_count) {
        (None, None, None) => None,
        (Some(global_parameter_count), Some(point_count), Some(latent_coordinate_count)) => {
            Some(StructuredLayout {
                global_parameter_count,
                point_count,
                latent_coordinate_count,
            })
        }
        _ => {
            return Err(python_fit_validation_error(
                "all structured least-squares dimensions must be supplied together",
            ));
        }
    };
    let global_parameter_count = layout.map_or(initial.len(), |value| value.global_parameter_count);
    let result = if let Some(layout) = layout {
        least_squares_structured(&initial, &lower, &upper, options, layout, evaluator)
    } else {
        least_squares(&initial, &lower, &upper, options, evaluator)
    }
    .map_err(to_python_fit_error)?;
    Ok(PyLeastSquaresResult {
        result,
        global_parameter_count,
        structured_layout: layout,
        predicted_pressure: None,
        parameter_covariance: None,
    })
}

fn to_python_fit_error(error: FitError) -> PyErr {
    match error {
        FitError::EosEvaluation(error) => {
            let class_name = if error.is_validation() {
                "FitEosValidationError"
            } else {
                "FitEosNumericalError"
            };
            peritheos_python_error(
                class_name,
                error.to_string(),
                error.code(),
                Some("fit"),
                error.field(),
            )
        }
        FitError::InvalidInput(message) => peritheos_python_error(
            "FitValidationError",
            message,
            "fit.invalid_input",
            Some("fit"),
            None,
        ),
        FitError::Evaluation(message) => peritheos_python_error(
            "FitError",
            message,
            "fit.evaluation_failed",
            Some("fit"),
            None,
        ),
        FitError::SingularSystem => peritheos_python_error(
            "FitNumericalError",
            "linear system is singular".to_owned(),
            "fit.singular_system",
            Some("fit"),
            None,
        ),
        other => peritheos_python_error(
            "FitError",
            other.to_string(),
            other.code(),
            Some("fit"),
            None,
        ),
    }
}

/// Private result of native delta-method covariance propagation.
#[pyclass(
    name = "LinearPropagation",
    frozen,
    module = "peritheos._rust",
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
struct PyLinearPropagation {
    result: LinearPropagation,
    output_count: usize,
}

#[pymethods]
impl PyLinearPropagation {
    #[getter]
    fn variance<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        vector_array(py, self.result.variance.clone())
    }

    #[getter]
    fn covariance<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArrayDyn<f64>>>> {
        let Some(values) = &self.result.covariance else {
            return Ok(None);
        };
        let output = ArrayD::from_shape_vec(
            numpy::ndarray::IxDyn(&[self.output_count, self.output_count]),
            values.clone(),
        )
        .map_err(|error| python_validation_error(error.to_string()))?;
        Ok(Some(output.into_pyarray(py)))
    }
}

#[pyfunction]
#[pyo3(signature = (jacobian, covariance, state_variance, full_covariance=false))]
fn linear_uncertainty(
    jacobian: PyReadonlyArrayDyn<'_, f64>,
    covariance: PyReadonlyArrayDyn<'_, f64>,
    state_variance: PyReadonlyArrayDyn<'_, f64>,
    full_covariance: bool,
) -> PyResult<PyLinearPropagation> {
    let jacobian_view = jacobian.as_array();
    if jacobian_view.ndim() != 2 {
        return Err(python_validation_error("jacobian must be two-dimensional"));
    }
    let output_count = jacobian_view.shape()[0];
    let parameter_count = jacobian_view.shape()[1];
    let jacobian_values = jacobian_view.iter().copied().collect::<Vec<_>>();
    let covariance_values = covariance.as_array().iter().copied().collect::<Vec<_>>();
    let state_values = state_variance
        .as_array()
        .iter()
        .copied()
        .collect::<Vec<_>>();
    let result = propagate_linear_uncertainty(
        &jacobian_values,
        output_count,
        parameter_count,
        &covariance_values,
        &state_values,
        full_covariance,
    )
    .map_err(to_python_fit_error)?;
    Ok(PyLinearPropagation {
        result,
        output_count,
    })
}

/// Private result of native Monte Carlo summary statistics.
#[pyclass(
    name = "MonteCarloSummary",
    frozen,
    module = "peritheos._rust",
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
struct PyMonteCarloSummary {
    result: MonteCarloSummary,
    output_count: usize,
}

#[pymethods]
impl PyMonteCarloSummary {
    #[getter]
    fn standard_error<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        vector_array(py, self.result.standard_error.clone())
    }

    #[getter]
    fn lower<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        vector_array(py, self.result.lower.clone())
    }

    #[getter]
    fn upper<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        vector_array(py, self.result.upper.clone())
    }

    #[getter]
    fn covariance<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArrayDyn<f64>>>> {
        let Some(values) = &self.result.covariance else {
            return Ok(None);
        };
        let output = ArrayD::from_shape_vec(
            numpy::ndarray::IxDyn(&[self.output_count, self.output_count]),
            values.clone(),
        )
        .map_err(|error| python_validation_error(error.to_string()))?;
        Ok(Some(output.into_pyarray(py)))
    }
}

#[pyfunction]
#[pyo3(signature = (samples, confidence, full_covariance=false))]
fn monte_carlo_summary(
    samples: PyReadonlyArrayDyn<'_, f64>,
    confidence: f64,
    full_covariance: bool,
) -> PyResult<PyMonteCarloSummary> {
    let view = samples.as_array();
    if view.ndim() != 2 {
        return Err(python_validation_error("samples must be two-dimensional"));
    }
    let sample_count = view.shape()[0];
    let output_count = view.shape()[1];
    let values = view.iter().copied().collect::<Vec<_>>();
    let result = summarize_monte_carlo(
        &values,
        sample_count,
        output_count,
        confidence,
        full_covariance,
    )
    .map_err(to_python_fit_error)?;
    Ok(PyMonteCarloSummary {
        result,
        output_count,
    })
}

#[pyfunction]
fn holzapfel_derivative_analytical(
    v0: f64,
    volume: f64,
    bulk_modulus: f64,
    k0: f64,
    c0: f64,
    c2: f64,
) -> PyResult<f64> {
    holzapfel_bulk_modulus_derivative_analytical(v0, volume, bulk_modulus, k0, c0, c2)
        .map_err(to_python_error)
}

#[pymodule]
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyRtEos>()?;
    module.add_class::<PyHugoniotEos>()?;
    module.add_class::<PyThermalEos>()?;
    module.add_class::<PyLeastSquaresResult>()?;
    module.add_class::<PyLinearPropagation>()?;
    module.add_class::<PyMonteCarloSummary>()?;
    module.add_function(wrap_pyfunction!(fit_least_squares, module)?)?;
    module.add_function(wrap_pyfunction!(native_fit::fit_rt_eos_native, module)?)?;
    module.add_function(wrap_pyfunction!(
        native_fit::fit_thermal_eos_native,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(linear_uncertainty, module)?)?;
    module.add_function(wrap_pyfunction!(monte_carlo_summary, module)?)?;
    module.add_function(wrap_pyfunction!(holzapfel_derivative_analytical, module)?)?;
    Ok(())
}
