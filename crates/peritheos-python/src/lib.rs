#![allow(clippy::needless_pass_by_value, clippy::similar_names)]

mod native_fit;

use numpy::ndarray::ArrayD;
use numpy::{IntoPyArray, PyArrayDyn, PyReadonlyArrayDyn};
use peritheos_core::isothermal::{
    Holzapfel, ModifiedTait, Murnaghan, NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2,
    BM3, BM4,
};
use peritheos_core::thermal::{
    MieGruneisenDebye, MieGruneisenEinstein, Sokolova2016, SokolovaParameters, ThermalModifiedTait,
};
use peritheos_core::{CaloricEos, EosError, EosResult, IsothermalEos, ThermalEos};
use peritheos_fit::{
    least_squares, least_squares_structured, propagate_linear_uncertainty, summarize_monte_carlo,
    FitError, LinearPropagation, Loss, MonteCarloSummary, SolverOptions, StructuredLayout,
};
use pyo3::exceptions::{
    PyArithmeticError, PyNotImplementedError, PyRuntimeError, PyTypeError, PyValueError,
};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyModule};

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

    fn volume_scalar(&self, pressure: f64) -> PyResult<f64> {
        self.model.volume(pressure).map_err(to_python_error)
    }

    fn pressure_array<'py>(
        &self,
        py: Python<'py>,
        volumes: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        map_array(py, volumes, |volume| self.model.pressure(volume))
    }

    fn bulk_modulus_array<'py>(
        &self,
        py: Python<'py>,
        volumes: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        map_array(py, volumes, |volume| self.model.bulk_modulus(volume))
    }

    fn volume_array<'py>(
        &self,
        py: Python<'py>,
        pressures: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        map_array(py, pressures, |pressure| self.model.volume(pressure))
    }

    fn __repr__(&self) -> String {
        format!("RtEos(model_name='{}')", self.model.name())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum ThermalModel {
    MieGruneisenDebye(MieGruneisenDebye<RtModel>),
    MieGruneisenEinstein(MieGruneisenEinstein<RtModel>),
    ThermalModifiedTait(ThermalModifiedTait),
    Sokolova2016(Sokolova2016),
}

impl ThermalModel {
    fn name(self) -> &'static str {
        match self {
            Self::MieGruneisenDebye(_) => "MieGruneisenDebye",
            Self::MieGruneisenEinstein(_) => "MieGruneisenEinstein",
            Self::ThermalModifiedTait(_) => "ThermalModifiedTait",
            Self::Sokolova2016(_) => "Sokolova2016",
        }
    }

    fn evaluate(self, quantity: &str, first: f64, second: f64) -> PyResult<f64> {
        match self {
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
        }
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
    fn mie_gruneisen_debye(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        theta0: f64,
        gamma0: f64,
        q: f64,
        n: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            model: ThermalModel::MieGruneisenDebye(
                MieGruneisenDebye::new(rt_eos.model, tr, theta0, gamma0, q, n)
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
    fn thermal_modified_tait(
        rt_eos: PyRef<'_, PyRtEos>,
        tr: f64,
        theta: f64,
        alpha0: f64,
        n: f64,
    ) -> PyResult<Self> {
        let RtModel::ModifiedTait(reference) = rt_eos.model else {
            return Err(PyTypeError::new_err(
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
    #[pyo3(signature = (
        rt_eos, tr, qe1o, me1, qe2o, me2, delta, t, a_0, m, g, e_0,
        beta=0.0, qbo=1.0, d=1.0, mb=0.0, qb1o=1.0, d1=1.0, mb1=0.0
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
    ) -> PyResult<Self> {
        let RtModel::Holzapfel(reference) = rt_eos.model else {
            return Err(PyTypeError::new_err(
                "Sokolova2016 requires a Holzapfel room-temperature EOS",
            ));
        };
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
        Ok(Self {
            model: ThermalModel::Sokolova2016(
                Sokolova2016::new(reference, parameters).map_err(to_python_error)?,
            ),
        })
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
        map_array2(py, first, second, |left, right| {
            self.model.evaluate(quantity, left, right)
        })
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
        "pressure" => model.pressure(first, second),
        "volume" => model.volume(first, second),
        "temperature" => model.temperature(first, second),
        "bulk_modulus" => model.bulk_modulus(first, second, 1.0e-6),
        "isothermal_compressibility" => model.isothermal_compressibility(first, second),
        "thermal_expansivity" => model.thermal_expansivity(first, second, 1.0e-5),
        _ => {
            return Err(PyNotImplementedError::new_err(format!(
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
    model: &peritheos_core::thermal::MieGruneisen<R, DEBYE>,
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

#[allow(clippy::needless_pass_by_value)]
fn map_array2<'py, F>(
    py: Python<'py>,
    first: PyReadonlyArrayDyn<'py, f64>,
    second: PyReadonlyArrayDyn<'py, f64>,
    mut evaluator: F,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>>
where
    F: FnMut(f64, f64) -> PyResult<f64>,
{
    let first_view = first.as_array();
    let second_view = second.as_array();
    if first_view.shape() != second_view.shape() {
        return Err(PyValueError::new_err(
            "native array inputs must have matching broadcasted shapes",
        ));
    }
    let values = first_view
        .iter()
        .zip(second_view.iter())
        .map(|(&left, &right)| evaluator(left, right))
        .collect::<PyResult<Vec<_>>>()?;
    let output = ArrayD::from_shape_vec(first_view.raw_dim(), values)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(output.into_pyarray(py))
}

#[allow(clippy::needless_pass_by_value)]
fn map_array<'py, F>(
    py: Python<'py>,
    input: PyReadonlyArrayDyn<'py, f64>,
    mut evaluator: F,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    let view = input.as_array();
    let values = view
        .iter()
        .map(|&value| evaluator(value).map_err(to_python_error))
        .collect::<PyResult<Vec<_>>>()?;
    let output = ArrayD::from_shape_vec(view.raw_dim(), values)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(output.into_pyarray(py))
}

#[allow(clippy::needless_pass_by_value)]
fn to_python_error(error: EosError) -> PyErr {
    match error {
        EosError::OutsideInvertibleRange => {
            PyValueError::new_err("pressure is outside the invertible expansion range")
        }
        EosError::InvalidParameter { .. }
        | EosError::InvalidState { .. }
        | EosError::BracketingFailed => PyValueError::new_err(error.to_string()),
        EosError::ConvergenceFailed | EosError::NonFiniteResult => {
            PyArithmeticError::new_err(error.to_string())
        }
        _ => PyRuntimeError::new_err(error.to_string()),
    }
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
    result: peritheos_fit::SolverResult,
    parameter_count: usize,
    predicted_pressure: Option<Vec<f64>>,
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
        let output = ArrayD::from_shape_vec(
            numpy::ndarray::IxDyn(&[self.result.residual_count, self.parameter_count]),
            self.result.jacobian.clone(),
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
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
            return Err(PyValueError::new_err(
                "all structured least-squares dimensions must be supplied together",
            ));
        }
    };
    let result = if let Some(layout) = layout {
        least_squares_structured(&initial, &lower, &upper, options, layout, evaluator)
    } else {
        least_squares(&initial, &lower, &upper, options, evaluator)
    }
    .map_err(to_python_fit_error)?;
    Ok(PyLeastSquaresResult {
        result,
        parameter_count: initial.len(),
        predicted_pressure: None,
    })
}

fn to_python_fit_error(error: FitError) -> PyErr {
    match error {
        FitError::InvalidInput(_) => PyValueError::new_err(error.to_string()),
        FitError::Evaluation(_) | FitError::SingularSystem => {
            PyRuntimeError::new_err(error.to_string())
        }
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
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
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
        return Err(PyValueError::new_err("jacobian must be two-dimensional"));
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
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
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
        return Err(PyValueError::new_err("samples must be two-dimensional"));
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

#[pymodule]
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyRtEos>()?;
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
    Ok(())
}
