use numpy::ndarray::ArrayD;
use numpy::{IntoPyArray, PyArrayDyn, PyReadonlyArrayDyn};
use peritheos_core::isothermal::{
    Holzapfel, ModifiedTait, Murnaghan, NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2,
    BM3, BM4,
};
use peritheos_core::{EosError, EosResult, IsothermalEos};
use pyo3::exceptions::{PyArithmeticError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyModule;

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
        EosError::InvalidParameter { .. }
        | EosError::InvalidState { .. }
        | EosError::OutsideInvertibleRange
        | EosError::BracketingFailed => PyValueError::new_err(error.to_string()),
        EosError::ConvergenceFailed | EosError::NonFiniteResult => {
            PyArithmeticError::new_err(error.to_string())
        }
        _ => PyRuntimeError::new_err(error.to_string()),
    }
}

#[pymodule]
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyRtEos>()?;
    Ok(())
}
