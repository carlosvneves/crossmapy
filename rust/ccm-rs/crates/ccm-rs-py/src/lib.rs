use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;

#[pyfunction]
fn ssr_pred_boot(_a: Vec<f64>, _b: Option<Vec<f64>>, _e: usize, _tau: usize, _predstep: usize) -> PyResult<PyObject> {
    Err(PyNotImplementedError::new_err("ccm-rs kernel not wired yet"))
}

#[pyfunction]
fn ccm_boot(_a: Vec<f64>, _b: Vec<f64>, _e: usize, _tau: usize, _iterations: usize) -> PyResult<PyObject> {
    Err(PyNotImplementedError::new_err("ccm-rs kernel not wired yet"))
}

#[pymodule]
fn _ccm_rs(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ssr_pred_boot, m)?)?;
    m.add_function(wrap_pyfunction!(ccm_boot, m)?)?;
    Ok(())
}
