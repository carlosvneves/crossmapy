use ccm_rs::{ccm, simplex};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyfunction]
#[pyo3(signature = (A, B=None, E=2, tau=1, predstep=1))]
fn ssr_pred_boot(
    py: Python<'_>,
    A: Vec<f64>,
    B: Option<Vec<f64>>,
    E: usize,
    tau: usize,
    predstep: usize,
) -> PyResult<PyObject> {
    let result = simplex::ssr_pred_boot(&A, B.as_deref(), E, tau, predstep)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let out = PyDict::new_bound(py);
    out.set_item("A", A.clone())?;
    out.set_item("Aest", result.aest)?;
    out.set_item("B", B.unwrap_or_else(|| A.clone()))?;
    out.set_item("E", E)?;
    out.set_item("tau", tau)?;
    out.set_item("pBlength", A.len())?;
    out.set_item("pAlength", A.len())?;
    out.set_item("predstep", predstep)?;
    out.set_item("rho", result.rho)?;
    out.set_item("acceptablelib", result.acceptablelib)?;
    out.set_item("plengthacceptablelib", result.plengthacceptablelib)?;
    Ok(out.into())
}

#[pyfunction]
#[pyo3(signature = (A, B, E, tau=1, DesiredL=None, iterations=100))]
fn ccm_boot(
    py: Python<'_>,
    A: Vec<f64>,
    B: Vec<f64>,
    E: usize,
    tau: usize,
    DesiredL: Option<Vec<usize>>,
    iterations: usize,
) -> PyResult<PyObject> {
    let result = ccm::ccm_boot(&A, &B, E, tau, DesiredL.as_deref(), iterations)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let out = PyDict::new_bound(py);
    out.set_item("A", A)?;
    out.set_item("Aest", result.aest)?;
    out.set_item("B", B)?;
    out.set_item("rho", result.rho)?;
    out.set_item("sdevrho", result.sdevrho)?;
    out.set_item("Lobs", result.lobs)?;
    out.set_item("E", E)?;
    out.set_item("tau", tau)?;
    out.set_item("FULLinfo", result.fullinfo)?;
    Ok(out.into())
}

#[pymodule]
fn _ccm_rs(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ssr_pred_boot, m)?)?;
    m.add_function(wrap_pyfunction!(ccm_boot, m)?)?;
    Ok(())
}
