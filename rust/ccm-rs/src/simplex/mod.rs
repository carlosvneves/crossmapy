use crate::core::error::CcmError;

pub fn ssr_pred_boot(_a: &[f64], _b: Option<&[f64]>, _e: usize, _tau: usize, _predstep: usize) -> Result<f64, CcmError> {
    Err(CcmError::NotImplemented("ssr_pred_boot kernel"))
}
