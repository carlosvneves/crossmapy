use crate::core::error::CcmError;

pub fn ccm_boot(_a: &[f64], _b: &[f64], _e: usize, _tau: usize, _iterations: usize) -> Result<(), CcmError> {
    Err(CcmError::NotImplemented("ccm_boot kernel"))
}
