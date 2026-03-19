use thiserror::Error;

#[derive(Debug, Error)]
pub enum CcmError {
    #[error("invalid input: {0}")]
    InvalidInput(&'static str),
}
