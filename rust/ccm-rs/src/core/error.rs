use thiserror::Error;

#[derive(Debug, Error)]
pub enum CcmError {
    #[error("not implemented: {0}")]
    NotImplemented(&'static str),
}
