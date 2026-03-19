"""Public Python API for crossmapy."""

from .api import CCM_boot, SSR_check_signal, SSR_pred_boot, ccmtest, load_ccm_data, make_ccm_data

__all__ = [
    "make_ccm_data",
    "load_ccm_data",
    "SSR_pred_boot",
    "SSR_check_signal",
    "CCM_boot",
    "ccmtest",
]

__version__ = "0.1.0"
