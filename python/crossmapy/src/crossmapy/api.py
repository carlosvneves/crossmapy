"""Stable public API facade for crossmapy."""

from .adapters.rust_backend import backend_ccm_boot, backend_ssr_pred_boot
from .ccm import ccmtest
from .data import load_ccm_data, make_ccm_data
from .signal import SSR_check_signal


def SSR_pred_boot(A, B=None, E=2, tau=1, predstep=1):
    return backend_ssr_pred_boot(A=A, B=B, E=E, tau=tau, predstep=predstep)


def CCM_boot(A, B, E, tau=1, iterations=100):
    return backend_ccm_boot(A=A, B=B, E=E, tau=tau, iterations=iterations)


__all__ = [
    "make_ccm_data",
    "load_ccm_data",
    "SSR_pred_boot",
    "SSR_check_signal",
    "CCM_boot",
    "ccmtest",
]
