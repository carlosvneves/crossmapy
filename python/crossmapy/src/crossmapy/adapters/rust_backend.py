"""Adapter for ccm-rs native backend with safe Python fallback."""

from .python_fallback import fallback_ccm_boot, fallback_ssr_pred_boot

try:
    from _ccm_rs import ccm_boot as rust_ccm_boot
    from _ccm_rs import ssr_pred_boot as rust_ssr_pred_boot
except Exception:
    rust_ccm_boot = None
    rust_ssr_pred_boot = None


def backend_ssr_pred_boot(*, A, B=None, E=2, tau=1, predstep=1):
    if rust_ssr_pred_boot is None:
        return fallback_ssr_pred_boot(A=A, B=B, E=E, tau=tau, predstep=predstep)
    return rust_ssr_pred_boot(A=A, B=B, E=E, tau=tau, predstep=predstep)


def backend_ccm_boot(*, A, B, E, tau=1, iterations=100):
    if rust_ccm_boot is None:
        return fallback_ccm_boot(A=A, B=B, E=E, tau=tau, iterations=iterations)
    return rust_ccm_boot(A=A, B=B, E=E, tau=tau, iterations=iterations)
