"""Adapter for ccm-rs native backend with safe Python fallback."""

import numpy as np

from .python_fallback import fallback_ccm_boot, fallback_ssr_pred_boot

try:
    from _ccm_rs import ccm_boot as rust_ccm_boot
    from _ccm_rs import ssr_pred_boot as rust_ssr_pred_boot
except Exception:
    rust_ccm_boot = None
    rust_ssr_pred_boot = None


def _normalize_ssr_result(result):
    result["A"] = np.asarray(result["A"], dtype=float)
    result["Aest"] = np.asarray(result["Aest"], dtype=float)
    result["B"] = np.asarray(result["B"], dtype=float)
    result["acceptablelib"] = np.asarray(result["acceptablelib"], dtype=int)
    return result


def _normalize_ccm_result(result):
    result["A"] = np.asarray(result["A"], dtype=float)
    result["Aest"] = np.asarray(result["Aest"], dtype=float)
    result["B"] = np.asarray(result["B"], dtype=float)
    result["rho"] = np.asarray(result["rho"], dtype=float)
    result["sdevrho"] = np.asarray(result["sdevrho"], dtype=float)
    result["Lobs"] = np.asarray(result["Lobs"], dtype=int)
    result["FULLinfo"] = np.asarray(result["FULLinfo"], dtype=float)
    return result


def backend_ssr_pred_boot(*, A, B=None, E=2, tau=1, predstep=1):
    if rust_ssr_pred_boot is None:
        return fallback_ssr_pred_boot(A=A, B=B, E=E, tau=tau, predstep=predstep)
    try:
        return _normalize_ssr_result(
            rust_ssr_pred_boot(A=A, B=B, E=E, tau=tau, predstep=predstep)
        )
    except Exception:
        return fallback_ssr_pred_boot(A=A, B=B, E=E, tau=tau, predstep=predstep)


def backend_ccm_boot(*, A, B, E, tau=1, iterations=100):
    if rust_ccm_boot is None:
        return fallback_ccm_boot(A=A, B=B, E=E, tau=tau, iterations=iterations)
    try:
        return _normalize_ccm_result(
            rust_ccm_boot(A=A, B=B, E=E, tau=tau, DesiredL=None, iterations=iterations)
        )
    except Exception:
        return fallback_ccm_boot(A=A, B=B, E=E, tau=tau, iterations=iterations)
