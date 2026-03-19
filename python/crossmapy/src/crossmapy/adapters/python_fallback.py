"""Pure Python fallback backend."""

from ..ccm import CCM_boot as py_CCM_boot
from ..simplex import SSR_pred_boot as py_SSR_pred_boot


def fallback_ssr_pred_boot(*, A, B=None, E=2, tau=1, predstep=1):
    return py_SSR_pred_boot(A=A, B=B, E=E, tau=tau, predstep=predstep)


def fallback_ccm_boot(*, A, B, E, tau=1, iterations=100):
    return py_CCM_boot(A=A, B=B, E=E, tau=tau, iterations=iterations)
