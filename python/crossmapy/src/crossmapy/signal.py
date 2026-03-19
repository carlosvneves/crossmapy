import numpy as np
import math
from .simplex import SSR_pred_boot


def SSR_check_signal(A, E, tau=1, predsteplist=None, matchSugi=0):
    """
    Check for nonlinear signal in the data.

    Tests whether predictive ability declines with increasing time distance.
    A declining signal indicates the data has nonlinear structure.

    Parameters
    ----------
    A : array-like
        Time series to check
    E : int
        Embedding dimension
    tau : int
        Time lag
    predsteplist : array-like, optional
        List of prediction steps to test (default: 1:10)
    matchSugi : int
        Cross-validation method (0 or 1)

    Returns
    -------
    dict
        Dictionary with 'predatout', 'rho_pre_slope', 'rho_predmaxCI'
    """
    if predsteplist is None:
        predsteplist = np.arange(1, 11)

    A = np.asarray(A, dtype=float).copy()

    B = A.copy()

    predatout = []

    for predstep in predsteplist:
        result = SSR_pred_boot(
            A, B, E=E, tau=tau, predstep=predstep, matchSugi=matchSugi
        )
        rho = result["rho"]
        predatout.append({"predstep": predstep, "rho": rho})

    predatout = np.array([(p["predstep"], p["rho"]) for p in predatout])
    predsteps = predatout[:, 0]
    rhos = predatout[:, 1]

    valid = np.isfinite(rhos)
    if np.sum(valid) > 1:
        x = predsteps[valid]
        y = rhos[valid]

        x_mean = np.mean(x)
        y_mean = np.mean(y)

        slope = np.sum((x - x_mean) * (y - y_mean)) / np.sum((x - x_mean) ** 2)

        sse = np.sum((y - (slope * x + (y_mean - slope * x_mean))) ** 2)
        dof = len(x) - 2
        mse = sse / dof if dof > 0 else np.inf
        se = (
            np.sqrt(mse / np.sum((x - x_mean) ** 2))
            if np.sum((x - x_mean) ** 2) > 0
            else np.inf
        )
        t_stat = slope / se if se > 0 else 0
        p_value = 2 * (1 - _t_cdf(np.abs(t_stat), dof)) if dof > 0 else np.nan

        rho_pre_slope = [slope, p_value]
    else:
        rho_pre_slope = [np.nan, np.nan]

    last_result = SSR_pred_boot(
        A, B, E=E, tau=tau, predstep=int(predsteplist[-1]), matchSugi=matchSugi
    )
    Aest = last_result["Aest"]
    A_clean = last_result["A"]

    valid_idx = np.isfinite(Aest) & np.isfinite(A_clean)

    if np.sum(valid_idx) > 2:
        n = np.sum(valid_idx)
        r, _ = _pearsonr(A_clean[valid_idx], Aest[valid_idx])

        z = 0.5 * np.log((1 + r) / (1 - r))
        se_z = 1 / np.sqrt(n - 3)

        z_critical = 1.96
        ci_low = (np.exp(2 * (z - z_critical * se_z)) - 1) / (
            np.exp(2 * (z - z_critical * se_z)) + 1
        )
        ci_high = (np.exp(2 * (z + z_critical * se_z)) - 1) / (
            np.exp(2 * (z + z_critical * se_z)) + 1
        )

        rho_predmaxCI = [ci_low, ci_high]
    else:
        rho_predmaxCI = [np.nan, np.nan]

    return {
        "predatout": predatout,
        "rho_pre_slope": rho_pre_slope,
        "rho_predmaxCI": rho_predmaxCI,
    }


def _pearsonr(x, y):
    """Calculate Pearson correlation coefficient."""
    x = np.asarray(x)
    y = np.asarray(y)

    xm = x - np.mean(x)
    ym = y - np.mean(y)

    r = np.sum(xm * ym) / np.sqrt(np.sum(xm**2) * np.sum(ym**2))

    return r, None


def _t_cdf(t, df):
    """Approximation of Student's t CDF using error function."""
    import math

    x = df / (df + t * t)
    prob = 0.5 * _betainc(df / 2, 0.5, x)
    if t > 0:
        return 1 - prob
    return prob


def _betainc(a, b, x):
    """Incomplete beta function approximation."""
    if x == 0:
        return 0
    if x == 1:
        return 1

    bt = math.exp(
        _gammaln(a + b)
        - _gammaln(a)
        - _gammaln(b)
        + a * math.log(x)
        + b * math.log(1 - x)
    )

    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b


def _betacf(a, b, x):
    """Continued fraction for incomplete beta function."""
    max_iter = 100
    eps = 3.0e-7
    qab = a + b
    qap = a + 1
    qam = a - 1
    c = 1
    d = 1 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1 / d
    h = d

    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1 / d
        delt = d * c
        h *= delt
        if abs(delt - 1) < eps:
            break

    return h


def _gammaln(x):
    """Logarithm of gamma function."""
    import math

    cof = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ]
    y = x
    tmp = x + 5.5
    tmp -= (x + 0.5) * math.log(tmp)
    ser = 1.000000000190015
    for j in range(6):
        ser += cof[j] / (y + j + 1)
    return -tmp + math.log(2.5066282746310005 * ser / x)
