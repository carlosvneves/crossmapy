import numpy as np


def _get_acceptable_lib_ccm(A, E, tau, plengtht):
    """
    Find acceptable starting positions for CCM library.
    Matches R: acceptablelib <= (plengtht-1), acceptablelib2 < (plengtht-1-tau)
    """
    gapdist = tau * (E - 1)

    acceptablelib = np.isfinite(A).astype(float)

    for i in range(1, gapdist + 1):
        shifted = np.concatenate([np.full(i, np.nan), A[:-i]])
        acceptablelib = acceptablelib * np.isfinite(shifted)

    acceptablelib = np.where(acceptablelib > 0)[0]
    acceptablelib = acceptablelib[acceptablelib <= (plengtht - 1)]
    acceptablelib2 = acceptablelib[acceptablelib < ((plengtht - 1) - tau)]

    return acceptablelib, acceptablelib2


def _build_lagged_array(B, E, tau):
    """Build time-delayed embedding of B."""
    from_idx = tau * (E - 1)
    n_points = len(B) - from_idx

    lagged = np.zeros((n_points, E))
    for k in range(E):
        lagged[:, k] = B[from_idx - tau * k : from_idx - tau * k + n_points]

    return lagged, from_idx


def _get_rho(A, Aest, acceptablelib):
    """Calculate Pearson correlation between A and Aest."""
    xbar = np.mean(A[acceptablelib])
    ybar = np.mean(Aest[acceptablelib])

    xyxybar = np.sum((A[acceptablelib] - xbar) * (Aest[acceptablelib] - ybar))
    xxbarsq = np.sum((A[acceptablelib] - xbar) ** 2)
    yybarsq = np.sum((Aest[acceptablelib] - ybar) ** 2)

    denom = np.sqrt(xxbarsq) * np.sqrt(yybarsq)

    if denom == 0:
        return 0

    rhocalc = xyxybar / denom

    if -1 <= rhocalc <= 1:
        return rhocalc
    return 0


def _ccm_single_iteration(A, B, E, tau, acceptablelib, DesiredL, from_idx):
    """
    Single bootstrap iteration of CCM - optimized version.
    """
    lengtht = len(A[~np.isnan(A)])
    LibLength = lengtht

    Aest = np.zeros(len(A))
    rho_out = np.zeros(len(DesiredL))

    B_lagged, _ = _build_lagged_array(B, E, tau)

    n_neighbors = E + 1
    min_weight = 0.000001

    for lidx, l in enumerate(DesiredL):
        if l < (from_idx + E + 1):
            l = from_idx + E + 1
        if l >= lengtht:
            l = lengtht - 1

        to = l

        lib_indices = np.random.choice(
            acceptablelib, size=to - from_idx + 1, replace=True
        )
        LibUse = np.zeros(len(A), dtype=int)
        LibUse[from_idx : to + 1] = lib_indices

        lib_lagged = np.zeros((to - from_idx + 1, E))
        for k in range(E):
            lib_lagged[:, k] = B[LibUse[from_idx : to + 1] - tau * k]

        Aest.fill(0)

        for ii, i in enumerate(acceptablelib):
            if i < from_idx or i >= len(B):
                continue

            point_lagged = np.zeros(E)
            for k in range(E):
                point_lagged[k] = B[i - tau * k]

            diff = lib_lagged - point_lagged
            distances = np.sqrt(np.sum(diff**2, axis=1))

            # Match C getorder: exclude self-reference (LibUse[ii] != i)
            self_mask = lib_indices == i
            distances_no_self = distances.copy()
            distances_no_self[self_mask] = np.inf

            sorted_idx = np.argsort(distances_no_self)
            neighbor_idx = sorted_idx[:n_neighbors]
            neighbors = lib_indices[neighbor_idx]
            neighbor_dists = distances_no_self[neighbor_idx]

            distsv = neighbor_dists[0]

            if distsv != 0:
                u = np.exp(-neighbor_dists / distsv)
                sumu = np.sum(u)

                w = u / sumu
                w = np.maximum(w, min_weight)
                sumw = np.sum(w)

                w = w / sumw
                Aest[i] = np.sum(A[neighbors] * w)
            else:
                w = np.ones(n_neighbors) * min_weight
                w[neighbor_dists == 0] = 1
                sumw = np.sum(w)
                w = w / sumw
                Aest[i] = np.sum(A[neighbors] * w)

        rho_out[lidx] = _get_rho(A, Aest, acceptablelib)

    return rho_out, Aest


def CCM_boot(A, B, E, tau=1, DesiredL=None, iterations=100):
    """
    Convergent Cross Mapping (CCM) with Bootstrap.
    """
    A = np.asarray(A, dtype=float).copy()
    B = np.asarray(B, dtype=float).copy()

    length_A = len(A)
    finite_A = A[np.isfinite(A)]
    plengtht = len(finite_A)

    if plengtht > len(A):
        plengtht = len(A)

    acceptablelib, acceptablelib2 = _get_acceptable_lib_ccm(A, E, tau, plengtht)
    lengthacceptablelib = len(acceptablelib)

    from_idx = tau * (E - 1)

    if DesiredL is None:
        DesiredL = np.arange(from_idx + E + 1, length_A - E + 2)
    else:
        DesiredL = np.asarray(DesiredL, dtype=int) + E - 2

    valid_L = []
    for dl in DesiredL:
        diff = np.abs(acceptablelib2 - dl)
        if len(diff) > 0:
            valid_L.append(acceptablelib2[np.argmin(diff)])
    DesiredL = np.unique(valid_L)

    if tau * (E + 1) > lengthacceptablelib:
        print(f"Error - too few records to test E = {E} and tau = {tau}")
        return {
            "A": A,
            "Aest": np.full_like(A, np.nan),
            "B": B,
            "rho": np.nan,
            "sdevrho": np.nan,
            "Lobs": np.nan,
            "E": E,
            "tau": tau,
            "FULLinfo": np.nan,
        }

    A[~np.isfinite(A)] = 0
    B[~np.isfinite(B)] = 0

    lpos = set()
    rho_results = []
    Aest_results = []

    np.random.seed(42)

    for it in range(iterations):
        rho_iter, Aest_iter = _ccm_single_iteration(
            A, B, E, tau, acceptablelib, DesiredL, from_idx
        )

        rho_results.append(rho_iter)   # keep full-length array (len == len(DesiredL))
        Aest_results.append(Aest_iter)

        for lob in (DesiredL - E + 1):
            lpos.add(lob)

        if it % 10 == 0:
            print(f"   Iteration {it + 1}/{iterations}...")

    lpos = sorted(list(lpos))
    lpos = np.array(lpos)

    rho_mat = np.full((len(DesiredL), iterations), np.nan)

    for itlst in range(iterations):
        rho_mat[:, itlst] = rho_results[itlst]

    rho_means = np.nanmean(rho_mat, axis=1)
    rho_sdev = np.nanstd(rho_mat, axis=1)

    Aest_avg = np.mean(Aest_results, axis=0)
    Aest_avg[Aest_avg == 0] = np.nan

    return {
        "A": A,
        "Aest": Aest_avg,
        "B": B,
        "rho": rho_means,
        "sdevrho": rho_sdev,
        "Lobs": lpos,
        "E": E,
        "tau": tau,
        "FULLinfo": rho_mat,
    }


def ccmtest(CCM_boot_A, CCM_boot_B):
    """
    CCM significance test.
    """
    FULLinfo_A = CCM_boot_A["FULLinfo"]
    FULLinfo_B = CCM_boot_B["FULLinfo"]

    nrows_A = FULLinfo_A.shape[0]
    nrows_B = FULLinfo_B.shape[0]
    ncols_A = FULLinfo_A.shape[1]
    ncols_B = FULLinfo_B.shape[1]

    if nrows_A > 1:
        pval_a_cause_b = (
            1 - np.sum(FULLinfo_A[0, :] < FULLinfo_A[nrows_A - 1, :]) / ncols_A
        )
    else:
        pval_a_cause_b = np.nan

    if nrows_B > 1:
        pval_b_cause_a = (
            1 - np.sum(FULLinfo_B[0, :] < FULLinfo_B[nrows_B - 1, :]) / ncols_B
        )
    else:
        pval_b_cause_a = np.nan

    return {"pval_a_cause_b": pval_a_cause_b, "pval_b_cause_a": pval_b_cause_a}
