import numpy as np


def _get_acceptable_lib(A, E, tau, predstep=1):
    """
    Find acceptable starting positions for the library.

    These are indices that have enough lagged values for reconstruction
    and don't jump over NA gaps.
    """
    gapdist = tau * (E - 1) + predstep
    length_A = len(A)

    acceptablelib = np.isfinite(A).astype(float)

    for i in range(1, gapdist + 1):
        shifted = np.concatenate([np.full(i, np.nan), A[:-i]])
        acceptablelib = acceptablelib * np.isfinite(shifted)

    acceptablelib = np.where(acceptablelib > 0)[0]
    acceptablelib = acceptablelib[acceptablelib <= length_A - 1]

    return acceptablelib


def _getorder_ssr(distances, E, acceptablelib, i, predstep):
    """
    Exact Python port of C's getorder_ssr.

    Critical: due to the j>0 condition in the inner loop,
    neighbors[0] is always fixed to acceptablelib[0] (or [1] if self-ref)
    and is never replaced. This matches R's C implementation exactly.
    """
    nneigh = 1
    n = 0
    length = len(acceptablelib)

    if acceptablelib[0] == i:
        n = 1
    if n >= length:
        n = length - 1

    neighbors = [0] * (E + 1)
    neighbors[0] = acceptablelib[n]

    for iii in range(n, length - predstep):
        ii = acceptablelib[iii]
        trip = False

        for j in range(nneigh):
            # j>0: neighbors[0] is NEVER replaced (matches C exactly)
            if (distances[ii] < distances[neighbors[j]]) and (ii != i) and (j > 0):
                for k in range(nneigh, j, -1):
                    if k < E + 1:
                        neighbors[k] = neighbors[k - 1]
                neighbors[j] = ii
                # nneigh NOT incremented: C checks (trip!=0) before setting trip=1
                trip = True
                break

        if (not trip) and (nneigh < E + 1) and (ii != i) and (neighbors[nneigh - 1] != ii):
            neighbors[nneigh] = ii
            if nneigh < E + 1:
                nneigh += 1

    return neighbors, nneigh


def SSR_pred_boot(A, B=None, E=2, tau=1, predstep=1, matchSugi=0):
    """
    Simplex Projection / Single Series Prediction with Bootstrap.

    Uses information in B to predict future values of A.
    If A=B, uses leave-one-out cross-validation.

    Parameters
    ----------
    A : array-like
        Target time series to predict
    B : array-like, optional
        Library time series (default: A)
    E : int
        Embedding dimension
    tau : int
        Time lag
    predstep : int
        Number of steps ahead to predict
    matchSugi : int
        If 1, removes only point i in cross-validation (Sugihara's method)
        If 0, removes all points within X(t-(E-1)):X(t+1)

    Returns
    -------
    dict
        Dictionary with 'A', 'Aest', 'B', 'rho', 'E', 'tau', 'predstep'
    """
    A = np.asarray(A, dtype=float).copy()

    if B is None:
        B = A.copy()
    else:
        B = np.asarray(B, dtype=float).copy()

    repvec = (
        1
        if (
            np.sum(np.isfinite(A) & np.isfinite(B) & (A == B)) == np.sum(np.isfinite(A))
        )
        and (len(A) == len(B))
        else 0
    )

    acceptablelib = _get_acceptable_lib(A, E, tau, predstep)
    lengthacceptablelib = len(acceptablelib)

    if tau * (E + 1) + predstep >= lengthacceptablelib:
        print(
            f"Error - too few records to test E = {E}, tau = {tau}, and predstep = {predstep}"
        )
        return {
            "A": A,
            "Aest": np.full_like(A, np.nan),
            "B": B,
            "E": E,
            "tau": tau,
            "pBlength": len(B),
            "pAlength": len(A),
            "predstep": predstep,
            "rho": np.nan,
            "acceptablelib": acceptablelib,
            "plengthacceptablelib": lengthacceptablelib,
        }

    A[~np.isfinite(A)] = 0
    B[~np.isfinite(B)] = 0

    from_idx = tau * (E - 1)
    nneigh = E + 1

    Aest = np.zeros(len(A))
    maxdist = 0

    for ii in range(lengthacceptablelib):
        i = acceptablelib[ii]

        distances = np.full(len(B), np.inf)

        if repvec == 1:
            for jj in range(lengthacceptablelib - predstep):
                j = acceptablelib[jj]

                if matchSugi == 1:
                    if i != j:
                        dist = 0.0
                        for k in range(E):
                            dist += (A[i - tau * k] - B[j - tau * k]) ** 2
                        distances[j] = np.sqrt(dist)
                        if distances[j] > maxdist:
                            maxdist = 999999999 * distances[j]
                    else:
                        distances[j] = maxdist
                else:
                    if (j > i + predstep) or (j <= i - E):
                        dist = 0.0
                        for k in range(E):
                            dist += (A[i - tau * k] - B[j - tau * k]) ** 2
                        distances[j] = np.sqrt(dist)
                        if distances[j] > maxdist:
                            maxdist = 999999999 * distances[j]
                    else:
                        distances[j] = maxdist
        else:
            for jj in range(lengthacceptablelib):
                j = acceptablelib[jj]

                dist = 0.0
                for k in range(E):
                    dist += (A[i - tau * k] - B[j - tau * k]) ** 2
                distances[j] = np.sqrt(dist)

        neighbors, found_nneigh = _getorder_ssr(distances, E, acceptablelib, i, predstep)
        distsv = distances[neighbors[0]]

        if found_nneigh < nneigh:
            Aest[i] = 0
            continue

        sumaest = 0.0

        if distsv != 0:
            u = np.zeros(nneigh)
            w = np.zeros(nneigh)

            for j in range(nneigh):
                u[j] = np.exp(-distances[neighbors[j]] / distsv)

            sumu = np.sum(u)

            for j in range(nneigh):
                w[j] = u[j] / sumu
                if w[j] < 0.000001:
                    w[j] = 0.000001

            sumw = np.sum(w)

            for j in range(nneigh):
                w[j] = w[j] / sumw
                sumaest += B[neighbors[j] + predstep] * w[j]
        else:
            w = np.zeros(nneigh)
            sumw = 0.0

            for j in range(nneigh):
                if distances[neighbors[j]] == 0:
                    w[j] = 1
                else:
                    w[j] = 0.000001
                sumw += w[j]

            for j in range(nneigh):
                w[j] = w[j] / sumw
                sumaest += A[neighbors[j]] * w[j]

        Aest[i] = sumaest

    # Match R: out$Aest[(1+predstep):N] <- out$Aest[1:(N-predstep)]
    # Then out$Aest[out$Aest==0] <- NA
    Aest_out = np.zeros(len(A))
    Aest_out[predstep:] = Aest[: len(Aest) - predstep]
    Aest_out[Aest_out == 0] = np.nan

    valid = np.isfinite(A) & np.isfinite(Aest_out)
    if np.sum(valid) > 1:
        rho = np.corrcoef(A[valid], Aest_out[valid])[0, 1]
    else:
        rho = np.nan

    return {
        "A": A,
        "Aest": Aest_out,
        "B": B,
        "E": E,
        "tau": tau,
        "pBlength": len(B),
        "pAlength": len(A),
        "predstep": predstep,
        "rho": rho,
        "acceptablelib": acceptablelib,
        "plengthacceptablelib": lengthacceptablelib,
    }
