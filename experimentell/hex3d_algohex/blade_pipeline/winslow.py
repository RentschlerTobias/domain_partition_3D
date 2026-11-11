"""Elliptic (Winslow) smoothing of one structured block, boundary fixed.

    winslow(X, iters=200, omega=0.5) -> X'

X is (ni, nj, nk, 3); its six faces stay where they are. Discretised in index
space (unit spacing, central differences):

    sum_i G^ii (x_ii + P_i x_i) + 2 sum_{i<j} G^ij x_ij = 0,

G^ij the contravariant metric (inverse of g_ij = x_i . x_j). The control
functions P_i = -(x_i . x_ii) / (x_i . x_i) are taken once from the INITIAL
grid along each grid line (Thomas-Middlecoff in its 1D form), so the point
distribution along every line -- the hub / shroud BL clustering -- is kept,
while the cross terms relax the angles. Point-Jacobi with relaxation omega,
vectorised over the block interior; every step that would fold a cell (any
corner Jacobian <= 0 in the block) is undone and omega halved.
"""
import numpy as np


def _d(X, ax):
    """First / second central differences along axis ax on the interior."""
    sl = lambda a, b: tuple(slice(a, X.shape[k] - b) if k == ax else slice(1, X.shape[k] - 1) for k in range(3))
    xp, xm, x0 = X[sl(2, 0)], X[sl(0, 2)], X[sl(1, 1)]
    return 0.5 * (xp - xm), xp - 2 * x0 + xm, xp, xm


def _cross(X, a, b):
    def sh(da, db):
        s = []
        for k in range(3):
            o = da if k == a else db if k == b else 0
            s.append(slice(1 + o, X.shape[k] - 1 + o))
        return X[tuple(s)]
    return 0.25 * (sh(1, 1) - sh(1, -1) - sh(-1, 1) + sh(-1, -1))


def _min_jac(X):
    """min corner Jacobian determinant sign test over all cells of the block."""
    a = X[1:, :-1, :-1] - X[:-1, :-1, :-1]
    b = X[:-1, 1:, :-1] - X[:-1, :-1, :-1]
    c = X[:-1, :-1, 1:] - X[:-1, :-1, :-1]
    j0 = np.einsum("...i,...i", a, np.cross(b, c))
    # the opposite corner too (catches most folds cheaply)
    a2 = X[:-1, 1:, 1:] - X[1:, 1:, 1:]
    b2 = X[1:, :-1, 1:] - X[1:, 1:, 1:]
    c2 = X[1:, 1:, :-1] - X[1:, 1:, 1:]
    j1 = -np.einsum("...i,...i", a2, np.cross(b2, c2))
    s = np.sign(np.median(j0))
    return min(float((s * j0).min()), float((s * j1).min()))


def winslow(X, iters=200, omega=0.5, tol=1e-9):
    X = np.array(X, float)
    if min(X.shape[:3]) < 3:
        return X, {"iters": 0}
    P = []
    for ax in range(3):
        x1, x2, _, _ = _d(X, ax)
        P.append(-np.einsum("...i,...i", x1, x2) / np.maximum(np.einsum("...i,...i", x1, x1), 1e-30))
    ok_j = _min_jac(X)
    it = 0
    for it in range(iters):
        D = [_d(X, ax) for ax in range(3)]
        g = np.empty(D[0][0].shape[:3] + (3, 3))
        for i in range(3):
            for j in range(3):
                g[..., i, j] = np.einsum("...k,...k", D[i][0], D[j][0])
        try:
            Gi = np.linalg.inv(g)
        except np.linalg.LinAlgError:
            break
        num = np.zeros_like(D[0][0])
        den = np.zeros(D[0][0].shape[:3])
        for i in range(3):
            x1, _, xp, xm = D[i]
            num += Gi[..., i, i, None] * (xp + xm + P[i][..., None] * x1)
            den += 2 * Gi[..., i, i]
        for a in range(3):
            for b in range(a + 1, 3):
                num += 2 * Gi[..., a, b, None] * _cross(X, a, b)
        new = num / den[..., None]
        inner = X[1:-1, 1:-1, 1:-1]
        step = omega * (new - inner)
        Xn = X.copy()
        Xn[1:-1, 1:-1, 1:-1] = inner + step
        if _min_jac(Xn) <= 0 < ok_j:
            omega *= 0.5
            if omega < 1e-3:
                break
            continue
        X = Xn
        if np.abs(step).max() < tol:
            break
    return X, {"iters": it + 1, "omega": omega}
