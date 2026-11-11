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


def winslow2d(X, project=None, iters=300, omega=0.5, every=10, tol=1e-10):
    """2D Winslow on a structured surface patch X (n0, n1, 3), boundary rows fixed.

    Same scheme as `winslow` in two index directions, written in the 3D
    embedding (the metric g_ij = x_i . x_j works on a curved patch);
    `project(points) -> points` puts the interior back onto the exact surface
    every `every` steps and at the end. Control functions from the initial
    patch keep the distribution along each grid line (the BL along span edges).
    A step that flips a quad against the patch's mean normal is undone."""
    X = np.array(X, float)
    if min(X.shape[:2]) < 3:
        return X, {"iters": 0}

    def d(Y, ax):
        if ax == 0:
            xp, xm, x0 = Y[2:, 1:-1], Y[:-2, 1:-1], Y[1:-1, 1:-1]
        else:
            xp, xm, x0 = Y[1:-1, 2:], Y[1:-1, :-2], Y[1:-1, 1:-1]
        return 0.5 * (xp - xm), xp - 2 * x0 + xm, xp, xm

    def flips(Y):
        a = Y[1:, :-1] - Y[:-1, :-1]
        b = Y[:-1, 1:] - Y[:-1, :-1]
        c = Y[1:, 1:] - Y[1:, :-1]
        e = Y[1:, 1:] - Y[:-1, 1:]
        n1 = np.cross(a, b)
        n2 = np.cross(c, -e)
        nm = (n1 + n2).reshape(-1, 3).sum(0)
        return int(((n1 @ nm) <= 0).sum() + ((n2 @ nm) <= 0).sum())

    P = []
    for ax in range(2):
        x1, x2, _, _ = d(X, ax)
        P.append(-np.einsum("...i,...i", x1, x2) / np.maximum(np.einsum("...i,...i", x1, x1), 1e-30))
    f0 = flips(X)
    it = 0
    for it in range(iters):
        D = [d(X, 0), d(X, 1)]
        g11 = np.einsum("...i,...i", D[0][0], D[0][0])
        g22 = np.einsum("...i,...i", D[1][0], D[1][0])
        g12 = np.einsum("...i,...i", D[0][0], D[1][0])
        # contravariant: G11 = g22 / J2, G22 = g11 / J2, G12 = -g12 / J2 (J2 cancels)
        xc = 0.25 * (X[2:, 2:] - X[2:, :-2] - X[:-2, 2:] + X[:-2, :-2])
        num = (g22[..., None] * (D[0][2] + D[0][3] + P[0][..., None] * D[0][0])
               + g11[..., None] * (D[1][2] + D[1][3] + P[1][..., None] * D[1][0])
               - 2 * g12[..., None] * xc)
        new = num / (2 * (g11 + g22))[..., None]
        step = omega * (new - X[1:-1, 1:-1])
        Xn = X.copy()
        Xn[1:-1, 1:-1] += step
        if project is not None and (it % every == every - 1):
            Xn[1:-1, 1:-1] = project(Xn[1:-1, 1:-1].reshape(-1, 3)).reshape(Xn[1:-1, 1:-1].shape)
        if flips(Xn) > f0:
            omega *= 0.5
            if omega < 1e-3:
                break
            continue
        X = Xn
        if np.abs(step).max() < tol:
            break
    if project is not None:
        X[1:-1, 1:-1] = project(X[1:-1, 1:-1].reshape(-1, 3)).reshape(X[1:-1, 1:-1].shape)
    return X, {"iters": it + 1, "omega": omega}
