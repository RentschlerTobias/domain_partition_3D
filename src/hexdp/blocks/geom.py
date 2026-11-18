"""Small geometry helpers shared by the building blocks."""

from __future__ import annotations

import numpy as np


def rot_z(X, a):
    c, s = np.cos(a), np.sin(a)
    Y = np.array(X, float, copy=True)
    Y[..., 0] = c * X[..., 0] - s * X[..., 1]
    Y[..., 1] = s * X[..., 0] + c * X[..., 1]
    return Y


def arcfrac(Q):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
    return s / s[-1]


def at_fraction(Q, t):
    s = arcfrac(Q)
    return np.array([np.interp(t, s, Q[:, k]) for k in range(Q.shape[1])])


def area2d(P):
    return 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])


def winslow_on_cylinder(L, iters=2000, tol=1e-10):
    """Interior of a structured section grid lying on r = const (approximately): unroll to (r phi, z), Winslow
    with the boundary fixed, roll back (exactly onto r when the boundary is)."""
    r = np.hypot(L[..., 0], L[..., 1]).mean()
    ph = np.unwrap(np.unwrap(np.arctan2(L[..., 1], L[..., 0]), axis=0), axis=1)
    X = np.dstack([r * ph, L[..., 2]])
    for _ in range(iters):
        E, W, N, S = X[2:, 1:-1], X[:-2, 1:-1], X[1:-1, 2:], X[1:-1, :-2]
        NE, NW, SE, SW = X[2:, 2:], X[:-2, 2:], X[2:, :-2], X[:-2, :-2]
        xi, eta = 0.5 * (E - W), 0.5 * (N - S)
        a = (eta ** 2).sum(-1)[..., None]
        g = (xi ** 2).sum(-1)[..., None]
        b = (xi * eta).sum(-1)[..., None]
        new = (a * (E + W) + g * (N + S) - 0.5 * b * (NE - NW - SE + SW)) / (2 * (a + g) + 1e-300)
        d = np.abs(new - X[1:-1, 1:-1]).max() if new.size else 0.0
        X[1:-1, 1:-1] = new
        if d < tol:
            break
    rb = np.r_[np.hypot(L[0, :, 0], L[0, :, 1]), np.hypot(L[-1, :, 0], L[-1, :, 1])]
    out = L.copy()
    rr = (np.full(L[1:-1, 1:-1, 0].shape, r) if np.ptp(rb) < 1e-6
          else np.hypot(L[1:-1, 1:-1, 0], L[1:-1, 1:-1, 1]))
    ph2 = X[1:-1, 1:-1, 0] / r
    out[1:-1, 1:-1, 0] = rr * np.cos(ph2)
    out[1:-1, 1:-1, 1] = rr * np.sin(ph2)
    out[1:-1, 1:-1, 2] = X[1:-1, 1:-1, 1]
    return out


def tfi(L):
    ni, nj = L.shape[:2]
    for i in range(1, ni - 1):
        for j in range(1, nj - 1):
            u, v = i / (ni - 1), j / (nj - 1)
            L[i, j] = ((1 - v) * L[i, 0] + v * L[i, -1] + (1 - u) * L[0, j] + u * L[-1, j]
                       - ((1 - u) * (1 - v) * L[0, 0] + u * (1 - v) * L[-1, 0] + u * v * L[-1, -1]
                          + (1 - u) * v * L[0, -1]))
    return L
