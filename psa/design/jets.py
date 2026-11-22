"""Closed-form Taylor jets of SU(2) exponentials (fast replacement of block-matrix exponentials).

For a constant step U(eps) = exp(-i (v0 + eps d) . sigma / 2) write x(eps) = |v0 + eps d|^2 and
    U = F(x) I - i G(x) (v0 + eps d) . sigma,   F(x) = cos(sqrt(x)/2),   G(x) = sin(sqrt(x)/2) / sqrt(x).
F and G are entire functions of x, so the eps-Taylor coefficients of U follow from their x-derivatives
(closed forms for large rotations, power series for small ones) composed with the quadratic x(eps).
Jets of products are obtained with a truncated Cauchy product in a tree reduction (time order kept).
"""
from __future__ import annotations

import math

import numpy as np

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

_SERIES_X = 4.0          # use the power series for x = r^2 below this (r < 2)
_N_TERMS = 16


def _series(coeffs_fn, x, order):
    """Values of the 0..order-th x-derivatives of sum_n a_n x^n."""
    out = np.zeros((order + 1,) + x.shape)
    for n in range(_N_TERMS):
        a = coeffs_fn(n)
        for k in range(min(order, n) + 1):
            out[k] += a * (math.factorial(n) / math.factorial(n - k)) * x ** (n - k)
    return out


def _f_coeff(n):
    return (-1) ** n / (4.0 ** n * math.factorial(2 * n))


def _g_coeff(n):
    return (-1) ** n / (2.0 * 4.0 ** n * math.factorial(2 * n + 1))


def fg_derivatives(x, order):
    """F, G and their x-derivatives up to ``order`` (<= 4) at x >= 0; arrays of shape (order+1, ...)."""
    x = np.asarray(x, dtype=np.float64)
    small = x < _SERIES_X
    F = np.empty((order + 1,) + x.shape)
    G = np.empty((order + 1,) + x.shape)
    if np.any(small):
        F[:, small] = _series(_f_coeff, x[small], order)
        G[:, small] = _series(_g_coeff, x[small], order)
    big = ~small
    if np.any(big):
        r = np.sqrt(x[big]); s = np.sin(r / 2); c = np.cos(r / 2)
        g = [s / r, (c * r - 2 * s) / (4 * r ** 3), (-6 * c * r - r ** 2 * s + 12 * s) / (16 * r ** 5),
             (-c * r ** 3 + 60 * c * r + 12 * r ** 2 * s - 120 * s) / (64 * r ** 7),
             (20 * c * r ** 3 - 840 * c * r + r ** 4 * s - 180 * r ** 2 * s + 1680 * s) / (256 * r ** 9)]
        f = [c] + [-gk / 4 for gk in g[:4]]                     # F^(k+1) = -G^(k) / 4
        for k in range(order + 1):
            F[k, big] = f[k]
            G[k, big] = g[k]
    return F, G


def _compose_series(D, x1, x2, order):
    """Taylor coefficients in eps of f(x0 + x1 eps + x2 eps^2) from derivatives D[k] = f^(k)(x0)."""
    out = [D[0]]
    if order >= 1:
        out.append(D[1] * x1)
    if order >= 2:
        out.append(D[1] * x2 + D[2] / 2 * x1 ** 2)
    if order >= 3:
        out.append(D[2] * x1 * x2 + D[3] / 6 * x1 ** 3)
    return out


def step_jets(v0, d, order=3):
    """Jets (N, order+1, 2, 2) of exp(-i (v0 + eps d) . sigma / 2) for rotation vectors v0, d of shape (N, 3)."""
    v0 = np.asarray(v0, dtype=np.float64).reshape(-1, 3)
    d = np.broadcast_to(np.asarray(d, dtype=np.float64), v0.shape)
    x0 = np.einsum("ij,ij->i", v0, v0)
    x1 = 2 * np.einsum("ij,ij->i", v0, d)
    x2 = np.einsum("ij,ij->i", d, d)
    F, G = fg_derivatives(x0, order)
    w = _compose_series(F, x1, x2, order)
    g = _compose_series(G, x1, x2, order)
    out = np.zeros((len(v0), order + 1, 2, 2), dtype=complex)
    for k in range(order + 1):
        vec = g[k][:, None] * v0 + (g[k - 1][:, None] * d if k else 0.0)
        out[:, k] = (w[k][:, None, None] * np.eye(2)
                     - 1j * (vec[:, 0, None, None] * SX + vec[:, 1, None, None] * SY + vec[:, 2, None, None] * SZ))
    return out


def compose(later, earlier):
    """Jets of later @ earlier (truncated Cauchy product), arrays (..., K, 2, 2)."""
    K = earlier.shape[-3]
    out = np.zeros(np.broadcast_shapes(later.shape, earlier.shape), dtype=complex)
    for k in range(K):
        for i in range(k + 1):
            out[..., k, :, :] += later[..., i, :, :] @ earlier[..., k - i, :, :]
    return out


def identity_jet(K):
    out = np.zeros((K, 2, 2), dtype=complex)
    out[0] = np.eye(2)
    return out


def product_jets(jets):
    """Jets of U_N ... U_2 U_1 for step jets (N, K, 2, 2) (first step first)."""
    K = jets.shape[1]
    if len(jets) == 0:
        return identity_jet(K)
    while len(jets) > 1:
        if len(jets) % 2:
            jets = np.concatenate([jets, identity_jet(K)[None]], axis=0)
        jets = compose(jets[1::2], jets[0::2])
    return jets[0]


def to_vector(M, k):
    """A with M = -i A . sigma / 2^k."""
    return np.real([1j * 2 ** (k - 1) * np.trace(M @ s) for s in (SX, SY, SZ)])


def magnus_vectors(X):
    """A1, A2, A3 from the jets X0..X3 of U(T; eps) (M_k of log(U0^dag U))."""
    U0d = X[0].conj().T
    D1, D2, D3 = U0d @ X[1], U0d @ X[2], U0d @ X[3]
    M = (D1, D2 - D1 @ D1 / 2, D3 - (D1 @ D2 + D2 @ D1) / 2 + D1 @ D1 @ D1 / 3)
    return [to_vector(m, k) for k, m in enumerate(M, start=1)]


def quaternion_of(U):
    """(w, x, y, z) of U = w I - i (x sx + y sy + z sz)."""
    return np.array([U[0, 0].real, -U[1, 0].imag, U[1, 0].real, -U[0, 0].imag])


# --------------------------------------------- derivatives with respect to controls ---

def _pmul(a, b, K):
    """Product of polynomials in (delta, eta), truncated at delta^K and eta^1; arrays (K+1, 2, ...)."""
    out = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=np.result_type(a, b))
    for k1 in range(K + 1):
        for j1 in range(2):
            for k2 in range(K + 1 - k1):
                for j2 in range(2 - j1):
                    out[k1 + k2, j1 + j2] += a[k1, j1] * b[k2, j2]
    return out


def step_jets_with_derivative(v0, d, e, f=None, order=3):
    """Jets of exp(-i w . sigma / 2) with w = v0 + delta d + eta e + delta eta f (all (N, 3)):
    returns (J, dJ), both (N, order+1, 2, 2): the delta-Taylor coefficients at eta = 0 and their
    derivative with respect to eta (e.g. a control direction)."""
    K = order
    v0 = np.asarray(v0, dtype=np.float64).reshape(-1, 3)
    N = len(v0)
    d = np.broadcast_to(np.asarray(d, dtype=np.float64), (N, 3))
    e = np.broadcast_to(np.asarray(e, dtype=np.float64), (N, 3))
    f = np.zeros((N, 3)) if f is None else np.broadcast_to(np.asarray(f, dtype=np.float64), (N, 3))
    W = np.zeros((K + 1, 2, N, 3))                   # w as polynomial
    W[0, 0] = v0
    if K >= 1:
        W[1, 0] = d
        W[1, 1] = f
    W[0, 1] = e
    dot = lambda a, b: np.einsum("...i,...i->...", a, b)
    X = np.zeros((K + 1, 2, N))                      # x = w . w (without the constant)
    if K >= 1:
        X[1, 0] = 2 * dot(v0, d)
        X[1, 1] = 2 * (dot(d, e) + dot(v0, f))
    if K >= 2:
        X[2, 0] = dot(d, d)
        X[2, 1] = 2 * dot(d, f)
    X[0, 1] = 2 * dot(v0, e)
    F, G = fg_derivatives(dot(v0, v0), K + 1)
    Fp = np.zeros((K + 1, 2, N)); Gp = np.zeros((K + 1, 2, N))
    power = np.zeros((K + 1, 2, N)); power[0, 0] = 1.0
    for j in range(K + 2):
        Fp += F[j] / math.factorial(j) * power
        Gp += G[j] / math.factorial(j) * power
        power = _pmul(power, X, K)
    V = _pmul(Gp[..., None], W, K)                   # (K+1, 2, N, 3)
    eye = np.eye(2)
    def mats(j):
        m = (Fp[:, j, :, None, None] * eye - 1j * (V[:, j, :, 0, None, None] * SX + V[:, j, :, 1, None, None] * SY
                                                     + V[:, j, :, 2, None, None] * SZ))
        return np.moveaxis(m, 0, 1)                  # (N, K+1, 2, 2)
    return mats(0), mats(1)


def prefix_products(J):
    """P[k] = J[k] o ... o J[0] for jets (N, K, 2, 2) (Hillis-Steele scan, log2 N vectorized steps)."""
    P = J.copy()
    s = 1
    while s < len(P):
        P[s:] = compose(P[s:], P[:-s])
        s *= 2
    return P


def suffix_products(J):
    """S[k] = J[N-1] o ... o J[k] for jets (N, K, 2, 2)."""
    S = J.copy()
    s = 1
    while s < len(S):
        S[:-s] = compose(S[s:], S[:-s])
        s *= 2
    return S
