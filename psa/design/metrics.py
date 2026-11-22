"""Gate, robustness and fidelity of a pulse sequence (on resonance, detuning noise along z).

Robustness vectors use the convention verified against the Magnus expansion:
    V(T; delta) = U_0(T)^dag U(T; delta) = exp(sum_k delta^k M_k),  M_k = -i A_k . sigma / 2^k,
    A_1 = r(T),   A_2 = -int r x r' dt,   A_3 = 2/3 int [S x (r' x P) + P x (r' x S)] dt,
with P(t) = r(t), S(t) = r(T) - r(t) (chords). Units: A_k in s^k (delta in rad/s).
``robustness`` computes them exactly for piecewise-constant pulses (Taylor jets of the propagator);
``robustness_geometric`` evaluates the curve formulas (used as an independent cross-check).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import expm

from psa.rotations import propagate_batch, quaternion_to_matrix
from psa.design.pulse import PulseSequence
from psa.design.conversion import ps_to_curve, ps_to_tangents

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


# ------------------------------------------------------------------- gate ---

def final_quaternion(ps: PulseSequence, offset_hz=0.0) -> np.ndarray:
    return propagate_batch(ps.as_array(), ps.dt, ps.umax, offset_hz)


def gate(ps: PulseSequence) -> dict:
    """Bloch rotation R of the whole pulse: axis (unit vector), angle in [0, 180] deg, matrix,
    and SU(2) quaternion (w, x, y, z) with U = w I - i (x sx + y sy + z sz)."""
    q = final_quaternion(ps)
    if q[0] < 0:
        q = -q
    w = min(1.0, max(-1.0, float(q[0])))
    angle = 2 * math.degrees(math.acos(w))
    s = np.linalg.norm(q[1:])
    axis = q[1:] / s if s > 1e-12 else np.array([0.0, 0.0, 1.0])
    return {"axis": axis, "angle_deg": angle, "matrix": quaternion_to_matrix(q), "quaternion": q}


# ------------------------------------------------------------- robustness ---

def _segment_jets(ps: PulseSequence) -> np.ndarray:
    """Taylor jets (X0..X3) in delta of every sub-pulse propagator, shape (N, 4, 2, 2)."""
    ux = 2 * np.pi * ps.umax * ps.amplitude / 100.0 * np.cos(np.radians(ps.phase))
    uy = 2 * np.pi * ps.umax * ps.amplitude / 100.0 * np.sin(np.radians(ps.phase))
    A0 = -0.5j * (ux[:, None, None] * SX + uy[:, None, None] * SY)
    A1 = np.broadcast_to(-0.5j * SZ, A0.shape)
    Z2 = np.zeros_like(A0)
    big = np.block([[A0, Z2, Z2, Z2], [A1, A0, Z2, Z2], [Z2, A1, A0, Z2], [Z2, Z2, A1, A0]]) * ps.dt
    E = expm(big)                                  # (N, 8, 8); first block column = jets of exp
    return np.stack([E[:, 2 * k:2 * k + 2, 0:2] for k in range(4)], axis=1)


def _compose(later, earlier):
    """Jets of the product later @ earlier (truncated Cauchy product), arrays (..., 4, 2, 2)."""
    out = np.zeros_like(earlier)
    for k in range(4):
        for i in range(k + 1):
            out[..., k, :, :] += later[..., i, :, :] @ earlier[..., k - i, :, :]
    return out


def propagator_jets(ps: PulseSequence) -> np.ndarray:
    """Jets X0..X3 of U(T; delta) = sum_k delta^k X_k (tree reduction, time order preserved)."""
    jets = _segment_jets(ps)
    if len(jets) == 0:
        out = np.zeros((4, 2, 2), dtype=complex); out[0] = np.eye(2)
        return out
    while len(jets) > 1:
        if len(jets) % 2:
            jets = np.concatenate([jets, _identity_jet()[None]], axis=0)
        jets = _compose(jets[1::2], jets[0::2])
    return jets[0]


def _identity_jet():
    out = np.zeros((4, 2, 2), dtype=complex)
    out[0] = np.eye(2)
    return out


def _to_vector(M, k):
    # M_k = -i A_k . sigma / 2^k  ->  A_k,j = i 2^(k-1) tr(M_k sigma_j)
    return np.real([1j * 2 ** (k - 1) * np.trace(M @ s) for s in (SX, SY, SZ)])


def robustness(ps: PulseSequence) -> dict:
    """Exact robustness vectors A1, A2, A3 (s, s^2, s^3) and normalized magnitudes |A_k| / T^k."""
    X = propagator_jets(ps)
    U0d = X[0].conj().T
    D1, D2, D3 = U0d @ X[1], U0d @ X[2], U0d @ X[3]
    M1 = D1
    M2 = D2 - D1 @ D1 / 2
    M3 = D3 - (D1 @ D2 + D2 @ D1) / 2 + D1 @ D1 @ D1 / 3
    A = [_to_vector(M, k) for k, M in enumerate((M1, M2, M3), start=1)]
    T = ps.duration
    return {"A1": A[0], "A2": A[1], "A3": A[2],
            "A1_norm": np.linalg.norm(A[0]) / T if T else 0.0,
            "A2_norm": np.linalg.norm(A[1]) / T ** 2 if T else 0.0,
            "A3_norm": np.linalg.norm(A[2]) / T ** 3 if T else 0.0}


def robustness_geometric(points, tangents=None) -> dict:
    """A1..A3 from the curve formulas (midpoint rule on the curve points)."""
    r = np.asarray(points, dtype=np.float64)
    dr = np.diff(r, axis=0)
    P = (r[:-1] + r[1:]) / 2
    S = r[-1] - P
    A1 = r[-1] - r[0]
    A2 = -np.sum(np.cross(P, dr), axis=0)
    A3 = 2.0 / 3.0 * np.sum(np.cross(S, np.cross(dr, P)) + np.cross(P, np.cross(dr, S)), axis=0)
    return {"A1": A1, "A2": A2, "A3": A3}


def amplitude_robustness(ps: PulseSequence) -> dict:
    """First-order amplitude-scaling (B1) robustness: the vector area of the tangent curve,
    -int t x t' dt (zero for first-order B1-robust UR pulses; see the paper)."""
    t = ps_to_tangents(ps)
    dt_ = np.diff(t, axis=0)
    mid = (t[:-1] + t[1:]) / 2
    B = -np.sum(np.cross(mid, dt_), axis=0) * ps.dt
    T = ps.duration
    return {"B1": B, "B1_norm": np.linalg.norm(B) / T if T else 0.0}


# --------------------------------------------------------------- fidelity ---

def fidelity_vs_offset(ps: PulseSequence, offsets_hz) -> np.ndarray:
    """Gate fidelity (|tr U(0)^dag U(offset)|^2 + 2) / 6 of the pulse at each offset relative to its
    own on-resonance gate."""
    q0 = final_quaternion(ps)
    q = propagate_batch(ps.as_array(), ps.dt, ps.umax, np.asarray(offsets_hz, dtype=np.float64))
    overlap = q @ q0                      # tr(U0^dag U) = 2 <q0, q>
    return (4 * overlap ** 2 + 2) / 6


def summary(ps: PulseSequence) -> dict:
    """Everything the editor's shared panel shows."""
    g = gate(ps)
    rb = robustness(ps)
    curve = ps_to_curve(ps)
    return {"gate": g, "robustness": rb, "curve": curve, "duration": ps.duration,
            "max_amplitude_percent": float(np.max(ps.amplitude)) if ps.N else 0.0}
