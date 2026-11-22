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

from psa.rotations import propagate_batch, quaternion_to_matrix
from psa.design.pulse import PulseSequence
from psa.design.conversion import ps_to_curve, ps_to_tangents
from psa.design import jets as J

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


# ------------------------------------------------------------------- gate ---

def final_quaternion(ps: PulseSequence, offset_hz=0.0) -> np.ndarray:
    return propagate_batch(ps.as_array(), ps.dt, ps.umax, offset_hz)


def gate(ps: PulseSequence) -> dict:
    """Bloch rotation R of the whole pulse: axis (unit vector), angle in [0, 180] deg, matrix,
    and SU(2) quaternion (w, x, y, z) with U = w I - i (x sx + y sy + z sz)."""
    return _gate_from_U(J.product_jets(_segment_jets(ps, order=0))[0])


def _gate_from_U(U):
    q = J.quaternion_of(U)
    q = q / np.linalg.norm(q)
    if q[0] < 0:
        q = -q
    s = np.linalg.norm(q[1:])
    angle = 2 * math.degrees(math.atan2(s, float(q[0])))          # accurate also for tiny angles
    axis = q[1:] / s if s > 1e-12 else np.array([0.0, 0.0, 1.0])
    return {"axis": axis, "angle_deg": angle, "matrix": quaternion_to_matrix(q), "quaternion": q}


def _segment_jets(ps: PulseSequence, perturbation="detuning", order=3) -> np.ndarray:
    """Taylor jets (N, order+1, 2, 2) of every sub-pulse propagator in the detuning delta
    (H -> H + delta sz / 2) or in the relative amplitude error epsilon (H -> (1 + epsilon) H)."""
    w = 2 * np.pi * ps.umax * ps.amplitude / 100.0 * ps.dt
    ph = np.radians(ps.phase)
    v0 = np.column_stack([w * np.cos(ph), w * np.sin(ph), np.zeros(ps.N)])
    d = v0 if perturbation == "amplitude" else np.array([0.0, 0.0, ps.dt])
    return J.step_jets(v0, d, order)


def propagator_jets(ps: PulseSequence, perturbation="detuning", order=3) -> np.ndarray:
    """Jets X0..X_order of U(T; delta) = sum_k delta^k X_k (tree reduction, time order preserved)."""
    return J.product_jets(_segment_jets(ps, perturbation, order))


# kept for older imports
_compose, _to_vector = J.compose, J.to_vector


def _identity_jet():
    return J.identity_jet(4)


def robustness(ps: PulseSequence) -> dict:
    """Exact robustness vectors A1, A2, A3 (s, s^2, s^3) and normalized magnitudes |A_k| / T^k."""
    return _robustness_from_jets(propagator_jets(ps), ps.duration)


def _robustness_from_jets(X, T):
    A = J.magnus_vectors(X)
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
    """First-order robustness against amplitude scaling Omega -> (1 + epsilon) Omega (exact):
    M1 = -i B . sigma / 2 with B = int Omega R^T n dt = -int t x t' dt (radians; the vector area of the
    tangent curve, see the paper). B1_norm = |B| / int Omega dt (0 = first-order B1 robust)."""
    if ps.N == 0:
        return {"B1": np.zeros(3), "B1_norm": 0.0}
    X = propagator_jets(ps, "amplitude", order=1)
    B = J.to_vector(X[0].conj().T @ X[1], 1)
    nutation = _nutation(ps)
    return {"B1": B, "B1_norm": np.linalg.norm(B) / nutation if nutation else 0.0}


def _nutation(ps):
    return float(np.sum(2 * np.pi * ps.umax * np.abs(ps.amplitude) / 100.0) * ps.dt)


def analyze(ps: PulseSequence, b1=True) -> dict:
    """Gate, A1..A3 and (optionally) B1 in one pass over the pulse (the gate comes with the jets)."""
    X = propagator_jets(ps)
    out = {"gate": _gate_from_U(X[0]), "robustness": _robustness_from_jets(X, ps.duration)}
    if b1:
        out["B1"] = amplitude_robustness(ps)
    return out


def amplitude_robustness_geometric(ps: PulseSequence) -> np.ndarray:
    """-int t x t' dt from the sampled tangents (independent cross-check of ``amplitude_robustness``)."""
    t = ps_to_tangents(ps)
    return -np.sum(np.cross((t[:-1] + t[1:]) / 2, np.diff(t, axis=0)), axis=0)


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
