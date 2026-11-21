"""Vectorized NumPy rotation backend (replaces the former Rust extensions).

Conventions are identical to the original implementation:

* Sub-pulse k (row ``[B1 %, phase deg]`` of ``PS``) has the rotation vector
  ``n_k = (B1/100 * Umax * inpoFact * cos(phase), ... * sin(phase), offset)`` and the
  rotation angle ``2*pi * T/inpoFact * |n_k|`` (``Umax`` and ``offset`` in Hz).
* ``R_k`` is the forward Rodrigues rotation of sub-pulse k, accumulated as
  ``R_(k) = R_k @ R_(k-1)`` with ``R_(0) = I``.
* Error-curve tangents are ``R_(k)^T (l * initialVector)``, k = 0..N.
* The "total rotation" used by the fidelity calculations is ``R_(N)^T``.

Rotations are composed as unit quaternions (w, x, y, z), which is cheaper than 3x3
matrix products. Batched functions take ``Umax`` and ``offset`` as arrays and propagate
all grid points simultaneously, looping only over the sub-pulses.
"""

import math
import numpy as np


def _rotation_vectors(PS, Umax, offset, inpoFact):
    PS = np.asarray(PS, dtype=np.float64)
    amplitude = PS[:, 0] * (Umax * inpoFact / 100.0)
    phase = np.radians(PS[:, 1])
    return np.stack([amplitude * np.cos(phase), amplitude * np.sin(phase), np.full(len(PS), float(offset))], axis=1)


def step_quaternions(PS, T, Umax, offset, inpoFact=1):
    """Quaternions (N, 4) of the forward rotations of all sub-pulses."""
    n = _rotation_vectors(PS, Umax, offset, inpoFact)
    norm = np.linalg.norm(n, axis=1)
    half_angle = (np.pi * T / inpoFact) * norm
    safe = np.where(norm != 0, norm, 1.0)
    scale = np.where(norm != 0, np.sin(half_angle) / safe, 0.0)
    return np.column_stack([np.cos(half_angle), n * scale[:, None]])


def accumulate_quaternions(steps):
    """Cumulative products q_(k) = step_k * q_(k-1), q_(0) = 1. Returns (N + 1, 4)."""
    q = np.empty((len(steps) + 1, 4))
    w, x, y, z = 1.0, 0.0, 0.0, 0.0
    q[0] = (w, x, y, z)
    for k, (c, a, b, d) in enumerate(steps.tolist(), start=1):
        w, x, y, z = (c * w - a * x - b * y - d * z,
                      c * x + a * w + b * z - d * y,
                      c * y + b * w + d * x - a * z,
                      c * z + d * w + a * y - b * x)
        q[k] = (w, x, y, z)
    return q


def quaternion_to_matrix(q):
    """Rotation matrices (..., 3, 3) of unit quaternions (..., 4)."""
    q = np.asarray(q, dtype=np.float64)
    w, x, y, z = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], axis=-1),
        np.stack([2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], axis=-1),
        np.stack([2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)], axis=-1),
    ], axis=-2)


def rotate(q, v):
    """R(q) v for quaternions (..., 4) and a vector v (3,) or vectors (..., 3).

    With q the final quaternion of a pulse, this is the Bloch vector after the pulse
    for initial Bloch vector v (forward, time-ordered evolution).
    """
    q = np.asarray(q, dtype=np.float64)
    w = q[..., 0:1]
    u = q[..., 1:4]
    v = np.broadcast_to(np.asarray(v, dtype=np.float64), u.shape)
    t = 2.0 * np.cross(u, v)
    return v + w * t + np.cross(u, t)


def inverse_rotate(q, v):
    """R(q)^T v for quaternions (..., 4) and a vector v (3,) or vectors (..., 3).

    With q_(k) this gives the error-curve tangent (toggling-frame image of v).
    """
    q = np.asarray(q, dtype=np.float64)
    return rotate(np.concatenate([q[..., 0:1], -q[..., 1:4]], axis=-1), v)


def tangent_vectors(PS, T, l, Umax, offset, inpoFact, initialVector):
    """Error-curve tangent vectors (N + 1, 3) including the initial one (former createVectors_Matrix)."""
    q = accumulate_quaternions(step_quaternions(PS, T, Umax, offset, inpoFact))
    return inverse_rotate(q, l * np.asarray(initialVector, dtype=np.float64))


def inverse_total_rotation(PS, T, Umax, offset, inpoFact=1):
    """R_(N)^T as a 3x3 matrix (former create_Rtot)."""
    q = accumulate_quaternions(step_quaternions(PS, T, Umax, offset, inpoFact))
    return quaternion_to_matrix(q[-1]).T


def propagate_batch(PS, T, Umax, offset, inpoFact=1, initialVector=None, l=1.0):
    """Propagate many (Umax, offset) points at once.

    Umax, offset: scalars or arrays broadcastable to a common shape S.
    Returns the final quaternions (S + (4,)). If ``initialVector`` is given, additionally
    returns the end point of the error curve with unit step length (S + (3,)), integrated with the
    trapezoid rule over the N + 1 tangents t_k = R_(k)^T (l * initialVector):
    sum_k t_k - (t_0 + t_N) / 2.
    """
    PS = np.asarray(PS, dtype=np.float64)
    Umax, offset = np.broadcast_arrays(np.asarray(Umax, dtype=np.float64), np.asarray(offset, dtype=np.float64))
    shape = Umax.shape
    Umax = Umax.ravel()
    offset = offset.ravel()
    G = Umax.size

    w = np.ones(G); x = np.zeros(G); y = np.zeros(G); z = np.zeros(G)
    angle_factor = np.pi * T / inpoFact
    amplitude_factor = Umax * inpoFact / 100.0
    oz2 = offset * offset

    track = initialVector is not None
    if track:
        v = l * np.asarray(initialVector, dtype=np.float64)
        tangent_sum = np.tile(v, (G, 1))

    for b1, phase in PS:
        a = b1 * amplitude_factor
        nx = a * math.cos(math.radians(phase))
        ny = a * math.sin(math.radians(phase))
        norm = np.sqrt(nx * nx + ny * ny + oz2)
        half = angle_factor * norm
        nonzero = norm != 0
        s = np.where(nonzero, np.sin(half) / np.where(nonzero, norm, 1.0), 0.0)
        c = np.cos(half)
        bx, by, bz = nx * s, ny * s, offset * s
        w, x, y, z = (c * w - bx * x - by * y - bz * z,
                      c * x + bx * w + by * z - bz * y,
                      c * y + by * w + bz * x - bx * z,
                      c * z + bz * w + bx * y - by * x)
        if track:
            tangent_sum += inverse_rotate(np.stack([w, x, y, z], axis=-1), v)

    q = np.stack([w, x, y, z], axis=-1).reshape(shape + (4,))
    if track:
        tangent_sum -= 0.5 * (v + inverse_rotate(q.reshape(-1, 4), v))
        return q, tangent_sum.reshape(shape + (3,))
    return q
