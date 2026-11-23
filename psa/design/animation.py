"""Spin dynamics during a pulse for the animation tab (rotating frame): an ensemble of spins (offset x B1
scaling) and one selected spin, at any time t of the pulse.

Conventions as psa.rotations: sub-pulse k rotates about n_k = (s A_k umax cos(phi_k), s A_k umax sin(phi_k),
offset) [Hz] by 2 pi dt |n_k| (forward, time ordered); m(t) = R(t) m0. Within a sub-pulse the axis is
constant, so a fraction a of it is the rotation by a times its angle: the motion is exact at every time,
not only at the sub-pulse boundaries.

Times are given as tau in [0, N] (sub-pulses elapsed, fractional).
Droplets (DROPS / SpinDrops style, one spin 1/2):
    propagator U = w 1 - i (x sx + y sy + z sz):  rank 0: f0 = w (sphere),  rank 1: f1(r) = -i (x, y, z) . r
    magnetization m:                               rank 1: f(r) = m . r (two lobes; the Q-bead colours a sphere)
Colours of complex values: phase 0 red, pi/2 yellow, pi green, 3 pi/2 blue (cyclic).
"""
from __future__ import annotations

import numpy as np

from psa.design.pulse import PulseSequence

PHASE_ANCHORS = np.array([[0.90, 0.15, 0.15], [0.95, 0.80, 0.15], [0.20, 0.75, 0.25], [0.20, 0.40, 0.95]])


def _qmul(a, b):
    """Quaternion products a * b (..., 4)."""
    w1, x1, y1, z1 = np.moveaxis(a, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(b, -1, 0)
    return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                     w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
                     w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2], axis=-1)


def rotate(q, v):
    """R(q) v for quaternions (..., 4) and vectors (..., 3) (or one vector (3,))."""
    q = np.asarray(q, float)
    w, u = q[..., :1], q[..., 1:]
    v = np.broadcast_to(np.asarray(v, float), u.shape)
    t = 2.0 * np.cross(u, v)
    return v + w * t + np.cross(u, t)


class Dynamics:
    """Precomputed rotations of an ensemble (every offset with every scaling) plus a selected spin (last)."""

    def __init__(self, ps: PulseSequence, offsets_hz, scales, m0=(0, 0, 1), selected=(0.0, 1.0)):
        self.ps = ps
        offs, scs = np.asarray(offsets_hz, float), np.asarray(scales, float)
        if offs.size == 0 or scs.size == 0 or np.any(scs <= 0):
            raise ValueError("the ensemble needs offsets and positive B1 values")
        grid_o, grid_s = np.meshgrid(offs, scs, indexing="ij")
        self.offsets = np.r_[grid_o.ravel(), float(selected[0])]          # (E + 1,), the selected spin last
        self.scales = np.r_[grid_s.ravel(), float(selected[1])]
        self.n_ensemble = self.offsets.size - 1
        m0 = np.asarray(m0, float)
        if np.linalg.norm(m0) == 0:
            raise ValueError("m0 must not be zero")
        self.m0 = m0 / np.linalg.norm(m0)
        a = np.asarray(ps.amplitude, float) * ps.umax / 100.0
        ph = np.radians(np.asarray(ps.phase, float))
        n = np.stack([a[None] * self.scales[:, None] * np.cos(ph)[None],
                      a[None] * self.scales[:, None] * np.sin(ph)[None],
                      np.broadcast_to(self.offsets[:, None], (self.offsets.size, ps.N))], axis=-1)   # (S, N, 3) Hz
        self.field = n
        norm = np.linalg.norm(n, axis=-1)
        self.half = np.pi * ps.dt * norm                                    # half rotation angle per sub-pulse
        self.axis = np.where(norm[..., None] > 0, n / np.where(norm > 0, norm, 1.0)[..., None], [0.0, 0.0, 1.0])
        steps = np.concatenate([np.cos(self.half)[..., None], self.axis * np.sin(self.half)[..., None]], axis=-1)
        q = np.empty((self.offsets.size, ps.N + 1, 4))
        q[:, 0] = (1.0, 0.0, 0.0, 0.0)
        for k in range(ps.N):
            q[:, k + 1] = _qmul(steps[:, k], q[:, k])
        self.q = q                                                          # (S, N + 1, 4)
        self.m = rotate(q, self.m0)                                         # (S, N + 1, 3)
        # error curve of the selected spin: tangents R^T z, trapezoid rule, in units of the duration
        qs = q[-1].copy(); qs[:, 1:] *= -1
        t = rotate(qs, [0.0, 0.0, 1.0])
        seg = (t[1:] + t[:-1]) / 2 * (ps.dt / ps.duration if ps.duration else 1.0)
        self.curve = np.vstack([np.zeros(3), np.cumsum(seg, axis=0)])        # (N + 1, 3), length 1
        self.tangents = t

    @property
    def N(self):
        return self.ps.N

    def quaternions(self, tau):
        """Quaternions (S, 4) of all spins at the (fractional) time tau in [0, N]."""
        tau = float(np.clip(tau, 0.0, self.N))
        k = min(int(np.floor(tau)), self.N - 1) if self.N else 0
        a = tau - k
        if self.N == 0 or a <= 0:
            return self.q[:, int(round(tau))]
        h = self.half[:, k] * a
        part = np.concatenate([np.cos(h)[:, None], self.axis[:, k] * np.sin(h)[:, None]], axis=-1)
        return _qmul(part, self.q[:, k])

    def magnetization(self, tau):
        return rotate(self.quaternions(tau), self.m0)                       # (S, 3)

    def trail(self, tau, length):
        """Magnetization paths (S, L, 3) over the last ``length`` sub-pulses up to tau (boundaries + tau)."""
        k1 = int(np.floor(np.clip(tau, 0, self.N)))
        k0 = max(0, k1 - int(length))
        path = self.m[:, k0:k1 + 1]
        return np.concatenate([path, self.magnetization(tau)[:, None]], axis=1)

    def sub_pulse(self, tau):
        """Index of the sub-pulse running at tau."""
        return int(min(max(np.floor(tau), 0), max(self.N - 1, 0)))

    def effective_field(self, tau, spin=-1):
        """Effective field (Hz vector) of a spin during the sub-pulse at tau (rotating frame)."""
        return self.field[spin, self.sub_pulse(tau)]

    def curve_until(self, tau):
        """Error-curve points of the selected spin up to tau (with the exact end point at tau)."""
        k = int(np.floor(np.clip(tau, 0, self.N)))
        pts = self.curve[:k + 1]
        if k < self.N:
            a = tau - k
            q = self.quaternions(tau)[-1].copy(); q[1:] *= -1
            t_tau = rotate(q, [0.0, 0.0, 1.0])
            frac = a * self.ps.dt / self.ps.duration
            pts = np.vstack([pts, self.curve[k] + (self.tangents[k] + t_tau) / 2 * frac])
        return pts


def unit_sphere(rows=24, cols=48):
    """Unit vectors on a (rows + 1) x cols grid (poles included) and triangle faces of the sphere mesh."""
    th = np.linspace(0, np.pi, rows + 1)
    ph = np.linspace(0, 2 * np.pi, cols, endpoint=False)
    T, P = np.meshgrid(th, ph, indexing="ij")
    dirs = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], axis=-1).reshape(-1, 3)
    faces = []
    for i in range(rows):
        for j in range(cols):
            a, b = i * cols + j, i * cols + (j + 1) % cols
            c, d = a + cols, b + cols
            faces += [[a, c, b], [b, c, d]]
    return dirs, np.array(faces, dtype=np.int32)


def phase_colors(values, alpha=1.0):
    """RGBA (n, 4) of complex values by their phase (0 red, pi/2 yellow, pi green, 3 pi/2 blue)."""
    phase = np.mod(np.angle(np.asarray(values, complex)), 2 * np.pi) / (np.pi / 2)
    i0 = np.floor(phase).astype(int) % 4
    f = (phase - np.floor(phase))[:, None]
    rgb = PHASE_ANCHORS[i0] * (1 - f) + PHASE_ANCHORS[(i0 + 1) % 4] * f
    return np.column_stack([rgb, np.full(len(rgb), alpha)])


def droplet(values, dirs):
    """Vertices of a droplet: every direction scaled by |f| (radius = |f(r)|)."""
    return dirs * np.abs(np.asarray(values))[:, None]


def propagator_droplets(q, dirs):
    """(rank-0 values, rank-1 values) of the propagator with quaternion q at the directions ``dirs``."""
    q = np.asarray(q, float)
    f0 = np.full(len(dirs), complex(q[0]))
    f1 = -1j * (dirs @ q[1:])
    return f0, f1


def magnetization_droplet(m, dirs):
    """Rank-1 values m . r of the magnetization (real: red positive, green negative)."""
    return (dirs @ np.asarray(m, float)).astype(complex)
