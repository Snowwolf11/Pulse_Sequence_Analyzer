"""Pulse sequence <-> error curve.

Forward (``ps_to_curve``): tangents t_k = R_(k)^T z at the sub-pulse boundaries, curve points by the
trapezoid rule (identical to the analyzer).

Inverse (``tangents_to_ps`` / ``curve_to_ps``): exact discrete inversion. With the accumulated rotation
R_(k) known, sub-pulse k+1 must satisfy R_(k+1)^T z = R_(k) t_(k+1) =: u. On resonance its rotation
axis is transverse, n = (cos phi, sin phi, 0), and its angle theta >= 0; then
    u = (-sin(theta) sin(phi), sin(theta) cos(phi), cos(theta)),
so theta = arccos(u_z), phi = atan2(-u_x, u_y), amplitude = theta / (2 pi dt umax) * 100 %.
The resulting pulse reproduces every target tangent exactly (up to round-off). Steps with theta ~ 0
(straight parts of the curve) get amplitude 0 and keep the previous phase; inflection points appear
as 180 degree phase jumps.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import CubicSpline

from psa.createCoordinates_Matrix import trapezoid_curve
from psa.rotations import tangent_vectors
from psa.design.pulse import PulseSequence

Z = np.array([0.0, 0.0, 1.0])


# ---------------------------------------------------------------- forward ---

def ps_to_tangents(ps: PulseSequence, initial=Z) -> np.ndarray:
    """Unit tangents t_0 ... t_N (N + 1 rows) at the sub-pulse boundaries."""
    return tangent_vectors(ps.as_array(), ps.dt, 1.0, ps.umax, 0.0, 1, np.asarray(initial, dtype=np.float64))


def ps_to_curve(ps: PulseSequence, initial=Z) -> np.ndarray:
    """Error-curve points r_0 ... r_N in seconds (trapezoid rule, unit speed)."""
    return trapezoid_curve(ps_to_tangents(ps, initial) * ps.dt)


# ---------------------------------------------------------------- inverse ---

def _rotation_from_axis_angle(axis, angle):
    x, y, z = axis
    c, s, C = math.cos(angle), math.sin(angle), 1.0 - math.cos(angle)
    return np.array([[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
                     [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
                     [z * x * C - y * s, z * y * C + x * s, c + z * z * C]])


def tangents_to_ps(tangents, dt=0.5e-6, umax=1e4, phase_start=0.0, zero_tol=1e-12) -> PulseSequence:
    """Exact pulse sequence (N = len(tangents) - 1 sub-pulses) whose error-curve tangents at the
    sub-pulse boundaries are ``tangents``. ``tangents[0]`` must be +z (see ``align_curve``).
    Phases are unwrapped (each phase is chosen closest to the previous one)."""
    t = np.asarray(tangents, dtype=np.float64)
    t = t / np.linalg.norm(t, axis=1, keepdims=True)
    if np.linalg.norm(t[0] - Z) > 1e-6:
        raise ValueError("the first tangent must be +z; align the curve first")
    N = len(t) - 1
    amp = np.zeros(N)
    phase = np.zeros(N)
    R = np.eye(3)                       # accumulated forward rotation R_(k)
    prev_phase = float(phase_start)
    for k in range(N):
        u = R @ t[k + 1]
        cos_theta = min(1.0, max(-1.0, u[2]))
        sin_theta = math.hypot(u[0], u[1])
        theta = math.atan2(sin_theta, cos_theta)
        if sin_theta < zero_tol:
            phi = prev_phase            # no transverse turn: phase undefined, keep previous
            if theta < 1e-6:
                theta = 0.0
        else:
            phi = math.degrees(math.atan2(-u[0], u[1]))
            phi += 360.0 * round((prev_phase - phi) / 360.0)
        amp[k] = theta / (2 * math.pi * dt * umax) * 100.0
        phase[k] = phi
        prev_phase = phi
        n = (math.cos(math.radians(phi)), math.sin(math.radians(phi)), 0.0)
        R = _rotation_from_axis_angle(n, theta) @ R
    return PulseSequence(amp, phase, dt, umax)


# --------------------------------------------------------- curve handling ---

def _clean(points):
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3 or len(pts) < 2:
        raise ValueError("points must be an (n, 3) array with n >= 2")
    step = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    scale = max(np.max(np.abs(pts)), 1e-300)
    keep = np.r_[True, step > 1e-13 * scale]
    return pts[keep]


class CurveSpline:
    """Cubic spline through curve points (chord-length parameter) with an arc-length table."""

    def __init__(self, points, oversample=None, max_fine=20000):
        pts = _clean(points)
        u = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
        self.spline = CubicSpline(u, pts, axis=0, bc_type="not-a-knot" if len(pts) > 3 else "natural")
        if oversample is None:          # dense input needs little oversampling: cap the table size
            oversample = int(min(20, max(2, math.ceil(max_fine / max(1, len(pts) - 1)))))
        self.u_fine = np.linspace(0.0, u[-1], (len(pts) - 1) * oversample + 1)
        d1 = self.spline(self.u_fine, 1)
        speed = np.linalg.norm(d1, axis=1)
        self.s_fine = np.r_[0.0, np.cumsum((speed[1:] + speed[:-1]) / 2 * np.diff(self.u_fine))]
        self.length = float(self.s_fine[-1])
        d2 = self.spline(self.u_fine, 2)
        self.curvature_fine = np.linalg.norm(np.cross(d1, d2), axis=1) / np.maximum(speed, 1e-300) ** 3

    def at_arclength(self, s):
        """Points and unit tangents at arc lengths ``s``."""
        u = np.interp(s, self.s_fine, self.u_fine)
        d1 = self.spline(u, 1)
        return self.spline(u), d1 / np.linalg.norm(d1, axis=1, keepdims=True)

    @property
    def max_curvature(self) -> float:
        return float(np.max(self.curvature_fine))


def rotation_to_z(t0) -> np.ndarray:
    """Smallest rotation Q with Q t0 = +z."""
    t0 = np.asarray(t0, dtype=np.float64) / np.linalg.norm(t0)
    axis = np.cross(t0, Z)
    s, c = np.linalg.norm(axis), float(t0 @ Z)
    if s < 1e-14:
        return np.eye(3) if c > 0 else _rotation_from_axis_angle((1.0, 0.0, 0.0), math.pi)
    return _rotation_from_axis_angle(axis / s, math.atan2(s, c))


def align_curve(points, initial_tangent=None):
    """Translate the curve to start at the origin and rotate it so that its initial tangent is +z.
    Returns (aligned points, rotation matrix Q)."""
    pts = _clean(points)
    if initial_tangent is None:
        initial_tangent = CurveSpline(pts).at_arclength(np.array([0.0]))[1][0]
    Q = rotation_to_z(initial_tangent)
    return (pts - pts[0]) @ Q.T, Q


def curve_to_ps(points, dt=0.5e-6, umax=1e4, duration=None, phase_start=0.0, align=True):
    """Pulse sequence whose error curve is ``points`` (any units; the shape is what matters).

    The curve is aligned (start at the origin, initial tangent +z), scaled so that its length equals
    the pulse duration, sampled at M + 1 equally spaced arc lengths and inverted exactly.
    ``duration`` (s) defaults to the shortest multiple of ``dt`` allowed by ``umax`` and the maximal
    curvature. Returns (PulseSequence, info dict)."""
    pts = _clean(points)
    spline = CurveSpline(pts)                      # one spline; alignment is a rigid motion applied afterwards
    if align:
        Q = rotation_to_z(spline.at_arclength(np.array([0.0]))[1][0])
    else:
        Q = np.eye(3)
    origin = pts[0]
    L, kmax = spline.length, spline.max_curvature
    def invert(duration):
        M = max(1, int(round(duration / dt)))
        s = np.linspace(0.0, L, M + 1)
        target_points, tangents = spline.at_arclength(s)
        target_points = (target_points - origin) @ Q.T
        tangents = tangents @ Q.T
        if align:
            tangents[0] = Z
        return M, target_points, tangents_to_ps(tangents, dt, umax, phase_start)

    # estimate from the spline curvature (overestimates at curvature jumps), then rescale by the
    # actually needed peak amplitude of the discrete pulse
    min_duration = L * kmax / (2 * math.pi * umax)
    if duration is None:
        for _ in range(3):
            duration = max(dt, math.ceil(min_duration / dt * (1 + 1e-9)) * dt)
            M, target_points, ps = invert(duration)
            peak = float(np.max(ps.amplitude)) if ps.N else 0.0
            if peak <= 100.0 + 1e-9 and peak > 99.0 or peak == 0.0:
                break
            min_duration = duration * peak / 100.0
        if ps.N and np.max(ps.amplitude) > 100.0 + 1e-9:   # rounding: one more step
            duration += dt
            M, target_points, ps = invert(duration)
        min_duration = duration * (float(np.max(ps.amplitude)) / 100.0 if ps.N else 0.0)
    else:
        M, target_points, ps = invert(duration)
    duration = M * dt
    scale = duration / L                     # seconds per curve unit
    info = {
        "duration": duration,
        "M": M,
        "length": L,
        "scale": scale,
        "min_duration": min_duration,
        "max_curvature": kmax,
        "max_amplitude_percent": float(np.max(ps.amplitude)) if ps.N else 0.0,
        "amplitude_ok": bool(ps.N == 0 or np.max(ps.amplitude) <= 100.0 + 1e-6),
        "target_points": target_points * scale,   # seconds, aligned
        "rotation": Q,
    }
    return ps, info
