"""Segment-based curve builder ("turtle" geometry).

A design is a list of segments; each starts where the previous one ended, in the same direction
(tangent continuity), rotated about the tangent by its ``twist`` (= phase jump of the pulse).
All lengths are in µs: the error curve is unit speed, so its length is the pulse duration.

Segment kinds (on resonance, R_min = 1 / (2 pi umax) is the tightest radius at full amplitude):
    line    length L                          -> amplitude 0 for L
    arc     radius R >= R_min, angle          -> constant amplitude R_min/R, constant phase
    helix   radius rho, pitch h (per turn),   -> curvature rho/(rho²+c²) = Omega, torsion c/(rho²+c²) = phase
            angle (about the helix axis)         rate (c = h / 2 pi); duration = angle * sqrt(rho² + c²)
    bezier  control points in the local frame -> exact discrete inversion (``conversion``) at the sub-pulse
            (x' = turning direction of the         duration
            current phase, y' = t x x', t); P1
            lies on the tangent (G1 continuity)

Continuous model: every segment is a list of "pieces" with constant Hamiltonian in a frame rotating
about z (amplitude Omega, start phase phi0, phase rate phidot, duration T). Gate and robustness vectors
A1..A3 are computed exactly for these pieces (z rotations commute with the detuning), which keeps the
"close / cancel" solver smooth. ``to_pulse`` discretizes the design into sub-pulses of duration dt.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

from psa.design.conversion import _rotation_from_axis_angle
from psa.design import jets as J
from psa.design.jets import SX, SY, SZ, to_vector as _to_vector
from psa.design.pulse import PulseSequence

KINDS = ("line", "arc", "helix", "bezier")
PARAMS = {
    "line": (("length", 20.0),),
    "arc": (("radius", 20.0), ("angle", 180.0)),
    "helix": (("radius", 20.0), ("pitch", 20.0), ("angle", 360.0)),
    "bezier": (("p1", 15.0), ("p2x", 15.0), ("p2y", 0.0), ("p2z", 25.0), ("p3x", 30.0), ("p3y", 5.0), ("p3z", 30.0)),
}
UNITS = {"length": "µs", "radius": "µs", "pitch": "µs/turn", "angle": "°", "twist": "°", "p1": "µs",
         "p2x": "µs", "p2y": "µs", "p2z": "µs", "p3x": "µs", "p3y": "µs", "p3z": "µs"}
BOUNDS = {"length": (0.0, 1e6), "radius": (1e-3, 1e6), "pitch": (-1e6, 1e6), "angle": (0.0, 36000.0),
          "twist": (-1e5, 1e5), "p1": (0.0, 1e6)}


@dataclass
class Segment:
    kind: str
    params: Dict[str, float] = field(default_factory=dict)
    twist: float = 0.0                       # degrees, phase jump at the start of the segment
    free: List[str] = field(default_factory=list)   # parameter names the solver may change

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"unknown segment kind {self.kind}")
        full = dict(PARAMS[self.kind])
        full.update(self.params)
        self.params = {k: float(v) for k, v in full.items()}

    def get(self, name):
        return self.twist if name == "twist" else self.params[name]

    def set(self, name, value):
        if name == "twist":
            self.twist = float(value)
        else:
            self.params[name] = float(value)

    def to_dict(self):
        return {"kind": self.kind, "params": dict(self.params), "twist": self.twist, "free": list(self.free)}

    @classmethod
    def from_dict(cls, d):
        return cls(d["kind"], d.get("params", {}), d.get("twist", 0.0), list(d.get("free", [])))


@dataclass
class Piece:
    Omega: float          # rad/s
    phi0: float           # rad
    phidot: float         # rad/s
    T: float              # s
    segment: int


def _bloch_rz(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def piece_rotation(p: Piece) -> np.ndarray:
    """Forward Bloch rotation of a piece: Rz(phi0 + phidot T) R_axis Rz(-phi0), R_axis about (Omega, 0, -phidot)."""
    w = math.hypot(p.Omega, p.phidot)
    if w * p.T == 0:
        core = np.eye(3)
    else:
        core = _rotation_from_axis_angle((p.Omega / w, 0.0, -p.phidot / w), w * p.T)
    return _bloch_rz(p.phi0 + p.phidot * p.T) @ core @ _bloch_rz(-p.phi0)


class CurveDesign:
    def __init__(self, segments=None, umax=1e4, dt=0.5e-6, phase0=0.0):
        self.segments: List[Segment] = list(segments or [])
        self.umax, self.dt = float(umax), float(dt)
        self.phase0 = float(phase0)          # absolute start phase in degrees

    # ---------------------------------------------------------------- basics ---
    @property
    def r_min_us(self):
        return 1e6 / (2 * math.pi * self.umax)

    def copy(self):
        return copy.deepcopy(self)

    def to_dict(self):
        return {"type": "builder", "umax": self.umax, "dt": self.dt, "phase0": self.phase0,
                "segments": [s.to_dict() for s in self.segments]}

    @classmethod
    def from_dict(cls, d):
        return cls([Segment.from_dict(s) for s in d["segments"]], d.get("umax", 1e4), d.get("dt", 0.5e-6), d.get("phase0", 0.0))

    # ---------------------------------------------------------------- pieces ---
    def pieces(self) -> Tuple[List[Piece], List[str]]:
        """Continuous model of the design and a list of warnings (e.g. amplitude above the limit)."""
        pieces, warnings = [], []
        self.segment_frames = []             # (R, phase) at every segment start (after its twist)
        Omega_max = 2 * math.pi * self.umax
        phase = math.radians(self.phase0)
        R = np.eye(3)                        # accumulated forward rotation
        for i, seg in enumerate(self.segments):
            phase += math.radians(seg.twist)
            self.segment_frames.append((R.copy(), phase))
            p = seg.params
            new = []
            if seg.kind == "line":
                if p["length"] > 0:
                    new = [Piece(0.0, phase, 0.0, p["length"] * 1e-6, i)]
            elif seg.kind == "arc":
                r = max(p["radius"], 1e-9) * 1e-6
                T = r * math.radians(p["angle"])
                if T > 0:
                    new = [Piece(1.0 / r, phase, 0.0, T, i)]
            elif seg.kind == "helix":
                rho = max(p["radius"], 1e-9) * 1e-6
                c = p["pitch"] * 1e-6 / (2 * math.pi)
                denom = rho ** 2 + c ** 2
                T = math.radians(p["angle"]) * math.sqrt(denom)
                if T > 0:
                    new = [Piece(rho / denom, phase, c / denom, T, i)]
            elif seg.kind == "bezier":
                new, R_end = self._bezier_pieces(seg, i, R, phase)
            for q in new:
                if q.Omega > Omega_max * (1 + 1e-9):
                    warnings.append(f"segment {i + 1} ({seg.kind}): amplitude {q.Omega / Omega_max * 100:.1f} % > 100 %")
                    break
            for q in new:
                R = piece_rotation(q) @ R
            if new:
                last = new[-1]
                phase = last.phi0 + last.phidot * last.T
            pieces += new
        return pieces, warnings

    @staticmethod
    def local_axes(R, phase):
        """Rows x', y', t of the local frame in lab coordinates: x' is the direction in which an arc with
        the current phase turns the tangent, t the tangent."""
        B = np.array([[-math.sin(phase), math.cos(phase), 0.0], [-math.cos(phase), -math.sin(phase), 0.0], [0.0, 0.0, 1.0]])
        return B @ R

    def _bezier_pieces(self, seg, index, R, phase):
        """Exact inversion of a cubic Bezier given in the local frame (x', y', t) (see local_axes)."""
        p = seg.params
        P = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, p["p1"]], [p["p2x"], p["p2y"], p["p2z"]],
                      [p["p3x"], p["p3y"], p["p3z"]]]) * 1e-6
        s = np.linspace(0.0, 1.0, 4001)
        d1 = (3 * (1 - s) ** 2)[:, None] * (P[1] - P[0]) + (6 * (1 - s) * s)[:, None] * (P[2] - P[1]) + \
             (3 * s ** 2)[:, None] * (P[3] - P[2])
        speed = np.linalg.norm(d1, axis=1)
        arclen = np.r_[0.0, np.cumsum((speed[1:] + speed[:-1]) / 2 * np.diff(s))]
        L = arclen[-1]
        n = int(round(L / self.dt))
        if n < 1 or L <= 0:
            return [], R
        target = np.interp(np.linspace(0, L, n + 1), arclen, s)
        dd = (3 * (1 - target) ** 2)[:, None] * (P[1] - P[0]) + (6 * (1 - target) * target)[:, None] * (P[2] - P[1]) + \
             (3 * target ** 2)[:, None] * (P[3] - P[2])
        norms = np.linalg.norm(dd, axis=1)
        if np.any(norms == 0):
            dd[norms == 0] = np.array([0.0, 0.0, 1.0])
            norms = np.linalg.norm(dd, axis=1)
        local_t = dd / norms[:, None]
        local_t[0] = (0.0, 0.0, 1.0)                     # G1: starts along the current tangent
        lab_t = local_t @ self.local_axes(R, phase)      # local -> lab
        step = L / n                                     # unit speed: step duration
        pieces, prev = [], phase
        for k in range(n):
            u = R @ lab_t[k + 1]
            sin_t = math.hypot(u[0], u[1])
            theta = math.atan2(sin_t, max(-1.0, min(1.0, u[2])))
            if sin_t < 1e-12:
                phi = prev
            else:
                phi = math.atan2(-u[0], u[1])
                phi += 2 * math.pi * round((prev - phi) / (2 * math.pi))
            piece = Piece(theta / step, phi, 0.0, step, index)
            pieces.append(piece)
            R = piece_rotation(piece) @ R
            prev = phi
        return pieces, R

    # ------------------------------------------------------- continuous model ---
    def duration(self, pieces=None):
        pieces = self.pieces()[0] if pieces is None else pieces
        return sum(p.T for p in pieces)

    def jets(self, pieces=None, perturbation="detuning", order=3):
        """Exact Taylor jets of the propagator of the continuous model in the detuning delta, or in the
        relative amplitude error epsilon (Omega -> (1 + epsilon) Omega) for perturbation="amplitude".
        Piece = Rz(phi0 + phidot T) exp(-i T (Omega/2 sx - phidot/2 sz + delta/2 sz)) Rz(phi0)^dag."""
        pieces = self.pieces()[0] if pieces is None else pieces
        if not pieces:
            return J.identity_jet(order + 1)
        Omega = np.array([p.Omega for p in pieces]); phid = np.array([p.phidot for p in pieces])
        T = np.array([p.T for p in pieces]); phi0 = np.array([p.phi0 for p in pieces])
        v0 = np.column_stack([Omega * T, np.zeros_like(T), -phid * T])
        d = np.column_stack([Omega * T, 0 * T, 0 * T]) if perturbation == "amplitude" else np.column_stack([0 * T, 0 * T, T])
        jets = J.step_jets(v0, d, order)                                   # rotating frame
        a = phi0 + phid * T                                                # Rz(a) X Rz(phi0)^dag, diagonal factors
        left = np.stack([np.exp(-0.5j * a), np.exp(0.5j * a)], axis=-1)
        right = np.stack([np.exp(0.5j * phi0), np.exp(-0.5j * phi0)], axis=-1)
        jets = left[:, None, :, None] * jets * right[:, None, None, :]
        return J.product_jets(jets)

    def robustness(self, pieces=None):
        """Exact A1, A2, A3 (s^k) of the continuous model and |A_k| / T^k."""
        pieces = self.pieces()[0] if pieces is None else pieces
        X = self.jets(pieces)
        U0d = X[0].conj().T
        D1, D2, D3 = U0d @ X[1], U0d @ X[2], U0d @ X[3]
        M = (D1, D2 - D1 @ D1 / 2, D3 - (D1 @ D2 + D2 @ D1) / 2 + D1 @ D1 @ D1 / 3)
        A = [_to_vector(m, k) for k, m in enumerate(M, start=1)]
        T = self.duration(pieces)
        return {"A1": A[0], "A2": A[1], "A3": A[2],
                "A1_norm": np.linalg.norm(A[0]) / T if T else 0.0,
                "A2_norm": np.linalg.norm(A[1]) / T ** 2 if T else 0.0,
                "A3_norm": np.linalg.norm(A[2]) / T ** 3 if T else 0.0}

    def curve(self, points_per_piece_us=0.25, pieces=None):
        """Points of the continuous error curve (s) and the indices where segments start.
        Each piece is subdivided into sub-steps of at most ``points_per_piece_us`` µs."""
        pieces = self.pieces()[0] if pieces is None else pieces
        pts = [np.zeros(3)]
        R = np.eye(3)
        starts = []
        last_seg = None
        for p in pieces:
            if p.segment != last_seg:
                starts.append((p.segment, len(pts) - 1))
                last_seg = p.segment
            m = max(2, int(math.ceil(p.T / (points_per_piece_us * 1e-6))))
            h = p.T / m
            for j in range(m):
                sub = Piece(p.Omega, p.phi0 + p.phidot * j * h, p.phidot, h, p.segment)
                t0 = R.T @ np.array([0.0, 0.0, 1.0])
                R = piece_rotation(sub) @ R
                t1 = R.T @ np.array([0.0, 0.0, 1.0])
                # chord of the sub-step: exact for arcs/lines up to O(h^3); midpoint tangent weighting
                pts.append(pts[-1] + (t0 + t1) / 2 * h)
        return np.array(pts), starts

    def frame_at_segment_starts(self):
        """For every segment: start position (s) and local axes (rows x', y', t in lab coordinates)."""
        pieces, _ = self.pieces()
        frames = list(self.segment_frames)
        pts, starts = self.curve(pieces=pieces)
        start_index = dict(starts)
        out, last = [], np.zeros(3)
        for i, (R, phase) in enumerate(frames):
            pos = pts[start_index[i]] if i in start_index else last
            out.append((pos.copy(), self.local_axes(R, phase)))
            last = pos
        return out

    def metrics(self, b1=False):
        """Metrics dict of the continuous model in the format of psa.design.optimize."""
        pieces, _ = self.pieces()
        rb = self.robustness(pieces)
        out = {"A": {1: rb["A1"], 2: rb["A2"], 3: rb["A3"]}, "T": self.duration(pieces), "R": self.gate_matrix(pieces)}
        if b1:
            X = self.jets(pieces, "amplitude")
            out["B1"] = _to_vector(X[0].conj().T @ X[1], 1)
            out["nutation"] = float(sum(abs(p.Omega) * p.T for p in pieces))
        return out

    def gate_matrix(self, pieces=None):
        pieces = self.pieces()[0] if pieces is None else pieces
        R = np.eye(3)
        for p in pieces:
            R = piece_rotation(p) @ R
        return R

    # ------------------------------------------------------------- discretize ---
    def to_pulse(self, dt=None) -> PulseSequence:
        """Sub-pulses of duration dt: each piece gets round(T/dt) sub-pulses (at least one for Bezier
        steps), constant amplitude, phases sampled at the sub-pulse midpoints."""
        dt = self.dt if dt is None else dt
        pieces, _ = self.pieces()
        amp, phase = [], []
        Omega_max = 2 * math.pi * self.umax
        for p in pieces:
            n = int(round(p.T / dt))
            if self.segments[p.segment].kind == "bezier":
                n = max(n, 1)
            if n == 0:
                continue
            amp += [p.Omega * p.T / (n * dt) / Omega_max * 100.0] * n
            h = p.T / n
            phase += list(np.degrees(p.phi0 + p.phidot * h * (np.arange(n) + 0.5)))
        return PulseSequence(np.array(amp), np.array(phase), dt, self.umax)

    # ----------------------------------------------------------------- solver ---
    def free_parameters(self):
        return [(i, name) for i, s in enumerate(self.segments) for name in s.free]

    def solve(self, orders=(1,), target=None, b1=False, max_nfev=None, callback=None, settings=None):
        """Adjust the free parameters for the goals (see psa.design.optimize): cancel the detuning orders,
        optionally reach a target gate (axis, angle_deg) and first-order B1 robustness.
        ``callback(nfev, best_cost)`` returning True stops the solver (best parameters so far are kept).
        Returns (new design, result dict)."""
        from psa.design.optimize import Goals, solve_problem
        solution, info = solve_problem("builder", self.to_dict(), None, Goals(tuple(orders), target, b1),
                                       callback, max_nfev, settings)
        design = CurveDesign.from_dict(solution)
        rb = design.robustness()
        info["norms"] = {k: rb[f"A{k}_norm"] for k in (1, 2, 3)}
        info["values"] = {f"{i + 1}:{n}": design.segments[i].get(n) for i, n in design.free_parameters()}
        return design, info
