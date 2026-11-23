"""GRAPE-like optimization of all sub-pulses with exact gradients.

Controls: Cartesian u_x, u_y of every sub-pulse in units of the maximum amplitude, x = (u_x[0..N-1], u_y[0..N-1]).
The amplitude limit |u| < 1 is built in: the solver works on unconstrained v with u = v tanh(|v|) / |v|
(direction kept, amplitude squashed smoothly into [0, 1)). Box bounds or penalties stall the trust region
solver at the limit, which is where good pulses usually sit.
Optionally the controls are smooth: cubic B-splines with M control points per channel (u = squash(B c)).

Smoothness (optional, all soft so that gradients stay useful): a Cartesian smoothness term
w N sum |u_(k+1) - u_k|^2, and one-sided limits of the amplitude rate (%/µs) and of the phase rate (kHz,
weighted by the amplitudes so that phase jumps at zero amplitude are free) that cost nothing below the
limit (squared hinge). Soft terms (smoothness, stay close) are lowered stage by stage (x0.1) until the goals
are reached exactly; each stage starts from the previous, smooth result.
Coarse to fine (optional): first stages on few spline points (large-scale shape only), then the chosen controls.

Exploration (Search): a population of starts (perturbations of the start pulse, smooth or white, or random
pulses) refined by GRAPE in rounds; successive halving keeps the best fraction and doubles the budget each
round; optional crossover children of survivors (blends in the control variables, experimental); then kicks
of the best result (perturb, re-optimize, keep if better, or with a temperature sometimes if worse).

Robustness: goals at a set of points (offset, amplitude scaling, weight). The gate goal is required at
every point ("order 0"); local mode adds the Taylor orders of the detuning error (A1..A3, SCQC) and
optionally the first-order B1 error at every point.
    local      gate + orders (+ B1) at the points; default points: the offsets of a list at scaling 1
    broadband  classical GRAPE: only the gate at the points (mean infidelity); default points: a grid of
               offsets x amplitude scalings
An explicit point list (``points``) replaces the default points in both modes; its weights (normalized to
mean 1) scale the rows of every point by sqrt(weight) and multiply the Lawson weights.
Gate goal: none, UR (rotation angle about an axis) or PP (initial magnetization m0 -> target).
For PP the local robustness terms are the components perpendicular to m0 (A_k x m0): an error rotation
about m0 does not change the transfer.

Exact gradients: every sub-pulse jet J_k and its derivative with respect to u_x[k], u_y[k] in closed form
(psa.design.jets); dX/du_k = S_(k+1) o dJ_k o P_(k-1) with prefix / suffix products; then the chain rule
through the Magnus terms, the gate and the final magnetization. Cost O(N) per gradient.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import List

import numpy as np
from scipy import sparse

from psa.design import jets as Jt
from psa.design.optimize import SolverSettings, rotation_angle_deg, rotation_from_axis_angle
from psa.design.pulse import PulseSequence

SIG = (Jt.SX, Jt.SY, Jt.SZ)


@dataclass
class GrapeConfig:
    mode: str = "local"                          # "local" or "broadband"
    offsets_hz: List[float] = field(default_factory=lambda: [0.0])
    orders: List[int] = field(default_factory=lambda: [1])
    b1: bool = False
    bb_offset_min_hz: float = -5000.0
    bb_offset_max_hz: float = 5000.0
    bb_offset_n: int = 11
    bb_scale_min: float = 1.0                    # amplitude scaling (1 = nominal)
    bb_scale_max: float = 1.0
    bb_scale_n: int = 1
    gate: str = "UR"                             # "none", "UR", "PP"
    target_axis: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0])
    target_angle_deg: float = 180.0
    m0: List[float] = field(default_factory=lambda: [0.0, 0.0, 1.0])
    m_target: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0])
    regularization: str = "none"                 # "none" or "stay_close"
    regularization_weight: float = 1e-3
    start_noise: float = 1e-2                    # tiny seeded perturbation of the start (breaks symmetries)
    seed: int = 0
    controls: str = "full"                       # "full" (every sub-pulse) or "spline" (smooth cubic B-splines)
    n_spline: int = 20                           # spline control points per channel (u_x, u_y)
    smooth_weight: float = 0.0                   # soft smoothness term  w N sum |u_(k+1) - u_k|^2  (0 = off)
    max_amp_rate: float = 0.0                    # soft limit of |d amplitude / dt| in %/µs (0 = off)
    max_phase_rate_khz: float = 0.0              # soft limit of |d phase / dt| / 2pi in kHz (0 = off)
    limit_weight: float = 1.0                    # weight of the rate-limit rows
    limits_late: bool = False                    # rate limits only in a last stage (after the free optimization)
    n_starts: int = 1                            # multi-start: number of starts
    multistart_noise: float = 0.5                # perturbation of starts 2, 3, ... (in the squashed variables)
    start_pool: str = "perturb"                  # starts 2, 3, ...: "perturb" the start pulse, "random" pulses, "mix"
    noise_correlation_us: float = 0.0            # correlation time of start / kick noise (0 = white noise)
    screen_rounds: int = 0                       # successive-halving rounds (0 = every start gets the full budget)
    screen_nfev: int = 40                        # evaluations per candidate in round 1 (then / keep_fraction per round)
    keep_fraction: float = 0.5                   # fraction of the candidates kept after each round
    crossover_rate: float = 0.0                  # dropped candidates replaced by crossover children (fraction, experimental)
    n_kicks: int = 0                             # kicks of the best result (perturb + re-optimize)
    kick_size: float = 0.3                       # kick perturbation (same units as the spread)
    kick_temperature: float = 0.0                # accept a worse kick with exp(-Δlog10(residual) / T) (0 = never)
    coarse_levels: List[int] = field(default_factory=list)   # spline points of coarse stages before the controls
    reweight_rounds: int = 0                     # Lawson rounds: raise the weights of the worst points (0 = off)
    reweight_power: float = 0.5                  # weight update w <- w (e / mean e)^power (1: classic Lawson, may oscillate)
    points: List[list] = field(default_factory=list)  # explicit points [offset_hz, scale, weight, goal] (empty: defaults)
    freedom: str = "none"                        # gate goals up to a z-rotation after the pulse: "none", "phase"
                                                 # (one free angle) or "linear" (angle + offset x delay: linear phase)
    max_delay_us: float = 0.0                    # "linear": |delay| at most this (0 = the pulse duration)
    start_kind: str = "perturb"                  # set by variant(): "perturb" or "random"

    def to_dict(self):
        return {k: (list(v) if isinstance(v, (list, tuple)) else v) for k, v in self.__dict__.items()}

    @classmethod
    def from_dict(cls, d):
        names = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in (d or {}).items() if k in names})

    def check(self):
        if self.mode not in ("local", "broadband"):
            raise ValueError("mode must be 'local' or 'broadband'")
        if self.gate not in ("none", "UR", "PP"):
            raise ValueError("gate must be none, UR or PP")
        if self.points:
            for pt in self.points:
                nums = list(pt[:3])
                if not 1 <= len(pt) <= 4 or not all(np.isfinite(float(v)) for v in nums):
                    raise ValueError("a point is [offset_hz, scale, weight, goal] (scale, weight and goal optional)")
                if len(pt) > 1 and float(pt[1]) <= 0:
                    raise ValueError("the amplitude scaling of a point must be positive")
                if len(pt) > 2 and float(pt[2]) <= 0:
                    raise ValueError("the weight of a point must be positive")
                if len(pt) > 3 and pt[3] is not None:
                    check_goal(pt[3])
        if self.freedom not in FREEDOMS:
            raise ValueError(f"freedom must be one of {FREEDOMS}")
        if self.max_delay_us < 0:
            raise ValueError("the maximal delay must not be negative")
        goals = self.point_goals()
        every_point_has_gate = self.gate != "none" or (goals and all(g is not None for g in goals))
        if self.mode == "local":
            if not self.offsets_hz and not self.points:
                raise ValueError("local robustness needs at least one offset")
            if not self.orders and not self.b1 and self.gate == "none" and not any(g is not None for g in goals):
                raise ValueError("choose at least one goal (an order, B1 or a gate)")
        else:
            if not every_point_has_gate:
                raise ValueError("broadband optimization needs a UR or PP gate goal (for every point)")
            if not self.points and (self.bb_offset_n < 1 or self.bb_scale_n < 1):
                raise ValueError("the broadband grids need at least one point")
        if self.start_noise < 0 or self.multistart_noise < 0:
            raise ValueError("the start noise must not be negative")
        if self.controls not in ("full", "spline"):
            raise ValueError("controls must be 'full' or 'spline'")
        if self.controls == "spline" and int(self.n_spline) < 4:
            raise ValueError("the spline needs at least 4 control points")
        if min(self.smooth_weight, self.max_amp_rate, self.max_phase_rate_khz, self.limit_weight,
               self.regularization_weight) < 0:
            raise ValueError("weights and rate limits must not be negative")
        if int(self.n_starts) < 1:
            raise ValueError("at least one start is needed")
        if self.start_pool not in ("perturb", "random", "mix"):
            raise ValueError("the start pool must be 'perturb', 'random' or 'mix'")
        if min(self.noise_correlation_us, self.kick_size, self.kick_temperature, self.crossover_rate) < 0:
            raise ValueError("noise correlation, kick size, temperature and crossover rate must not be negative")
        if int(self.screen_rounds) < 0 or int(self.n_kicks) < 0:
            raise ValueError("rounds and kicks must not be negative")
        if int(self.screen_nfev) < 2:
            raise ValueError("the first screening round needs at least 2 evaluations")
        if not 0 < self.keep_fraction < 1:
            raise ValueError("the kept fraction must be between 0 and 1")
        if self.crossover_rate > 1:
            raise ValueError("the crossover rate must be at most 1")
        if int(self.reweight_rounds) < 0 or self.reweight_power <= 0:
            raise ValueError("reweighting needs rounds >= 0 and a positive power")
        if any(int(m) < 4 for m in self.coarse_levels):
            raise ValueError("coarse levels need at least 4 spline points")
        if self.gate == "PP" and (np.linalg.norm(self.m0) == 0 or np.linalg.norm(self.m_target) == 0):
            raise ValueError("PP vectors must not be zero")
        if self.gate == "UR" and np.linalg.norm(self.target_axis) == 0:
            raise ValueError("the target rotation axis must not be zero")

    def variant(self, i):
        """Configuration of start i of a multi-start (other seed; larger perturbation for i > 0)."""
        if i == 0:
            return self
        d = self.to_dict()
        rand = self.start_pool == "random" or (self.start_pool == "mix" and i % 2 == 1)
        d.update(seed=int(self.seed) + i, start_noise=self.multistart_noise, start_kind="random" if rand else "perturb")
        return GrapeConfig.from_dict(d)

    def replace(self, **kw):
        d = self.to_dict()
        d.update(kw)
        return GrapeConfig.from_dict(d)

    def n_points(self):
        """Number of points that reweighting balances."""
        return len(self.point_list())

    def point_list(self):
        """(offset_hz, scale, weight) of every point: the explicit points, else the default points of the mode
        (local: the offsets at scaling 1; broadband: the offset x scaling grid), weights 1."""
        if self.points:
            return [(float(p[0]), float(p[1]) if len(p) > 1 else 1.0, float(p[2]) if len(p) > 2 else 1.0)
                    for p in self.points]
        if self.mode == "local":
            return [(float(o), 1.0, 1.0) for o in self.offsets_hz]
        offs = np.linspace(self.bb_offset_min_hz, self.bb_offset_max_hz, int(self.bb_offset_n))
        scales = np.linspace(self.bb_scale_min, self.bb_scale_max, int(self.bb_scale_n))
        return [(float(o), float(s), 1.0) for o in offs for s in scales]

    def point_goals(self):
        """Own goal of every point (a goal dict, see check_goal) or None (= the gate goal of the config)."""
        if self.points:
            return [p[3] if len(p) > 3 else None for p in self.points]
        return [None] * len(self.point_list())

    def point_targets(self):
        """Effective goal of every point: its own goal, else the config's gate goal (None if "none")."""
        default = None
        if self.gate == "UR":
            default = {"gate": "UR", "axis": list(self.target_axis), "angle": float(self.target_angle_deg)}
        elif self.gate == "PP":
            default = {"gate": "PP", "m0": list(self.m0), "mt": list(self.m_target)}
        return [g if g is not None else default for g in self.point_goals()]

    def ensemble(self):
        """(offset_hz, scale) of every point."""
        return [(o, s) for o, s, _ in self.point_list()]

    def user_weights(self):
        """Point weights normalized to mean 1, or None if all are equal."""
        w = np.array([p[2] for p in self.point_list()], float)
        return None if np.all(w == w[0]) else w / np.mean(w)


FREEDOMS = ("none", "phase", "linear")


def rz(theta):
    """Rotations about z by the angles theta: (..., 3, 3)."""
    c, s_ = np.cos(theta), np.sin(theta)
    out = np.zeros(np.shape(theta) + (3, 3))
    out[..., 0, 0], out[..., 0, 1], out[..., 1, 0], out[..., 1, 1], out[..., 2, 2] = c, -s_, s_, c, 1.0
    return out


def _overlaps(R, target):
    """(a, b) with  overlap(theta) = a cos(theta) + b sin(theta) + const  for the target rotated by Rz(theta)
    after the pulse (UR: tr((Rz Rt)^T R); PP: (R m0) . (Rz mt)). R: (..., 3, 3)."""
    kind, Rt, m0, mt = target
    if kind == "UR":
        M = R @ Rt.T
        return M[..., 0, 0] + M[..., 1, 1], M[..., 1, 0] - M[..., 0, 1]
    m = R @ m0
    return m[..., 0] * mt[0] + m[..., 1] * mt[1], m[..., 1] * mt[0] - m[..., 0] * mt[1]


def fit_freedom(a, b, offsets_hz, weights, freedom, tau_max):
    """Best free z-rotation after the pulse, theta_i = phi + 2 pi offset_i tau (tau = 0 for "phase"): maximizes
    sum w_i (a_i cos theta_i + b_i sin theta_i), i.e. minimizes the weighted squared gate residual. phi in closed
    form (atan2), tau by a grid search refined by golden section within |tau| <= tau_max. Returns (phi, tau)."""
    a, b, f, w = (np.asarray(v, float) for v in (a, b, offsets_hz, weights))

    def best_phi(tau):
        d = 2 * math.pi * f * tau
        A = np.sum(w * (a * np.cos(d) + b * np.sin(d)))
        B = np.sum(w * (b * np.cos(d) - a * np.sin(d)))
        return math.atan2(B, A), math.hypot(A, B)

    if freedom != "linear" or tau_max <= 0 or np.ptp(f) == 0:
        return best_phi(0.0)[0], 0.0
    taus = np.linspace(-tau_max, tau_max, 401)
    d = 2 * math.pi * np.outer(taus, f)
    A = (np.cos(d) * a + np.sin(d) * b) @ w
    B = (np.cos(d) * b - np.sin(d) * a) @ w
    k = int(np.argmax(np.hypot(A, B)))
    lo, hi = taus[max(k - 1, 0)], taus[min(k + 1, len(taus) - 1)]
    g = (math.sqrt(5) - 1) / 2
    x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo)
    f1, f2 = best_phi(x1)[1], best_phi(x2)[1]
    for _ in range(40):
        if f1 > f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - g * (hi - lo); f1 = best_phi(x1)[1]
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + g * (hi - lo); f2 = best_phi(x2)[1]
    tau = 0.5 * (lo + hi)
    return best_phi(tau)[0], tau


def check_goal(goal):
    """A per-point goal: {"gate": "UR", "axis": [x, y, z], "angle": degrees} or {"gate": "PP", "m0": [...],
    "mt": [...]}. Raises ValueError if invalid."""
    if not isinstance(goal, dict) or goal.get("gate") not in ("UR", "PP"):
        raise ValueError("a point goal is {'gate': 'UR', 'axis', 'angle'} or {'gate': 'PP', 'm0', 'mt'}")
    vecs = ("axis",) if goal["gate"] == "UR" else ("m0", "mt")
    for k in vecs:
        v = np.asarray(goal.get(k, ()), float)
        if v.shape != (3,) or not np.all(np.isfinite(v)) or np.linalg.norm(v) == 0:
            raise ValueError(f"point goal: '{k}' must be a non-zero 3-vector")
    if goal["gate"] == "UR" and not np.isfinite(float(goal.get("angle", np.nan))):
        raise ValueError("point goal: the UR angle must be a number")


def goal_target(goal):
    """("UR", Rt, None, None) or ("PP", None, m0, mt) of a goal dict; None for None."""
    if goal is None:
        return None
    if goal["gate"] == "UR":
        return ("UR", rotation_from_axis_angle(_unit(goal["axis"]), float(goal["angle"])), None, None)
    return ("PP", None, _unit(goal["m0"]), _unit(goal["mt"]))


def goal_key(goal):
    """Hashable key of a goal (to group equal goals)."""
    if goal is None:
        return None
    if goal["gate"] == "UR":
        return ("UR", tuple(np.round(_unit(goal["axis"]), 12)), round(float(goal["angle"]), 9))
    return ("PP", tuple(np.round(_unit(goal["m0"]), 12)), tuple(np.round(_unit(goal["mt"]), 12)))


# ------------------------------------------------------------------ engine ---

def propagator_and_gradient(ux, uy, dt, umax, offset_hz=0.0, scale=1.0, order=3, perturbation="detuning"):
    """Jets X (order+1, 2, 2) of U(T) expanded around ``offset_hz`` in the detuning (or in the relative
    amplitude error for perturbation="amplitude"), and their derivatives dX (2N, order+1, 2, 2) with
    respect to u_x[0..N-1] then u_y[0..N-1]."""
    N = len(ux)
    W = 2 * math.pi * umax * scale * dt
    v0 = np.column_stack([W * ux, W * uy, np.full(N, 2 * math.pi * offset_hz * dt)])
    ex = np.array([W, 0.0, 0.0]); ey = np.array([0.0, W, 0.0])
    if perturbation == "amplitude":
        d = np.column_stack([W * ux, W * uy, np.zeros(N)])
        fx, fy = ex, ey
    else:
        d = np.array([0.0, 0.0, dt]); fx = fy = None
    J, dJx = Jt.step_jets_with_derivative(v0, d, ex, fx, order)
    _, dJy = Jt.step_jets_with_derivative(v0, d, ey, fy, order)
    P = Jt.prefix_products(J)
    S = Jt.suffix_products(J)
    ident = Jt.identity_jet(order + 1)[None]
    Pm1 = np.concatenate([ident, P[:-1]], axis=0)       # J[k-1] o ... o J[0]
    Sp1 = np.concatenate([S[1:], ident], axis=0)        # J[N-1] o ... o J[k+1]
    dX = np.concatenate([Jt.compose(Sp1, Jt.compose(dJx, Pm1)), Jt.compose(Sp1, Jt.compose(dJy, Pm1))], axis=0)
    return P[-1], dX


def _vec(M, k):
    """A (..., 3) with M = -i A . sigma / 2^k, for M (..., 2, 2)."""
    return np.real(np.stack([1j * 2 ** (k - 1) * np.einsum("...ij,ji->...", M, s) for s in SIG], axis=-1))


def magnus_with_gradient(X, dX):
    """A1..A3 and their derivatives (M, 3) from jets X and derivative jets dX."""
    K = len(X) - 1                                   # available orders
    U0d = X[0].conj().T
    dU0d = np.conj(np.swapaxes(dX[:, 0], -1, -2))
    D = [None] + [U0d @ X[k] for k in range(1, K + 1)]
    dD = [None] + [dU0d @ X[k] + U0d @ dX[:, k] for k in range(1, K + 1)]
    M, dM = [], []
    if K >= 1:
        M.append(D[1]); dM.append(dD[1])
    if K >= 2:
        D1, D2, d1, d2 = D[1], D[2], dD[1], dD[2]
        M.append(D2 - D1 @ D1 / 2); dM.append(d2 - (d1 @ D1 + D1 @ d1) / 2)
    if K >= 3:
        D3, d3 = D[3], dD[3]
        M.append(D3 - (D1 @ D2 + D2 @ D1) / 2 + D1 @ D1 @ D1 / 3)
        dM.append(d3 - (d1 @ D2 + D1 @ d2 + d2 @ D1 + D2 @ d1) / 2 + (d1 @ D1 @ D1 + D1 @ d1 @ D1 + D1 @ D1 @ d1) / 3)
    return [(_vec(M[k], k + 1), _vec(dM[k], k + 1)) for k in range(K)]


def rotation_with_gradient(U, dU):
    """Bloch rotation R_ij = tr(s_i U s_j U^dag) / 2 and dR (M, 3, 3)."""
    Ud = U.conj().T
    R = np.real(np.array([[np.trace(si @ U @ sj @ Ud) / 2 for sj in SIG] for si in SIG]))
    dR = np.real(np.stack([np.stack([np.einsum("ij,mjk,kl,li->m", si, dU, sj, Ud) for sj in SIG], -1) for si in SIG], -2))
    return R, dR


LIMIT_SHAPE = "smooth"


def _limit_shape(h):
    """Residual of a relative limit excess h >= 0 and its derivative: h^2 / (1 + h) is C1 at the limit
    (a kink there spoils the Gauss-Newton model) and grows linearly far beyond it."""
    if LIMIT_SHAPE == "linear":
        return h, np.ones_like(h)
    return h * h / (1 + h), h * (2 + h) / (1 + h) ** 2


def _angle_between_deg(a, b):
    return math.degrees(math.atan2(float(np.linalg.norm(np.cross(a, b))), float(np.dot(a, b))))


def _unit(v):
    v = np.asarray(v, dtype=np.float64)
    return v / np.linalg.norm(v)


# ------------------------------------------------------------- smooth basis ---

def spline_basis(N, M):
    """(N, M) matrix of M clamped cubic B-splines on [0, 1], sampled at the sub-pulse centres."""
    from scipy.interpolate import BSpline
    if M < 4:
        raise ValueError("the spline needs at least 4 control points")
    if M > N:
        raise ValueError(f"more spline control points ({M}) than sub-pulses ({N})")
    t = np.r_[np.zeros(4), np.linspace(0, 1, M - 2)[1:-1], np.ones(4)]
    x = (np.arange(N) + 0.5) / N
    B = BSpline(t, np.eye(M), 3, extrapolate=False)(x)
    return np.nan_to_num(B)


def control_noise(N, dt, std, correlation_us, rng):
    """Random controls (2N: x then y) with standard deviation ``std``: white, or smooth with a Gaussian
    correlation time ``correlation_us`` (then the starts differ in their large-scale shape, not by spikes)."""
    z = rng.normal(0.0, 1.0, (2, N))
    sigma = correlation_us * 1e-6 / dt
    if sigma >= 0.5 and N > 1:
        from scipy.ndimage import gaussian_filter1d
        z = gaussian_filter1d(z, sigma, axis=1, mode="reflect")
        z /= np.sqrt(np.mean(z * z, axis=1, keepdims=True)) + 1e-300
    return std * np.r_[z[0], z[1]]


def has_soft(cfg):
    return (cfg.regularization == "stay_close" and cfg.regularization_weight > 0) or cfg.smooth_weight > 0


def stage_plan(cfg, N):
    """Stages (controls, spline points, soft-term factor, kind): "coarse" spline levels (soft terms at full
    weight), then the chosen controls with the soft terms lowered stage by stage (SOFT_SCHEDULE, kind "main"),
    then "reweight" stages (Lawson: before each, the weights of the points grow with their errors)."""
    soft = SOFT_SCHEDULE if has_soft(cfg) else (1.0,)
    final_m = int(cfg.n_spline) if cfg.controls == "spline" else N
    levels = sorted({int(m) for m in cfg.coarse_levels if 4 <= int(m) < final_m and int(m) <= N})
    plan = [("spline", m, soft[0], "coarse") for m in levels] + [(cfg.controls, int(cfg.n_spline), f, "main") for f in soft]
    if cfg.n_points() > 1:
        plan += [(cfg.controls, int(cfg.n_spline), soft[-1], "reweight")] * int(cfg.reweight_rounds)
    if cfg.limits_late and (cfg.max_amp_rate > 0 or cfg.max_phase_rate_khz > 0):
        # limits switched on late: the free optimization first (limits from the start can block the way to
        # good solutions, which often pass through rough intermediate pulses), then a stage with the limits
        i = next(k for k, st in enumerate(plan) if st[3] == "reweight") if plan[-1][3] == "reweight" else len(plan)
        plan.insert(i, (cfg.controls, int(cfg.n_spline), soft[-1], "limits"))
    return plan


def lawson_update(weights, errors, power=1.0):
    """Lawson step towards the minimax solution: w <- w (e / mean e)^power, normalized to mean 1 and kept
    within [1e-4, 1e4]."""
    w = np.asarray(weights, float)
    e = np.maximum(np.asarray(errors, float), 1e-300)
    w = w * (e / np.mean(e)) ** power
    w = np.clip(w / np.mean(w), 1e-4, 1e4)
    return w / np.mean(w)


# ----------------------------------------------------------------- problem ---

INFO_KEYS = ("mean infidelity", "gate error [°]", "free phase [°]", "free delay [µs]", "max amplitude [%]", "max amplitude rate [%/µs]", "max phase rate [kHz]", "roughness",
             "point error spread")
SOFT_SCHEDULE = (1.0, 0.1, 0.01, 1e-3, 1e-4, 0.0)      # weights of the soft terms in successive stages
LIMIT_EXCESS_TOL = 1e-3
SOFT_STAGE_TOL = 1e-6                                   # xtol / ftol / gtol of the soft and coarse stages


class GrapeProblem:
    """Residuals (goal blocks + limits + soft terms) of the controls u and their exact Jacobian, plus the map
    from the optimization parameters p to u (amplitude squash, optionally a spline basis)."""

    def __init__(self, ps: PulseSequence, config: GrapeConfig, settings: SolverSettings = None):
        config.check()
        self.cfg, self.st = config, settings or SolverSettings()
        self.dt, self.umax, self.N = ps.dt, ps.umax, ps.N
        if self.N == 0:
            raise ValueError("the start pulse is empty")
        a = np.minimum(ps.amplitude / 100.0, 1.0)
        ph = np.radians(ps.phase)
        self.x0 = np.r_[a * np.cos(ph), a * np.sin(ph)]
        self.T = ps.duration
        self.m0, self.mt = _unit(config.m0), _unit(config.m_target)
        self.Rt = rotation_from_axis_angle(config.target_axis, config.target_angle_deg) if config.gate == "UR" else None
        self.targets = [goal_target(g) for g in config.point_targets()]      # per point (None: no gate goal)
        self.freedom_fit = (0.0, 0.0)        # fitted free z-rotation after the pulse (phi [rad], tau [s])
        self.soft_scale = 1.0
        self.limit_scale = 1.0               # 0 switches the rate limits off (stages before "limits")
        self.point_weights = None            # Lawson weights of the points (offsets / grid points), None = all 1
        self.point_errors = None             # goal residual of every point, user-weighted (set by goals)
        self.user_w = config.user_weights()  # user weights of the points (mean 1), None = all equal
        D = sparse.diags([-np.ones(self.N - 1), np.ones(self.N - 1)], [0, 1], shape=(self.N - 1, self.N))
        self.D2 = sparse.block_diag([D, D], format="csr")             # first differences of u_x and u_y
        w0 = self.unsquash(self.x0)
        if config.controls == "spline":
            B = spline_basis(self.N, int(config.n_spline))
            self.basis = sparse.block_diag([sparse.csr_matrix(B), sparse.csr_matrix(B)], format="csr")
            self._pinv = np.linalg.pinv(B)
        else:
            self.basis = None
        self.p0 = self.params(w0)

    def controls(self, x):
        return x[:self.N], x[self.N:]

    def pulse(self, x) -> PulseSequence:
        ux, uy = self.controls(x)
        amp = np.hypot(ux, uy)
        return PulseSequence(np.minimum(amp, 1.0) * 100.0, np.degrees(np.unwrap(np.arctan2(uy, ux))), self.dt, self.umax)

    def has_soft_terms(self):
        return has_soft(self.cfg)

    # ------------------------------------------------- parameters -> controls ---
    def wvec(self, p):
        """Unsquashed controls w (2N) of the parameters p."""
        return p if self.basis is None else self.basis @ p

    def params(self, w):
        """Parameters of the unsquashed controls w (least-squares fit for spline controls)."""
        w = np.asarray(w, dtype=np.float64)
        if self.basis is None:
            return w.copy()
        return np.r_[self._pinv @ w[:self.N], self._pinv @ w[self.N:]]

    def to_controls(self, p):
        """Controls u(p) and the sparse Jacobian du/dp."""
        w = p if self.basis is None else self.basis @ p
        u, M = self.squash(w)
        return u, (M if self.basis is None else M @ self.basis)

    def squash(self, v):
        """Controls x = u(v) with |u| < 1 and the sparse Jacobian du/dv."""
        N = self.N
        vx, vy = v[:N], v[N:]
        rho = np.hypot(vx, vy)
        small = rho < 1e-6
        safe = np.where(small, 1.0, rho)
        th = np.tanh(rho)
        s = np.where(small, 1.0 - rho ** 2 / 3, th / safe)
        q = np.where(small, -2.0 / 3, (1.0 - th ** 2 - th / safe) / safe ** 2)     # (ds/drho) / rho
        jxx, jyy, jxy = s + q * vx * vx, s + q * vy * vy, q * vx * vy
        M = sparse.bmat([[sparse.diags(jxx), sparse.diags(jxy)], [sparse.diags(jxy), sparse.diags(jyy)]], format="csr")
        return np.r_[s * vx, s * vy], M

    def unsquash(self, x, max_amplitude=0.99):
        """Inverse of squash (amplitudes capped at ``max_amplitude`` so the map is not saturated)."""
        ux, uy = self.controls(x)
        a = np.hypot(ux, uy)
        b = np.minimum(a, max_amplitude)
        f = np.where(a > 1e-12, np.arctanh(b) / np.where(a > 1e-12, a, 1.0), 1.0)
        return np.r_[f * ux, f * uy]

    # ------------------------------------------------------------- goals ---
    def _gate_rows(self, R, dR, weight, target=None):
        kind, Rt, m0, mt = target if target is not None else (self.cfg.gate, self.Rt, self.m0, self.mt)
        if kind == "UR":
            return weight * (R - Rt).ravel() / 2, weight * dR.reshape(len(dR), 9).T / 2
        m = R @ m0
        return weight * (m - mt) / 2, weight * (dR @ m0).T / 2

    def goals(self, x, with_jac=True):
        """Goal residuals (dense) and Jacobian (rows x 2N), plus a summary dict.

        One builder for every point (offset, scaling): the gate rows (order 0), and in local mode the
        Taylor orders of the detuning error and the first-order B1 error. All rows are divided by
        sqrt(number of points) (with user weights: sqrt(w_i / sum w)), so the squared residual is a
        (weighted) mean over the points and its balance with the soft terms does not depend on how many
        points there are. Summary: worst gate infidelity ("max infidelity", the solved criterion), mean
        infidelity and the largest gate error in degrees (info), and the largest |A_k|/T^k, |B1|/(Ω T)."""
        ux, uy = self.controls(x)
        cfg, st = self.cfg, self.st
        rows, jac, summ = [], [], {}
        marks = []                                   # index in rows where each point starts
        pts = cfg.ensemble()
        local = cfg.mode == "local"
        orders = list(cfg.orders) if local else []
        K = max(orders) if orders else 0             # Taylor jets only if an order is wanted
        norm = 1.0 / math.sqrt(len(pts))
        inf, gerr = [], []
        props = [propagator_and_gradient(ux, uy, self.dt, self.umax, off, sc, K) for off, sc in pts]
        rots = [rotation_with_gradient(X[0], dX[:, 0]) if tgt is not None else None
                for (X, dX), tgt in zip(props, self.targets)]
        targets = list(self.targets)
        if cfg.freedom != "none" and any(t is not None for t in targets):
            # free z-rotation after the pulse, fitted for the current pulse (variable projection: the Jacobian
            # is taken at the fitted rotation)
            gi = [i for i, t in enumerate(targets) if t is not None]
            ab = [_overlaps(rots[i][0], targets[i]) for i in gi]
            w = np.ones(len(pts))
            if self.user_w is not None:
                w = w * self.user_w
            if self.point_weights is not None and len(self.point_weights) == len(pts):
                w = w * np.asarray(self.point_weights, float)
            tau_max = (cfg.max_delay_us * 1e-6) if cfg.max_delay_us > 0 else self.T
            phi, tau = fit_freedom([v[0] for v in ab], [v[1] for v in ab], [pts[i][0] for i in gi], w[gi],
                                   cfg.freedom, tau_max)
            self.freedom_fit = (phi, tau)
            for i in gi:
                kind, Rt, m0, mt = targets[i]
                Z = rz(phi + 2 * math.pi * pts[i][0] * tau)
                targets[i] = (kind, Z @ Rt, m0, None) if kind == "UR" else (kind, None, m0, Z @ mt)
            summ["free phase [°]"] = float((math.degrees(phi) + 180) % 360 - 180)
            if cfg.freedom == "linear":
                summ["free delay [µs]"] = float(tau * 1e6)
        for (off, sc), tgt, (X, dX), rot in zip(pts, targets, props, rots):
            marks.append(len(rows))
            pp_m0 = tgt[2] if tgt is not None and tgt[0] == "PP" else None     # PP: only the error that moves m0
            if orders:
                mg = magnus_with_gradient(X, dX)
                for k in orders:
                    A, dA = mg[k - 1]
                    if pp_m0 is not None:
                        A, dA = np.cross(A, pp_m0), np.cross(dA, pp_m0)
                    w = st.weight(k) / self.T ** k * norm
                    rows.append(w * A); jac.append(w * dA.T)
                    key = f"|A{k}|/T^{k}" + (" ⊥m0" if cfg.gate == "PP" else "")
                    summ[key] = max(summ.get(key, 0.0), float(np.linalg.norm(A) / self.T ** k))
            if local and cfg.b1:
                Xa, dXa = propagator_and_gradient(ux, uy, self.dt, self.umax, off, sc, 1, "amplitude")
                U0d = Xa[0].conj().T
                B = _vec(U0d @ Xa[1], 1)
                dB = _vec(np.conj(np.swapaxes(dXa[:, 0], -1, -2)) @ Xa[1] + U0d @ dXa[:, 1], 1)
                if pp_m0 is not None:
                    B, dB = np.cross(B, pp_m0), np.cross(dB, pp_m0)
                nut = 2 * math.pi * self.umax * self.T              # maximal nutation angle
                w = st.weight("b1") / nut * norm
                rows.append(w * B); jac.append(w * dB.T)
                summ["|B1|/(Ω T)"] = max(summ.get("|B1|/(Ω T)", 0.0), float(np.linalg.norm(B) / nut))
            if tgt is not None:
                R, dR = rot
                r, j = self._gate_rows(R, dR, st.weight("gate") * norm if local else st.weight("gate") / math.sqrt(len(pts)),
                                       tgt)
                rows.append(r); jac.append(j)
                kind, Rt, m0, mt = tgt
                if kind == "UR":                     # 1 - (3 + tr(Rt^T R)) / 6 = |R - Rt|_F^2 / 12 (no cancellation)
                    inf.append(float(np.sum((R - Rt) ** 2)) / 12)
                    gerr.append(rotation_angle_deg(Rt.T @ R))
                else:                                # 1 - (1 + m . mt) / 2 = |m - mt|^2 / 4
                    m = R @ m0
                    inf.append(float(np.sum((m - mt) ** 2)) / 4)
                    gerr.append(_angle_between_deg(m, mt))
        if inf:
            summ["max infidelity"] = float(np.max(inf))
            summ["mean infidelity"] = float(np.mean(inf))
            summ["gate error [°]"] = float(np.max(gerr))
        sizes = np.cumsum([0] + [len(v) for v in rows])
        bounds = [(sizes[a], sizes[b]) for a, b in zip(marks, marks[1:] + [len(rows)])]
        r = np.concatenate(rows) if rows else np.zeros(0)
        Jm = np.vstack(jac) if jac else np.zeros((0, 2 * self.N))
        uw = self.user_w
        if uw is not None:                           # user weights: rows of point i x sqrt(w_i)
            scale = np.ones(len(r))
            for (a, b), wt in zip(bounds, uw):
                scale[a:b] = math.sqrt(wt)
            r, Jm = r * scale, Jm * scale[:, None]
        self.point_errors = np.array([np.linalg.norm(r[a:b]) for a, b in bounds])
        if self.point_weights is not None and len(self.point_weights) == len(bounds):
            scale = np.ones(len(r))
            for (a, b), wt in zip(bounds, self.point_weights):
                scale[a:b] = math.sqrt(wt)
            r, Jm = r * scale, Jm * scale[:, None]
        if len(bounds) > 1:
            summ["point error spread"] = float(np.max(self.point_errors) / max(np.mean(self.point_errors), 1e-300))
        return r, Jm, summ

    def reached(self, summ):
        st = self.st
        for k, v in summ.items():
            if k in INFO_KEYS:
                continue
            if k == "max infidelity":
                if v >= st.gate_infidelity_tolerance:
                    return False
            elif k == "rate-limit excess":
                if v >= LIMIT_EXCESS_TOL:
                    return False
            elif v >= st.solved_tolerance:
                return False
        return True

    # ------------------------------------------------ limits and soft terms ---
    def rates(self, x):
        """Amplitude rate [%/µs] and phase rate [kHz] between neighbouring sub-pulses (phase rate only where
        both amplitudes exceed 1 %)."""
        ux, uy = self.controls(x)
        a = np.hypot(ux, uy)
        amp_rate = np.diff(a) * 100 / (self.dt * 1e6)
        dphi = np.arctan2(ux[:-1] * uy[1:] - uy[:-1] * ux[1:], ux[:-1] * ux[1:] + uy[:-1] * uy[1:])
        f = np.abs(dphi) / (2 * math.pi * self.dt) / 1e3
        f = np.where((a[:-1] > 0.01) & (a[1:] > 0.01), f, 0.0)
        return amp_rate, f

    def extra(self, x):
        """Rows of the rate limits (one-sided, squared-hinge cost) and of the soft terms (stay close,
        smoothness; scaled by ``soft_scale``), their sparse Jacobian with respect to x, and info values."""
        cfg, N = self.cfg, self.N
        ux, uy = self.controls(x)
        rows, mats, info = [], [], {}
        excess = 0.0
        if cfg.max_amp_rate > 0 or cfg.max_phase_rate_khz > 0:
            eps = 1e-9
            a = np.sqrt(ux * ux + uy * uy + eps * eps)
            x1, y1, x2, y2, a1, a2 = ux[:-1], uy[:-1], ux[1:], uy[1:], a[:-1], a[1:]
            i = np.arange(N - 1)
            cols = [i, i + N, i + 1, i + 1 + N]                         # x1, y1, x2, y2
            w = cfg.limit_weight * self.limit_scale
            if cfg.max_amp_rate > 0:
                X = cfg.max_amp_rate
                k = 100 / (self.dt * 1e6)
                rate = (a2 - a1) * k
                h = np.maximum(np.abs(rate) - X, 0.0) / X
                fh, dfh = _limit_shape(h)
                g = np.where(h > 0, w * dfh * np.sign(rate) * k / X, 0.0)
                vals = [-g * x1 / a1, -g * y1 / a1, g * x2 / a2, g * y2 / a2]
                rows.append(w * fh); mats.append(sparse.csr_matrix((np.concatenate(vals), (np.tile(i, 4), np.concatenate(cols))),
                                                                  shape=(N - 1, 2 * N)))
                excess = max(excess, float(np.max(h, initial=0.0)))
            if cfg.max_phase_rate_khz > 0:
                Y = cfg.max_phase_rate_khz * 1e3
                c = x1 * y2 - y1 * x2
                d = x1 * x2 + y1 * y2
                den = np.maximum(c * c + d * d, 1e-300)
                dphi = np.arctan2(c, d)
                f = dphi / (2 * math.pi * self.dt)
                h = np.maximum(np.abs(f) - Y, 0.0) / Y
                gw = a1 * a2                                             # phase jumps at zero amplitude are free
                fh, dfh = _limit_shape(h)
                r = w * gw * fh
                sg = np.where(h > 0, dfh * np.sign(f) / (2 * math.pi * self.dt * Y), 0.0)
                dphi_d = [(d * y2 - c * x2) / den, (-d * x2 - c * y2) / den, (-d * y1 - c * x1) / den, (d * x1 - c * y1) / den]
                dg = [x1 / a1 * a2, y1 / a1 * a2, x2 / a2 * a1, y2 / a2 * a1]
                vals = [w * (dg[j] * fh + gw * sg * dphi_d[j]) for j in range(4)]
                rows.append(r); mats.append(sparse.csr_matrix((np.concatenate(vals), (np.tile(i, 4), np.concatenate(cols))),
                                                              shape=(N - 1, 2 * N)))
                excess = max(excess, float(np.max(gw * h, initial=0.0)))
            info["rate-limit excess"] = excess
        if self.soft_scale > 0:
            if cfg.regularization == "stay_close" and cfg.regularization_weight > 0:
                s = math.sqrt(cfg.regularization_weight * self.soft_scale)
                rows.append(s * (x - self.x0)); mats.append(sparse.identity(2 * N, format="csr") * s)
            if cfg.smooth_weight > 0 and N > 1:
                s = math.sqrt(cfg.smooth_weight * self.soft_scale * N)
                rows.append(s * (self.D2 @ x)); mats.append(self.D2 * s)
        info["roughness"] = float(N * np.sum((self.D2 @ x) ** 2))
        amp_rate, ph_rate = self.rates(x)
        info["max amplitude rate [%/µs]"] = float(np.max(np.abs(amp_rate), initial=0.0))
        info["max phase rate [kHz]"] = float(np.max(ph_rate, initial=0.0))
        if not rows:
            return np.zeros(0), sparse.csr_matrix((0, 2 * N)), info
        return np.concatenate(rows), sparse.vstack(mats, format="csr"), info

    def evaluate_blocks(self, x):
        """Goal residuals with dense Jacobian, extra residuals (limits, soft terms) with sparse Jacobian (both
        with respect to the controls x), the summary and the norm of the goal residuals."""
        r, Jg, summ = self.goals(x)
        re, Je, info = self.extra(x)
        summ.update(info)
        summ["max amplitude [%]"] = float(100 * np.max(np.hypot(*self.controls(x))))
        goal = float(np.linalg.norm(self.point_errors)) if self.point_errors is not None else float(np.linalg.norm(r))
        return r, Jg, re, Je, summ, goal                     # goal residual without the Lawson weights

    def evaluate(self, x):
        """All residuals, their sparse Jacobian (with respect to the controls x), summary and goal residual."""
        r, Jg, re, Je, summ, goal = self.evaluate_blocks(x)
        return np.r_[r, re], sparse.vstack([sparse.csr_matrix(Jg), Je], format="csr"), summ, goal


# ------------------------------------------------------------------ solver ---

class _Stop(Exception):
    def __init__(self, reason):
        self.reason = reason


def lm_step(G, rg, S, rs, mu):
    """Levenberg-Marquardt step d minimizing |rg + G d|^2 + |rs + S d|^2 + mu |d|^2 for few dense rows G and
    sparse rows S (smoothness, limits, stay close: banded): A = S^T S + mu I is sparse, the dense rows enter
    through the Woodbury identity (an m x m system, m = number of goal rows). Without S this is the
    minimum-norm step -G^T (G G^T + mu I)^-1 rg of an underdetermined problem."""
    from scipy.sparse.linalg import splu
    n = G.shape[1]
    b = -(G.T @ rg)
    if S is not None and S.shape[0]:
        b = b - S.T @ rs
        A = (S.T @ S + mu * sparse.identity(n)).tocsc()
        lu = splu(A)
        solve = lambda v: lu.solve(np.asarray(v, dtype=float))
        Ab = solve(b)
        if G.shape[0] == 0:
            return Ab
        AG = solve(np.ascontiguousarray(G.T))                     # (n, m)
        K = np.eye(G.shape[0]) + G @ AG
        return Ab - AG @ np.linalg.solve(K, G @ Ab)
    # A = mu I
    K = mu * np.eye(G.shape[0]) + G @ G.T
    return b / mu - G.T @ np.linalg.solve(K, G @ b) / mu


def levenberg_marquardt(fun, p, max_nfev, xtol, ftol, gtol, lam_box=None):
    """Minimize |r(p)|^2 with r, Jacobians from fun(p) -> (rg, G dense, rs, S sparse). Every call of fun is
    one evaluation (residuals and exact Jacobian together). ``lam_box`` ([damping] or [None]) carries the
    damping in and out, so that a resumed run continues with it. Returns (p, status message)."""
    rg, G, rs, S = fun(p)
    cost = rg @ rg + rs @ rs
    scale = (np.sum(G * G) + (S.multiply(S).sum() if S is not None else 0.0)) / max(len(p), 1)
    lam = lam_box[0] if lam_box is not None and lam_box[0] is not None else 1e-2
    nfev = 1
    while nfev < max_nfev:
        grad = G.T @ rg + (S.T @ rs if S is not None and S.shape[0] else 0.0)
        if gtol and np.max(np.abs(grad), initial=0.0) <= gtol * max(1.0, cost):
            return p, "gradient tolerance reached"
        mu = lam * max(scale, 1e-300)
        d = lm_step(G, rg, S, rs, mu)
        pred = cost - (np.sum((rg + G @ d) ** 2) + (np.sum((rs + S @ d) ** 2) if S is not None and S.shape[0] else 0.0))
        rg2, G2, rs2, S2 = fun(p + d)
        nfev += 1
        cost2 = rg2 @ rg2 + rs2 @ rs2
        if cost2 < cost:
            small_step = xtol and np.linalg.norm(d) <= xtol * (xtol + np.linalg.norm(p))
            small_gain = ftol and (cost - cost2) <= ftol * cost and pred <= ftol * cost
            p, rg, G, rs, S, cost = p + d, rg2, G2, rs2, S2, cost2
            scale = (np.sum(G * G) + (S.multiply(S).sum() if S is not None else 0.0)) / max(len(p), 1)
            lam = max(lam / 3.0, 1e-15)
            if lam_box is not None:
                lam_box[0] = lam
            if small_step or small_gain:
                return p, "step / improvement tolerance reached"
        else:
            lam *= 4.0
            if lam_box is not None:
                lam_box[0] = lam
            if lam > 1e16:
                return p, "no further improvement possible"
    return p, "evaluation limit reached"


def _start_parameters(prob, cfg):
    """Parameters of the first stage of a new candidate: the start pulse (tiny symmetry-breaking noise, or the
    spread for starts 2, 3, ...: white or smooth) or, for start_kind "random", a random pulse."""
    rng = np.random.default_rng(cfg.seed)
    if cfg.start_kind == "random":
        return prob.params(control_noise(prob.N, prob.dt, 1.0, cfg.noise_correlation_us, rng))
    if cfg.start_noise <= 0:
        return prob.p0.copy()
    if cfg.noise_correlation_us > 0:
        w = prob.unsquash(prob.x0) + control_noise(prob.N, prob.dt, cfg.start_noise, cfg.noise_correlation_us, rng)
        return prob.params(w)
    # an exactly symmetric start (e.g. a constant-phase pulse) can keep the solver in a symmetric subspace
    return prob.p0 + rng.normal(0.0, cfg.start_noise, prob.p0.shape)


STOP_MESSAGES = {"solved": "all goals reached", "limit": "evaluation limit reached", "time": "time limit reached",
                 "user": "stopped by the user", "paused": "evaluation budget of this round used"}


def optimize(ps: PulseSequence, config: GrapeConfig, settings: SolverSettings = None, callback=None,
             snapshot=None, snapshot_interval_s=1.0, start_index=0, state=None, budget=None, start_summary=True):
    """Run the optimization from the start pulse ``ps`` (start ``start_index`` of a multi-start: other
    seed and, for index > 0, the spread or a random start). ``callback(nfev, goal_residual)`` returning True
    stops it; ``snapshot(pulse, goal_residual)`` receives the best pulse every ``snapshot_interval_s``.

    Solver: Levenberg-Marquardt (lm_step) on the parameters, in the stages of stage_plan: coarse spline
    levels, then the controls with the soft terms (stay close, smoothness) lowered stage by stage until the
    goals are reached. Non-final stages get a share of the budget and relaxed tolerances.

    Resumable: ``budget`` limits the evaluations of this call (the run then pauses, stop_reason "paused");
    ``state`` (info["state"] of a paused run, or {"stage", "w"} for a new candidate given by its unsquashed
    controls w) continues it. settings.max_nfev is the total of all calls. Returns (pulse, info dict)."""
    st = settings or SolverSettings()
    cfg = config.variant(start_index)
    cfg.check()
    plan = stage_plan(cfg, ps.N)
    problems = {}

    def problem(s):
        ctrl, m = plan[s][:2]
        key = (ctrl, m if ctrl == "spline" else 0)
        if key not in problems:
            problems[key] = GrapeProblem(ps, cfg.replace(controls=ctrl, n_spline=m), st)
        return problems[key]

    t0 = time.time()
    last_snap = [t0]
    weights = None
    lam_box = [None]                                     # LM damping of the current stage (kept when resumed)
    minimax = {"worst": np.inf, "w": None, "weights": None}
    if state is None:
        s, used, nfev0 = 0, 0, 0
        p = _start_parameters(problem(0), cfg)
    else:
        s, used, nfev0 = min(int(state.get("stage", 0)), len(plan) - 1), int(state.get("stage_nfev", 0)), int(state.get("nfev", 0))
        p = (np.asarray(state["p"], float) if state.get("p") is not None
             else problem(s).params(np.asarray(state["w"], float)))
        if state.get("weights") is not None:
            weights = np.asarray(state["weights"], float)
        lam_box[0] = state.get("lam")
    count = {"nfev": nfev0, "history": []}
    limit = st.max_nfev if budget is None else min(st.max_nfev, nfev0 + max(1, int(budget)))
    tols = [t if t and t > 0 else 0.0 for t in (st.xtol, st.ftol, st.gtol)]
    n_pre = sum(1 for stage in plan if stage[3] in ("coarse", "main")) - 1      # stages before the last main stage
    plan_kinds = [stage[3] for stage in plan]
    late_limits = "limits" in plan_kinds
    share = max(10, st.max_nfev // (2 * max(n_pre, 1)))
    reason, message, cost = "converged", "", float("nan")

    def advance():
        nonlocal s, p, used, weights
        w = problem(s).wvec(p)
        s, used = s + 1, 0
        lam_box[0] = None
        p = problem(s).params(w)
        if plan[s][3] == "reweight":                     # Lawson: weights grow with the errors of the points
            nxt = problem(s)
            nxt.point_weights = weights
            nxt.goals(nxt.to_controls(p)[0])
            worst = float(np.max(nxt.point_errors))
            if worst < minimax["worst"]:                  # keep the round with the smallest worst-point error
                minimax.update(worst=worst, w=nxt.wvec(p).copy(), weights=weights)
            weights = lawson_update(np.ones(len(nxt.point_errors)) if weights is None else weights,
                                    nxt.point_errors, cfg.reweight_power)

    while True:
        prob = problem(s)
        prob.soft_scale = plan[s][2]
        prob.point_weights = weights
        final = s == len(plan) - 1
        kind = plan[s][3]
        last_main = kind == "main" and not final and plan[s + 1][3] != "main"
        prob.limit_scale = 0.0 if late_limits and plan_kinds.index("limits") > s else 1.0
        if last_main and plan[s + 1][3] == "limits":      # leave half of the rest for the stage with the limits
            stage_left = (st.max_nfev - count["nfev"] + used) // 2 - used
        elif final or last_main or (kind == "limits" and plan[s + 1][3] == "reweight"):
            # the main stage runs as without reweighting; the rounds share the rest
            stage_left = st.max_nfev - count["nfev"]
        elif kind == "reweight":
            stage_left = max(10, (st.max_nfev - count["nfev"] + used) // (len(plan) - s)) - used
        else:
            stage_left = share - used
        call_left = limit - count["nfev"]
        if count["nfev"] >= st.max_nfev:
            reason = "limit"
            break
        if stage_left <= 0:
            if final:
                reason = "limit"
                break
            advance()
            continue
        if call_left <= 0:
            reason = "paused"
            break
        best = {"p": p.copy(), "cost": np.inf, "goal": np.inf, "summ": None}

        def fun(q):
            if st.time_limit_s > 0 and time.time() - t0 > st.time_limit_s:
                raise _Stop("time")
            u, Mp = prob.to_controls(q)
            rg, Jg, rs, Js, summ, goal = prob.evaluate_blocks(u)
            G = np.asarray(Jg @ Mp) if sparse.issparse(Mp) else Jg @ Mp
            S = sparse.csr_matrix(Js @ Mp) if Js.shape[0] else None
            c = float(np.sqrt(rg @ rg + rs @ rs))
            count["nfev"] += 1
            if c < best["cost"]:
                best.update(p=np.array(q), cost=c, goal=goal, summ=summ)
            count["history"].append(best["goal"])
            if final and st.stop_when_solved and prob.reached(summ):
                best.update(p=np.array(q), cost=c, goal=goal, summ=summ)
                raise _Stop("solved")
            if callback is not None and callback(count["nfev"], best["goal"]):
                raise _Stop("user")
            if snapshot is not None and time.time() - last_snap[0] > snapshot_interval_s:
                last_snap[0] = time.time()
                snapshot(prob.pulse(prob.to_controls(best["p"])[0]), best["goal"])
            return rg, G, rs, (S if S is not None else sparse.csr_matrix((0, len(q))))

        stage_tols = tols if final or last_main or kind in ("reweight", "limits") else [max(t, SOFT_STAGE_TOL) for t in tols]
        n_before = count["nfev"]
        stopped = None
        try:
            _, message = levenberg_marquardt(fun, p, min(stage_left, call_left), *stage_tols, lam_box=lam_box)
        except _Stop as e:
            stopped = e.reason
        used += count["nfev"] - n_before
        if best["summ"] is not None:
            p, cost = best["p"], best["cost"]
        if stopped is not None:
            reason = stopped
            break
        if best["summ"] is not None and prob.reached(best["summ"]) and (st.stop_when_solved or final):
            reason = "solved"
            break
        converged = message != "evaluation limit reached"
        stage_done = converged or count["nfev"] - n_before >= stage_left        # else: the call budget ran out
        if not stage_done:
            reason = "paused"
            break
        if final:
            reason = "converged" if converged else "limit"
            break
        advance()
    message = STOP_MESSAGES.get(reason, message)
    prob = problem(s)
    if minimax["w"] is not None and reason in ("converged", "limit") and plan[s][3] == "reweight":
        prob.point_weights = None
        prob.goals(prob.to_controls(p)[0])
        if float(np.max(prob.point_errors)) > minimax["worst"]:
            p, weights = prob.params(minimax["w"]), minimax["weights"]
    u = prob.to_controls(p)[0]
    prob.soft_scale, prob.point_weights = plan[-1][2], None
    _, _, summ, goal = prob.evaluate(u)
    errors = prob.point_errors
    start = prob.evaluate(prob.x0)[2] if start_summary else {}
    prob.soft_scale = plan[s][2]
    return prob.pulse(u), {"stop_reason": reason, "message": message, "reached": prob.reached(summ),
                           "values": summ, "start_summary": start, "nfev": count["nfev"], "goal_residual": goal,
                           "cost": cost if np.isfinite(cost) else goal, "elapsed_s": time.time() - t0,
                           "history": count["history"], "start_index": start_index, "cancelled": reason == "user",
                           "finished": reason != "paused", "stage": s, "n_stages": len(plan),
                           "state": {"stage": s, "p": p.tolist(), "stage_nfev": used, "nfev": count["nfev"],
                                     "weights": None if weights is None else weights.tolist(), "lam": lam_box[0]},
                           "point_weights": None if weights is None else weights.tolist(),
                           "point_errors": None if errors is None else errors.tolist(),
                           "w": prob.wvec(p).tolist()}


def band_infidelity(ps: PulseSequence, config: GrapeConfig, offsets_hz, scale=1.0):
    """Infidelity of the pulse at every offset of a (fine) test grid, without gradients: UR 1 - (3 + tr(Rt^T R)) / 6,
    PP 1 - (1 + (R m0) . m_target) / 2 (the same measures as the broadband goal)."""
    from psa.rotations import propagate_batch, quaternion_to_matrix
    R = quaternion_to_matrix(propagate_batch(ps.as_array(), ps.dt, ps.umax * scale, np.asarray(offsets_hz, float)))
    if config.gate == "PP":
        m0, mt = _unit(config.m0), _unit(config.m_target)
        return 1 - (1 + (R @ m0) @ mt) / 2
    Rt = rotation_from_axis_angle(config.target_axis, config.target_angle_deg)
    return 1 - (3 + np.einsum("ij,kij->k", Rt, R)) / 6


def band_quality(ps: PulseSequence, config: GrapeConfig, lo_hz, hi_hz, n=401, scale_lo=1.0, scale_hi=1.0, n_scales=1,
                 thresholds=(1e-3, 1e-4, 1e-5)):
    """Worst and mean infidelity over [lo, hi] on n test points (not the optimization points: shows dips
    between them), and the fraction of the band below each threshold. With n_scales > 1: over the area
    [lo, hi] x [scale_lo, scale_hi] (B1 scaling), plus the position of the worst value."""
    offs = np.linspace(lo_hz, hi_hz, int(n))
    if int(n_scales) <= 1:
        inf = np.maximum(band_infidelity(ps, config, offs), 0.0)
        out = {"worst infidelity": float(np.max(inf)), "mean infidelity": float(np.mean(inf))}
    else:
        scales = np.linspace(scale_lo, scale_hi, int(n_scales))
        area = np.maximum(np.array([band_infidelity(ps, config, offs, s) for s in scales]), 0.0)
        i, j = np.unravel_index(int(np.argmax(area)), area.shape)
        inf = area.ravel()
        out = {"worst infidelity": float(area[i, j]), "mean infidelity": float(np.mean(area)),
               "worst at offset [kHz]": float(offs[j] / 1e3), "worst at B1 [%]": float(scales[i] * 100)}
    for t in thresholds:
        out[f"fraction below {t:g}"] = float(np.mean(inf <= t))
    return out


def nearest_point_index(config: GrapeConfig, offsets_hz, scales):
    """Index of the nearest point (offset and B1 distances relative to the spans of the points) for every
    (scale, offset) of a grid: array (len(scales), len(offsets_hz))."""
    pts = np.array(config.ensemble(), float)
    so = max(float(np.ptp(pts[:, 0])), 1.0)
    ss = max(float(np.ptp(pts[:, 1])), 1e-3)
    S, O = np.meshgrid(np.asarray(scales, float), np.asarray(offsets_hz, float), indexing="ij")
    d = ((O[..., None] - pts[:, 0]) / so) ** 2 + ((S[..., None] - pts[:, 1]) / ss) ** 2
    return np.argmin(d, axis=-1)


def fitted_freedom(ps: PulseSequence, config: GrapeConfig):
    """(phi, tau) of the config's free z-rotation after the pulse, fitted at its points for this pulse (as in the
    optimization, with the user weights); (0, 0) without freedom."""
    if config.freedom == "none":
        return 0.0, 0.0
    from psa.design.fidelity import matrices
    pts = config.ensemble()
    targets = [goal_target(g) for g in config.point_targets()]
    gi = [i for i, t in enumerate(targets) if t is not None]
    if not gi:
        return 0.0, 0.0
    R = [matrices(ps, [pts[i][0]], [pts[i][1]])[0, 0] for i in gi]
    ab = [_overlaps(r, targets[i]) for r, i in zip(R, gi)]
    uw = config.user_weights()
    w = np.ones(len(pts)) if uw is None else uw
    tau_max = config.max_delay_us * 1e-6 if config.max_delay_us > 0 else ps.duration
    return fit_freedom([v[0] for v in ab], [v[1] for v in ab], [pts[i][0] for i in gi], w[gi], config.freedom, tau_max)


def goal_infidelity_grid(ps: PulseSequence, config: GrapeConfig, offsets_hz, scales):
    """1 - F on a (scale, offset) grid where every location is judged by the goal of its nearest point (the
    config's gate goal or the point's own goal); NaN where that point has no gate goal. With freedom, the targets
    are rotated by the z-rotation fitted at the points (phi + 2 pi offset tau, see fitted_freedom)."""
    from psa.design.fidelity import infidelity_from_matrices, matrices
    R = matrices(ps, offsets_hz, scales)
    targets = config.point_targets()
    if config.freedom != "none":
        phi, tau = fitted_freedom(ps, config)
        Z = rz(phi + 2 * math.pi * np.asarray(offsets_hz, float) * tau)             # (n_off, 3, 3), per column
        idx = (nearest_point_index(config, offsets_hz, scales) if len({goal_key(t) for t in targets}) > 1
               else np.zeros((len(scales), len(offsets_hz)), int))
        out = np.full(idx.shape, np.nan)
        for i_key, t in {goal_key(t): t for t in targets}.items():
            if t is None:
                continue
            members = [i for i, u in enumerate(targets) if goal_key(u) == i_key]
            mask = np.isin(idx, members)
            if not np.any(mask):
                continue
            kind, Rt, m0, mt = goal_target(t)
            if kind == "UR":
                full = np.sum((R - Z @ Rt) ** 2, axis=(-2, -1)) / 12
            else:
                full = np.sum((R @ m0 - Z @ mt) ** 2, axis=-1) / 4
            out[mask] = full[mask]
        return out
    if len({goal_key(t) for t in targets}) == 1:                   # one goal everywhere
        idx = np.zeros((len(scales), len(offsets_hz)), int)
    else:
        idx = nearest_point_index(config, offsets_hz, scales)
    out = np.full(idx.shape, np.nan)
    groups = {}
    for i, t in enumerate(targets):
        groups.setdefault(goal_key(t), (t, []))[1].append(i)
    for key, (t, members) in groups.items():
        if t is None:
            continue
        mask = np.isin(idx, members)
        if np.any(mask):
            out[mask] = infidelity_from_matrices(R[mask], goal_reference(t))
    return out


def goal_reference(goal):
    """The fidelity-module reference dict of a goal dict."""
    if goal["gate"] == "UR":
        return {"kind": "UR", "axis": goal["axis"], "angle": goal["angle"]}
    return {"kind": "PP", "m0": goal["m0"], "mt": goal["mt"]}


def area_quality(ps: PulseSequence, config: GrapeConfig, n_offsets=201, n_b1=11, thresholds=(1e-3, 1e-4, 1e-5)):
    """Gate infidelity on a fine test grid over the area of the point selection (offset range x B1 range;
    one B1 value if all points have the same scaling). Shows what happens between the points. Returns
    None without a gate goal or if all points have the same offset; else worst / mean infidelity, the fraction
    of the area below each threshold, the position of the worst value and the grid size."""
    targets = config.point_targets()
    if all(t is None for t in targets):
        return None
    pts = config.ensemble()
    offs, scales = [p[0] for p in pts], [p[1] for p in pts]
    if max(offs) - min(offs) <= 0:
        return None
    o = np.linspace(min(offs), max(offs), int(n_offsets))
    sc = np.linspace(min(scales), max(scales), int(n_b1)) if max(scales) > min(scales) else np.array([scales[0]])
    inf = goal_infidelity_grid(ps, config, o, sc)
    if np.all(np.isnan(inf)):
        return None
    inf = np.where(np.isnan(inf), 0.0, np.maximum(inf, 0.0))
    i, j = np.unravel_index(int(np.argmax(inf)), inf.shape)
    out = {"worst infidelity": float(inf[i, j]), "mean infidelity": float(np.mean(inf)),
           "worst at offset [kHz]": float(o[j] / 1e3), "worst at B1 [%]": float(sc[i] * 100), "grid": inf.shape}
    for t in thresholds:
        out[f"fraction below {t:g}"] = float(np.mean(inf <= t))
    return out


def pulse_rates(ps: PulseSequence, min_amp=1.0):
    """Largest amplitude change rate [%/µs] and phase change rate [kHz] of the sampled pulse (phase steps
    only where both neighbouring sub-pulses have more than ``min_amp`` % amplitude; wrapped to ±180°)."""
    a = np.asarray(ps.amplitude, float)
    dt_us = ps.dt * 1e6
    amp_rate = float(np.max(np.abs(np.diff(a)) / dt_us, initial=0.0))
    dphi = np.angle(np.exp(1j * np.radians(np.diff(np.asarray(ps.phase, float)))))
    f = np.abs(dphi) / (2 * math.pi * ps.dt) / 1e3
    f = np.where((a[:-1] > min_amp) & (a[1:] > min_amp), f, 0.0)
    return amp_rate, float(np.max(f, initial=0.0))


def selection_key(config: GrapeConfig, info):
    """Order of results: reached first; among those the smoothest (if smoothing is on), else the lowest goal
    residual."""
    second = info["values"].get("roughness", 0.0) if config.smooth_weight > 0 and info["reached"] else info["goal_residual"]
    return (not info["reached"], second)


# ------------------------------------------------------------- exploration ---

EXPLORATION_PRESETS = {
    # n_starts etc.; "correlation" is a fraction of the pulse duration (-> noise_correlation_us)
    "exploit": dict(n_starts=1, multistart_noise=0.5, start_pool="perturb", correlation=0.0, screen_rounds=0,
                    screen_nfev=40, keep_fraction=0.5, crossover_rate=0.0, n_kicks=0, kick_size=0.3,
                    kick_temperature=0.0, coarse_levels=[]),
    "balanced": dict(n_starts=8, multistart_noise=0.5, start_pool="perturb", correlation=0.05, screen_rounds=2,
                     screen_nfev=30, keep_fraction=0.5, crossover_rate=0.0, n_kicks=4, kick_size=0.3,
                     kick_temperature=0.0, coarse_levels=[]),
    "explore": dict(n_starts=24, multistart_noise=1.0, start_pool="mix", correlation=0.1, screen_rounds=3,
                    screen_nfev=25, keep_fraction=0.5, crossover_rate=0.0, n_kicks=8, kick_size=0.5,
                    kick_temperature=0.3, coarse_levels=[8, 16]),
}


def exploration_preset(name, duration_s):
    """GrapeConfig fields of an exploration preset for a pulse of the given duration."""
    d = dict(EXPLORATION_PRESETS[name])
    d["noise_correlation_us"] = round(d.pop("correlation") * duration_s * 1e6, 3)
    return d


class Search:
    """Exploration around GRAPE with a population of candidates, independent of where the work runs.

    Round r: every live candidate is refined (resumably) up to screen_nfev / keep_fraction^r evaluations;
    then the best keep_fraction survive (reached goals first, then the lowest cost), and a fraction of the
    dropped places is refilled with crossover children of two survivors. After screen_rounds rounds the
    survivors run to the end (settings.max_nfev per candidate). Then n_kicks kicks of the best result (in
    batches of ``parallel``): smooth perturbation of size kick_size, full re-optimization, accepted if better
    (Metropolis with kick_temperature on log10 of the goal residual, or of the roughness once reached).

    Drive it with tasks() (tasks ready now; run each with run_task, anywhere) and submit(result) for every
    finished one; ``done`` tells when it is over, best() gives the result."""

    def __init__(self, ps: PulseSequence, config: GrapeConfig, settings: SolverSettings = None, parallel=1):
        config.check()
        self.ps, self.cfg, self.st = ps, config, settings or SolverSettings()
        self.parallel = max(1, int(parallel))
        self.rng = np.random.default_rng(int(config.seed) + 104729)
        self.plan = stage_plan(config, ps.N)
        self.n_coarse = sum(1 for stage in self.plan if stage[3] == "coarse")
        self._pulse = [np.minimum(ps.amplitude, 100.0).tolist(), np.asarray(ps.phase, float).tolist(), ps.dt, ps.umax]
        self.cands = {}
        self._ready, self._out = [], {}
        self._tid = 0
        self.round = 0
        self.phase = "rounds"
        self.kicks_left = int(config.n_kicks)
        self.kicks_accepted = 0
        self.kicks_done = 0
        self.current = None
        self.stopped = False
        self.log = []
        self.errors = []
        for i in range(int(config.n_starts)):
            self._add("start", start_index=i)
        self._issue_round()

    # ------------------------------------------------------------ public ---
    @property
    def done(self):
        return self.phase == "done"

    def tasks(self):
        out, self._ready = self._ready, []
        return out

    @property
    def running(self):
        return len(self._out)

    def submit(self, result):
        cid = self._out.pop(result["tid"], None)
        if cid is None:
            return
        c = self.cands[cid]
        if "error" in result:
            c["alive"], c["finished"] = False, True
            self.errors.append(result["error"])
        else:
            info = result["info"]
            c.update(info=info, amp=result["amp"], phase=result["phase"], state=info["state"], finished=info["finished"])
            if info.get("cancelled"):
                self.stopped = True
        if not self._out and not self._ready:
            self._advance()

    def stop(self):
        """No new tasks; the search ends when the running ones have returned."""
        self.stopped = True
        for t in self._ready:
            self._out.pop(t["tid"], None)
        self._ready = []
        if not self._out:
            self.phase = "done"

    def best(self):
        """(candidate id, pulse, info) of the best result; info["search"] summarizes the search."""
        have = [cid for cid, c in self.cands.items() if c["info"] is not None]
        if not have:
            raise RuntimeError(self.errors[0] if self.errors else "no result")
        cid = min(have, key=lambda k: selection_key(self.cfg, self.cands[k]["info"]))
        c = self.cands[cid]
        info = dict(c["info"])
        info["search"] = {"best": cid, "kind": c["kind"], "candidates": len(self.cands), "log": list(self.log),
                          "reached": sum(1 for k in have if self.cands[k]["info"]["reached"]),
                          "evaluations": sum(self.cands[k]["info"]["nfev"] for k in have),
                          "kicks": self.kicks_done, "kicks_accepted": self.kicks_accepted}
        return cid, PulseSequence(np.asarray(c["amp"]), np.asarray(c["phase"]), self.ps.dt, self.ps.umax), info

    def infos(self):
        return [c["info"] for c in self.cands.values() if c["info"] is not None]

    # ----------------------------------------------------------- internal ---
    def _add(self, kind, start_index=0, state=None, parents=()):
        cid = len(self.cands)
        self.cands[cid] = {"kind": kind, "start_index": start_index, "state": state, "parents": tuple(parents),
                           "info": None, "amp": None, "phase": None, "alive": True, "finished": False}
        return cid

    def _task(self, cid, budget):
        c = self.cands[cid]
        task = {"tid": self._tid, "cid": cid, "pulse": self._pulse, "config": self.cfg.to_dict(),
                "settings": self.st.to_dict(), "start_index": c["start_index"], "state": c["state"], "budget": budget}
        self._out[self._tid] = cid
        self._ready.append(task)
        self._tid += 1

    def _round_budget(self):
        if self.round >= int(self.cfg.screen_rounds):
            return None
        return int(math.ceil(self.cfg.screen_nfev / self.cfg.keep_fraction ** self.round))

    def _issue_round(self):
        total = self._round_budget()
        for cid, c in self.cands.items():
            if c["alive"] and not c["finished"]:
                nfev = c["state"]["nfev"] if c["state"] and "nfev" in c["state"] else 0
                self._task(cid, None if total is None else max(1, total - nfev))

    def _screen_key(self, cid):
        info = self.cands[cid]["info"]
        if info["reached"]:
            return (0, info["values"].get("roughness", 0.0) if self.cfg.smooth_weight > 0 else info["goal_residual"])
        return (1, info["cost"])

    def _mask(self, N):
        if self.rng.random() < 0.5:                                    # blend of the two parents
            return np.full(N, self.rng.uniform(0.25, 0.75))
        t0 = self.rng.uniform(0.2, 0.8) * N                            # first part of one, rest of the other
        width = max(1.0, 0.05 * N)
        return 0.5 * (1 + np.tanh((np.arange(N) + 0.5 - t0) / width))

    def _crossover(self, a, b):
        ca, cb = self.cands[a], self.cands[b]
        m = self._mask(self.ps.N)
        m2 = np.r_[m, m]
        w = m2 * np.asarray(ca["info"]["w"]) + (1 - m2) * np.asarray(cb["info"]["w"])
        stage = min(ca["state"]["stage"], cb["state"]["stage"])
        return self._add("child", state={"stage": stage, "w": w.tolist(), "stage_nfev": 0, "nfev": 0}, parents=(a, b))

    def _kick_score(self, info):
        if info["reached"]:
            return (0, math.log10(max(info["values"].get("roughness", 0.0), 1e-300)) if self.cfg.smooth_weight > 0 else 0.0)
        return (1, math.log10(max(info["goal_residual"], 1e-300)))

    def _accept(self, new, cur):
        a, b = self._kick_score(new), self._kick_score(cur)
        if a[0] != b[0]:
            return a[0] < b[0]
        d = a[1] - b[1]
        if d < 0:
            return True
        T = self.cfg.kick_temperature
        return T > 0 and self.rng.random() < math.exp(-d / T)

    def _worth_kicking(self, info):
        return not (info["reached"] and self.cfg.smooth_weight == 0)

    def _advance(self):
        while not self._ready and not self._out and self.phase != "done":
            if self.stopped:
                self.phase = "done"
                return
            if self.phase == "rounds":
                self._next_round()
            else:
                self._next_kicks()

    def _next_round(self):
        cfg = self.cfg
        alive = [cid for cid, c in self.cands.items() if c["alive"] and c["info"] is not None]
        if not alive:
            self.phase = "done"
            return
        ranked = sorted(alive, key=self._screen_key)
        lead = self.cands[ranked[0]]["info"]
        if lead["reached"] and cfg.smooth_weight == 0 and self.st.stop_when_solved:
            self.log.append(f"round {self.round + 1}: goals reached by {sum(self.cands[k]['info']['reached'] for k in alive)} "
                            f"of {len(alive)} candidates")
            self._enter_kicks(ranked[0])
            return
        if self.round < int(cfg.screen_rounds) and len(alive) > 1:
            k = max(1, int(math.ceil(cfg.keep_fraction * len(alive))))
            for cid in ranked[k:]:
                self.cands[cid]["alive"] = False
            survivors = ranked[:k]
            n_child = int(round(cfg.crossover_rate * (len(alive) - k))) if k >= 2 else 0
            for _ in range(n_child):
                i, j = self.rng.choice(len(survivors), 2, replace=False)
                self._crossover(survivors[int(i)], survivors[int(j)])
            self.log.append(f"round {self.round + 1}: {len(alive)} candidates, {self._round_budget()} evaluations each "
                            f"→ kept {k}" + (f" + {n_child} crossover children" if n_child else ""))
            self.round += 1
            if k == 1 and not n_child:
                self.round = int(cfg.screen_rounds)
            self._issue_round()
            return
        if self.round < int(cfg.screen_rounds) or any(not self.cands[c]["finished"] for c in alive):
            self.round = int(cfg.screen_rounds)           # final round: the rest of the budget
            self._issue_round()
            if self._ready:
                return
        n_reached = sum(1 for k in alive if self.cands[k]["info"]["reached"])
        self.log.append(f"final: {len(alive)} candidate{'s' if len(alive) > 1 else ''}, {n_reached} reached the goals")
        self._enter_kicks(ranked[0] if alive else None)

    def _enter_kicks(self, cid):
        self.phase = "kicks"
        have = [k for k, c in self.cands.items() if c["info"] is not None]
        self.current = min(have, key=lambda k: selection_key(self.cfg, self.cands[k]["info"])) if have else cid
        self._batch = []
        if self.kicks_left > 0 and not self._worth_kicking(self.cands[self.current]["info"]):
            self.log.append("kicks skipped: the goals are already reached")
            self.kicks_left = 0
        self._issue_kicks()

    def _issue_kicks(self):
        if self.kicks_left <= 0 or self.current is None:
            self.phase = "done"
            return
        cur = self.cands[self.current]
        w0 = np.asarray(cur["info"]["w"])
        self._batch = []
        for _ in range(min(self.kicks_left, self.parallel)):
            w = w0 + control_noise(self.ps.N, self.ps.dt, self.cfg.kick_size, self.cfg.noise_correlation_us, self.rng)
            cid = self._add("kick", state={"stage": self.n_coarse, "w": w.tolist(), "stage_nfev": 0, "nfev": 0},
                            parents=(self.current,))
            self._batch.append(cid)
            self._task(cid, None)
        self.kicks_left -= len(self._batch)
        self.kicks_done += len(self._batch)

    def _next_kicks(self):
        accepted = 0
        for cid in self._batch:
            info = self.cands[cid]["info"]
            if info is not None and self._accept(info, self.cands[self.current]["info"]):
                self.current = cid
                accepted += 1
        self.kicks_accepted += accepted
        done = self.kicks_done
        self.log.append(f"kicks {done - len(self._batch) + 1}–{done}: {accepted} accepted")
        if not self._worth_kicking(self.cands[self.current]["info"]):
            self.kicks_left = 0
        self._issue_kicks()


def run_task(task, queue=None, stop_event=None):
    """Run one Search task (in any process). With a queue it reports ("progress", cid, nfev, goal_residual)
    and ("snapshot", cid, amp, phase, goal_residual); a set stop_event stops it. Returns the result dict
    for Search.submit."""
    tid, cid = task["tid"], task["cid"]
    try:
        amp, phase, dt, umax = task["pulse"]
        ps = PulseSequence(np.asarray(amp), np.asarray(phase), dt, umax)
        last = [-10]

        def callback(nfev, cost):
            if queue is not None and nfev - last[0] >= 5:
                last[0] = nfev
                queue.put(("progress", cid, nfev, cost))
            return stop_event is not None and stop_event.is_set()

        def snapshot(p, cost):
            queue.put(("snapshot", cid, p.amplitude.tolist(), p.phase.tolist(), cost))

        cfg = GrapeConfig.from_dict(task["config"])
        result, info = optimize(ps, cfg, SolverSettings.from_dict(task.get("settings")),
                                callback, snapshot if queue is not None and task.get("snapshots", True) else None,
                                start_index=task["start_index"], state=task.get("state"), budget=task.get("budget"))
        info["history"] = info["history"][-5000:]
        out = {"tid": tid, "cid": cid, "amp": result.amplitude.tolist(), "phase": result.phase.tolist(), "info": info}
        if task.get("band") is not None:                 # score on a fine test grid (quick mode)
            out["band"] = band_quality(result, cfg, *task["band"])
            out["rates"] = pulse_rates(result)
        return out
    except Exception as e:
        return {"tid": tid, "cid": cid, "error": f"{type(e).__name__}: {e}"}


def grape_worker(tasks, results, stop_event):
    """Worker process: runs tasks from the ``tasks`` queue until it gets None; puts ("done", result) (and the
    progress messages of run_task) into ``results``."""
    while True:
        task = tasks.get()
        if task is None:
            return
        results.put(("done", run_task(task, results, stop_event)))


def search(ps: PulseSequence, config: GrapeConfig, settings: SolverSettings = None, processes=1):
    """Run the whole Search (tasks in parallel processes if ``processes`` > 1). Returns (best pulse, its info,
    infos of all candidates)."""
    s = Search(ps, config, settings, parallel=processes)
    pool = None
    try:
        if processes > 1:
            import multiprocessing as mp
            pool = mp.get_context("spawn").Pool(processes)
        while not s.done:
            tasks = s.tasks()
            results = pool.map(run_task, tasks) if pool is not None and len(tasks) > 1 else [run_task(t) for t in tasks]
            for r in results:
                s.submit(r)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    _, pulse, info = s.best()
    return pulse, info, s.infos()


optimize_multistart = search
