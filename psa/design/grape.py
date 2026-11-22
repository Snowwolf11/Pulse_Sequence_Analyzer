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
Exploration: several starts (seeded perturbations) in parallel processes; the best one is kept.

Robustness (choose one):
    local      Taylor orders of the detuning error (A1..A3, SCQC) and optionally first-order B1 scaling,
               expanded around each offset of a list (not only 0)
    broadband  classical GRAPE: the gate goal is required on a grid of offsets x amplitude scalings
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
    n_starts: int = 1                            # multi-start: number of starts
    multistart_noise: float = 0.5                # perturbation of starts 2, 3, ... (in the squashed variables)

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
        if self.mode == "local":
            if not self.offsets_hz:
                raise ValueError("local robustness needs at least one offset")
            if not self.orders and not self.b1 and self.gate == "none":
                raise ValueError("choose at least one goal (an order, B1 or a gate)")
        else:
            if self.gate == "none":
                raise ValueError("broadband optimization needs a UR or PP gate goal")
            if self.bb_offset_n < 1 or self.bb_scale_n < 1:
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
        if self.gate == "PP" and (np.linalg.norm(self.m0) == 0 or np.linalg.norm(self.m_target) == 0):
            raise ValueError("PP vectors must not be zero")
        if self.gate == "UR" and np.linalg.norm(self.target_axis) == 0:
            raise ValueError("the target rotation axis must not be zero")

    def variant(self, i):
        """Configuration of start i of a multi-start (other seed; larger perturbation for i > 0)."""
        if i == 0:
            return self
        d = self.to_dict()
        d.update(seed=int(self.seed) + i, start_noise=self.multistart_noise)
        return GrapeConfig.from_dict(d)

    def ensemble(self):
        """(offset_hz, scale) pairs of the broadband grid."""
        offs = np.linspace(self.bb_offset_min_hz, self.bb_offset_max_hz, int(self.bb_offset_n))
        scales = np.linspace(self.bb_scale_min, self.bb_scale_max, int(self.bb_scale_n))
        return [(float(o), float(s)) for o in offs for s in scales]


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


# ----------------------------------------------------------------- problem ---

INFO_KEYS = ("mean infidelity", "max amplitude [%]", "max amplitude rate [%/µs]", "max phase rate [kHz]", "roughness")
SOFT_SCHEDULE = (1.0, 0.1, 0.01, 1e-3, 1e-4, 0.0)      # weights of the soft terms in successive stages
LIMIT_EXCESS_TOL = 1e-3
SOFT_STAGE_TOL = 1e-6                                   # xtol / ftol / gtol of the soft stages                                 # rate limits count as kept below this relative excess


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
        self.soft_scale = 1.0
        D = sparse.diags([-np.ones(self.N - 1), np.ones(self.N - 1)], [0, 1], shape=(self.N - 1, self.N))
        self.D2 = sparse.block_diag([D, D], format="csr")             # first differences of u_x and u_y
        w0 = self.unsquash(self.x0)
        if config.controls == "spline":
            B = spline_basis(self.N, int(config.n_spline))
            self.basis = sparse.block_diag([sparse.csr_matrix(B), sparse.csr_matrix(B)], format="csr")
            c = np.linalg.lstsq(B, np.column_stack([w0[:self.N], w0[self.N:]]), rcond=None)[0]
            self.p0 = np.r_[c[:, 0], c[:, 1]]
        else:
            self.basis = None
            self.p0 = w0

    def controls(self, x):
        return x[:self.N], x[self.N:]

    def pulse(self, x) -> PulseSequence:
        ux, uy = self.controls(x)
        amp = np.hypot(ux, uy)
        return PulseSequence(np.minimum(amp, 1.0) * 100.0, np.degrees(np.unwrap(np.arctan2(uy, ux))), self.dt, self.umax)

    def has_soft_terms(self):
        c = self.cfg
        return (c.regularization == "stay_close" and c.regularization_weight > 0) or c.smooth_weight > 0

    # ------------------------------------------------- parameters -> controls ---
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
    def _gate_rows(self, R, dR, weight):
        if self.cfg.gate == "UR":
            return weight * (R - self.Rt).ravel() / 2, weight * dR.reshape(len(dR), 9).T / 2
        m = R @ self.m0
        return weight * (m - self.mt) / 2, weight * (dR @ self.m0).T / 2

    def goals(self, x, with_jac=True):
        """Goal residuals (dense) and Jacobian (rows x 2N), plus a summary dict."""
        ux, uy = self.controls(x)
        cfg, st = self.cfg, self.st
        rows, jac, summ = [], [], {}
        if cfg.mode == "local":
            K = max(cfg.orders) if cfg.orders else 0
            for off in cfg.offsets_hz:
                X, dX = propagator_and_gradient(ux, uy, self.dt, self.umax, off, 1.0, max(K, 1))
                if cfg.orders:
                    mg = magnus_with_gradient(X, dX)
                    for k in cfg.orders:
                        A, dA = mg[k - 1]
                        if cfg.gate == "PP":
                            A, dA = np.cross(A, self.m0), np.cross(dA, self.m0)
                        w = st.weight(k) / self.T ** k
                        rows.append(w * A); jac.append(w * dA.T)
                        key = f"|A{k}|/T^{k}" + (" ⊥m0" if cfg.gate == "PP" else "")
                        summ[key] = max(summ.get(key, 0.0), float(np.linalg.norm(A) / self.T ** k))
                if cfg.b1:
                    Xa, dXa = propagator_and_gradient(ux, uy, self.dt, self.umax, off, 1.0, 1, "amplitude")
                    U0d = Xa[0].conj().T
                    B = _vec(U0d @ Xa[1], 1)
                    dB = _vec(np.conj(np.swapaxes(dXa[:, 0], -1, -2)) @ Xa[1] + U0d @ dXa[:, 1], 1)
                    if cfg.gate == "PP":
                        B, dB = np.cross(B, self.m0), np.cross(dB, self.m0)
                    norm = 2 * math.pi * self.umax * self.T             # maximal nutation angle
                    w = st.weight("b1") / norm
                    rows.append(w * B); jac.append(w * dB.T)
                    summ["|B1|/(Ω T)"] = max(summ.get("|B1|/(Ω T)", 0.0), float(np.linalg.norm(B) / norm))
                if cfg.gate != "none":
                    R, dR = rotation_with_gradient(X[0], dX[:, 0])
                    r, j = self._gate_rows(R, dR, st.weight("gate"))
                    rows.append(r); jac.append(j)
                    err = (rotation_angle_deg(self.Rt.T @ R) if cfg.gate == "UR"
                           else _angle_between_deg(R @ self.m0, self.mt))
                    summ["gate error [°]"] = max(summ.get("gate error [°]", 0.0), err)
        else:
            ens = cfg.ensemble()
            w = st.weight("gate") / math.sqrt(len(ens))
            inf = []
            for off, sc in ens:
                X, dX = propagator_and_gradient(ux, uy, self.dt, self.umax, off, sc, 0)
                R, dR = rotation_with_gradient(X[0], dX[:, 0])
                r, j = self._gate_rows(R, dR, w)
                rows.append(r); jac.append(j)
                if cfg.gate == "UR":
                    inf.append(1 - (3 + np.trace(self.Rt.T @ R)) / 6)
                else:
                    inf.append(1 - (1 + float(R @ self.m0 @ self.mt)) / 2)
            summ["mean infidelity"] = float(np.mean(inf))
            summ["max infidelity"] = float(np.max(inf))
        r = np.concatenate(rows) if rows else np.zeros(0)
        Jm = np.vstack(jac) if jac else np.zeros((0, 2 * self.N))
        return r, Jm, summ

    def reached(self, summ):
        st = self.st
        for k, v in summ.items():
            if k in INFO_KEYS:
                continue
            if k == "gate error [°]":
                if v >= st.gate_tolerance_deg:
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
            w = cfg.limit_weight
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
        return r, Jg, re, Je, summ, float(np.linalg.norm(r))

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


def levenberg_marquardt(fun, p, max_nfev, xtol, ftol, gtol):
    """Minimize |r(p)|^2 with r, Jacobians from fun(p) -> (rg, G dense, rs, S sparse). Every call of fun is
    one evaluation (residuals and exact Jacobian together). Returns (p, status message)."""
    rg, G, rs, S = fun(p)
    cost = rg @ rg + rs @ rs
    scale = (np.sum(G * G) + (S.multiply(S).sum() if S is not None else 0.0)) / max(len(p), 1)
    lam = 1e-2
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
            if small_step or small_gain:
                return p, "step / improvement tolerance reached"
        else:
            lam *= 4.0
            if lam > 1e16:
                return p, "no further improvement possible"
    return p, "evaluation limit reached"


def optimize(ps: PulseSequence, config: GrapeConfig, settings: SolverSettings = None, callback=None,
             snapshot=None, snapshot_interval_s=1.0, start_index=0):
    """Run the optimization from the start pulse ``ps`` (start ``start_index`` of a multi-start: other
    seed and, for index > 0, a larger start perturbation). ``callback(nfev, goal_residual)`` returning True
    stops it; ``snapshot(pulse, goal_residual)`` receives the best pulse every ``snapshot_interval_s``.
    Solver: Levenberg-Marquardt (lm_step) on the parameters. With soft terms (stay close, smoothness) it
    runs in stages with decreasing soft weights (SOFT_SCHEDULE), each started from the previous result,
    until the goals are reached. Returns (result pulse, info dict)."""
    st = settings or SolverSettings()
    cfg = config.variant(start_index)
    prob = GrapeProblem(ps, cfg, st)
    t0 = time.time()
    last_snap = [t0]
    p = prob.p0.copy()
    if cfg.start_noise > 0:
        # an exactly symmetric start (e.g. a constant-phase pulse) can keep the solver in a symmetric subspace
        p = p + np.random.default_rng(cfg.seed).normal(0.0, cfg.start_noise, p.shape)
    count = {"nfev": 0, "history": []}
    tols = [t if t and t > 0 else 0.0 for t in (st.xtol, st.ftol, st.gtol)]
    schedule = SOFT_SCHEDULE if prob.has_soft_terms() else (1.0,)
    reason, message = "converged", ""
    for factor in schedule:
        prob.soft_scale = factor
        final = factor == schedule[-1]
        best = {"p": p.copy(), "cost": np.inf, "goal": np.inf, "summ": None}

        def fun(q):
            if count["nfev"] >= st.max_nfev:
                raise _Stop("limit")
            if st.time_limit_s > 0 and time.time() - t0 > st.time_limit_s:
                raise _Stop("time")
            u, Mp = prob.to_controls(q)
            rg, Jg, rs, Js, summ, goal = prob.evaluate_blocks(u)
            G = np.asarray(Jg @ Mp) if sparse.issparse(Mp) else Jg @ Mp
            S = sparse.csr_matrix(Js @ Mp) if Js.shape[0] else None
            cost = float(np.sqrt(rg @ rg + rs @ rs))
            count["nfev"] += 1
            if cost < best["cost"]:
                best.update(p=np.array(q), cost=cost, goal=goal, summ=summ)
            count["history"].append(best["goal"])
            if final and st.stop_when_solved and prob.reached(summ):
                best.update(p=np.array(q), cost=cost, goal=goal, summ=summ)
                raise _Stop("solved")
            if callback is not None and callback(count["nfev"], best["goal"]):
                raise _Stop("user")
            if snapshot is not None and time.time() - last_snap[0] > snapshot_interval_s:
                last_snap[0] = time.time()
                snapshot(prob.pulse(prob.to_controls(best["p"])[0]), best["goal"])
            return rg, G, rs, (S if S is not None else sparse.csr_matrix((0, len(q))))

        remaining = st.max_nfev - count["nfev"]
        if final:
            stage_nfev, stage_tols = remaining, tols
        else:                     # soft stages: a share of the budget, end when the progress stalls
            stage_nfev = max(10, st.max_nfev // (2 * (len(schedule) - 1)))
            stage_tols = [max(t, SOFT_STAGE_TOL) for t in tols]
        try:
            _, message = levenberg_marquardt(fun, p, min(stage_nfev, remaining), *stage_tols)
            reason = "limit" if message == "evaluation limit reached" and (final or count["nfev"] >= st.max_nfev) else "converged"
        except _Stop as s:
            reason = s.reason
            message = {"solved": "all goals reached", "limit": "evaluation limit reached",
                       "time": "time limit reached", "user": "stopped by the user"}[s.reason]
        p = best["p"]
        if reason != "converged":
            break
        if best["summ"] is not None and prob.reached(best["summ"]) and (st.stop_when_solved or final):
            reason, message = "solved", "all goals reached"
            break
    u = prob.to_controls(p)[0]
    prob.soft_scale = schedule[-1]
    _, _, summ, goal = prob.evaluate(u)
    start = prob.evaluate(prob.to_controls(prob.p0)[0])[2]
    return prob.pulse(u), {"stop_reason": reason, "message": message, "reached": prob.reached(summ),
                           "values": summ, "start_summary": start, "nfev": count["nfev"], "goal_residual": goal,
                           "elapsed_s": time.time() - t0, "history": count["history"], "start_index": start_index,
                           "cancelled": reason == "user"}


def selection_key(config: GrapeConfig, info):
    """Order of multi-start results: reached first; among those the smoothest (if smoothing is on), else
    the lowest goal residual."""
    second = info["values"].get("roughness", 0.0) if config.smooth_weight > 0 and info["reached"] else info["goal_residual"]
    return (not info["reached"], second)


def _run_start(args):
    amp, phase, dt, umax, cfg, st, i = args
    ps, info = optimize(PulseSequence(np.asarray(amp), np.asarray(phase), dt, umax), GrapeConfig.from_dict(cfg),
                        SolverSettings.from_dict(st), start_index=i)
    return ps.amplitude, ps.phase, info


def optimize_multistart(ps: PulseSequence, config: GrapeConfig, settings: SolverSettings = None, processes=1):
    """All ``config.n_starts`` starts (in parallel processes if ``processes`` > 1); returns the best pulse
    (selection_key), its info and the infos of all starts."""
    st = settings or SolverSettings()
    jobs = [(ps.amplitude, ps.phase, ps.dt, ps.umax, config.to_dict(), st.to_dict(), i) for i in range(int(config.n_starts))]
    if processes > 1 and len(jobs) > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(min(processes, len(jobs))) as pool:
            results = pool.map(_run_start, jobs)
    else:
        results = [_run_start(j) for j in jobs]
    amp, phase, info = min(results, key=lambda r: selection_key(config, r[2]))
    return PulseSequence(amp, phase, ps.dt, ps.umax), info, [r[2] for r in results]


def run_grape_job(job, queue, stop_event):
    """Process entry point for one start: job = {"pulse": [amp, phase, dt, umax], "config", "settings",
    "start_index"}. Puts ("progress", i, nfev, goal_residual), ("snapshot", i, amp, phase, goal_residual),
    ("done", i, amp, phase, info) or ("error", i, msg)."""
    i = int(job.get("start_index", 0))
    try:
        amp, phase, dt, umax = job["pulse"]
        ps = PulseSequence(np.asarray(amp), np.asarray(phase), dt, umax)
        last = [0]

        def callback(nfev, cost):
            if nfev - last[0] >= 5 or nfev == 1:
                last[0] = nfev
                queue.put(("progress", i, nfev, cost))
            return stop_event.is_set()

        def snapshot(p, cost):
            queue.put(("snapshot", i, p.amplitude.tolist(), p.phase.tolist(), cost))

        result, info = optimize(ps, GrapeConfig.from_dict(job["config"]), SolverSettings.from_dict(job.get("settings")),
                                callback, snapshot, start_index=i)
        info["history"] = info["history"][-5000:]
        queue.put(("done", i, result.amplitude.tolist(), result.phase.tolist(), info))
    except Exception as e:
        queue.put(("error", i, f"{type(e).__name__}: {e}"))
