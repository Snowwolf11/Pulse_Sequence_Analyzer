"""Parametrized curve and pulse families for the editor.

Curve families return points of a curve (any units); the editor converts them with
``conversion.curve_to_ps``. Pulse families (2H, 4H, 6H) return the pulse sequence directly from the
closed-form construction (phase ramps), so no conversion is needed.
"""
from __future__ import annotations

import ast
import math
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
from scipy.optimize import brentq

from psa.design.pulse import PulseSequence


@dataclass
class Param:
    name: str
    label: str
    default: float
    minimum: float
    maximum: float
    integer: bool = False
    unit: str = ""


class Family:
    name = ""
    description = ""
    kind = "curve"                      # "curve" or "pulse"
    params: List[Param] = []

    def defaults(self) -> Dict[str, float]:
        return {p.name: p.default for p in self.params}

    def check(self, values):
        out = self.defaults()
        out.update(values or {})
        for p in self.params:
            v = min(p.maximum, max(p.minimum, float(out[p.name])))
            out[p.name] = int(round(v)) if p.integer else v
        return out

    def points(self, values, n=4000) -> np.ndarray:
        raise NotImplementedError

    def pulse(self, values, dt, umax) -> PulseSequence:
        raise NotImplementedError


# ---------------------------------------------------------- curve families ---

class Helix(Family):
    name = "Helix"
    description = "x = R cos(2πu), y = R sin(2πu), z = h·u,  u ∈ [0, turns].  Constant curvature and torsion."
    params = [Param("R", "radius R", 1.0, 0.05, 10.0), Param("h", "pitch h (per turn)", 1.0, -20.0, 20.0),
              Param("turns", "turns", 2.0, 0.05, 50.0)]

    def points(self, values, n=4000):
        v = self.check(values)
        u = np.linspace(0, v["turns"], n)
        return np.c_[v["R"] * np.cos(2 * np.pi * u), v["R"] * np.sin(2 * np.pi * u), v["h"] * u]


class LissajousKnot(Family):
    name = "Lissajous knot"
    description = ("x = cos(nx·t + φx), y = cos(ny·t + φy), z = cos(nz·t + φz),  t ∈ [0, 2π].  "
                   "Closed for integer frequencies; a knot if nx, ny, nz are pairwise coprime.")
    params = [Param("nx", "nx", 3, 1, 20, True), Param("ny", "ny", 2, 1, 20, True), Param("nz", "nz", 7, 1, 20, True),
              Param("phx", "φx", 0.7, -math.pi, math.pi, unit="rad"), Param("phy", "φy", 0.2, -math.pi, math.pi, unit="rad"),
              Param("phz", "φz", 0.0, -math.pi, math.pi, unit="rad")]

    def points(self, values, n=6000):
        v = self.check(values)
        t = np.linspace(0, 2 * np.pi, n)
        return np.c_[np.cos(v["nx"] * t + v["phx"]), np.cos(v["ny"] * t + v["phy"]), np.cos(v["nz"] * t + v["phz"])]


class TorusKnot(Family):
    name = "Torus knot"
    description = ("(p, q) torus knot: x = (R + r cos qt) cos pt, y = (R + r cos qt) sin pt, z = r sin qt,  "
                   "t ∈ [0, 2π].  A knot if p and q are coprime.")
    params = [Param("p", "p", 2, 1, 20, True), Param("q", "q", 3, 1, 20, True),
              Param("R", "torus radius R", 2.0, 0.1, 10.0), Param("r", "tube radius r", 0.8, 0.01, 10.0)]

    def points(self, values, n=6000):
        v = self.check(values)
        t = np.linspace(0, 2 * np.pi, n)
        rad = v["R"] + v["r"] * np.cos(v["q"] * t)
        return np.c_[rad * np.cos(v["p"] * t), rad * np.sin(v["p"] * t), v["r"] * np.sin(v["q"] * t)]


class CleliaCurve(Family):
    name = "Clelia curve"
    description = "Curve on the unit sphere: x = sin(n t) cos(m t), y = sin(n t) sin(m t), z = cos(n t),  t ∈ [0, T]."
    params = [Param("m", "m", 7.0, 0.1, 30.0), Param("n", "n", 3.342, 0.1, 30.0), Param("T", "T", 8.0, 0.1, 50.0)]

    def points(self, values, n=6000):
        v = self.check(values)
        t = np.linspace(0, v["T"], n)
        return np.c_[np.sin(v["n"] * t) * np.cos(v["m"] * t), np.sin(v["n"] * t) * np.sin(v["m"] * t), np.cos(v["n"] * t)]


# ----------------------------------------------------------- custom family ---

SAFE_FUNCTIONS = {name: getattr(np, name) for name in
                  ("sin", "cos", "tan", "arcsin", "arccos", "arctan", "arctan2", "sinh", "cosh", "tanh",
                   "exp", "log", "log10", "sqrt", "abs", "sign", "floor", "ceil", "minimum", "maximum")}
SAFE_FUNCTIONS.update({"asin": np.arcsin, "acos": np.arccos, "atan": np.arctan, "atan2": np.arctan2})
SAFE_CONSTANTS = {"pi": np.pi, "e": np.e}
_ALLOWED_NODES = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Call, ast.Name, ast.Load, ast.Constant,
                  ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod, ast.USub, ast.UAdd, ast.FloorDiv)


def compile_expression(text, variables):
    """Compile a formula like ``a*cos(u) + u**2`` that may only use + - * / ** %, numbers, the
    functions in SAFE_FUNCTIONS, pi, e and the given variable names. Raises ValueError otherwise."""
    text = text.strip().replace("^", "**")
    if not text:
        raise ValueError("empty expression")
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as e:
        raise ValueError(f"syntax error in '{text}': {e.msg}") from None
    allowed_names = set(SAFE_FUNCTIONS) | set(SAFE_CONSTANTS) | set(variables)
    called = {id(n.func) for n in ast.walk(tree) if isinstance(n, ast.Call)}
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and node.id in SAFE_FUNCTIONS and id(node) not in called and node.id not in variables:
            raise ValueError(f"function '{node.id}' must be called, e.g. {node.id}(u)")
        if not isinstance(node, _ALLOWED_NODES):
            raise ValueError(f"'{type(node).__name__}' is not allowed in '{text}'")
        if isinstance(node, ast.Name) and node.id not in allowed_names:
            raise ValueError(f"unknown name '{node.id}' in '{text}'")
        if isinstance(node, ast.Call) and not (isinstance(node.func, ast.Name) and node.func.id in SAFE_FUNCTIONS):
            raise ValueError(f"only the functions {sorted(SAFE_FUNCTIONS)} may be called")
        if isinstance(node, ast.Constant) and not isinstance(node.value, (int, float)):
            raise ValueError(f"only numbers are allowed as constants in '{text}'")
    code = compile(tree, "<expression>", "eval")

    def evaluate(**values):
        env = {"__builtins__": {}}
        env.update(SAFE_FUNCTIONS); env.update(SAFE_CONSTANTS); env.update(values)
        return eval(code, env)
    return evaluate


class CustomCurve(Family):
    name = "Custom (formula)"
    description = ("Type x(u), y(u), z(u) with parameters a, b, c, d. Allowed: + - * / ** (or ^), numbers, pi, e, "
                   "sin cos tan asin acos atan atan2 sinh cosh tanh exp log log10 sqrt abs sign minimum maximum.")
    params = [Param("u0", "u start", 0.0, -1000.0, 1000.0), Param("u1", "u end", 2 * math.pi, -1000.0, 1000.0),
              Param("a", "a", 1.0, -20.0, 20.0), Param("b", "b", 2.0, -20.0, 20.0),
              Param("c", "c", 3.0, -20.0, 20.0), Param("d", "d", 0.5, -20.0, 20.0)]
    default_expressions = ("cos(a*u)", "sin(b*u)", "d*sin(c*u)")

    def points(self, values, n=4000, expressions=None):
        v = self.check(values)
        exprs = expressions or self.default_expressions
        u = np.linspace(v["u0"], v["u1"], n)
        cols = []
        for e in exprs:
            val = compile_expression(e, ("u", "a", "b", "c", "d"))(u=u, a=v["a"], b=v["b"], c=v["c"], d=v["d"])
            val = np.broadcast_to(np.asarray(val, dtype=np.float64), u.shape)
            if not np.all(np.isfinite(val)):
                raise ValueError(f"'{e}' gives non-finite values on [u start, u end]")
            cols.append(val)
        return np.column_stack(cols)


# ----------------------------------------------------------- pulse families ---

def phase_ramp(phase_start_deg, rate_deg_per_s, n, dt):
    """Phases of a linear ramp sampled at the sub-pulse midpoints."""
    return phase_start_deg + rate_deg_per_s * dt * (np.arange(n) + 0.5)


class TwoHelix(Family):
    name = "2H (closed two-helix loop)"
    kind = "pulse"
    description = ("Two phase ramps of duration τ at full amplitude whose error curve closes (first-order robust). "
                   "Closed form: T from Ω τ = T / sqrt(1 − sin T / T), T ∈ (π, 2π); ν = sqrt(−sin T / T); "
                   "ramp rate −Ω ν; phase jump Δφ(T) between the ramps.  τ ∈ (π/Ω, 2π/Ω).")
    params = [Param("tau", "ramp duration τ", 75.0, 50.5, 99.5, unit="µs (at 10 kHz)"),
              Param("phi0", "start phase φ0", 0.0, -180.0, 180.0, unit="°")]

    @staticmethod
    def parameters(Omega, tau):
        """Closed-form 2H parameters (Omega in rad/s, tau in s): ramp rate (deg/s), phase jump (deg), T."""
        f = lambda T: T / math.sqrt(1 - math.sin(T) / T) - Omega * tau
        T = brentq(f, math.pi + 1e-12, 2 * math.pi - 1e-12)
        nu = math.sqrt(-math.sin(T) / T)
        jump = math.degrees(math.atan2(-2 * nu * math.sqrt(1 + nu ** 2) * math.sin(T), -(1 + (1 + 2 * nu ** 2) * math.cos(T))))
        rate = math.degrees(-Omega * nu)
        return rate, jump, T

    def pulse(self, values, dt, umax):
        v = self.check(values)
        Omega = 2 * math.pi * umax
        tau_lo, tau_hi = math.pi / Omega, 2 * math.pi / Omega
        n = max(1, int(round(v["tau"] * 1e-6 / dt)))
        tau = min(max(n * dt, tau_lo * (1 + 1e-6)), tau_hi * (1 - 1e-6))
        rate, jump, _ = self.parameters(Omega, tau)
        ph1 = phase_ramp(v["phi0"], rate, n, tau / n)
        ph2 = phase_ramp(v["phi0"] + rate * tau + jump, rate, n, tau / n)
        return PulseSequence(np.full(2 * n, 100.0), np.r_[ph1, ph2], dt, umax)

    def tau_range_us(self, umax):
        return 0.5e6 / umax, 1e6 / umax


class FourHelix(Family):
    name = "4H (0° UR, 2nd-order robust)"
    kind = "pulse"
    description = ("Four phase ramps PR(φ0 → φ0+φ̇τ), PR(φ0+φ̇τ → φ0), and the same shifted by 180°; φ̇ is the first "
                   "root of Ω² sin(τ√(Ω²+φ̇²)) + τ φ̇² √(Ω²+φ̇²) = 0 (curve closure). Identity gate, "
                   "first- and second-order robust against detuning. Solutions exist for Ωτ ∈ (π, 2π) + 2πk, "
                   "i.e. τ ∈ (50, 100) µs, (150, 200) µs, … at 10 kHz.")
    params = [Param("tau", "ramp duration τ", 60.0, 50.5, 999.5, unit="µs"),
              Param("phi0", "start phase φ0", 0.0, -180.0, 180.0, unit="°")]

    @staticmethod
    def phase_rate(Omega, tau):
        """First positive root φ̇ (rad/s) of the closure condition."""
        def f(p):
            w = math.sqrt(Omega ** 2 + p ** 2)
            return Omega ** 2 * math.sin(tau * w) + tau * p ** 2 * w
        grid = np.linspace(0.0, 50 * Omega, 20001)[1:]
        vals = np.array([f(p) for p in grid])
        sign = np.flatnonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))
        if not len(sign):
            k = Omega * tau / (2 * math.pi)
            raise ValueError(f"no closing phase-ramp rate for τ = {tau * 1e6:.2f} µs (Ωτ/2π = {k:.3f}); "
                             "4H solutions need Ωτ/2π in (k + 1/2, k + 1)")
        i = sign[0]
        return brentq(f, grid[i], grid[i + 1])

    def pulse(self, values, dt, umax):
        v = self.check(values)
        Omega = 2 * math.pi * umax
        n = max(1, int(round(v["tau"] * 1e-6 / dt)))
        tau = n * dt
        rate = math.degrees(self.phase_rate(Omega, tau))
        p0 = v["phi0"]
        ramps = [phase_ramp(p0, rate, n, dt), phase_ramp(p0 + rate * tau, -rate, n, dt),
                 phase_ramp(p0 + 180, rate, n, dt), phase_ramp(p0 + 180 + rate * tau, -rate, n, dt)]
        return PulseSequence(np.full(4 * n, 100.0), np.concatenate(ramps), dt, umax)


# Solutions of the continuous 6H problem (loop B a closed 2H loop): (T_A, Omega_B / Omega_A, T_B, twist A -> B
# in rad). Found by numerical continuation from the 6H pulse library (6H_Full_Derivation, T_A = 3.3 ... 4.6);
# the branch ends in a fold at T_A ≈ 4.65. Used as starting points of the solve.
SIX_H_TABLE = (
    (3.18, 0.31677394, 5.88125475, -0.14832762),
    (3.19, 0.34689011, 5.88134094, -0.15917521),
    (3.20, 0.37304846, 5.88143718, -0.16754823),
    (3.21, 0.39630722, 5.88154272, -0.17417659),
    (3.22, 0.41733371, 5.88165700, -0.17950560),
    (3.23, 0.43657914, 5.88177959, -0.18382710),
    (3.24, 0.45436457, 5.88191014, -0.18734266),
    (3.25, 0.47092756, 5.88204838, -0.19019713),
    (3.26, 0.48644939, 5.88219405, -0.19249793),
    (3.27, 0.50107194, 5.88234697, -0.19432674),
    (3.28, 0.51490853, 5.88250696, -0.19574694),
    (3.29, 0.52805131, 5.88267388, -0.19680853),
    (3.30, 0.54057634, 5.88284759, -0.19755154),
    (3.35, 0.59585520, 5.88381476, -0.19743864),
    (3.40, 0.64218619, 5.88494014, -0.19241365),
    (3.45, 0.68232469, 5.88621876, -0.18348138),
    (3.50, 0.71788730, 5.88764825, -0.17100934),
    (3.55, 0.74991511, 5.88922795, -0.15504575),
    (3.60, 0.77912452, 5.89095829, -0.13545713),
    (3.65, 0.80603402, 5.89284055, -0.11199904),
    (3.70, 0.83103429, 5.89487660, -0.08435894),
    (3.75, 0.85442979, 5.89706883, -0.05218534),
    (3.80, 0.87646487, 5.89941999, -0.01510846),
    (3.85, 0.89734086, 5.90193315, 0.02724669),
    (3.90, 0.91722786, 5.90461169, 0.07525549),
    (3.95, 0.93627300, 5.90745916, 0.12929835),
    (4.00, 0.95460662, 5.91047934, 0.18977522),
    (4.05, 0.97234698, 5.91367618, 0.25713068),
    (4.10, 0.98960401, 5.91705376, 0.33188968),
    (4.15, 1.00648258, 5.92061626, 0.41470571),
    (4.20, 1.02308530, 5.92436801, 0.50642571),
    (4.25, 1.03951542, 5.92831336, 0.60818069),
    (4.30, 1.05587964, 5.93245676, 0.72152077),
    (4.35, 1.07229151, 5.93680269, 0.84863174),
    (4.40, 1.08887528, 5.94135563, 0.99271322),
    (4.45, 1.10577104, 5.94612008, 1.15870875),
    (4.50, 1.12314147, 5.95110052, 1.35491599),
    (4.55, 1.14118152, 5.95630134, 1.59731920),
    (4.60, 1.16013261, 5.96172692, 1.92628693),
    (4.64, 1.17615466, 5.96623207, 2.37905204),
)


def two_helix_jump(T):
    """Phase jump (rad) between the two ramps of a closed 2H loop with helix angle T (6H_curves.pdf)."""
    nu = math.sqrt(-math.sin(T) / T)
    return math.atan2(-2 * nu * math.sqrt(1 + nu ** 2) * math.sin(T), -(1 + (1 + 2 * nu ** 2) * math.cos(T)))


class SixHelix(Family):
    name = "6H (0° UR, 3rd-order robust)"
    kind = "pulse"
    description = ("Three closed 2H loops in the order A₁ B₁ B₂ A₂ C₁ C₂ as six phase ramps; the second half is the "
                   "first shifted by 180° (C = B). Loop A: amplitude Ω_A, helix angle T_A, rate −Ω_A ν(T_A), "
                   "ν = sqrt(−sin T / T); loop B: rate +Ω_B ν(T_B) with the 2H phase jump. Ω_B/Ω_A, T_B and the twist "
                   "A → B are solved for A2 = A3 = 0 and the identity gate (A1 = 0 by closure). The sampled pulse is "
                   "then polished on the sub-pulse grid, so it is exact as shown (~1 s). T_A ∈ (π, 4.65); "
                   "shortest at T_A ≈ 4.15 (≈ 490 µs at 10 kHz).")
    params = [Param("TA", "loop A angle T_A (shape)", 4.15, 3.18, 4.62, unit="rad"),
              Param("amp", "max. amplitude", 100.0, 1.0, 100.0, unit="%"),
              Param("phi0", "start phase φ0", 0.0, -180.0, 180.0, unit="°")]

    def __init__(self):
        self._cache = {}

    @staticmethod
    def _nu(T):
        return math.sqrt(-math.sin(T) / T)

    def _loops(self, TA, z):
        """(B rate, B amplitude, B ramp duration, jump inside B, twist A -> B) from (Omega_B, T_B, twist), Omega_A = 1."""
        oB, TB, j1 = z
        nu = self._nu(TB)
        return oB * nu, oB, TB / (oB * math.sqrt(1 + nu ** 2)), -two_helix_jump(TB), j1

    def _pieces(self, TA, x):
        """Continuous model (Omega_A = 1) as builder pieces."""
        from psa.design.builder import Piece
        rB, oB, tB, dB, j1 = x
        nu = self._nu(TA)
        rA, tA = -nu, TA / math.sqrt(1 + nu ** 2)
        out = []
        for half in (0, 1):
            ph = -math.pi * half
            for k, (O, r, t, jump) in enumerate(((1.0, rA, tA, j1), (oB, rB, tB, dB), (oB, rB, tB, 0.0))):
                out.append(Piece(O, ph, r, t, 3 * half + k))
                ph += r * t + jump
        return out

    def shape(self, TA):
        """Continuous solution (B rate, B amplitude, B ramp duration, jump inside B, twist A -> B), Omega_A = 1."""
        from scipy.optimize import least_squares
        from psa.design.builder import CurveDesign
        d = CurveDesign()

        def residual(z):
            p = self._pieces(TA, self._loops(TA, z))
            rb, T = d.robustness(p), d.duration(p)
            R = d.gate_matrix(p)
            return np.r_[rb["A1"] / T, rb["A2"] / T ** 2, rb["A3"] / T ** 3, R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]

        if not math.pi < TA < 2 * math.pi:
            raise ValueError("T_A must lie in (π, 2π)")
        table = np.array(SIX_H_TABLE)
        for i in np.argsort(np.abs(table[:, 0] - TA))[:3]:          # nearest seeds first
            sol = least_squares(residual, table[i, 1:], xtol=1e-15, ftol=1e-15, gtol=1e-15)
            if np.linalg.norm(sol.fun) < 1e-11 and sol.x[0] > 0 and math.pi < sol.x[1] < 2 * math.pi:
                return np.array(self._loops(TA, sol.x))
        raise ValueError(f"no 6H solution found for T_A = {TA:.4f} (the family ends at T_A ≈ 4.65)")

    @staticmethod
    def _sample(y, nA, nB, dt, umax, phi0_rad):
        OA, rA, OB, rB, dB, j1 = y
        Om = 2 * math.pi * umax
        amp, ph = [], []
        for half in (0, 1):
            p = phi0_rad - math.pi * half
            for O, r, n, jump in ((OA, rA, nA, j1), (OB, rB, nB, dB), (OB, rB, nB, 0.0)):
                amp.append(np.full(n, O / Om * 100.0))
                ph.append(p + r * dt * (np.arange(n) + 0.5))
                p += r * n * dt + jump
        return PulseSequence(np.concatenate(amp), np.degrees(np.concatenate(ph)), dt, umax)

    def pulse(self, values, dt, umax):
        from scipy.optimize import least_squares
        from psa.design.metrics import analyze
        v = self.check(values)
        key = (round(v["TA"], 12), round(v["amp"], 12), dt, umax)
        if key not in self._cache:
            rB, oB, tB, dB, j1 = self.shape(v["TA"])
            OA = 2 * math.pi * umax * v["amp"] / 100.0 / max(1.0, oB)      # the larger loop gets the max. amplitude
            nu = self._nu(v["TA"])
            tA, tBs = v["TA"] / math.sqrt(1 + nu ** 2) / OA, tB / OA
            nA, nB = max(2, math.ceil(tA / dt - 1e-9)), max(2, math.ceil(tBs / dt - 1e-9))   # rounding up lowers the amplitudes
            base = np.array([OA * tA / (nA * dt), -OA * nu * tA / (nA * dt),
                             oB * OA * tBs / (nB * dt), rB * OA * tBs / (nB * dt), 1.0, 1.0])
            y0 = np.array([1.0, 1.0, 1.0, 1.0, dB, j1])

            def residual(y):                      # exact metrics of the sampled pulse
                ps = self._sample(y * base, nA, nB, dt, umax, 0.0)
                a = analyze(ps, b1=False)
                rb, T = a["robustness"], ps.duration
                return np.r_[rb["A1"] / T, rb["A2"] / T ** 2, rb["A3"] / T ** 3, a["gate"]["quaternion"][1:]]

            sol = least_squares(residual, y0, method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14)
            y = sol.x if np.linalg.norm(sol.fun) < np.linalg.norm(residual(y0)) else y0
            if len(self._cache) > 64:
                self._cache.clear()
            self._cache[key] = (y * base, nA, nB)
        y, nA, nB = self._cache[key]
        return self._sample(y, nA, nB, dt, umax, math.radians(v["phi0"]))


FAMILIES = [Helix(), LissajousKnot(), TorusKnot(), CleliaCurve(), CustomCurve(), TwoHelix(), FourHelix(), SixHelix()]


def family_by_name(name) -> Family:
    for f in FAMILIES:
        if f.name == name:
            return f
    raise KeyError(name)
