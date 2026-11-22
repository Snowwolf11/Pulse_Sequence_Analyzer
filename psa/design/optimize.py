"""Shared solver of the editor: adjust free variables (keyframe values, family parameters or builder
segment parameters) so that the chosen goals are met.

Goals (any combination):
    orders   detuning robustness: |A_k| / T^k -> 0 for k in orders (1 closed curve, 2 areas, 3 third order)
    target   optional target gate (rotation axis, angle in degrees): R -> R_target
    b1       first-order robustness against amplitude scaling: |B1| / int Omega dt -> 0

The residual is the concatenation of A_k / T^k, (R - R_target) and B1 / nutation angle; it is
minimized with scipy's least_squares inside bounds. ``callback(nfev, best_cost)`` may stop the run
(the best variables found so far are returned).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np
from scipy.optimize import least_squares

from psa.design import metrics
from psa.design.pulse import PulseSequence

TOLERANCE = 1e-8          # default: a goal counts as reached below this (normalized) value


@dataclass
class SolverSettings:
    """User-adjustable solver settings (editor: Settings → Solver)."""
    max_nfev: int = 2000                 # limit of ALL residual evaluations (incl. those for the derivatives)
    time_limit_s: float = 0.0            # stop after this many seconds (0 = no limit)
    stop_when_solved: bool = True        # stop as soon as every goal is below its "solved" threshold
    xtol: float = 1e-12                  # stop when the step is this small (relative)
    ftol: float = 1e-12                  # stop when the residual improves less than this (relative)
    gtol: float = 1e-12                  # stop when the gradient is this small
    diff_step: float = 0.0               # finite-difference relative step (0 = scipy default ~1.5e-8)
    solved_tolerance: float = TOLERANCE  # "solved": every robustness goal below this
    gate_tolerance_deg: float = 1e-5     # "solved": gate error below this
    weight_order1: float = 1.0
    weight_order2: float = 1.0
    weight_order3: float = 1.0
    weight_gate: float = 1.0
    weight_b1: float = 1.0

    def to_dict(self):
        return dict(self.__dict__)

    @classmethod
    def from_dict(cls, d):
        names = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in (d or {}).items() if k in names})

    def weight(self, key):
        return {1: self.weight_order1, 2: self.weight_order2, 3: self.weight_order3,
                "gate": self.weight_gate, "b1": self.weight_b1}[key]


def rotation_from_axis_angle(axis, angle_deg):
    a = np.asarray(axis, dtype=np.float64)
    n = np.linalg.norm(a)
    if n == 0:
        raise ValueError("the target rotation axis must not be zero")
    x, y, z = a / n
    th = math.radians(angle_deg)
    c, s, C = math.cos(th), math.sin(th), 1 - math.cos(th)
    return np.array([[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
                     [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
                     [z * x * C - y * s, z * y * C + x * s, c + z * z * C]])


def rotation_angle_deg(R):
    """Rotation angle of R; atan2 of the skew and trace parts is accurate also for tiny angles (arccos of
    the trace alone loses half of the digits there)."""
    w = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / 2
    return math.degrees(math.atan2(float(np.linalg.norm(w)), (float(np.trace(R)) - 1) / 2))


@dataclass
class Goals:
    orders: Tuple[int, ...] = (1,)
    target: Optional[Tuple[Tuple[float, float, float], float]] = None   # (axis, angle_deg)
    b1: bool = False

    def __post_init__(self):
        self.orders = tuple(int(k) for k in self.orders)
        if not self.orders and self.target is None and not self.b1:
            raise ValueError("choose at least one goal (an order, the target gate or B1 robustness)")

    def to_dict(self):
        return {"orders": list(self.orders), "target": None if self.target is None else [list(self.target[0]), self.target[1]],
                "b1": self.b1}

    @classmethod
    def from_dict(cls, d):
        t = d.get("target")
        return cls(tuple(d.get("orders", ())), None if t is None else (tuple(t[0]), float(t[1])), bool(d.get("b1", False)))


# ----------------------------------------------------------------- metrics ---

def pulse_metrics(ps: PulseSequence, b1=False) -> dict:
    """A_k, duration, Bloch rotation and (optionally) B1 vector of a pulse sequence (exact)."""
    a = metrics.analyze(ps, b1)
    rb = a["robustness"]
    out = {"A": {1: rb["A1"], 2: rb["A2"], 3: rb["A3"]}, "T": ps.duration, "R": a["gate"]["matrix"]}
    if b1:
        out.update(B1=a["B1"]["B1"], nutation=metrics._nutation(ps))
    return out


def residual(m: dict, goals: Goals, settings: SolverSettings = None) -> np.ndarray:
    s = settings or SolverSettings()
    T = m["T"]
    parts = []
    for k in goals.orders:
        parts.append(s.weight(k) * (m["A"][k] / T ** k if T > 0 else np.ones(3)))
    if goals.target is not None:
        parts.append(s.weight("gate") * (m["R"] - rotation_from_axis_angle(*goals.target)).ravel() / 2.0)
    if goals.b1:
        parts.append(s.weight("b1") * (m["B1"] / m["nutation"] if m.get("nutation") else np.ones(3)))
    return np.concatenate(parts) if parts else np.zeros(1)


def summary(m: dict, goals: Goals, settings: SolverSettings = None) -> dict:
    """Normalized values of all goals (for display) and whether each is reached."""
    T = m["T"]
    out = {}
    for k in goals.orders:
        out[f"|A{k}|/T^{k}"] = float(np.linalg.norm(m["A"][k]) / T ** k) if T > 0 else float("inf")
    if goals.target is not None:
        out["gate error [°]"] = rotation_angle_deg(rotation_from_axis_angle(*goals.target).T @ m["R"])
    if goals.b1:
        out["|B1|/∫Ω"] = float(np.linalg.norm(m["B1"]) / m["nutation"]) if m.get("nutation") else float("inf")
    st = settings or SolverSettings()
    reached = all((v < st.gate_tolerance_deg if k == "gate error [°]" else v < st.solved_tolerance) for k, v in out.items())
    return {"values": out, "reached": reached}


# ------------------------------------------------------------------ solver ---

class _Stop(Exception):
    def __init__(self, reason):
        super().__init__(reason)
        self.reason = reason


STOP_MESSAGES = {"solved": "all goals reached", "limit": "evaluation limit reached",
                 "time": "time limit reached", "user": "stopped by the user"}


def least_squares_goals(evaluate, x0, lo, hi, goals: Goals, callback=None, max_nfev=None, settings=None):
    """Minimize the goal residual over x in [lo, hi]. ``evaluate(x)`` returns a metrics dict.

    Stops (returning the best x seen) when scipy converges (xtol / ftol / gtol), when all goals are
    reached (``stop_when_solved``), after ``max_nfev`` evaluations in total (derivative evaluations
    included), after ``time_limit_s`` seconds, or when ``callback(nfev, best_cost)`` returns True.
    Returns (best x, info dict with "stop_reason" in solved / limit / time / user / converged)."""
    import time
    st = settings or SolverSettings()
    max_nfev = int(max_nfev or st.max_nfev)
    t0 = time.time()
    x0 = np.clip(np.asarray(x0, dtype=np.float64), lo, hi)
    best = {"x": x0.copy(), "cost": np.inf, "nfev": 0, "solved_x": None}

    def fun(x):
        if best["nfev"] >= max_nfev:
            raise _Stop("limit")
        if st.time_limit_s > 0 and time.time() - t0 > st.time_limit_s:
            raise _Stop("time")
        m = evaluate(x)
        r = residual(m, goals, st)
        cost = float(np.linalg.norm(r))
        best["nfev"] += 1
        if cost < best["cost"]:
            best.update(x=np.array(x, dtype=np.float64), cost=cost)
        if st.stop_when_solved and summary(m, goals, st)["reached"]:
            best.update(solved_x=np.array(x, dtype=np.float64), cost=cost)
            raise _Stop("solved")
        if callback is not None and callback(best["nfev"], best["cost"]):
            raise _Stop("user")
        return r

    reason, message, success = "converged", "", False
    try:
        # a tolerance of 0 switches that criterion off (scipy: None); scipy needs at least one of them
        tols = [t if t and t > 0 else None for t in (st.xtol, st.ftol, st.gtol)]
        if all(t is None for t in tols):
            tols[2] = 1e-15
        res = least_squares(fun, x0, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev,
                            xtol=tols[0], ftol=tols[1], gtol=tols[2], diff_step=st.diff_step or None)
        message, success = res.message, bool(res.success)
        if res.status == 0:                                   # scipy's own evaluation limit
            reason = "limit"
    except _Stop as stop:
        reason = stop.reason
        message = STOP_MESSAGES[reason]
        success = reason == "solved"
    x = best["solved_x"] if best["solved_x"] is not None else best["x"]
    return x, {"success": success, "cancelled": reason == "user", "stop_reason": reason, "message": message,
               "cost": best["cost"], "nfev": best["nfev"], "elapsed_s": time.time() - t0}


# ---------------------------------------------------------------- problems ---
# Each problem knows its variables (names, start values, bounds), evaluates metrics for a variable
# vector and writes a solution back into its design dictionary.

class WaveformProblem:
    """Free variables: values of selected keyframes ("amplitude"/"phase", keyframe sub-pulse index)."""

    def __init__(self, design: dict, free):
        from psa.design.waveform import Waveform
        self.design = design
        self.free = [(str(c), int(i)) for c, i in free]
        self.waveform = Waveform.from_dict(design)
        self.lookup = []
        for ch, idx in self.free:
            channel = getattr(self.waveform, ch)
            pos = channel.find(idx)
            if pos is None:
                raise ValueError(f"no {ch} keyframe at sub-pulse {idx}")
            self.lookup.append((channel, pos))

    def start(self):
        x0 = np.array([ch.keyframes[p].value for ch, p in self.lookup])
        lo = np.array([max(ch.minimum, -1e7) for ch, _ in self.lookup])
        hi = np.array([min(ch.maximum, 1e7) for ch, _ in self.lookup])
        return x0, lo, hi

    def apply(self, x):
        for (ch, p), v in zip(self.lookup, x):
            ch.keyframes[p].value = float(v)
        return self.waveform

    def evaluate(self, x, goals):
        return pulse_metrics(self.apply(x).to_pulse(), goals.b1)

    def solution(self, x):
        return self.apply(x).to_dict()


class FamilyProblem:
    """Free variables: selected (non-integer) family parameters. Curve families use a fixed duration
    so that the number of sub-pulses does not change during the solve."""

    def __init__(self, design: dict, free):
        from psa.design.families import family_by_name
        self.design = design
        self.family = family_by_name(design["family"])
        self.values = dict(self.family.check(design.get("values", {})))
        self.free = list(free)
        self.params = {p.name: p for p in self.family.params}
        for n in self.free:
            if n not in self.params or self.params[n].integer:
                raise ValueError(f"'{n}' cannot be a free variable")

    def start(self):
        x0 = np.array([self.values[n] for n in self.free])
        lo = np.array([self.params[n].minimum for n in self.free])
        hi = np.array([self.params[n].maximum for n in self.free])
        return x0, lo, hi

    def pulse(self, x):
        from psa.design.conversion import curve_to_ps
        values = dict(self.values)
        values.update({n: float(v) for n, v in zip(self.free, x)})
        d = self.design
        if self.family.kind == "pulse":
            return self.family.pulse(values, d["dt"], d["umax"])
        kwargs = {"expressions": d["expressions"]} if d["family"].startswith("Custom") else {}
        ps, _ = curve_to_ps(self.family.points(values, **kwargs), d["dt"], d["umax"], duration=d["duration"])
        return ps

    def evaluate(self, x, goals):
        return pulse_metrics(self.pulse(x), goals.b1)

    def solution(self, x):
        out = dict(self.design)
        out["values"] = dict(self.values)
        out["values"].update({n: float(v) for n, v in zip(self.free, x)})
        return out


class BuilderProblem:
    """Free variables: the builder segment parameters marked 'free' (continuous, exact model)."""

    def __init__(self, design: dict, free=None):
        from psa.design.builder import BOUNDS, CurveDesign
        self.design = CurveDesign.from_dict(design)
        self.free = self.design.free_parameters()
        if not self.free:
            raise ValueError("no free parameters: tick 'free' for some segment parameters")
        self.bounds = BOUNDS

    def start(self):
        d = self.design
        x0 = np.array([d.segments[i].get(n) for i, n in self.free])
        lo = np.array([self.bounds.get(n, (-1e6, 1e6))[0] for _, n in self.free])
        hi = np.array([self.bounds.get(n, (-1e6, 1e6))[1] for _, n in self.free])
        for k, (i, n) in enumerate(self.free):
            if n == "radius" and d.segments[i].kind == "arc":
                lo[k] = max(lo[k], d.r_min_us)
        return x0, lo, hi

    def apply(self, x):
        for (i, n), v in zip(self.free, x):
            self.design.segments[i].set(n, v)
        return self.design

    def evaluate(self, x, goals):
        return self.apply(x).metrics(goals.b1)

    def solution(self, x):
        return self.apply(x).to_dict()


PROBLEMS = {"waveform": WaveformProblem, "family": FamilyProblem, "builder": BuilderProblem}


def solve_problem(kind, design, free, goals: Goals, callback=None, max_nfev=None, settings=None):
    """Run the solver for a design dictionary of the given kind. Returns (solution dict, result dict)."""
    st = settings or SolverSettings()
    problem = PROBLEMS[kind](design, free)
    x0, lo, hi = problem.start()
    x, info = least_squares_goals(lambda x: problem.evaluate(x, goals), x0, lo, hi, goals, callback, max_nfev, st)
    m = problem.evaluate(x, goals)
    info.update(summary(m, goals, st))
    info["start_summary"] = summary(problem.evaluate(x0, goals), goals, st)["values"]
    return problem.solution(x), info


def run_job(job, queue, stop_event):
    """Entry point of a solver process. ``job`` = {"kind", "design", "free", "goals", "settings"}.
    Puts ("progress", nfev, cost), then ("done", solution, info) or ("error", message) on ``queue``."""
    try:
        last = [0]

        def callback(nfev, cost):
            if nfev - last[0] >= 5 or nfev == 1:
                last[0] = nfev
                queue.put(("progress", nfev, cost))
            return stop_event.is_set()

        settings = SolverSettings.from_dict(job.get("settings"))
        solution, info = solve_problem(job["kind"], job["design"], job.get("free"), Goals.from_dict(job["goals"]),
                                       callback, job.get("max_nfev"), settings)
        queue.put(("done", solution, info))
    except Exception as e:                                       # reported to the GUI
        queue.put(("error", f"{type(e).__name__}: {e}"))
