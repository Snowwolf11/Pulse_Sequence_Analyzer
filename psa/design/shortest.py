"""Shortest robust pulse: the shortest pulse found (on the sub-pulse grid dt) that implements a gate exactly and
cancels the detuning error through order p (A1..Ap = 0) and optionally the first-order B1 error, with the rf
amplitude at most umax.

Method (direct, "descent with continuation"): GRAPE in local mode (psa.design.grape, exact gradients, the
amplitude limit built in) solves the goals as equations at a fixed number of sub-pulses N. The search
    1. finds a solution at the start duration (more time if none of the starts converges),
    2. shortens it step by step: the best solution, squeezed in time (gate preserving: amplitude scaled up,
       clipped at 100 %), plus fresh starts (random pulses, 2H/4H/6H, CORPSE) are solved at the shorter N;
       success keeps the step, failure halves it, down to one sub-pulse,
    3. verifies every accepted pulse independently (psa.design.metrics: exact jets, not the optimizer's terms).
The result is the shortest pulse *found*: a local method, no proof of optimality. Reported with simple lower
bounds: gate only, T >= angle / (360° nu_max); with A1 = 0 (closed error curve, curvature <= 2 pi nu_max)
Fenchel's theorem gives T >= 1 / nu_max (100 µs at 10 kHz).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Callable, List, Optional

import numpy as np

from psa.design.grape import GrapeConfig, optimize
from psa.design.metrics import analyze
from psa.design.optimize import SolverSettings
from psa.design.pulse import PulseSequence

STALL_WINDOW = 150          # evaluations: a start whose goal residual improved less than STALL_GAIN (decades)
STALL_GAIN = 0.3            # over the last STALL_WINDOW evaluations while still above STALL_LEVEL is given up
STALL_LEVEL = 1e-7

# verification thresholds (independent of the optimizer's "solved")
GATE_TOL = 1e-10            # average gate infidelity (UR) / transfer infidelity (PP)
ORDER_TOL = 1e-8            # |A_k| / T^k
B1_TOL = 1e-8               # |B| / nutation
AMP_TOL = 100.0 + 1e-9      # %


@dataclass
class ShortestSpec:
    gate: str = "UR"                                       # "UR" or "PP"
    axis: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0])
    angle_deg: float = 180.0
    m0: List[float] = field(default_factory=lambda: [0.0, 0.0, 1.0])
    mt: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0])
    order: int = 1                                         # A1..A_order = 0 (0: the gate only)
    b1: bool = False                                       # first-order B1 robustness too
    dt: float = 0.5e-6
    umax: float = 1e4                                      # Hz (100 %)

    def check(self):
        if self.gate not in ("UR", "PP"):
            raise ValueError("the gate must be UR or PP")
        if not 0 <= int(self.order) <= 3:
            raise ValueError("the order must be 0, 1, 2 or 3")
        if self.dt <= 0 or self.umax <= 0:
            raise ValueError("dt and umax must be positive")
        if self.gate == "UR" and np.linalg.norm(self.axis) == 0:
            raise ValueError("the rotation axis must not be zero")

    def config(self, seed=0, random_start=False) -> GrapeConfig:
        cfg = GrapeConfig(mode="local", offsets_hz=[0.0], orders=list(range(1, int(self.order) + 1)), b1=bool(self.b1),
                          gate=self.gate, target_axis=list(map(float, self.axis)), target_angle_deg=float(self.angle_deg),
                          m0=list(map(float, self.m0)), m_target=list(map(float, self.mt)), seed=int(seed),
                          start_kind="random" if random_start else "perturb")
        cfg.check()
        return cfg

    def lower_bound_s(self):
        """A simple lower bound of the duration (see the module doc)."""
        if self.order >= 1:
            return 1.0 / self.umax
        if self.gate == "UR":
            return abs(self.angle_deg) / 360.0 / self.umax
        m0, mt = (np.asarray(v, float) / np.linalg.norm(v) for v in (self.m0, self.mt))
        return math.degrees(math.acos(np.clip(m0 @ mt, -1, 1))) / 360.0 / self.umax


def verify(ps: PulseSequence, spec: ShortestSpec) -> dict:
    """Independent check of a pulse (exact jets): gate error, |A_k| / T^k, B1, max. amplitude, ok."""
    a = analyze(ps, b1=True)
    q = a["gate"]["quaternion"]
    if spec.gate == "UR":
        n = np.asarray(spec.axis, float) / np.linalg.norm(spec.axis)
        h = math.radians(spec.angle_deg) / 2
        c = abs(q @ np.r_[math.cos(h), math.sin(h) * n])
        gate_err = 1 - (2 + 4 * c * c) / 6
    else:
        m0, mt = (np.asarray(v, float) / np.linalg.norm(v) for v in (spec.m0, spec.mt))
        gate_err = 1 - (1 + (a["gate"]["matrix"] @ m0) @ mt) / 2
    rb = a["robustness"]
    orders = {f"A{k}": float(rb[f"A{k}_norm"]) for k in (1, 2, 3)}
    if spec.gate == "PP":                       # only the error components perpendicular to m0 matter
        m0 = np.asarray(spec.m0, float) / np.linalg.norm(spec.m0)
        for k in (1, 2, 3):
            orders[f"A{k}"] = float(np.linalg.norm(np.cross(rb[f"A{k}"], m0)) / ps.duration ** k)
    b1 = float(a["B1"]["B1_norm"])
    ok = (gate_err < GATE_TOL and all(orders[f"A{k}"] < ORDER_TOL for k in range(1, spec.order + 1))
          and (not spec.b1 or b1 < B1_TOL) and float(np.max(ps.amplitude)) <= AMP_TOL)
    return {"ok": bool(ok), "gate error": float(max(gate_err, 0.0)), **orders, "B1": b1,
            "max amplitude [%]": float(np.max(ps.amplitude)), "duration [µs]": ps.duration * 1e6, "N": int(ps.N)}


# ---------------------------------------------------------------- starts ---

def fit_steps(ps: PulseSequence, n) -> PulseSequence:
    """The pulse on n sub-pulses of the same dt, gate preserving: time stretched / squeezed, amplitude scaled by
    the inverse factor (clipped at 100 %), unwrapped phase interpolated."""
    n = max(2, int(n))
    x_old, x_new = (np.arange(ps.N) + 0.5) / ps.N, (np.arange(n) + 0.5) / n
    amp = np.interp(x_new, x_old, np.asarray(ps.amplitude, float)) * ps.N / n
    phase = np.degrees(np.interp(x_new, x_old, np.unwrap(np.radians(np.asarray(ps.phase, float)))))
    return PulseSequence(np.clip(amp, 0.0, 100.0), (phase + 180.0) % 360.0 - 180.0, ps.dt, ps.umax)


def corpse(angle_deg, phase_deg, dt, umax) -> PulseSequence:
    """Short CORPSE (first-order detuning robust rotation about an axis in the xy-plane: rotations theta/2 - k,
    2 pi - 2k, theta/2 - k at phases phi, phi + 180°, phi, sin k = sin(theta/2) / 2; 116.7 µs for 180° at
    10 kHz) at full amplitude, sampled on dt (rounded: a start, not exact)."""
    t = math.radians(angle_deg)
    k = math.asin(math.sin(t / 2) / 2)
    angles = (t / 2 - k, 2 * math.pi - 2 * k, t / 2 - k)
    amp, ph = [], []
    for a, p in zip(angles, (phase_deg, phase_deg + 180.0, phase_deg)):
        n = max(1, int(round(a / (2 * math.pi * umax * dt))))
        amp.append(np.full(n, 100.0)); ph.append(np.full(n, p))
    return PulseSequence(np.concatenate(amp), (np.concatenate(ph) + 180.0) % 360.0 - 180.0, dt, umax)


def start_pulses(spec: ShortestSpec, n, n_random=3, seed=0):
    """Named start pulses on n sub-pulses: random pulses, CORPSE (UR about an axis in the xy-plane), 2H / 4H / 6H
    (identity gates with robust error curves; the solver has to change the gate)."""
    from psa.design.families import family_by_name
    rng = np.random.default_rng(seed)
    out = [(f"random {i + 1}", PulseSequence(rng.uniform(30, 90, n), rng.uniform(-180, 180, n), spec.dt, spec.umax))
           for i in range(n_random)]
    ax = np.asarray(spec.axis, float)
    if spec.gate == "UR" and abs(ax[2]) < 1e-9 * np.linalg.norm(ax) and spec.angle_deg % 360:
        out.append(("CORPSE", fit_steps(corpse(spec.angle_deg, math.degrees(math.atan2(ax[1], ax[0])), spec.dt,
                                               spec.umax), n)))
    for name, vals in (("2H (closed two-helix loop)", {"tau": 75}), ("4H (0° UR, 2nd-order robust)", {"tau": 99.5}),
                       ("6H (0° UR, 3rd-order robust)", {"TA": 4.15})):
        try:
            out.append((name.split(" ")[0], fit_steps(family_by_name(name).pulse(vals, spec.dt, spec.umax), n)))
        except ValueError:
            pass
    return out


# ---------------------------------------------------------------- search ---

def stall_stopper(window=STALL_WINDOW, gain=STALL_GAIN, level=STALL_LEVEL):
    """optimize() callback that gives a start up once it stalls: the goal residual improved less than ``gain``
    decades over the last ``window`` evaluations while still above ``level`` (a shorter duration that cannot
    be reached costs a few hundred evaluations instead of the full budget)."""
    hist = []

    def callback(nfev, goal):
        hist.append(math.log10(max(float(goal), 1e-300)))
        return (len(hist) > window and goal > level and hist[-window - 1] - hist[-1] < gain)
    return callback


def solve(ps: PulseSequence, spec: ShortestSpec, seed=0, max_nfev=3000, noise=None, stall=True):
    """Local-mode GRAPE from ps at its N (given up when it stalls, see stall_stopper). Returns (pulse,
    verification)."""
    cfg = spec.config(seed)
    if noise is not None:
        cfg.start_noise = float(noise)
    out, info = optimize(ps, cfg, SolverSettings(max_nfev=int(max_nfev)), start_summary=False,
                         callback=stall_stopper() if stall else None)
    return out, verify(out, spec)


def _task(args):
    ps_arr, spec_d, seed, max_nfev, noise, dt, umax = args
    spec = ShortestSpec(**spec_d)
    ps = PulseSequence(ps_arr[:, 0], ps_arr[:, 1], dt, umax)
    out, v = solve(ps, spec, seed, max_nfev, noise)
    return np.column_stack([out.amplitude, out.phase]), v


class SearchStopped(Exception):
    """Raised by a mapper (or should_stop) to end a search early; search() returns the best pulse so far."""


def search(spec: ShortestSpec, start_us, n_random=3, perturbed=2, max_nfev=3000, time_limit_s=600.0,
           max_grow=4, mapper: Optional[Callable] = None, log: Optional[Callable] = None, start: PulseSequence = None,
           seed=0, on_best: Optional[Callable] = None):
    """Shortest pulse found, starting at start_us (see the module doc).

    At every duration the candidates are: the best solution so far squeezed to N (plus ``perturbed`` noisy
    copies), ``n_random`` random pulses and the family / CORPSE starts. ``mapper(fn, iterable)`` runs them
    (default: sequential; e.g. a process pool's map); it may raise SearchStopped to end the search (the best
    pulse so far is returned, message "stopped"). ``on_best(pulse, verification)`` is called for every
    accepted (shorter) pulse. Returns a dict with the best pulse, its verification, the history [(N, ok, best
    gate error, best A, seconds)] and the lower bound. ``seed`` changes every random start and perturbation
    (repeat a search with another seed to check the result)."""
    spec.check()
    mapper = mapper or (lambda f, it: list(map(f, it)))
    log = log or (lambda s: None)
    on_best = on_best or (lambda p, v: None)
    t0 = time.time()
    spec_d = dict(spec.__dict__)
    history = []

    def batch(cands):
        args = [(np.column_stack([p.amplitude, p.phase]), spec_d, s + 1000 * int(seed), nfev, noise, spec.dt, spec.umax)
                for _, p, s, noise, nfev in cands]
        return list(zip(cands, mapper(_task, args)))

    def run(n, base):
        # stage 1: continuation of the best solution (cheap, usually enough); stage 2: fresh starts
        res = []
        if base is not None:
            sq = fit_steps(base, n)
            res += batch([("continuation", sq, 0, None, max_nfev)]
                         + [(f"continuation ~{i + 1}", sq, 100 + i, 0.05, max_nfev) for i in range(perturbed)])
        if not any(v["ok"] for _, (_, v) in res):
            fresh = [(name, p, 200 + i, None, max(200, max_nfev // 2)) for i, (name, p) in
                     enumerate(start_pulses(spec, n, n_random, seed=n + 7919 * int(seed)))]
            if start is not None:
                fresh.append(("given start", fit_steps(start, n), 300, None, max_nfev))
            res += batch(fresh)
        best, best_v, best_name = None, None, None
        for (name, *_), (arr, v) in res:
            key = (not v["ok"], v["gate error"] + sum(v[f"A{k}"] for k in range(1, spec.order + 1)) + (v["B1"] if spec.b1 else 0))
            if best_v is None or key < best_key:
                best, best_v, best_name, best_key = PulseSequence(arr[:, 0], arr[:, 1], spec.dt, spec.umax), v, name, key
        res = [v for _, (_, v) in res]
        history.append({"N": n, "duration [µs]": n * spec.dt * 1e6, "ok": best_v["ok"], "from": best_name,
                        "gate error": best_v["gate error"], "worst order": max([best_v[f"A{k}"] for k in range(1, spec.order + 1)] or [0.0]),
                        "B1": best_v["B1"], "seconds": time.time() - t0,
                        "solved": sum(v["ok"] for v in res), "candidates": len(res)})
        log(f"{n * spec.dt * 1e6:8.1f} µs: {'ok  ' if best_v['ok'] else 'fail'} ({history[-1]['solved']}/{len(res)} solved, "
            f"best from {best_name}) {time.time() - t0:6.1f} s")
        return (best, best_v) if best_v["ok"] else (None, best_v)

    # 1. a solution at the start duration (grow if needed)
    n = max(2, int(round(start_us * 1e-6 / spec.dt)))
    best, stopped = None, False
    try:
        for _ in range(max_grow + 1):
            best, v = run(n, None)
            if best is not None or time.time() - t0 > time_limit_s:
                break
            n = int(math.ceil(n * 1.25))
    except SearchStopped:
        stopped = True
    if best is None:
        return {"pulse": None, "verification": None, "history": history, "lower bound [µs]": spec.lower_bound_s() * 1e6,
                "message": "stopped" if stopped else "no solution found at the start duration (or longer)",
                "elapsed_s": time.time() - t0}
    best_v = v
    on_best(best, best_v)
    # 2. descent
    # steps of 10 % down; after the first failure: bisection between the longest failure and the best success
    n_floor = max(2, int(math.ceil(spec.lower_bound_s() / spec.dt - 1e-9)))
    step = max(1, int(round(0.1 * best.N)))
    failed = n_floor - 1                       # longest duration known to fail (n_floor - 1: below the bound)
    message = "no shorter pulse found (bisection down to one sub-pulse)"
    try:
        while not stopped and best.N - failed > 1:
            if time.time() - t0 > time_limit_s:
                message = "time limit"
                break
            n_try = max(best.N - step, failed + 1)
            if failed >= n_floor:                  # a failure is known: bisect
                n_try = (best.N + failed) // 2
            cand, v = run(n_try, best)
            if cand is not None:
                best, best_v = cand, v
                on_best(best, best_v)
            else:
                failed = n_try
        else:
            if not stopped and failed < n_floor:
                message = "reached the lower bound"
    except SearchStopped:
        stopped = True
    if stopped:
        message = "stopped"
    return {"pulse": best, "verification": best_v, "history": history, "lower bound [µs]": spec.lower_bound_s() * 1e6,
            "message": message, "elapsed_s": time.time() - t0}


# ------------------------------------------------------- background run ---

BLAS_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
             "NUMEXPR_NUM_THREADS")


def search_process(spec_d, start_us, kwargs, queue, stop_event, n_workers=1, start_arr=None):
    """Entry point of a background search (e.g. the GUI's shortest-pulse window; run it in a spawned process).
    The candidates of every duration run in a pool of n_workers processes (single-threaded BLAS). Messages on
    queue: ("log", text), ("best", amp, phase, verification) for every accepted pulse, and at the end ("done",
    result dict without the pulse, amp, phase) or ("error", text). stop_event ends the search (the pool is
    terminated, the best pulse so far is returned)."""
    import multiprocessing as mp
    import os
    import traceback
    try:
        spec = ShortestSpec(**spec_d)
        start = None if start_arr is None else PulseSequence(start_arr[:, 0], start_arr[:, 1], spec.dt, spec.umax)
        for k in BLAS_VARS:
            os.environ[k] = "1"
        with mp.get_context("spawn").Pool(max(1, int(n_workers))) as pool:
            def mapper(f, it):
                if stop_event.is_set():
                    raise SearchStopped()
                r = pool.map_async(f, list(it))
                while not r.ready():
                    if stop_event.is_set():
                        raise SearchStopped()
                    r.wait(0.2)
                return r.get()

            def on_best(p, v):
                queue.put(("best", list(map(float, p.amplitude)), list(map(float, p.phase)), v))

            res = search(spec, start_us, mapper=mapper, log=lambda s: queue.put(("log", s)), start=start,
                         on_best=on_best, **kwargs)
            pool.terminate()
        p = res.pop("pulse")
        queue.put(("done", res, None if p is None else list(map(float, p.amplitude)),
                   None if p is None else list(map(float, p.phase))))
    except Exception as e:                                      # noqa: BLE001 - reported to the GUI
        queue.put(("error", f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=3)}"))
