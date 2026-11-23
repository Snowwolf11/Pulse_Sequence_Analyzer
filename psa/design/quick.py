"""Quick GRAPE mode: target gate, duration, bandwidth and a time budget in, a smooth broadband pulse out.

The recipe is the winner of the quick-settings study (studies/quick_settings_study.py, run 2026-11-22/23,
268 runs; report: studies/quick_settings_report.py):

    robustness   broadband grid over the band, 2 points per 1/T (T = duration), no Lawson reweighting
    controls     cubic splines, one control point per 6 µs, smoothness weight 1e-6
    limits       amplitude rate / phase rate as soft limits, switched on only in a last stage
    exploration  22 starts (mix of perturbed and random, spread 1.0, noise correlation 5 % of T),
                 one screening round of 30 evaluations, the best third continues (3000 evaluations each), no kicks

Driver (QuickRun, like the study's emulated quick mode): the Search tasks run in slices of a few
evaluations, so the run can follow a wall-clock budget. Every slice is scored on a fine test grid (401
offsets over the band). After 80 % of the budget (or when the search is over) the 3 best pulses without
the limits are re-optimized with the limits (the "limits" stage) for the rest of the time. The result is
the pulse with the smallest worst-case infidelity on the test grid among those that keep the rate limits
(5 % tolerance); if none does, the best one overall (flagged).

QuickRun is executor-agnostic like Search: tasks() gives the tasks ready now (run each with
grape.run_task, anywhere), submit(result) takes them back; ``done`` and best() tell the outcome.
"""
from __future__ import annotations

import math
import time

import numpy as np

from psa.design.grape import GrapeConfig, Search, stage_plan
from psa.design.optimize import SolverSettings
from psa.design.pulse import PulseSequence

RECIPE = dict(density=2.0, spline_us=6.0, smooth_weight=1e-6, n_starts=22, start_pool="mix", spread=1.0,
              correlation=0.05, screen_rounds=1, screen_nfev=30, keep_fraction=0.33, n_kicks=0, max_nfev=3000)
DEFAULT_LIMITS = dict(max_amp_rate=20.0, max_phase_rate_khz=50.0)
ACCEPT = 1.05                    # a result may exceed the rate limits by 5 %
CHUNK = 20                       # evaluations per slice
POLISH = 0.2                     # last fraction of the budget: the best free pulses get the limits
K_POLISH = 3
TEST_POINTS = 401


def quick_config(duration_s, band_hz=(-1e4, 1e4), gate="UR", target_axis=(1, 0, 0), target_angle_deg=180.0,
                 m0=(0, 0, 1), m_target=(1, 0, 0), max_amp_rate=20.0, max_phase_rate_khz=50.0, seed=0, dt=None):
    """(GrapeConfig, SolverSettings) of the quick recipe for a pulse of ``duration_s`` over ``band_hz``."""
    lo, hi = float(band_hz[0]), float(band_hz[1])
    if not hi > lo:
        raise ValueError("the band must have upper > lower offset")
    if gate not in ("UR", "PP"):
        raise ValueError("quick mode needs a gate goal (UR or PP)")
    us = duration_s * 1e6
    n_pts = int(round(RECIPE["density"] * (hi - lo) * duration_s)) + 1
    n_spline = max(4, int(round(us / RECIPE["spline_us"])))
    if dt is not None:
        n_spline = min(n_spline, max(4, int(round(duration_s / dt))))
    cfg = GrapeConfig(mode="broadband", bb_offset_min_hz=lo, bb_offset_max_hz=hi, bb_offset_n=max(2, n_pts),
                      gate=gate, target_axis=list(map(float, target_axis)), target_angle_deg=float(target_angle_deg),
                      m0=list(map(float, m0)), m_target=list(map(float, m_target)), seed=int(seed),
                      controls="spline", n_spline=n_spline, smooth_weight=RECIPE["smooth_weight"],
                      max_amp_rate=float(max_amp_rate), max_phase_rate_khz=float(max_phase_rate_khz), limits_late=True,
                      n_starts=RECIPE["n_starts"], start_pool=RECIPE["start_pool"], multistart_noise=RECIPE["spread"],
                      noise_correlation_us=RECIPE["correlation"] * us, screen_rounds=RECIPE["screen_rounds"],
                      screen_nfev=RECIPE["screen_nfev"], keep_fraction=RECIPE["keep_fraction"],
                      n_kicks=RECIPE["n_kicks"], crossover_rate=0.0, coarse_levels=[])
    cfg.check()
    # "solved" at gate infidelity 1e-8, as in the study (the GRAPE default is stricter)
    return cfg, SolverSettings(max_nfev=RECIPE["max_nfev"], gate_infidelity_tolerance=1e-8)


def rectangular_start(duration_s, dt, umax):
    n = max(4, int(round(duration_s / dt)))
    return PulseSequence(np.full(n, 100.0), np.zeros(n), dt, umax)


class QuickRun:
    """Time-budgeted driver of the quick recipe (see the module doc)."""

    def __init__(self, ps: PulseSequence, config: GrapeConfig, settings: SolverSettings, budget_s, parallel=1,
                 clock=time.time, chunk=CHUNK, polish=POLISH, k_polish=K_POLISH):
        self.ps, self.cfg, self.st = ps, config, settings
        self.band = (float(config.bb_offset_min_hz), float(config.bb_offset_max_hz), TEST_POINTS)
        self.budget = float(budget_s)
        self.clock, self.chunk, self.k_polish = clock, int(chunk), int(k_polish)
        self.t0 = clock()
        self.search_end = self.budget * (1 - polish)
        self.search = Search(ps, config, settings, parallel=parallel)
        self.limits = (float(config.max_amp_rate), float(config.max_phase_rate_khz))
        kinds = [k for *_, k in stage_plan(config, ps.N)]
        self.limits_stage = kinds.index("limits") if "limits" in kinds else None
        self.phase = "search"                      # search -> polish -> done
        self.incumbent = None                      # best pulse within the limits: dict(score, band, rates, amp, phase)
        self.best_any = None                       # best pulse regardless of the limits
        self.pool = {}                             # cid -> (score, w): best free pulse per candidate
        self.improved = False                      # set by submit when the incumbent changed
        self.stopped = False
        self.log = []
        self.slices = 0
        self._ready, self._out, self._tid = [], {}, 0
        self._polished = 0
        self.consider(ps.amplitude, ps.phase, None, None)          # the start pulse is the baseline
        self._absorb()

    # ------------------------------------------------------------ public ---
    @property
    def done(self):
        return self.phase == "done"

    @property
    def running(self):
        return len(self._out)

    @property
    def round(self):
        return self.search.round

    @property
    def cands(self):
        return self.search.cands

    def elapsed(self):
        return self.clock() - self.t0

    def tasks(self):
        out, self._ready = self._ready, []
        return out

    def stop(self):
        """No new slices; done when the running ones have returned (the best pulse so far is kept)."""
        self.stopped = True
        self.phase = "stopping" if self._out else "done"
        self._ready = []

    def acceptable(self, rates):
        ar, pr = rates
        la, lp = self.limits
        return (la <= 0 or ar <= la * ACCEPT) and (lp <= 0 or pr <= lp * ACCEPT)

    def consider(self, amp, phase, band, rates, cid=None, w=None):
        """Record a pulse (band = band_quality on the test grid, rates = pulse_rates); returns True if it is
        the new incumbent."""
        if band is None:
            from psa.design.grape import band_quality, pulse_rates
            p = PulseSequence(np.asarray(amp), np.asarray(phase), self.ps.dt, self.ps.umax)
            band, rates = band_quality(p, self.cfg, *self.band), pulse_rates(p)
        score = math.log10(max(band["worst infidelity"], 1e-16))
        rec = dict(score=score, band=band, rates=tuple(rates), amp=list(amp), phase=list(phase))
        if self.best_any is None or score < self.best_any["score"]:
            self.best_any = rec
        if cid is not None and w is not None and score < self.pool.get(cid, (math.inf,))[0]:
            self.pool[cid] = (score, w)
        if self.acceptable(rates) and (self.incumbent is None or score < self.incumbent["score"]):
            self.incumbent = rec
            return True
        return False

    def submit(self, result):
        rec = self._out.pop(result["tid"], None)
        if rec is None:
            return
        t = rec["task"]
        if "error" in result:
            self.log.append(f"slice error: {result['error']}")
            if not rec["polish"]:
                self.search.submit({"tid": t["tid"], "cid": t["cid"], "error": result["error"]})
        else:
            info = result["info"]
            if self.consider(result["amp"], result["phase"], result["band"], result["rates"],
                             None if rec["polish"] else t["cid"], info.get("w")):
                self.improved = True
            prev = rec["state"]["nfev"] if rec["state"] and "nfev" in rec["state"] else 0
            rec["state"] = info["state"]
            if rec["left"] is not None:
                rec["left"] -= max(info["nfev"] - prev, 1)
            finished = info["finished"] or (rec["left"] is not None and rec["left"] <= 0)
            end = self.search_end if self.phase == "search" else self.budget
            if finished:
                if not rec["polish"]:
                    info["history"] = []
                    self.search.submit({"tid": t["tid"], "cid": t["cid"], "amp": result["amp"],
                                        "phase": result["phase"], "info": info})
                    if self.phase == "search":
                        self._absorb()
            elif self.phase in ("search", "polish") and self.elapsed() < end:
                self._issue(rec)
        self._advance()

    def best(self):
        """(pulse, summary) of the result: the incumbent (within the limits) or, if there is none, the best
        pulse overall (summary["within_limits"] False)."""
        rec = self.incumbent if self.incumbent is not None else self.best_any
        pulse = PulseSequence(np.asarray(rec["amp"]), np.asarray(rec["phase"]), self.ps.dt, self.ps.umax)
        return pulse, {"within_limits": self.incumbent is not None, "band": rec["band"], "rates": rec["rates"],
                       "slices": self.slices, "candidates": len(self.search.cands), "log": list(self.log),
                       "elapsed_s": self.elapsed(), "polished": self._polished}

    # ----------------------------------------------------------- internal ---
    def _issue(self, rec):
        t = rec["task"]
        n = self.chunk if rec["left"] is None else max(1, min(self.chunk, rec["left"]))
        task = dict(t, tid=self._tid, state=rec["state"], budget=n, band=self.band, snapshots=False)
        self._out[self._tid] = rec
        self._ready.append(task)
        self._tid += 1
        self.slices += 1

    def _absorb(self):
        for t in self.search.tasks():
            self._issue({"task": t, "state": t["state"], "left": t["budget"], "polish": False})

    def _advance(self):
        if self._out or self._ready:
            return
        if self.phase == "search" and (self.search.done or self.elapsed() >= self.search_end):
            if not self.search.done:
                self.search.stop()
            self.log.append(f"search: {len(self.search.cands)} candidates, {self.slices} slices, "
                            f"{self.elapsed():.0f} s")
            self._start_polish()
        elif self.phase in ("search", "polish", "stopping"):
            if self.phase == "polish":
                self.log.append(f"polish: {self._polished} pulses re-optimized with the rate limits, {self.elapsed():.0f} s")
            self.phase = "done"

    def _start_polish(self):
        self.phase = "polish"
        if self.limits_stage is None or not self.pool or self.elapsed() >= self.budget:
            self.phase = "done"
            return
        picks = sorted(self.pool.values(), key=lambda v: v[0])[:self.k_polish]
        base = {"tid": None, "cid": -1, "pulse": self.search._pulse, "config": self.cfg.to_dict(),
                "settings": self.st.to_dict(), "start_index": 0}
        for i, (_, w) in enumerate(picks):
            state = {"stage": self.limits_stage, "w": list(w), "stage_nfev": 0, "nfev": 0}
            self._issue({"task": dict(base, cid=-1 - i), "state": state, "left": None, "polish": True})
        self._polished = len(picks)
