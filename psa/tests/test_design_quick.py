"""Tests of the quick GRAPE mode (psa.design.quick)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.grape import band_quality, pulse_rates, run_task
from psa.design.pulse import PulseSequence
from psa.design.quick import QuickRun, quick_config, rectangular_start

DT, UMAX = 0.5e-6, 1e4


def test_quick_config_is_the_study_recipe():
    cfg, st = quick_config(450e-6, (-1e4, 1e4), "UR", (0, 1, 0), 90.0, dt=DT)
    assert (cfg.mode, cfg.bb_offset_n, cfg.bb_offset_min_hz, cfg.bb_offset_max_hz) == ("broadband", 19, -1e4, 1e4)
    assert (cfg.controls, cfg.n_spline, cfg.smooth_weight, cfg.limits_late) == ("spline", 75, 1e-6, True)
    assert (cfg.n_starts, cfg.start_pool, cfg.screen_rounds, cfg.screen_nfev, cfg.n_kicks) == (22, "mix", 1, 30, 0)
    assert abs(cfg.noise_correlation_us - 22.5) < 1e-9 and st.max_nfev == 3000
    assert (cfg.max_amp_rate, cfg.max_phase_rate_khz, cfg.target_angle_deg) == (20.0, 50.0, 90.0)
    assert quick_config(10e-6, dt=DT)[0].n_spline == 4                 # never more control points than sub-pulses
    for bad in (dict(band_hz=(1e3, -1e3)), dict(gate="none")):
        try:
            quick_config(100e-6, **bad)
        except ValueError:
            continue
        raise AssertionError(f"accepted {bad}")


def test_pulse_rates():
    ps = PulseSequence(np.array([0.0, 50.0, 100.0, 100.0]), np.array([0.0, 0.0, 0.0, 36.0]), DT, UMAX)
    ar, pr = pulse_rates(ps)
    assert abs(ar - 100.0) < 1e-9                                    # 50 % per 0.5 µs
    assert abs(pr - 0.1 / DT / 1e3) < 1e-9                           # 36° = 0.1 turn per sub-pulse = 200 kHz
    ps = PulseSequence(np.array([0.0, 100.0]), np.array([0.0, 170.0]), DT, UMAX)
    assert pulse_rates(ps)[1] == 0.0                                 # phase jumps at zero amplitude are free


def _drive(q, clock, step):
    phases = []
    while not q.done:
        results = [run_task(t) for t in q.tasks()]
        clock[0] += step
        for r in results:
            q.submit(r)
        if not phases or phases[-1] != q.phase:
            phases.append(q.phase)
    return phases


def test_quick_run_searches_then_polishes_within_the_budget():
    T = 100e-6
    cfg, st = quick_config(T, (-4e3, 4e3), "UR", (1, 0, 0), 90.0, dt=DT, max_amp_rate=20.0, max_phase_rate_khz=50.0)
    ps0 = rectangular_start(T, DT, UMAX)
    clock = [0.0]
    q = QuickRun(ps0, cfg, st, budget_s=10.0, parallel=11, clock=lambda: clock[0], chunk=5)
    start_score = q.best_any["score"]
    assert len(q._ready) == 22 and q.phase == "search"
    phases = _drive(q, clock, 1.0)
    assert phases[-1] == "done" and "polish" in phases, phases
    assert clock[0] <= 10.0 + 1.0 and q._polished == 3                # ends within one slice of the budget
    pulse, s = q.best()
    assert s["within_limits"] and pulse.N == ps0.N and np.max(pulse.amplitude) <= 100 + 1e-9
    ar, pr = pulse_rates(pulse)
    assert ar <= 21.0 and pr <= 52.5
    b = band_quality(pulse, cfg, -4e3, 4e3)
    assert abs(b["worst infidelity"] - s["band"]["worst infidelity"]) < 1e-12
    assert np.log10(b["worst infidelity"]) < start_score - 1          # far better than the rectangular start
    assert any(line.startswith("search:") for line in s["log"]) and any(line.startswith("polish:") for line in s["log"])


def test_quick_run_stop_keeps_the_best_so_far():
    T = 100e-6
    cfg, st = quick_config(T, (-4e3, 4e3), dt=DT)
    clock = [0.0]
    q = QuickRun(rectangular_start(T, DT, UMAX), cfg, st, budget_s=100.0, parallel=11, clock=lambda: clock[0], chunk=5)
    results = [run_task(t) for t in q.tasks()]
    q.stop()
    assert q.phase == "stopping" and q.tasks() == []
    for r in results:
        q.submit(r)
    assert q.done and q.stopped and q.tasks() == []
    pulse, s = q.best()
    assert pulse.N == 200 and s["slices"] == 22


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t(); print(f"[PASS] {t.__name__}")
        except AssertionError as e:
            failed += 1; print(f"[FAIL] {t.__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
