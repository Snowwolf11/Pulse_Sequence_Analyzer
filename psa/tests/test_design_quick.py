"""Tests of the quick GRAPE mode (psa.design.quick)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.grape import band_quality, pulse_rates, run_task
from psa.design.pulse import PulseSequence
from psa.design.quick import QuickRun, b1_count, quick_config, rectangular_start, stretch_pulse

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
    start_score = q.start_rec["score"]
    assert q.incumbent is None and q.best_any is None             # the start pulse is only the baseline
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


def test_quick_b1_range():
    """B1 range: the offsets of the recipe times n_B1 = 2 (s_hi - s_lo) nu_rf T + 1 (>= 3) scalings; the run is
    scored on the area. Without a range the config is the plain recipe."""
    plain = quick_config(450e-6, (-1e4, 1e4), dt=DT)[0]
    same = quick_config(450e-6, (-1e4, 1e4), dt=DT, b1_range=(1.0, 1.0), umax_hz=UMAX)[0]
    assert plain.to_dict() == same.to_dict() and plain.bb_scale_n == 1
    assert b1_count((0.8, 1.2), UMAX, 450e-6) == 5 and b1_count((0.95, 1.05), UMAX, 100e-6) == 3
    assert b1_count((1.0, 1.0), UMAX, 450e-6) == 1
    cfg = quick_config(450e-6, (-1e4, 1e4), dt=DT, b1_range=(0.8, 1.2), umax_hz=UMAX)[0]
    assert (cfg.bb_offset_n, cfg.bb_scale_min, cfg.bb_scale_max, cfg.bb_scale_n) == (19, 0.8, 1.2, 5)
    assert len(cfg.ensemble()) == 95
    for bad in (dict(b1_range=(1.2, 0.8), umax_hz=UMAX), dict(b1_range=(0.8, 1.2)), dict(b1_range=(0.0, 1.0), umax_hz=UMAX)):
        try:
            quick_config(100e-6, dt=DT, **bad)
        except ValueError:
            continue
        raise AssertionError(f"accepted {bad}")
    T = 100e-6
    cfg, st = quick_config(T, (-4e3, 4e3), "UR", (1, 0, 0), 90.0, dt=DT, b1_range=(0.9, 1.1), umax_hz=UMAX)
    clock = [0.0]
    q = QuickRun(rectangular_start(T, DT, UMAX), cfg, st, budget_s=4.0, parallel=11, clock=lambda: clock[0], chunk=5)
    assert q.band == (-4e3, 4e3, 201, 0.9, 1.1, 9)
    _drive(q, clock, 1.0)
    pulse, s = q.best()
    b = band_quality(pulse, cfg, *q.band)
    assert abs(b["worst infidelity"] - s["band"]["worst infidelity"]) < 1e-12 and "worst at B1 [%]" in b
    nominal = band_quality(pulse, cfg, -4e3, 4e3, 201)
    assert b["worst infidelity"] >= nominal["worst infidelity"] - 1e-15    # the area contains the nominal line


def test_quick_run_reserves_time_for_the_polish():
    """Long slices (big problems): the search ends early enough for POLISH_SLICES polish slices, and the result
    is never the unoptimized start pulse when optimized pulses exist."""
    from psa.design.quick import POLISH_SLICES
    T = 100e-6
    cfg, st = quick_config(T, (-4e3, 4e3), "UR", (1, 0, 0), 90.0, dt=DT)
    clock = [0.0]
    q = QuickRun(rectangular_start(T, DT, UMAX), cfg, st, budget_s=100.0, parallel=11, clock=lambda: clock[0], chunk=5)
    phases = _drive(q, clock, 15.0)                      # every cycle of slices takes 15 s of a 100 s budget
    assert "polish" in phases and q._polished == 3 and q.slice_s == 15.0
    assert q.search_deadline == 100.0 - POLISH_SLICES * 15.0          # earlier than at 80 s
    pulse, s = q.best()
    assert s["optimized"] and s["within_limits"]
    # budget over before any polish: one polish round still runs, the result is an optimized pulse
    clock = [0.0]
    q = QuickRun(rectangular_start(T, DT, UMAX), cfg, st, budget_s=5.0, parallel=11, clock=lambda: clock[0], chunk=5)
    _drive(q, clock, 30.0)
    pulse, s = q.best()
    assert q._polished == 3 and s["optimized"] and s["band"]["worst infidelity"] < q.start_rec["band"]["worst infidelity"]


def test_stretch_pulse():
    ps = PulseSequence(np.linspace(10, 90, 100), np.linspace(-170, 530, 100) % 360 - 180, DT, UMAX)
    long = stretch_pulse(ps, 100e-6)
    assert long.N == 200 and long.dt == DT and long.umax == UMAX
    assert np.isclose(long.amplitude[0], ps.amplitude[0], atol=0.5) and np.isclose(long.amplitude[-1], 90, atol=0.5)
    unwrapped = np.degrees(np.unwrap(np.radians(long.phase)))
    assert np.all(np.diff(unwrapped) > 0) and np.max(np.abs(np.diff(unwrapped))) < 5    # no jumps from the wrap
    same = stretch_pulse(ps, 50e-6)
    assert np.allclose(same.amplitude, ps.amplitude) and np.allclose(np.cos(np.radians(same.phase - ps.phase)), 1)


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
