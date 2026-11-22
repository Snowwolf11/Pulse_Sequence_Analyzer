"""Tests of the keyframe model of the visual editor (psa.design.waveform)."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.pulse import PulseSequence
from psa.design.waveform import Waveform, Channel, Keyframe, HOLD, LINEAR, simplify


def test_linear_and_hold_interpolation():
    ch = Channel([Keyframe(0, 0.0), Keyframe(10, 100.0, HOLD), Keyframe(20, 50.0)], 0, 100)
    v = ch.values(25)
    assert np.allclose(v[:10], np.arange(10) * 10.0)      # linear 0 -> 100
    assert np.allclose(v[10:20], 100.0)                    # hold
    assert np.allclose(v[20:], 50.0)                       # after the last keyframe


def test_values_before_first_keyframe_and_clamping():
    ch = Channel([Keyframe(5, 30.0)], 0, 100)
    assert np.allclose(ch.values(8), 30.0)
    ch.set(2, 150.0)
    assert ch.keyframes[0].index == 2 and ch.keyframes[0].value == 100.0


def test_move_snaps_and_respects_neighbours():
    ch = Channel([Keyframe(0, 0.0), Keyframe(10, 1.0), Keyframe(20, 2.0)])
    ch.move(1, 13.6, 5.0, N=30)
    assert ch.keyframes[1].index == 14 and ch.keyframes[1].value == 5.0
    ch.move(1, 25.0, 5.0, N=30)                            # cannot pass the next keyframe
    assert ch.keyframes[1].index == 19
    ch.move(2, 99.0, 0.0, N=30)                            # clamped to the last sub-pulse
    assert ch.keyframes[2].index == 29


def test_phase_is_not_wrapped():
    w = Waveform(N=100)
    w.phase.set(99, 1080.0)
    ps = w.to_pulse()
    assert ps.phase[-1] == 1080.0 and np.all(np.diff(ps.phase) > 0)


def test_from_pulse_exact():
    rng = np.random.default_rng(0)
    amp = np.clip(np.r_[np.full(30, 100.0), np.linspace(100, 20, 40), rng.uniform(0, 100, 30)], 0, 100)
    phase = np.r_[np.zeros(30), np.full(40, 90.0), np.cumsum(rng.normal(0, 40, 30))]
    ps = PulseSequence(amp, phase)
    w = Waveform.from_pulse(ps)
    back = w.to_pulse()
    assert np.allclose(back.amplitude, ps.amplitude)
    assert np.allclose((back.phase - ps.phase + 180) % 360 - 180, 0, atol=1e-9)
    assert len(w.amplitude.keyframes) < ps.N                      # linear parts collapse
    assert len(w.amplitude.keyframes) <= 2 + 2 + 30 + 1


def test_from_pulse_simplified_within_tolerance():
    t = np.linspace(0, 1, 500)
    ps = PulseSequence(50 + 40 * np.sin(6 * t), 720 * t ** 2)
    w = Waveform.from_pulse(ps, amplitude_tolerance=0.5, phase_tolerance=1.0)
    back = w.to_pulse()
    assert np.max(np.abs(back.amplitude - ps.amplitude)) <= 0.5 + 1e-9
    assert np.max(np.abs(back.phase - ps.phase)) <= 1.0 + 1e-9
    assert len(w.amplitude.keyframes) < 40 and len(w.phase.keyframes) < 40


def test_set_duration_and_serialization():
    w = Waveform(N=200, dt=0.5e-6, umax=2e4)
    w.amplitude.set(150, 20.0)
    w.set_duration(50e-6)                                  # 100 sub-pulses: keyframe 150 removed
    assert w.N == 100 and all(k.index < 100 for k in w.amplitude.keyframes)
    w.phase.set(40, 45.0, HOLD)
    d = json.loads(json.dumps(w.to_dict()))
    w2 = Waveform.from_dict(d)
    assert np.allclose(w2.to_pulse().as_array(), w.to_pulse().as_array()) and w2.umax == 2e4


def test_set_dt_keeps_times():
    w = Waveform(N=200, dt=0.5e-6)
    w.amplitude.set(100, 20.0)                             # at 50 us
    w.set_dt(1e-6)
    assert w.N == 100 and w.amplitude.find(50) is not None and abs(w.duration - 100e-6) < 1e-15
    w.set_dt(0.25e-6)
    assert w.N == 400 and w.amplitude.find(200) is not None


def test_simplify_endpoints():
    assert simplify([1, 2, 3, 4], 0) == [0, 3]
    assert simplify([0, 0, 5, 0, 0], 0) == [0, 1, 2, 3, 4]


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
