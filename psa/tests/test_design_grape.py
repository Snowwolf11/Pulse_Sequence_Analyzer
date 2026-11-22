"""Tests of the GRAPE optimizer (psa.design.grape): exact gradients, consistency, optimizations."""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design import metrics
from psa.design.grape import (GrapeConfig, GrapeProblem, optimize, optimize_multistart, propagator_and_gradient,
                               rotation_with_gradient, spline_basis)
from psa.design.optimize import SolverSettings
from psa.design.pulse import PulseSequence
from psa.rotations import propagate_batch, quaternion_to_matrix

DT, UMAX = 0.5e-6, 1e4


def random_pulse(N=40, seed=0):
    rng = np.random.default_rng(seed)
    return PulseSequence(rng.uniform(20, 90, N), rng.uniform(-180, 180, N), DT, UMAX)


CONFIGS = {
    "local UR orders 1-3, two offsets": GrapeConfig(mode="local", offsets_hz=[0.0, 3000.0], orders=[1, 2, 3], gate="UR",
                                                    target_axis=[0, 1, 0], target_angle_deg=90),
    "local PP orders 1-2 + B1": GrapeConfig(mode="local", offsets_hz=[-2000.0], orders=[1, 2], b1=True, gate="PP",
                                            m0=[0, 0, 1], m_target=[1, 0, 0]),
    "local B1 only": GrapeConfig(mode="local", offsets_hz=[0.0], orders=[], b1=True, gate="none"),
    "broadband UR offsets x scales": GrapeConfig(mode="broadband", gate="UR", bb_offset_n=3, bb_scale_min=0.9,
                                                 bb_scale_max=1.1, bb_scale_n=2),
    "broadband PP": GrapeConfig(mode="broadband", gate="PP", m0=[0, 0, 1], m_target=[0, -1, 0], bb_offset_n=4),
    "stay close": GrapeConfig(mode="local", orders=[1], gate="UR", regularization="stay_close",
                              regularization_weight=0.3),
    "smoothness + rate limits": GrapeConfig(mode="local", orders=[1], gate="UR", smooth_weight=0.05, max_amp_rate=20,
                                            max_phase_rate_khz=50),
    "spline controls + limits": GrapeConfig(mode="local", offsets_hz=[1000.0], orders=[1, 2], gate="PP", controls="spline",
                                            n_spline=9, max_amp_rate=5, max_phase_rate_khz=5, smooth_weight=0.01),
}


def test_gradients_match_finite_differences():
    """Jacobian with respect to the optimization parameters (through squash map and spline basis)."""
    ps = random_pulse(24)
    for name, cfg in CONFIGS.items():
        prob = GrapeProblem(ps, cfg)
        rng = np.random.default_rng(1)
        p = prob.p0 + rng.normal(0, 0.3, prob.p0.shape)
        res = lambda q: prob.evaluate(prob.to_controls(q)[0])[0]
        u, Mp = prob.to_controls(p)
        J = (prob.evaluate(u)[1] @ Mp).toarray()
        h = 1e-6
        fd = np.column_stack([(res(p + h * e) - res(p - h * e)) / (2 * h) for e in np.eye(len(p))])
        scale = np.max(np.abs(fd))
        assert np.max(np.abs(J - fd)) < 1e-5 * scale, (name, np.max(np.abs(J - fd)) / scale)
        if "limits" in name:                                     # the limits are active at this point
            assert prob.evaluate(u)[2]["rate-limit excess"] > 0


def test_amplitude_limit_map():
    prob = GrapeProblem(random_pulse(30, seed=5), GrapeConfig())
    rng = np.random.default_rng(2)
    v = np.r_[rng.normal(0, 2, 60)]
    v[:3] = 0.0; v[30:33] = [1e-8, 0.0, -3e-7]               # tiny vectors: series branch
    u, M = prob.squash(v)
    assert np.all(np.hypot(u[:30], u[30:]) < 1)
    M = M.toarray(); h = 1e-7
    fd = np.column_stack([(prob.squash(v + h * e)[0] - prob.squash(v - h * e)[0]) / (2 * h) for e in np.eye(60)])
    assert np.max(np.abs(M - fd)) < 1e-6
    assert np.allclose(prob.squash(prob.unsquash(prob.x0))[0], prob.x0, atol=1e-12)     # random_pulse amplitudes < 99 %


def test_engine_consistent_with_metrics():
    ps = random_pulse(60, seed=3)
    ux = ps.amplitude / 100 * np.cos(np.radians(ps.phase)); uy = ps.amplitude / 100 * np.sin(np.radians(ps.phase))
    X, _ = propagator_and_gradient(ux, uy, DT, UMAX, 0.0, 1.0, 3)
    from psa.design.jets import magnus_vectors
    rb = metrics.robustness(ps)
    for k, A in enumerate(magnus_vectors(X), start=1):
        assert np.allclose(A, rb[f"A{k}"], rtol=1e-10, atol=1e-14 * ps.duration ** k)
    for off, sc in ((2500.0, 1.0), (-4000.0, 1.1)):
        X, dX = propagator_and_gradient(ux, uy, DT, UMAX, off, sc, 0)
        R, _ = rotation_with_gradient(X[0], dX[:, 0])
        q = propagate_batch(ps.as_array(), DT, UMAX * sc, off)
        assert np.allclose(R, quaternion_to_matrix(q), atol=1e-12)


def rectangular(amp, phase, us):
    n = int(round(us * 1e-6 / DT))
    return PulseSequence(np.full(n, amp), np.full(n, phase), DT, UMAX)


def test_local_ur_first_order_robust_inversion():
    """Robust 180°x (first-order detuning) from a rectangular start."""
    cfg = GrapeConfig(mode="local", offsets_hz=[0.0], orders=[1], gate="UR", target_axis=[1, 0, 0], target_angle_deg=180)
    out, info = optimize(rectangular(100, 0, 150), cfg)
    assert info["reached"] and info["stop_reason"] == "solved", info["values"]
    g = metrics.gate(out)
    assert abs(g["angle_deg"] - 180) < 1e-4 and np.max(out.amplitude) <= 100 + 1e-9
    assert metrics.robustness(out)["A1_norm"] < 1e-8


def test_local_two_offsets_and_pp():
    # (120 us is too short for this at 10 kHz: every solver variant ends at the same residual there)
    cfg = GrapeConfig(mode="local", offsets_hz=[-3000.0, 3000.0], orders=[1], gate="PP", m0=[0, 0, 1], m_target=[1, 0, 0])
    out, info = optimize(rectangular(60, 90, 200), cfg)
    assert info["reached"], info["values"]
    for off in (-3000.0, 3000.0):
        R = quaternion_to_matrix(propagate_batch(out.as_array(), DT, UMAX, off))
        assert np.allclose(R @ [0, 0, 1], [1, 0, 0], atol=1e-4)


def test_broadband_pp_improves():
    cfg = GrapeConfig(mode="broadband", gate="PP", m0=[0, 0, 1], m_target=[1, 0, 0], bb_offset_min_hz=-8000,
                      bb_offset_max_hz=8000, bb_offset_n=9)
    out, info = optimize(rectangular(100, 90, 100), cfg, SolverSettings(max_nfev=120))
    before, after = info["start_summary"]["mean infidelity"], info["values"]["mean infidelity"]
    assert after < before / 20, (before, after)
    assert np.max(out.amplitude) <= 100 + 1e-9


def test_stay_close_regularization():
    start = rectangular(70, 30, 80)
    base = GrapeConfig(mode="local", orders=[1], gate="none")
    free, _ = optimize(start, base, SolverSettings(stop_when_solved=False, max_nfev=200))
    reg = GrapeConfig(mode="local", orders=[1], gate="none", regularization="stay_close", regularization_weight=1.0)
    close, _ = optimize(start, reg, SolverSettings(stop_when_solved=False, max_nfev=200))
    dist = lambda p: np.linalg.norm(p.amplitude - start.amplitude) + np.linalg.norm(p.phase - start.phase)
    assert dist(close) < dist(free)


def test_stopping_rules():
    cfg = GrapeConfig(mode="broadband", gate="UR", bb_offset_n=7)
    _, info = optimize(rectangular(100, 0, 100), cfg, SolverSettings(max_nfev=7, stop_when_solved=False))
    assert info["stop_reason"] == "limit" and info["nfev"] == 7
    _, info = optimize(rectangular(100, 0, 100), cfg, callback=lambda n, c: n >= 4)
    assert info["stop_reason"] == "user" and info["nfev"] == 4


def test_spline_basis():
    B = spline_basis(200, 12)
    assert B.shape == (200, 12) and np.allclose(B.sum(axis=1), 1) and np.all(B >= 0)
    x = (np.arange(200) + 0.5) / 200                       # cubic polynomials are reproduced exactly
    c = np.linalg.lstsq(B, x ** 3 - x, rcond=None)[0]
    assert np.allclose(B @ c, x ** 3 - x, atol=1e-12)


def robust_inversion(**kw):
    """Robust 180°x (first order); 200 µs leaves room for smoothness and rate limits."""
    return GrapeConfig(mode="local", offsets_hz=[0.0], orders=[1], gate="UR", target_axis=[1, 0, 0], target_angle_deg=180, **kw)


def test_smoothness_and_rate_limits():
    """From a random start the plain optimization is spiky; smoothness and rate limits give smooth pulses
    that still reach all goals exactly."""
    start = random_pulse(400, seed=7)
    plain, i_plain = optimize(start, robust_inversion())
    smooth, i_smooth = optimize(start, robust_inversion(smooth_weight=1e-4))
    assert i_plain["reached"] and i_smooth["reached"], (i_plain["values"], i_smooth["values"])
    assert i_smooth["values"]["roughness"] < i_plain["values"]["roughness"] / 100, (i_smooth["values"], i_plain["values"])
    assert metrics.robustness(smooth)["A1_norm"] < 1e-8                       # (normalized: |A1|/T)
    limited, info = optimize(start, robust_inversion(max_amp_rate=40, max_phase_rate_khz=40))
    assert info["reached"], info["values"]
    v = info["values"]
    assert v["max amplitude rate [%/µs]"] < 40 * 1.02 and v["max phase rate [kHz]"] < 40 * 1.02, v
    assert i_plain["values"]["max phase rate [kHz]"] > 200                     # (the plain result is not limited)
    assert abs(metrics.gate(limited)["angle_deg"] - 180) < 1e-4


def test_spline_controls():
    out, info = optimize(random_pulse(400, seed=7), robust_inversion(controls="spline", n_spline=16))
    assert info["reached"], info["values"]
    assert GrapeProblem(random_pulse(400), robust_inversion(controls="spline", n_spline=16)).p0.shape == (32,)
    _, plain = optimize(random_pulse(400, seed=7), robust_inversion())
    assert info["values"]["roughness"] < plain["values"]["roughness"] / 20, (info["values"], plain["values"])
    assert abs(metrics.gate(out)["angle_deg"] - 180) < 1e-4


def test_multistart():
    cfg = robust_inversion(n_starts=3, smooth_weight=1e-4)
    start = rectangular(100, 0, 200)
    best, info, all_infos = optimize_multistart(start, cfg)
    assert len(all_infos) == 3 and info["reached"]
    assert len({round(i["values"]["roughness"], 6) for i in all_infos}) == 3          # different starts, different results
    assert info["values"]["roughness"] == min(i["values"]["roughness"] for i in all_infos if i["reached"])
    best2, info2, _ = optimize_multistart(start, cfg, processes=3)                     # parallel: same results
    assert np.allclose(best2.amplitude, best.amplitude) and info2["start_index"] == info["start_index"]


def test_config_validation():
    for bad in (GrapeConfig(mode="broadband", gate="none"), GrapeConfig(mode="local", orders=[], gate="none"),
                GrapeConfig(mode="local", offsets_hz=[]), GrapeConfig(controls="spline", n_spline=3),
                GrapeConfig(max_phase_rate_khz=-1), GrapeConfig(n_starts=0)):
        try:
            bad.check()
        except ValueError:
            continue
        raise AssertionError(bad)
    assert GrapeConfig.from_dict(GrapeConfig(offsets_hz=[1.0, 2.0]).to_dict()).offsets_hz == [1.0, 2.0]


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
