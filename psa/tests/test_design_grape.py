"""Tests of the GRAPE optimizer (psa.design.grape): exact gradients, consistency, optimizations."""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design import metrics
from psa.design.grape import (GrapeConfig, GrapeProblem, Search, control_noise, exploration_preset, optimize,
                               optimize_multistart, propagator_and_gradient, rotation_with_gradient, run_task, search,
                               spline_basis, stage_plan, lawson_update, band_infidelity,
                               band_quality, area_quality)
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
    "local weighted points at scalings": GrapeConfig(mode="local", points=[[-1e3, 0.9, 1.0], [2e3, 1.1, 3.0]],
                                                     orders=[1, 2], b1=True, gate="UR"),
    "broadband weighted points": GrapeConfig(mode="broadband", points=[[-1e3, 0.95, 2.0], [0.0, 1.0, 1.0], [4e3, 1.05, 0.5]],
                                             gate="PP"),
}


def test_gradients_match_finite_differences():
    """Jacobian with respect to the optimization parameters (through squash map and spline basis)."""
    ps = random_pulse(24)
    for name, cfg in CONFIGS.items():
        prob = GrapeProblem(ps, cfg)
        rng = np.random.default_rng(1)
        if cfg.n_points() > 1:                                   # Lawson weights of the points
            prob.point_weights = rng.uniform(0.2, 3.0, cfg.n_points())
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


def test_resumable_optimize():
    """A run in small budgets (paused, resumed from its state) reaches the goals like one long run; the
    evaluation count adds up; a state may also give the controls w of a new candidate."""
    cfg = robust_inversion(smooth_weight=1e-4)
    start = rectangular(100, 0, 200)
    st = SolverSettings(max_nfev=400)
    state, calls, total = None, 0, 0
    while True:
        out, info = optimize(start, cfg, st, state=state, budget=5)
        calls += 1
        assert info["nfev"] - (state["nfev"] if state else 0) <= 5
        total = info["nfev"]
        if info["finished"]:
            break
        assert info["stop_reason"] == "paused" and not info["cancelled"]
        state = info["state"]
    assert info["reached"] and calls > 2 and total <= 400 and info["stage"] == info["n_stages"] - 1
    w = np.asarray(info["w"])
    out2, info2 = optimize(start, cfg, st, state={"stage": info["n_stages"] - 1, "w": w.tolist()})
    assert info2["reached"] and info2["nfev"] <= 5                   # already a solution


def test_coarse_to_fine():
    cfg = robust_inversion(coarse_levels=[8, 16, 5000])
    plan = stage_plan(cfg, 300)
    assert [(c, m) for c, m, _, _ in plan] == [("spline", 8), ("spline", 16), ("full", 20)]
    assert [m for _, m, _, _ in stage_plan(robust_inversion(coarse_levels=[8, 16], controls="spline", n_spline=12), 300)] == [8, 12]
    out, info = optimize(rectangular(100, 0, 150), cfg, SolverSettings(max_nfev=400))
    assert info["reached"] and info["stage"] == 2 and out.N == 300
    assert metrics.robustness(out)["A1_norm"] < 1e-8


def test_control_noise():
    rng = np.random.default_rng(0)
    white = control_noise(1000, DT, 0.5, 0.0, rng)
    smooth = control_noise(1000, DT, 0.5, 20.0, rng)
    for z in (white, smooth):
        assert z.shape == (2000,) and abs(np.sqrt(np.mean(z[:1000] ** 2)) - 0.5) < 0.05       # RMS = size
    rough = lambda z: np.mean(np.diff(z[:1000]) ** 2)
    assert rough(smooth) < rough(white) / 100


def test_search_screening_crossover_kicks():
    """Successive halving (6 -> 3 -> 2 with crossover children), then kicks (smoothing on: smoother wins);
    serial and parallel runs give the same result."""
    cfg = robust_inversion(smooth_weight=1e-4, n_starts=6, start_pool="mix", noise_correlation_us=10.0,
                           screen_rounds=2, screen_nfev=8, keep_fraction=0.5, crossover_rate=1.0, n_kicks=2,
                           kick_size=0.2, kick_temperature=0.3)
    start = rectangular(100, 0, 150)
    st = SolverSettings(max_nfev=200)
    s = Search(start, cfg, st, parallel=2)
    kinds = []
    while not s.done:
        for t in s.tasks():
            kinds.append(s.cands[t["cid"]]["kind"])
            s.submit(run_task(t))
    cid, best, info = s.best()
    assert info["reached"] and kinds.count("kick") == 2 and "child" in kinds
    assert s.log[0].startswith("round 1: 6 candidates, 8 evaluations each → kept 3 + 3 crossover children"), s.log
    assert s.log[1].startswith("round 2: 6 candidates, 16 evaluations each → kept 3 + 3"), s.log
    assert info["search"]["candidates"] == 6 + 3 + 3 + 2 and info["search"]["kicks"] == 2
    reached = [i for i in s.infos() if i["reached"]]
    assert info["values"]["roughness"] == min(i["values"]["roughness"] for i in reached)
    best2, info2, _ = search(start, cfg, st, processes=3)
    assert np.allclose(best2.amplitude, best.amplitude) and info2["search"]["best"] == cid
    # reached without smoothing: no kicks needed, the search ends after the round that solved it
    _, info3, _ = search(start, robust_inversion(n_starts=2, n_kicks=3), st)
    assert info3["reached"] and info3["search"]["kicks"] == 0


def test_exploration_presets():
    for name in ("exploit", "balanced", "explore"):
        d = exploration_preset(name, 200e-6)
        robust_inversion(**d).check()
    assert exploration_preset("explore", 200e-6)["noise_correlation_us"] == 20.0


def test_lawson_reweighting():
    """Weights grow with the point errors (normalized to mean 1); reweight stages only with several points;
    the worst point gets better, the errors more even (local and broadband)."""
    w = lawson_update(np.ones(3), [1.0, 2.0, 4.0], 1.0)
    assert abs(np.mean(w) - 1) < 1e-12 and w[0] < w[1] < w[2] and abs(w[2] / w[0] - 4) < 1e-9
    assert [k for *_, k in stage_plan(robust_inversion(reweight_rounds=3), 100)] == ["main"]          # one offset
    assert [k for *_, k in stage_plan(robust_inversion(reweight_rounds=2).replace(offsets_hz=[-1e3, 1e3]), 100)] == \
        ["main", "reweight", "reweight"]
    st = SolverSettings(max_nfev=1500)
    bb = GrapeConfig(mode="broadband", gate="PP", m0=[0, 0, 1], m_target=[1, 0, 0], bb_offset_min_hz=-10000,
                     bb_offset_max_hz=10000, bb_offset_n=7)
    _, plain = optimize(rectangular(100, 0, 60), bb, st)
    _, even = optimize(rectangular(100, 0, 60), bb.replace(reweight_rounds=10), st)
    assert not plain["reached"] and plain["point_weights"] is None and len(even["point_weights"]) == 7
    assert max(even["point_errors"]) < 0.95 * max(plain["point_errors"])                  # smaller worst case
    assert even["values"]["point error spread"] < plain["values"]["point error spread"] - 0.1
    # local mode (offsets are the points): the rounds run and never make the worst offset worse
    loc = GrapeConfig(mode="local", offsets_hz=[-1000.0, 4000.0], orders=[1], gate="PP", m0=[0, 0, 1], m_target=[1, 0, 0])
    _, plain = optimize(rectangular(100, 0, 30), loc, st)
    _, even = optimize(rectangular(100, 0, 30), loc.replace(reweight_rounds=5), st)
    assert plain["stop_reason"] == "converged" and even["nfev"] > plain["nfev"] and len(even["point_weights"]) == 2
    assert max(even["point_errors"]) <= max(plain["point_errors"]) * (1 + 1e-9)

def test_band_quality_matches_broadband_goal():
    """The test-grid scorer (no gradients) gives the same infidelities as the broadband goal at its points."""
    ps = random_pulse(50, seed=4)
    for cfg in (GrapeConfig(mode="broadband", gate="UR", target_axis=[0, 1, 1], target_angle_deg=120, bb_offset_n=5),
                GrapeConfig(mode="broadband", gate="PP", m0=[0, 0, 1], m_target=[0, 1, 0], bb_offset_n=5)):
        summ = GrapeProblem(ps, cfg).goals(GrapeProblem(ps, cfg).x0)[2]
        inf = band_infidelity(ps, cfg, np.linspace(cfg.bb_offset_min_hz, cfg.bb_offset_max_hz, 5))
        assert abs(np.max(inf) - summ["max infidelity"]) < 1e-12 and abs(np.mean(inf) - summ["mean infidelity"]) < 1e-12
    q = band_quality(ps, cfg, -5000, 5000, n=101)
    assert q["worst infidelity"] >= q["mean infidelity"] and 0 <= q["fraction below 0.001"] <= 1


def test_limits_late():
    """Rate limits in a last stage only: the plan gets a "limits" stage (before the reweight rounds); the result
    keeps the limits and reaches the goals."""
    cfg = robust_inversion(max_amp_rate=5, max_phase_rate_khz=20, limits_late=True)
    assert [k for *_, k in stage_plan(cfg, 300)] == ["main", "limits"]
    two = cfg.replace(offsets_hz=[-1e3, 1e3], reweight_rounds=2)
    assert [k for *_, k in stage_plan(two, 300)] == ["main", "limits", "reweight", "reweight"]
    start = random_pulse(300, seed=2)
    out, info = optimize(start, cfg, SolverSettings(max_nfev=1500))
    assert info["reached"] and info["stage"] == 1, info["message"]
    assert info["values"]["rate-limit excess"] < 1e-3 and metrics.robustness(out)["A1_norm"] < 1e-8


def test_config_validation():
    for bad in (GrapeConfig(mode="broadband", gate="none"), GrapeConfig(mode="local", orders=[], gate="none"),
                GrapeConfig(mode="local", offsets_hz=[]), GrapeConfig(controls="spline", n_spline=3),
                GrapeConfig(max_phase_rate_khz=-1), GrapeConfig(n_starts=0),
                GrapeConfig(start_pool="x"), GrapeConfig(keep_fraction=1.0), GrapeConfig(crossover_rate=2),
                GrapeConfig(coarse_levels=[3]), GrapeConfig(reweight_power=0), GrapeConfig(reweight_rounds=-1), GrapeConfig(screen_nfev=1), GrapeConfig(n_kicks=-1)):
        try:
            bad.check()
        except ValueError:
            continue
        raise AssertionError(bad)
    assert GrapeConfig.from_dict(GrapeConfig(offsets_hz=[1.0, 2.0]).to_dict()).offsets_hz == [1.0, 2.0]


def _goals(ps, **kw):
    prob = GrapeProblem(ps, GrapeConfig(**kw))
    r, J, _ = prob.goals(prob.x0)
    return r, J, prob


def test_explicit_points_reproduce_the_default_points():
    """points = the default points (weights 1, or all equal) gives exactly the old residuals in both modes."""
    ps = random_pulse(50, 4)
    loc = dict(mode="local", offsets_hz=[-2e3, 0.0, 1e3], orders=[1, 2], b1=True, gate="UR")
    bb = dict(mode="broadband", bb_offset_min_hz=-4e3, bb_offset_max_hz=4e3, bb_offset_n=5, bb_scale_min=0.9,
              bb_scale_max=1.1, bb_scale_n=3, gate="PP")
    for d in (loc, bb):
        r0, J0, p0 = _goals(ps, **d)
        pts = [[o, s, 3.0] for o, s, _ in p0.cfg.point_list()]
        r1, J1, p1 = _goals(ps, **dict(d, points=pts, offsets_hz=[7.0], bb_offset_n=2))   # old fields ignored
        assert np.array_equal(r0, r1) and np.array_equal(J0, J1) and p1.user_w is None
        assert p1.cfg.n_points() == len(pts)


def test_point_weights_and_scalings():
    ps = random_pulse(50, 5)
    pts = [[-2e3, 1.0, 1.0], [0.0, 1.0, 1.0], [2e3, 1.0, 4.0]]          # normalized: 0.5, 0.5, 2
    r0, _, p0 = _goals(ps, mode="broadband", points=[p[:2] for p in pts], gate="UR")
    r1, _, p1 = _goals(ps, mode="broadband", points=pts, gate="UR")
    n = len(r0) // 3
    assert np.allclose(r1[:2 * n], r0[:2 * n] * math.sqrt(0.5)) and np.allclose(r1[2 * n:], r0[2 * n:] * math.sqrt(2))
    assert np.allclose(p1.point_errors, p0.point_errors * np.sqrt([0.5, 0.5, 2]))   # Lawson sees weighted errors
    # local goals at amplitude scaling s = the goals of the pulse with s x the amplitude at scaling 1
    s = 0.85
    scaled = PulseSequence(ps.amplitude * s, ps.phase, ps.dt, ps.umax)
    for kw in (dict(orders=[1, 2, 3], b1=False, gate="UR"), dict(orders=[1], b1=True, gate="PP")):
        ra, _, _ = _goals(ps, mode="local", points=[[1e3, s]], **kw)
        rb, _, _ = _goals(scaled, mode="local", offsets_hz=[1e3], **kw)
        assert np.allclose(ra, rb, atol=1e-12), kw
    for bad in ([[0.0, -1.0]], [[0.0, 1.0, 0.0]], [[0.0, 1.0, 1.0, 5.0]]):
        try:
            GrapeConfig(mode="broadband", points=bad).check()
        except ValueError:
            continue
        raise AssertionError(f"accepted points {bad}")


def test_weighted_points_optimize():
    """An optimization with weighted explicit points runs (gradients consistent) and favours the heavy point."""
    ps = PulseSequence(np.full(200, 100.0), np.zeros(200), DT, UMAX)
    st = SolverSettings(max_nfev=150)
    res = {}
    for w in (1.0, 50.0):
        cfg = GrapeConfig(mode="broadband", points=[[-3e3, 1.0, 1.0], [0.0, 1.0, 1.0], [3e3, 1.0, w]], gate="UR",
                          target_angle_deg=90, seed=1)
        out, info = optimize(ps, cfg, st)
        res[w] = band_infidelity(out, cfg, [3e3])[0]
    assert res[50.0] <= res[1.0] * 1.01 + 1e-12, res


def test_one_goal_builder_normalized_by_the_points():
    """Gate only gives the same residual in both modes; all goal rows are divided by sqrt(number of points),
    so duplicating every point leaves the squared residual unchanged (a mean, not a sum)."""
    ps = random_pulse(40, 7)
    pts = [[-2e3, 0.95, 1.0], [1e3, 1.05, 1.0]]
    for gate in ("UR", "PP"):
        rl, Jl, _ = _goals(ps, mode="local", points=pts, orders=[], gate=gate)
        rb, Jb, _ = _goals(ps, mode="broadband", points=pts, gate=gate)
        assert np.allclose(rl, rb, rtol=1e-14, atol=1e-16) and np.allclose(Jl, Jb, rtol=1e-14, atol=1e-16)
    for kw in (dict(mode="local", orders=[1, 2], b1=True, gate="UR"), dict(mode="broadband", gate="PP")):
        r1, _, _ = _goals(ps, points=pts, **kw)
        r2, _, _ = _goals(ps, points=pts + pts, **kw)
        assert abs(r1 @ r1 - r2 @ r2) < 1e-12 * (r1 @ r1), kw
    # weights: sqrt(w_i / sum w); doubling all weights changes nothing
    r1, _, _ = _goals(ps, mode="local", points=[p[:2] + [1.0] for p in pts], orders=[1], gate="UR")
    r2, _, _ = _goals(ps, mode="local", points=[p[:2] + [2.0] for p in pts], orders=[1], gate="UR")
    assert np.array_equal(r1, r2)


def test_unified_solved_rule():
    """Gate: worst infidelity at the points below gate_infidelity_tolerance (both modes); the gate error in
    degrees is only information."""
    ps = random_pulse(40, 8)
    for mode in ("local", "broadband"):
        prob = GrapeProblem(ps, GrapeConfig(mode=mode, points=[[0.0, 1.0], [2e3, 1.0]], orders=[], gate="UR"))
        _, _, summ = prob.goals(prob.x0)
        assert {"max infidelity", "mean infidelity", "gate error [°]"} <= set(summ)
        assert summ["max infidelity"] >= summ["mean infidelity"] > 0
        R = quaternion_to_matrix(propagate_batch(ps.as_array(), DT, UMAX, np.array([0.0, 2e3])))
        Rt = prob.Rt
        assert abs(summ["max infidelity"] - max(1 - (3 + np.trace(Rt.T @ r)) / 6 for r in R)) < 1e-14
        prob.st.gate_infidelity_tolerance = summ["max infidelity"] * 1.001
        assert prob.reached(summ)
        prob.st.gate_infidelity_tolerance = summ["max infidelity"] * 0.999
        prob.st.gate_tolerance_deg = 1e9                     # the degree tolerance is not used by GRAPE
        assert not prob.reached(summ)


def test_area_quality():
    ps = random_pulse(60, 9)
    cfg = GrapeConfig(mode="broadband", points=[[-3e3, 0.9], [2e3, 1.1], [0.0, 1.0]], gate="UR")
    q = area_quality(ps, cfg, n_offsets=31, n_b1=5)
    assert q["grid"] == (5, 31) and -3 <= q["worst at offset [kHz]"] <= 2 and 90 <= q["worst at B1 [%]"] <= 110
    inf = np.array([band_infidelity(ps, cfg, np.linspace(-3e3, 2e3, 31), scale=s) for s in np.linspace(0.9, 1.1, 5)])
    assert abs(q["worst infidelity"] - inf.max()) < 1e-14 and abs(q["mean infidelity"] - inf.mean()) < 1e-14
    q1 = area_quality(ps, cfg.replace(points=[[-3e3, 1.0], [2e3, 1.0]]), n_offsets=31)
    assert q1["grid"] == (1, 31)                                         # one B1 value
    assert area_quality(ps, cfg.replace(points=[[0.0, 0.9], [0.0, 1.1]])) is None        # no offset range
    assert area_quality(ps, cfg.replace(mode="local", gate="none", orders=[1])) is None  # no gate goal


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
