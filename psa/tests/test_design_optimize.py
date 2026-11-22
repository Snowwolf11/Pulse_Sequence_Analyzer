"""Tests of the shared solver (psa.design.optimize)."""
import math
import multiprocessing as mp
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design import metrics
from psa.design.builder import CurveDesign, Segment
from psa.design.optimize import Goals, pulse_metrics, rotation_angle_deg, rotation_from_axis_angle, run_job, solve_problem
from psa.design.waveform import HOLD, Waveform

UMAX, DT = 1e4, 0.5e-6
RMIN = 1e6 / (2 * math.pi * UMAX)


def test_goals_validation():
    try:
        Goals((), None, False)
    except ValueError:
        pass
    else:
        raise AssertionError("empty goals accepted")
    g = Goals((1, 2), ((1, 0, 0), 180.0), True)
    assert Goals.from_dict(g.to_dict()) == g


def test_settings_weights_and_tolerances():
    from psa.design.optimize import SolverSettings, residual, summary
    st = SolverSettings(weight_order1=10.0, weight_gate=0.5, solved_tolerance=1e-3, gate_tolerance_deg=2.0)
    assert SolverSettings.from_dict(st.to_dict()) == st and SolverSettings.from_dict({"unknown": 1}) == SolverSettings()
    m = {"A": {1: np.array([1e-4, 0, 0]), 2: np.zeros(3), 3: np.zeros(3)}, "T": 1.0,
         "R": rotation_from_axis_angle((1, 0, 0), 181.0)}
    g = Goals((1,), ((1, 0, 0), 180.0))
    r = residual(m, g, st)
    assert np.allclose(r[:3], [1e-3, 0, 0]) and abs(np.linalg.norm(r[3:]) - 0.5 * math.sqrt(2) * math.sin(math.radians(0.5))) < 1e-12
    assert summary(m, g, st)["reached"] and not summary(m, g)["reached"]


def _hard_design(n=7):
    return CurveDesign([Segment("arc", {"radius": 1.3 * RMIN, "angle": 100 + 37 * k}, twist=60 * k,
                                free=(["angle", "twist", "radius"] if k else [])) for k in range(n)], UMAX, DT)


def test_stop_when_solved():
    from psa.design.optimize import SolverSettings
    d = _hard_design(4)
    on = solve_problem("builder", d.to_dict(), None, Goals((1, 2)), settings=SolverSettings())[1]
    off = solve_problem("builder", d.to_dict(), None, Goals((1, 2)), settings=SolverSettings(stop_when_solved=False))[1]
    assert on["reached"] and on["stop_reason"] == "solved"
    assert off["reached"] and off["stop_reason"] == "converged"
    assert on["nfev"] < off["nfev"], (on["nfev"], off["nfev"])
    assert all(v < 1e-8 for v in on["values"].values())


def test_total_evaluation_limit_and_time_limit():
    import time
    from psa.design.optimize import SolverSettings
    d = _hard_design(7)
    for limit in (10, 37):
        info = solve_problem("builder", d.to_dict(), None, Goals((1, 2, 3)), settings=SolverSettings(max_nfev=limit))[1]
        assert info["nfev"] == limit and info["stop_reason"] == "limit", info
        assert info["values"]["|A1|/T^1"] <= info["start_summary"]["|A1|/T^1"] + 1e-15 or info["cost"] < np.inf
    t0 = time.time()
    info = solve_problem("builder", d.to_dict(), None, Goals((1, 2, 3)),
                         settings=SolverSettings(time_limit_s=0.3, max_nfev=10 ** 7, stop_when_solved=False,
                                                 xtol=0, ftol=0, gtol=0))[1]
    assert info["stop_reason"] == "time" and time.time() - t0 < 2.0, info
    user = solve_problem("builder", d.to_dict(), None, Goals((1, 2, 3)), callback=lambda n, c: n >= 12)[1]
    assert user["stop_reason"] == "user" and user["cancelled"] and user["nfev"] == 12


def test_default_tolerances():
    from psa.design.optimize import SolverSettings
    s = SolverSettings()
    assert (s.xtol, s.ftol, s.gtol) == (1e-12, 1e-12, 1e-12) and s.stop_when_solved and s.time_limit_s == 0


def test_builder_target_gate_with_closure():
    """Three arcs (CORPSE-like) with free angles and twists -> robust 180° x rotation."""
    d = CurveDesign([Segment("arc", {"radius": RMIN, "angle": 60}, free=["angle"]),
                     Segment("arc", {"radius": RMIN, "angle": 300}, twist=180, free=["angle", "twist"]),
                     Segment("arc", {"radius": RMIN, "angle": 60}, twist=180, free=["angle", "twist"])], UMAX, DT,
                    phase0=0.0)
    solved, res = d.solve(orders=(1,), target=((1, 0, 0), 180.0))
    assert res["reached"], res["values"]
    R = solved.gate_matrix()
    assert rotation_angle_deg(rotation_from_axis_angle((1, 0, 0), 180).T @ R) < 1e-5
    assert solved.robustness()["A1_norm"] < 1e-8


def test_builder_b1_matches_discretized_pulse():
    d = CurveDesign([Segment("arc", {"radius": 20, "angle": 140}),
                     Segment("helix", {"radius": 18, "pitch": 30, "angle": 300}, twist=40),
                     Segment("arc", {"radius": 35, "angle": 200}, twist=-75)], UMAX, 0.01e-6)
    m = d.metrics(b1=True)
    ab = metrics.amplitude_robustness(d.to_pulse())
    assert np.linalg.norm(m["B1"] - ab["B1"]) < 2e-3 * np.linalg.norm(ab["B1"]), (m["B1"], ab["B1"])


def test_builder_b1_goal():
    """Two loops of opposite orientation can be made B1 robust (like the 4H construction)."""
    d = CurveDesign([Segment("arc", {"radius": 1.2 * RMIN, "angle": 360}),
                     Segment("arc", {"radius": 1.2 * RMIN, "angle": 300}, twist=160, free=["angle", "twist"])], UMAX, DT)
    solved, res = d.solve(orders=(), b1=True)
    assert res["reached"], res["values"]


def waveform_design():
    w = Waveform(N=240, dt=DT, umax=UMAX)
    w.amplitude.set(0, 80.0, HOLD); w.amplitude.set(80, 90.0, HOLD); w.amplitude.set(160, 70.0, HOLD)
    w.phase.set(0, 0.0, HOLD); w.phase.set(80, 120.0, HOLD); w.phase.set(160, 30.0, HOLD)
    return w


def test_waveform_target_gate():
    w = waveform_design()
    free = [("amplitude", 0), ("amplitude", 80), ("amplitude", 160), ("phase", 80), ("phase", 160)]
    sol, info = solve_problem("waveform", w.to_dict(), free, Goals((), ((0, 1, 0), 90.0)))
    assert info["reached"], info["values"]
    ps = Waveform.from_dict(sol).to_pulse()
    g = metrics.gate(ps)
    assert abs(g["angle_deg"] - 90) < 1e-4 and np.allclose(g["axis"], [0, 1, 0], atol=1e-5)
    assert sol["phase"][0][1] == 0.0                                   # keyframe not listed as free: unchanged


def test_waveform_closure_improves():
    w = waveform_design()
    free = [(c, i) for c in ("amplitude", "phase") for i in (0, 80, 160)]
    sol, info = solve_problem("waveform", w.to_dict(), free, Goals((1,)))
    start = info["start_summary"]["|A1|/T^1"]
    assert info["values"]["|A1|/T^1"] < 1e-8 or info["values"]["|A1|/T^1"] < start / 100, (start, info["values"])
    ps = Waveform.from_dict(sol).to_pulse()
    assert np.all(ps.amplitude <= 100 + 1e-9) and np.all(ps.amplitude >= 0)


def test_family_parameter_closes_custom_curve():
    design = {"type": "family", "family": "Custom (formula)",
              "values": {"u0": 0.0, "u1": 2 * math.pi, "a": 0.3, "b": 1.0, "c": 2.0, "d": 0.4},
              "expressions": ["cos(u)", "sin(b*u) + a*u", "d*sin(c*u)"], "dt": DT, "umax": UMAX, "duration": 400e-6}
    sol, info = solve_problem("family", design, ["a"], Goals((1,)))
    assert abs(sol["values"]["a"]) < 1e-3, sol["values"]["a"]
    assert info["values"]["|A1|/T^1"] < 1e-4


def test_run_job_in_a_process():
    ctx = mp.get_context("spawn")
    q, stop = ctx.Queue(), ctx.Event()
    d = CurveDesign([Segment("arc", {"radius": 20, "angle": 90}),
                     Segment("arc", {"radius": 20, "angle": 200}, twist=150, free=["angle", "twist"]),
                     Segment("arc", {"radius": 20, "angle": 120}, twist=150, free=["angle", "twist"])], UMAX, DT)
    job = {"kind": "builder", "design": d.to_dict(), "free": None, "goals": Goals((1,)).to_dict()}
    p = ctx.Process(target=run_job, args=(job, q, stop))
    p.start()
    msgs = []
    while True:
        m = q.get(timeout=120)
        msgs.append(m)
        if m[0] in ("done", "error"):
            break
    p.join(10)
    assert msgs[-1][0] == "done", msgs[-1]
    assert any(m[0] == "progress" for m in msgs)
    assert msgs[-1][2]["reached"]


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
