"""Tests of the shortest robust pulse search (psa.design.shortest)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.metrics import gate  # noqa: E402
from psa.design.pulse import PulseSequence  # noqa: E402
from psa.design.shortest import ShortestSpec, corpse, fit_steps, search, solve, start_pulses, verify  # noqa: E402

DT, UMAX = 0.5e-6, 1e4


def test_verify_and_bounds():
    spec = ShortestSpec(angle_deg=90.0, order=0)
    v = verify(PulseSequence(np.full(50, 100.0), np.zeros(50), DT, UMAX), spec)      # rectangular 90° x: exact
    assert v["ok"] and v["gate error"] < 1e-14 and v["A1"] > 0.1
    assert not verify(PulseSequence(np.full(50, 100.0), np.zeros(50), DT, UMAX), ShortestSpec(angle_deg=90.0, order=1))["ok"]
    assert np.isclose(spec.lower_bound_s(), 25e-6) and np.isclose(ShortestSpec(order=1).lower_bound_s(), 100e-6)
    pp = ShortestSpec(gate="PP", m0=[0, 0, 1], mt=[0, -1, 0], order=0)              # z -> -y by 90° about x
    assert verify(PulseSequence(np.full(50, 100.0), np.zeros(50), DT, UMAX), pp)["ok"]


def test_fit_steps_keeps_the_gate():
    ps = PulseSequence(np.full(80, 50.0), np.full(80, 30.0), DT, UMAX)
    for n in (40, 100):                                                             # squeeze (amp 100 %) / stretch
        g0, g1 = gate(ps), gate(fit_steps(ps, n))
        assert fit_steps(ps, n).N == n and np.isclose(g0["angle_deg"], g1["angle_deg"]) and np.allclose(g0["axis"], g1["axis"])
    c = corpse(180.0, 0.0, DT / 100, UMAX)                                         # fine sampling: nearly exact
    v = verify(c, ShortestSpec(order=1))
    assert v["gate error"] < 1e-6 and v["A1"] < 1e-4 and np.isclose(c.duration, 7 / 6 / UMAX, rtol=1e-3)
    names = [n for n, _ in start_pulses(ShortestSpec(), 240)]
    assert "CORPSE" in names and "4H" in names and all(p.N == 240 for _, p in start_pulses(ShortestSpec(), 240))


def test_solve_first_order_with_b1():
    spec = ShortestSpec(angle_deg=180.0, order=1, b1=True)
    rng = np.random.default_rng(1)
    out, v = solve(PulseSequence(rng.uniform(30, 90, 400), rng.uniform(-180, 180, 400), DT, UMAX), spec, seed=1)
    assert v["ok"] and v["A1"] < 1e-8 and v["B1"] < 1e-8 and v["gate error"] < 1e-10, v


def test_search_gate_only_reaches_the_bound():
    res = search(ShortestSpec(angle_deg=90.0, order=0), 35.0, n_random=1, perturbed=0, time_limit_s=120)
    v = res["verification"]
    assert v["ok"] and v["N"] == 50 and res["message"] == "reached the lower bound"       # 25 µs = 90° at 10 kHz
    assert [h["ok"] for h in res["history"]][0] and np.isclose(res["lower bound [µs]"], 25.0)


def test_search_stop_keeps_the_best_and_reports_every_pulse():
    from psa.design.shortest import SearchStopped
    calls, accepted = [], []

    def mapper(f, it):                         # the third duration is stopped
        calls.append(1)
        if len(calls) == 3:
            raise SearchStopped()
        return list(map(f, it))
    res = search(ShortestSpec(angle_deg=90.0, order=0), 60.0, n_random=1, perturbed=0, time_limit_s=120,
                 mapper=mapper, on_best=lambda p, v: accepted.append(p.N))
    assert res["message"] == "stopped" and res["pulse"] is not None and res["verification"]["ok"]
    assert accepted == [120, 108] and res["pulse"].N == 108


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
