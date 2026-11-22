"""Tests of the curve / pulse families (psa.design.families)."""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design import metrics
from psa.design.conversion import curve_to_ps, ps_to_curve
from psa.design.families import FAMILIES, CustomCurve, FourHelix, TwoHelix, compile_expression, family_by_name

DT, UMAX = 0.5e-6, 1e4


def test_curve_families_convert_and_reproduce():
    for fam in FAMILIES:
        if fam.kind != "curve":
            continue
        pts = fam.points(fam.defaults())
        assert pts.shape[1] == 3 and np.all(np.isfinite(pts)), fam.name
        ps, info = curve_to_ps(pts, DT, UMAX)
        assert info["amplitude_ok"], fam.name
        err = np.max(np.linalg.norm(ps_to_curve(ps) - info["target_points"], axis=1)) / info["duration"]
        assert err < 2e-4, (fam.name, err)


def test_closed_curve_families_are_first_order_robust():
    """Closed families (Lissajous, torus knot) give A1 ~ 0 after conversion."""
    for name in ("Lissajous knot", "Torus knot"):
        fam = family_by_name(name)
        ps, _ = curve_to_ps(fam.points(fam.defaults()), DT, UMAX)
        assert metrics.robustness(ps)["A1_norm"] < 1e-4, name


def test_check_clamps_and_rounds():
    v = family_by_name("Torus knot").check({"p": 2.7, "q": 99, "R": -1})
    assert v["p"] == 3 and v["q"] == 20 and v["R"] == 0.1


def test_two_helix_closes():
    fam = TwoHelix()
    for tau in (55.0, 75.0, 95.0):
        ps = fam.pulse({"tau": tau}, DT, UMAX)
        assert ps.N == 2 * round(tau * 1e-6 / DT) and np.allclose(ps.amplitude, 100)
        A1 = metrics.robustness(ps)["A1_norm"]
        fine = metrics.robustness(fam.pulse({"tau": tau}, DT / 10, UMAX))["A1_norm"]
        assert A1 < 2e-3 and fine < A1 / 20, (tau, A1, fine)          # midpoint sampling: O(dt^2)


def test_four_helix_is_second_order_robust_identity():
    fam = FourHelix()
    try:
        fam.pulse({"tau": 125.0}, DT, UMAX)              # between the solution windows
    except ValueError:
        pass
    else:
        raise AssertionError("expected no 4H solution at 125 us")
    for tau in (55.0, 75.0, 95.0, 175.0):
        ps = fam.pulse({"tau": tau}, DT / 5, UMAX)
        rb = metrics.robustness(ps)
        g = metrics.gate(ps)
        assert rb["A1_norm"] < 1e-4 and rb["A2_norm"] < 1e-4, (tau, rb["A1_norm"], rb["A2_norm"])
        assert g["angle_deg"] < 0.1, (tau, g["angle_deg"])


def test_four_helix_closure_condition():
    """The returned rate is a root of the closure condition and gives a closed 2-helix loop."""
    Omega, tau = 2 * math.pi * UMAX, 60e-6
    p = FourHelix.phase_rate(Omega, tau)
    w = math.sqrt(Omega ** 2 + p ** 2)
    assert abs(Omega ** 2 * math.sin(tau * w) + tau * p ** 2 * w) < 1e-6 * Omega ** 2


def test_custom_expressions():
    fam = CustomCurve()
    pts = fam.points(fam.defaults())
    u = np.linspace(0, 2 * np.pi, len(pts))
    assert np.allclose(pts[:, 0], np.cos(u)) and np.allclose(pts[:, 2], 0.5 * np.sin(3 * u))
    pts = fam.points({"u0": 0, "u1": 1}, expressions=("u", "u^2", "3"))        # ^ means power, constants broadcast
    assert np.allclose(pts[:, 1], pts[:, 0] ** 2) and np.allclose(pts[:, 2], 3)
    f = compile_expression("a*sqrt(u) + pi", ("u", "a"))
    assert abs(f(u=4.0, a=2.0) - (4 + np.pi)) < 1e-12


def test_custom_expressions_reject_unsafe_input():
    bad = ["__import__('os')", "u.real", "lambda: 1", "open('x')", "x", "sin", "[u]", "'text'", "u if u else 1"]
    for text in bad:
        try:
            compile_expression(text, ("u",))
        except ValueError:
            continue
        raise AssertionError(f"accepted unsafe expression {text!r}")
    try:
        CustomCurve().points({"u0": -1, "u1": 1}, expressions=("u", "log(u)", "u"))
    except ValueError:
        pass
    else:
        raise AssertionError("non-finite values were not reported")


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
