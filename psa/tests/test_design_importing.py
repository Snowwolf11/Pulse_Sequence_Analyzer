"""Tests of curve import (psa.design.importing)."""
import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.conversion import align_curve, curve_to_ps, ps_to_curve
from psa.design.importing import curve_from_curvature_torsion, import_curve, load_table

DT, UMAX = 0.5e-6, 1e4


def write(text):
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


def test_load_table_separators_and_headers():
    for text in ("x,y,z\n1,2,3\n4,5,6\n", "# comment\n1;2;3\n4;5;6\n", "1 2 3\n4\t5\t6\n\n", "##TITLE=a\n1, 2, 3\n4,5,6\n"):
        t = load_table(write(text))
        assert t.shape == (2, 3) and np.allclose(t, [[1, 2, 3], [4, 5, 6]]), text
    t = load_table(write("1e-3,-2.5,+.5\n"))
    assert np.allclose(t, [[1e-3, -2.5, 0.5]])


def test_frenet_integration_helix():
    R, c = 2.0, 0.5
    kappa, tau = R / (R * R + c * c), c / (R * R + c * c)
    s = np.linspace(0, 4 * math.pi * math.sqrt(R * R + c * c), 4001)   # two turns
    pts = curve_from_curvature_torsion(s, np.full_like(s, kappa), np.full_like(s, tau))
    # helix property: distance between points one turn apart = 2 pi c, independent of the start
    n_turn = 2000
    d = np.linalg.norm(pts[n_turn:] - pts[:-n_turn], axis=1)
    assert np.allclose(d, 2 * math.pi * c, rtol=1e-5)
    # numerical curvature and torsion of the result
    d1 = np.gradient(pts, axis=0); d2 = np.gradient(d1, axis=0); d3 = np.gradient(d2, axis=0)
    cr = np.cross(d1, d2)
    k_num = np.linalg.norm(cr, axis=1) / np.linalg.norm(d1, axis=1) ** 3
    t_num = np.einsum("ij,ij->i", cr, d3) / np.einsum("ij,ij->i", cr, cr)
    assert np.allclose(np.median(k_num), kappa, rtol=1e-4) and np.allclose(np.median(t_num), tau, rtol=1e-3)


def test_frenet_with_zero_curvature_and_jumps():
    # straight - arc - straight: kappa piecewise, tau = 0
    s = np.linspace(0, 6, 6001)
    kappa = np.where((s > 2) & (s < 2 + math.pi), 1.0, 0.0)              # half circle of radius 1
    pts = curve_from_curvature_torsion(s, kappa, np.zeros_like(s))
    assert np.allclose(pts[2000], [0, 0, 2], atol=1e-9)                     # first straight part along z
    end = pts[-1]
    rest = 6 - 2 - math.pi                                                 # straight length after the U-turn
    assert abs(end[2] - (2.0 - rest)) < 5e-3 and abs(abs(end[0]) - 2.0) < 5e-3   # back down, 2 units aside
    ps, info = curve_to_ps(pts, DT, UMAX)
    assert np.max(ps.amplitude[: int(0.25 * ps.N)]) < 1.0                  # straight part: no RF


def test_import_points_and_tables_match():
    """The same helix given as points or as curvature/torsion gives the same pulse."""
    R, c = 1.0, 0.3
    u = np.linspace(0, 6 * math.pi, 6000)
    pts = np.c_[R * np.cos(u), R * np.sin(u), c * u]
    p_points = write("\n".join(f"{x},{y},{z}" for x, y, z in pts))
    L = 6 * math.pi * math.sqrt(R * R + c * c)
    s = np.linspace(0, L, 6000)
    k, t = R / (R * R + c * c), c / (R * R + c * c)
    p_table = write("s,kappa,tau\n" + "\n".join(f"{a},{k},{t}" for a in s))
    ps_a, info_a = curve_to_ps(import_curve(p_points, "points"), DT, UMAX)
    ps_b, info_b = curve_to_ps(import_curve(p_table, "s_kappa_tau"), DT, UMAX)
    assert abs(ps_a.N - ps_b.N) <= 1
    n = min(ps_a.N, ps_b.N)
    assert np.allclose(ps_a.amplitude[5:n - 5], ps_b.amplitude[5:n - 5], atol=0.5)
    # helix -> constant amplitude and constant phase rate
    assert np.std(ps_a.amplitude[10:-10]) < 0.2


def test_kappa_tau_equal_spacing():
    path = write("\n".join("0.5, 0.1" for _ in range(500)))
    pts = import_curve(path, "kappa_tau")
    assert pts.shape == (500, 3)


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
