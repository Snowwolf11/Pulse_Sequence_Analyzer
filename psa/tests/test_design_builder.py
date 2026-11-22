"""Tests of the segment-based curve builder (psa.design.builder)."""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design import metrics
from psa.design.builder import CurveDesign, Segment
from psa.design.conversion import ps_to_curve

UMAX, DT = 1e4, 0.5e-6
RMIN = 1e6 / (2 * math.pi * UMAX)        # 15.915 µs


def design(*segments, dt=DT):
    return CurveDesign(list(segments), UMAX, dt)


def test_full_circle_closes_and_uses_full_amplitude():
    d = design(Segment("arc", {"radius": RMIN, "angle": 360}))
    pts, _ = d.curve()
    assert np.linalg.norm(pts[-1]) < 1e-12
    assert d.robustness()["A1_norm"] < 1e-12
    ps = d.to_pulse()
    assert ps.N == 200 and np.allclose(ps.amplitude, 100.0, atol=0.01)   # 2 pi R_min = 100 µs
    assert metrics.robustness(ps)["A1_norm"] < 1e-3


def test_arc_geometry():
    for R, ang in ((20.0, 90.0), (50.0, 250.0)):
        d = design(Segment("arc", {"radius": R, "angle": ang}))
        pts, _ = d.curve()
        chord = 2 * R * math.sin(math.radians(ang) / 2) * 1e-6
        assert abs(np.linalg.norm(pts[-1]) - chord) < 1e-4 * chord        # drawn with 0.25 µs chords
        assert abs(d.duration() - R * math.radians(ang) * 1e-6) < 1e-15
        ps = d.to_pulse()
        assert np.allclose(ps.amplitude, RMIN / R * 100, rtol=1e-2)


def test_helix_curvature_torsion_and_pitch():
    for pitch in (40.0, -25.0):
        d = design(Segment("helix", {"radius": 20.0, "pitch": pitch, "angle": 360.0}))
        pts, _ = d.curve(points_per_piece_us=0.05)
        assert abs(np.linalg.norm(pts[-1]) - abs(pitch) * 1e-6) < 1e-3 * abs(pitch) * 1e-6     # one turn: displacement = pitch
        d1 = np.gradient(pts, axis=0); d2 = np.gradient(d1, axis=0); d3 = np.gradient(d2, axis=0)
        cr = np.cross(d1, d2)
        kappa = np.median(np.linalg.norm(cr, axis=1) / np.linalg.norm(d1, axis=1) ** 3)
        tau = np.median(np.einsum("ij,ij->i", cr, d3) / np.einsum("ij,ij->i", cr, cr))
        c = pitch / (2 * math.pi)
        assert abs(kappa - 20.0 / (400 + c * c) * 1e6) < 1e-3 * kappa
        assert abs(tau - c / (400 + c * c) * 1e6) < 2e-3 * abs(tau)


def test_twist_makes_curve_nonplanar():
    def planarity(twist):
        d = design(Segment("arc", {"radius": 20, "angle": 120}), Segment("arc", {"radius": 25, "angle": 120}, twist=twist))
        pts, _ = d.curve()
        centered = pts - pts.mean(axis=0)
        return np.linalg.svd(centered, compute_uv=False)[-1] / np.linalg.svd(centered, compute_uv=False)[0]
    assert planarity(0.0) < 1e-9 and planarity(180.0) < 1e-9          # 180° twist = inflection, still planar
    assert planarity(90.0) > 1e-2


def test_line_segment_is_a_pause():
    d = design(Segment("arc", {"radius": 30, "angle": 90}), Segment("line", {"length": 10}), Segment("arc", {"radius": 30, "angle": 90}))
    ps = d.to_pulse()
    n_arc = int(round(30 * math.pi / 2 / 0.5))
    assert np.allclose(ps.amplitude[n_arc:n_arc + 20], 0.0)
    _, starts = d.curve()
    assert [s for s, _ in starts] == [0, 1, 2]


def test_continuous_model_matches_discretized_pulse():
    d = design(Segment("arc", {"radius": 20, "angle": 140}),
               Segment("helix", {"radius": 18, "pitch": 30, "angle": 300}, twist=40),
               Segment("line", {"length": 7}),
               Segment("arc", {"radius": 35, "angle": 200}, twist=-75), dt=0.01e-6)
    rb = d.robustness()
    ps = d.to_pulse()
    rb_ps = metrics.robustness(ps)
    for k in (1, 2, 3):
        a, b = rb[f"A{k}"], rb_ps[f"A{k}"]
        assert np.linalg.norm(a - b) < 2e-3 * np.linalg.norm(a), (k, a, b)
    R = d.gate_matrix()
    assert np.allclose(R, metrics.gate(ps)["matrix"], atol=2e-3)
    pts, _ = d.curve(points_per_piece_us=0.01)
    assert np.linalg.norm(pts[-1] - ps_to_curve(ps)[-1]) < 1e-3 * d.duration()


def test_bezier_segment():
    # Bezier with all control points on the tangent: a straight piece, zero amplitude
    d = design(Segment("bezier", {"p1": 10, "p2x": 0, "p2y": 0, "p2z": 20, "p3x": 0, "p3y": 0, "p3z": 30}))
    assert np.max(d.to_pulse().amplitude) < 1e-6
    # curved Bezier after an arc: tangent continuous at the junction, curve matches the analytic Bezier
    arc = Segment("arc", {"radius": 25, "angle": 70})
    bez = Segment("bezier", {"p1": 12, "p2x": 10, "p2y": 4, "p2z": 25, "p3x": 25, "p3y": -6, "p3z": 35}, twist=30)
    d = design(arc, bez, dt=0.05e-6)
    pts, starts = d.curve(points_per_piece_us=0.05)
    j = dict(starts)[1]
    t_before = pts[j] - pts[j - 1]; t_after = pts[j + 1] - pts[j]
    cosang = t_before @ t_after / np.linalg.norm(t_before) / np.linalg.norm(t_after)
    assert cosang > 1 - 1e-4, cosang
    pos, axes = d.frame_at_segment_starts()[1]
    P = np.array([[0, 0, 0], [0, 0, 12], [10, 4, 25], [25, -6, 35]]) * 1e-6 @ axes + pos
    assert np.linalg.norm(pts[-1] - P[-1]) < 2e-3 * d.duration(), np.linalg.norm(pts[-1] - P[-1]) / d.duration()


def test_solver_closes_curve_first_and_second_order():
    d = design(Segment("arc", {"radius": 20, "angle": 90}),
               Segment("arc", {"radius": 20, "angle": 200}, twist=150, free=["angle", "twist"]),
               Segment("arc", {"radius": 20, "angle": 120}, twist=150, free=["angle", "twist"]))
    assert d.robustness()["A1_norm"] > 1e-2
    solved, res = d.solve(orders=(1,))
    assert res["reached"] and res["norms"][1] < 1e-8, res          # stops once below the "solved" threshold
    d2 = design(Segment("arc", {"radius": 20, "angle": 180}),
                Segment("arc", {"radius": 20, "angle": 180}, twist=170, free=["angle", "twist"]),
                Segment("arc", {"radius": 20, "angle": 180}, twist=100, free=["angle", "twist"]),
                Segment("arc", {"radius": 20, "angle": 180}, twist=170, free=["angle", "twist"]),
                Segment("arc", {"radius": 20, "angle": 180}, twist=100, free=["angle", "twist"]))
    solved2, res2 = d2.solve(orders=(1, 2))
    assert res2["norms"][1] < 1e-8 and res2["norms"][2] < 1e-8, res2
    # the solved design, discretized at a fine dt, is robust too (up to discretization)
    ps = solved2.to_pulse(0.01e-6)
    assert metrics.robustness(ps)["A1_norm"] < 1e-3


def test_solver_can_be_stopped_and_returns_best_so_far():
    d = design(Segment("arc", {"radius": 20, "angle": 90}),
               Segment("arc", {"radius": 20, "angle": 200}, twist=150, free=["angle", "twist"]),
               Segment("arc", {"radius": 20, "angle": 120}, twist=150, free=["angle", "twist"]))
    start = d.robustness()["A1_norm"]
    calls = []
    solved, res = d.solve(orders=(1, 2, 3), callback=lambda n, c: calls.append(c) or n >= 5)
    assert res["cancelled"] and res["nfev"] == 5 and len(calls) == 5
    assert res["cost"] <= calls[0] + 1e-15                          # best of the evaluations so far
    assert solved.robustness()["A1_norm"] <= start + 1e-12
    assert d.segments[1].params["angle"] == 200                     # the original design is untouched


def test_serialization():
    d = design(Segment("arc", {"radius": 22}, twist=10, free=["angle"]), Segment("bezier", {"p3y": 3.0}))
    d2 = CurveDesign.from_dict(json.loads(json.dumps(d.to_dict())))
    assert np.allclose(d2.curve()[0], d.curve()[0]) and d2.segments[0].free == ["angle"]


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
