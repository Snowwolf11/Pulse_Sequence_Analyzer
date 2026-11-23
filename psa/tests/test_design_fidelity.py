"""Tests of the fidelity maps (psa.design.fidelity)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design import metrics
from psa.design.fidelity import infidelity_grid, run_job, summary
from psa.design.pulse import PulseSequence
from psa.rotations import propagate_batch, quaternion_to_matrix

DT, UMAX = 0.5e-6, 1e4


def pulse(N=120, seed=0):
    rng = np.random.default_rng(seed)
    return PulseSequence(rng.uniform(30, 100, N), rng.uniform(-180, 180, N), DT, UMAX)


def test_own_reference_matches_the_panel():
    ps = pulse()
    offs = np.linspace(-8e3, 8e3, 33)
    inf = infidelity_grid(ps, offs, [0.9, 1.0, 1.1])
    assert inf.shape == (3, 33)
    assert np.allclose(inf[1], 1 - metrics.fidelity_vs_offset(ps, offs), atol=1e-13)
    assert infidelity_grid(ps, [0.0], [1.0])[0, 0] < 1e-28                   # exact at the reference point


def test_ur_and_pp_references():
    ps = pulse(seed=1)
    g = metrics.gate(ps)
    offs, scales = np.linspace(-5e3, 5e3, 11), [0.95, 1.05]
    ur = infidelity_grid(ps, offs, scales, {"kind": "UR", "axis": g["axis"], "angle": g["angle_deg"]})
    assert np.allclose(ur, infidelity_grid(ps, offs, scales), atol=1e-12)      # own gate = this rotation
    m0, mt = np.array([0, 0, 1.0]), np.array([0, 1.0, 0])
    pp = infidelity_grid(ps, offs, scales, {"kind": "PP", "m0": m0, "mt": mt})
    for i, s in enumerate(scales):
        R = quaternion_to_matrix(propagate_batch(ps.as_array(), DT, UMAX * s, offs))
        assert np.allclose(pp[i], 1 - (1 + (R @ m0) @ mt) / 2, atol=1e-13)
    # amplitude scaling s = the pulse with s x the amplitude
    scaled = PulseSequence(ps.amplitude * 0.95, ps.phase, DT, UMAX)
    assert np.allclose(pp[0], infidelity_grid(scaled, offs, [1.0], {"kind": "PP", "m0": m0, "mt": mt})[0], atol=1e-13)


def test_run_job_and_summary():
    ps = pulse(seed=2)
    job = {"id": 7, "pulse": [ps.amplitude.tolist(), ps.phase.tolist(), DT, UMAX], "reference": {"kind": "own"},
           "map": {"offsets_hz": [-1e3, 0.0, 1e3], "scales": [0.9, 1.0]},
           "cut_offset": {"offsets_hz": [-2e3, 2e3], "scale": 1.0}, "cut_amp": {"scales": [0.8, 1.0, 1.2], "offset_hz": 0.0}}
    r = run_job(job)
    assert r["id"] == 7 and r["map"].shape == (2, 3) and r["cut_offset"].shape == (2,) and r["cut_amp"].shape == (3,)
    assert r["map"][1, 1] == r["cut_amp"][1] and r["cut_amp"][1] < 1e-28         # 0 kHz, 100 %: the reference
    bad = run_job(dict(job, reference={"kind": "UR", "axis": [0, 0, 0], "angle": 90}))
    assert bad["id"] == 7 and "error" in bad
    s = summary([1e-6, 1e-4, 1e-2, 0.5])
    assert s["worst"] == 0.5 and s["below 0.0001"] == 0.5 and s["below 1e-05"] == 0.25


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
