"""Tests of the fidelity maps (psa.design.fidelity)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design import metrics
from psa.design.fidelity import infidelity_grid, phase_profile, run_job, summary
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


def test_freedom_fixed_and_fitted():
    """A pulse followed by a free delay tau0 is, up to a linear phase with delay tau0, the pulse itself."""
    a = pulse(80, seed=3)
    n0 = 200                                                    # 100 µs of free evolution after the pulse
    b = PulseSequence(np.r_[a.amplitude, np.zeros(n0)], np.r_[a.phase, np.zeros(n0)], DT, UMAX)
    offs, scales = np.linspace(-3e3, 3e3, 25), [0.95, 1.0, 1.05]
    ref = {"kind": "UR", **{k: v for k, v in zip(("axis", "angle"), (metrics.gate(a)["axis"], metrics.gate(a)["angle_deg"]))}}
    base = infidelity_grid(a, offs, scales, ref)
    job = {"id": 1, "pulse": [b.amplitude.tolist(), b.phase.tolist(), DT, UMAX],
           "map": {"offsets_hz": offs.tolist(), "scales": scales}, "cut_offset": None, "cut_amp": None}
    fixed = [run_job(dict(job, reference=dict(ref, freedom="linear", fit=False, phi_deg=0.0, tau_us=sgn * n0 * DT * 1e6)))
             for sgn in (1, -1)]
    err = [np.max(np.abs(r["map"] - base)) for r in fixed]
    assert min(err) < 1e-12, err                                # one sign = the convention of the free evolution
    assert fixed[0]["correction"]["source"] == "fixed"
    none = run_job(dict(job, reference=ref))
    assert none["correction"] is None and np.mean(none["map"]) > 2 * np.mean(base)
    fit = run_job(dict(job, reference=dict(ref, freedom="linear", fit=True)))
    c = fit["correction"]
    assert c["source"] == "fitted on the map area" and np.mean(fit["map"]) <= np.mean(base) + 1e-12
    for dphi, dtau in ((1.0, 0), (-1.0, 0), (0, 0.5), (0, -0.5)):    # the fit is the best of the neighbours
        r = run_job(dict(job, reference=dict(ref, freedom="linear", fit=False, phi_deg=c["phi_deg"] + dphi,
                                             tau_us=c["tau_us"] + dtau)))
        assert np.mean(r["map"]) >= np.mean(fit["map"]) - 1e-15
    ph = run_job(dict(job, reference=dict(ref, freedom="phase", fit=True)))
    assert ph["correction"]["tau_us"] == 0 and np.mean(ph["map"]) <= np.mean(none["map"]) + 1e-15
    # cuts use the same correction; without a map the fit needs an area
    cut = run_job(dict(job, map=None, fit_area=job["map"], cut_offset={"offsets_hz": offs.tolist(), "scale": 1.0},
                       reference=dict(ref, freedom="linear", fit=True)))
    assert np.allclose(cut["cut_offset"], fit["map"][1]) and np.isclose(cut["correction"]["tau_us"], c["tau_us"])
    assert "error" in run_job(dict(job, map=None, cut_offset=cut and {"offsets_hz": [0.0], "scale": 1.0},
                                   reference=dict(ref, freedom="phase", fit=True)))
    # own reference and PP with freedom
    for r in ({"kind": "own"}, {"kind": "PP", "m0": [0, 0, 1], "mt": [1, 0, 0]}):
        out = run_job(dict(job, reference=dict(r, freedom="linear", fit=True)))
        assert "error" not in out and np.mean(out["map"]) <= np.mean(run_job(dict(job, reference=r))["map"]) + 1e-15


def test_phase_profile():
    """A delay tau0 after the pulse adds a linear phase: the fitted delay grows by tau0 (in the convention of
    the "linear phase" freedom), the residual stays the same; the fit removes it."""
    a = pulse(60, seed=5)
    n0 = 60                                                     # 30 µs free evolution
    b = PulseSequence(np.r_[a.amplitude, np.zeros(n0)], np.r_[a.phase, np.zeros(n0)], DT, UMAX)
    offs = np.linspace(-4e3, 4e3, 161)
    pa, pb = phase_profile(a, offs, max_delay_us=200), phase_profile(b, offs, max_delay_us=200)
    assert np.allclose(pa["mxy"], pb["mxy"]) and np.isclose(abs(pb["tau_us"] - pa["tau_us"]), 30.0, atol=1e-3)
    assert np.allclose(pa["residual_deg"], pb["residual_deg"], atol=1e-6)
    # the line is the "linear phase" freedom: with it fixed, PP z -> transverse direction at phi is exact
    m0 = np.array([0, 0, 1.0])
    ref = {"kind": "PP", "m0": m0, "mt": [1, 0, 0], "freedom": "linear", "fit": False,
           "phi_deg": pb["phi_deg"], "tau_us": pb["tau_us"]}
    job = {"id": 1, "pulse": [b.amplitude.tolist(), b.phase.tolist(), DT, UMAX], "reference": ref,
           "map": None, "cut_offset": {"offsets_hz": offs.tolist(), "scale": 1.0}, "cut_amp": None,
           "phase": {"offsets_hz": offs.tolist(), "scale": 1.0, "m0": [0, 0, 1], "max_delay_us": 200}}
    r = run_job(job)
    keep_z = infidelity_grid(b, offs, [1.0], {"kind": "PP", "m0": m0, "mt": [0, 0, 1]})[0]   # = (1 - mz) / 2
    mxy, res = pb["mxy"], np.radians(pb["residual_deg"])
    mz = 1 - 2 * keep_z
    expected = ((mxy * np.cos(res) - 1) ** 2 + (mxy * np.sin(res)) ** 2 + mz ** 2) / 4
    assert np.allclose(r["cut_offset"], expected, atol=1e-12)
    assert np.allclose(r["phase"]["tau_us"], pb["tau_us"]) and np.allclose(r["phase"]["fit_deg"], pb["fit_deg"])
    # the shown phase is continuous next to the line (no 360° jumps)
    assert np.all(np.abs(np.diff(pb["phase_deg"])) < 180)


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
