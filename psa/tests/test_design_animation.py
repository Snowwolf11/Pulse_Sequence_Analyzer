"""Tests of the animation dynamics (psa.design.animation)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.animation import (Dynamics, droplet, magnetization_droplet, phase_colors, propagator_droplets,
                                  rotate, unit_sphere)
from psa.design.conversion import ps_to_curve
from psa.design.pulse import PulseSequence
from psa.rotations import propagate_batch, rotate as rot_ref

DT, UMAX = 0.5e-6, 1e4


def pulse(N=80, seed=0):
    rng = np.random.default_rng(seed)
    return PulseSequence(rng.uniform(20, 100, N), rng.uniform(-180, 180, N), DT, UMAX)


def test_dynamics_match_the_propagator():
    ps = pulse()
    offs, scs = np.array([-3e3, 0.0, 2e3]), np.array([0.9, 1.1])
    d = Dynamics(ps, offs, scs, m0=(1, 0, 0), selected=(1e3, 1.05))
    assert d.n_ensemble == 6 and d.q.shape == (7, 81, 4) and d.m.shape == (7, 81, 3)
    for i in range(7):
        qf = propagate_batch(ps.as_array(), DT, UMAX * d.scales[i], d.offsets[i])
        assert np.allclose(d.q[i, -1], qf, atol=1e-12)
        assert np.allclose(d.m[i, -1], rot_ref(qf, [1, 0, 0]), atol=1e-12)
    # every boundary: the propagator of the first k sub-pulses
    k = 37
    part = PulseSequence(ps.amplitude[:k], ps.phase[:k], DT, UMAX)
    assert np.allclose(d.q[2, k], propagate_batch(part.as_array(), DT, UMAX * d.scales[2], d.offsets[2]), atol=1e-12)


def test_fractional_times_are_exact():
    ps = pulse(10, seed=1)
    d = Dynamics(ps, [0.0, 4e3], [1.0])
    assert np.allclose(d.quaternions(3.0), d.q[:, 3]) and np.allclose(d.quaternions(10), d.q[:, 10])
    # half of sub-pulse 3 = a pulse whose 4th sub-pulse lasts dt / 2
    half = d.quaternions(3.5)
    two = PulseSequence(np.r_[ps.amplitude[:3]], np.r_[ps.phase[:3]], DT, UMAX)
    q3 = propagate_batch(two.as_array(), DT, UMAX, 4e3)
    qh = propagate_batch(np.array([[ps.amplitude[3], ps.phase[3]]]), DT / 2, UMAX, 4e3)
    from psa.design.animation import _qmul
    assert np.allclose(half[1], _qmul(qh, q3), atol=1e-12)
    assert np.allclose(np.linalg.norm(d.magnetization(2.3), axis=1), 1.0)
    tr = d.trail(5.4, 3)
    assert tr.shape == (3, 5, 3) and np.allclose(tr[:, -1], d.magnetization(5.4))


def test_error_curve_and_field():
    ps = pulse(60, seed=2)
    d = Dynamics(ps, [0.0], [1.0])                                   # selected spin on resonance, nominal
    ref = ps_to_curve(ps) / ps.duration
    assert np.allclose(d.curve[-1] - d.curve[0], ref[-1] - ref[0], atol=1e-9)
    assert np.isclose(np.sum(np.linalg.norm(np.diff(d.curve, axis=0), axis=1)), 1.0, atol=1e-3)   # length = T
    assert np.allclose(d.curve_until(60), d.curve) and len(d.curve_until(20.5)) == 22
    f = d.effective_field(7.2)
    assert np.allclose(f[:2], ps.amplitude[7] * UMAX / 100 * np.array([np.cos(np.radians(ps.phase[7])),
                                                                    np.sin(np.radians(ps.phase[7]))]))


def test_droplets():
    dirs, faces = unit_sphere(8, 16)
    assert dirs.shape == (9 * 16, 3) and faces.max() < len(dirs)
    f0, f1 = propagator_droplets([1.0, 0, 0, 0], dirs)               # identity: sphere, no lobes
    assert np.allclose(np.abs(f0), 1) and np.allclose(f1, 0)
    f0, f1 = propagator_droplets([0.0, 1.0, 0, 0], dirs)             # 180° about x: lobes along x
    assert np.allclose(f0, 0) and np.isclose(np.max(np.abs(f1)), 1, atol=0.01)
    v = droplet(magnetization_droplet([0, 0, 1], dirs), dirs)
    assert np.allclose(np.linalg.norm(v, axis=1), np.abs(dirs[:, 2]))
    c = phase_colors([1, -1, 1j])
    assert c[0, 0] > c[0, 1] and c[1, 1] > c[1, 0] and c.shape == (3, 4)              # red, green
    assert np.allclose(rotate([1.0, 0, 0, 0], [0, 0, 1]), [0, 0, 1])


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
