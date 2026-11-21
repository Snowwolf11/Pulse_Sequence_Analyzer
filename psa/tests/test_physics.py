"""Independent physics check: the rotation backend against exact SU(2) propagation.

Run with ``python -m pytest psa/tests`` or ``python psa/tests/test_physics.py``.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.getPulseSequence import getPulseSequence
from psa.rotations import tangent_vectors, inverse_total_rotation, propagate_batch, inverse_rotate, quaternion_to_matrix

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1, -1]).astype(complex)
PS_UR = getPulseSequence(os.path.join(HERE, "data", "PS_ur.bruker"))
T, UMAX = 0.5e-6, 1e4


def su2_reference(PS, T, Umax, offset):
    """Tangents U_k^dag sz U_k (k = 0..N) and the Bloch rotation matrix R of U_N."""
    U = np.eye(2, dtype=complex)
    def tangent(U):
        M = U.conj().T @ SZ @ U
        return np.real([np.trace(M @ s) / 2 for s in (SX, SY, SZ)])
    tangents = [tangent(U)]
    for b1, phase in PS:
        w = 2 * np.pi * np.array([b1 / 100 * Umax * np.cos(np.radians(phase)), b1 / 100 * Umax * np.sin(np.radians(phase)), offset])
        a = np.linalg.norm(w)
        n = w / a if a else np.array([0, 0, 1.0])
        U = (np.cos(a * T / 2) * np.eye(2) - 1j * np.sin(a * T / 2) * (n[0] * SX + n[1] * SY + n[2] * SZ)) @ U
        tangents.append(tangent(U))
    R = np.array([[np.real(np.trace(si @ U @ sj @ U.conj().T)) / 2 for sj in (SX, SY, SZ)] for si in (SX, SY, SZ)])
    return np.array(tangents), R


def test_tangents_match_su2():
    for offset in (0.0, 3000.0, -7000.0):
        ref, _ = su2_reference(PS_UR, T, UMAX, offset)
        assert np.max(np.abs(tangent_vectors(PS_UR, T, 1.0, UMAX, offset, 1, [0, 0, 1]) - ref)) < 1e-10


def test_inverse_total_rotation_is_RT():
    for offset in (0.0, 3000.0):
        _, R = su2_reference(PS_UR, T, UMAX, offset)
        assert np.max(np.abs(inverse_total_rotation(PS_UR, T, UMAX, offset) - R.T)) < 1e-10


def test_batch_matches_single():
    amps = np.array([9e3, 1e4, 1.1e4]); offs = np.array([-2e3, 0.0, 5e3])
    q, tsum = propagate_batch(PS_UR[:200], T, amps, offs, initialVector=[0, 0, 1])
    for i in range(3):
        t = tangent_vectors(PS_UR[:200], T, 1.0, amps[i], offs[i], 1, [0, 0, 1])
        assert np.allclose(inverse_rotate(q[i], [0, 0, 1]), t[-1], atol=1e-12)
        assert np.allclose(tsum[i], t.sum(axis=0) - (t[0] + t[-1]) / 2, atol=1e-10)
        assert np.allclose(quaternion_to_matrix(q[i]).T, inverse_total_rotation(PS_UR[:200], T, amps[i], offs[i]), atol=1e-12)


def test_pp_fidelity_is_physical_transfer():
    """PP fidelity = (1 + R(delta) m0 . R(0) m0) / 2 for the pulse as given (BUG_TRACKER PSA-01)."""
    from psa.calculate_pulse_sequence_quality import calculate_pulse_sequence_quality
    PS = getPulseSequence(os.path.join(HERE, "data", "quality_dir", "pulse_zz_ur.bruker"))[:400]
    for m0 in (np.array([0, 0, 1.0]), np.array([1.0, 0, 0]), np.array([0.3, -0.5, 0.8]) / np.linalg.norm([0.3, -0.5, 0.8])):
        _, R0 = su2_reference(PS, T, UMAX, 0.0)
        target = R0 @ m0
        for offset in (2000.0, -4500.0):
            _, Rd = su2_reference(PS, T, UMAX, offset)
            expected = (1 + (Rd @ m0) @ target) / 2
            got = calculate_pulse_sequence_quality(PS, T, 1.0, UMAX, offset, 1, [0, 0, 1], target, 1, initialMagnetization=m0)
            assert abs(got - expected) < 1e-10, (m0, offset, got, expected)


def test_angle_mode_is_flip_angle_of_magnetization():
    """Stability mode 3 = polar angle of R m0 (BUG_TRACKER PSA-05)."""
    from psa.calculate_curve_stability import calculate_curve_stability
    PS = getPulseSequence(os.path.join(HERE, "data", "quality_dir", "pulse035.bruker"))
    X, Y, Z, *_ = calculate_curve_stability(PS, T, 1.0, UMAX, 20, 20, 3, [0, 0, 1], [0, 0, 1])
    for j, i in ((0, 0), (77, 150), (199, 199)):
        _, R = su2_reference(PS, T, X[j, i], Y[j, i])
        m = R @ [0, 0, 1]
        assert abs(Z[j, i] - np.degrees(np.arctan2(np.hypot(m[0], m[1]), m[2]))) < 1e-8


def exact_endpoint(PS, offset, sub=40):
    """r(T) from finely subdivided exact propagation (trapezoid on a 40x finer grid), units of T."""
    ts, _ = su2_reference(np.repeat(PS, sub, axis=0), T / sub, UMAX, offset)
    return (ts[:-1] + ts[1:]).sum(axis=0) / 2 / sub


def test_curve_endpoint_trapezoid():
    """Rotation-matrix curve: M + 1 points, end point = r(T) (BUG_TRACKER PSA-04)."""
    from psa.createCurve import createCurve
    for offset, inpo, tol in ((3000.0, 1, 2e-3), (3000.0, 5, 1e-4), (-6000.0, 5, 1e-4)):
        CM, VM, PSm = createCurve(PS_UR.copy(), T=T, l=1, maximumAmplitude=UMAX, offset=offset, inpoFact=inpo,
                                  calculationMethod=1, initialVector=np.array([0, 0, 1]))
        assert len(CM) == len(PS_UR) * inpo + 1 and len(VM) == len(CM)
        assert np.linalg.norm(CM[-1] / T - exact_endpoint(PS_UR, offset)) < tol


def test_closed_curve_is_closed():
    from psa.createCurve import createCurve
    from psa.calculate_pulse_sequence_quality import calculate_pulse_sequence_quality
    PS = getPulseSequence(os.path.join(HERE, "data", "PS_4H_closed.bruker"))
    assert np.linalg.norm(exact_endpoint(PS, 0.0)) < 1e-6
    for inpo in (1, 5):
        CM, _, _ = createCurve(PS.copy(), T=T, l=1, maximumAmplitude=UMAX, offset=0.0, inpoFact=inpo,
                               calculationMethod=1, initialVector=np.array([0, 0, 1]))
        assert np.linalg.norm(CM[-1]) / T < 1e-9
    assert calculate_pulse_sequence_quality(PS, T, 1.0, UMAX, 0.0, 1, [0, 0, 1], [0, 0, 1], 3) < 1e-9


def test_orientation_quantities_use_true_tangents():
    """Rotation angle / angle to x-axis use t_0 = initial vector and t_M = R^T initial vector."""
    import psa.app_backend as ab
    for method in (1, 2):
        r = ab.analyze_pulse_sequence(ab.AnalyzerParams(os.path.join(HERE, "data", "PS_ur.bruker"), 10, 1, 0.5, 5, 0, 2, 100,
                                                        method, np.array([0, 0, 1])))
        _, R = su2_reference(PS_UR, T, UMAX, 2000.0)
        tM = R.T @ [0, 0, 1]
        assert abs(r.totalRot - np.degrees(np.arccos(np.clip(tM[2], -1, 1)))) < 1e-8
        assert abs(r.phi - np.degrees(np.arctan2(tM[1], tM[0]))) < 1e-8
        assert np.allclose(r.VM[0] / np.linalg.norm(r.VM[0]), [0, 0, 1])


def test_zero_amplitude_and_offset_is_identity():
    PS = np.array([[0.0, 0.0], [0.0, 90.0]])
    assert np.allclose(inverse_total_rotation(PS, T, UMAX, 0.0), np.eye(3))
    assert np.allclose(tangent_vectors(PS, T, 2.0, UMAX, 0.0, 1, [0, 0, 1]), [[0, 0, 2]] * 3)


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for t in tests:
        t(); print(f"[PASS] {t.__name__}")
    print(f"{len(tests)}/{len(tests)} passed")
