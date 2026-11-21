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
        assert np.allclose(tsum[i], t.sum(axis=0), atol=1e-10)
        assert np.allclose(quaternion_to_matrix(q[i]).T, inverse_total_rotation(PS_UR[:200], T, amps[i], offs[i]), atol=1e-12)


def test_zero_amplitude_and_offset_is_identity():
    PS = np.array([[0.0, 0.0], [0.0, 90.0]])
    assert np.allclose(inverse_total_rotation(PS, T, UMAX, 0.0), np.eye(3))
    assert np.allclose(tangent_vectors(PS, T, 2.0, UMAX, 0.0, 1, [0, 0, 1]), [[0, 0, 2]] * 3)


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for t in tests:
        t(); print(f"[PASS] {t.__name__}")
    print(f"{len(tests)}/{len(tests)} passed")
