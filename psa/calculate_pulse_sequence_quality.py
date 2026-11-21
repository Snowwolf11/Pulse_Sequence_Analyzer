import numpy as np

from psa.rotations import propagate_batch, inverse_rotate, quaternion_to_matrix


def rotation_matrix(axis, angle_deg):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    th = np.deg2rad(angle_deg)
    c, s, C = np.cos(th), np.sin(th), 1 - np.cos(th)
    return np.array([
        [c + x*x*C,     x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s,   c + y*y*C,   y*z*C - x*s],
        [z*x*C - y*s,   z*y*C + x*s, c + z*z*C]
    ])


def calculate_pulse_sequence_quality_batch(PS, T, l, Umax, off, initialVector=np.array([0, 0, 1]),
                                           shallVec=np.array([0, 0, 1]), calcType=1, target_Rotation=np.eye(3)):
    """Quality values for arrays of amplitudes ``Umax`` and offsets ``off`` (broadcast together).

    calcType 1: (1 + t_final . shallVec) / 2 with t_final = R^T (l * initialVector)
    calcType 2: (1 + tr(target^T R^T) / 3) / 2
    calcType 3: 2 * |curve end point| (closedness, unit step length l)
    """
    initialVector = np.asarray(initialVector, dtype=np.float64)
    if calcType == 1:
        q = propagate_batch(PS, T, Umax, off)
        Q = inverse_rotate(q, l * initialVector) @ np.asarray(shallVec, dtype=np.float64)
        return (1 + Q) / 2
    if calcType == 2:
        q = propagate_batch(PS, T, Umax, off)
        Rtot = np.swapaxes(quaternion_to_matrix(q), -1, -2)
        Q = np.einsum("ij,...ij->...", np.asarray(target_Rotation, dtype=np.float64), Rtot) / 3.0
        return (1 + Q) / 2
    if calcType == 3:
        _, end_point = propagate_batch(PS, T, Umax, off, initialVector=initialVector, l=l)
        return 2 * np.linalg.norm(end_point, axis=-1)
    raise ValueError("calcType must be 1, 2 or 3")


def calculate_pulse_sequence_quality(PS, T, l, Umax, off, inpoFact=1, initialVector=np.array([0, 0, 1]),
                                     shallVec=np.array([0, 0, 1]), calcType=1, target_Rotation=np.eye(3)):
    """Quality value for a single amplitude/offset point (inpoFact is fixed to 1 as before)."""
    return float(calculate_pulse_sequence_quality_batch(PS, T, l, Umax, off, initialVector, shallVec, calcType, target_Rotation))
