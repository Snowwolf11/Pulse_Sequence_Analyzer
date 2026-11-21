import time

import numpy as np

from psa.rotations import inverse_total_rotation, propagate_batch, rotate
from psa.calculate_pulse_sequence_quality import calculate_pulse_sequence_quality_batch


def calculate_curve_stability(PS, T, l, Umax, scalingRange_percent=20, offsetRange_kHz=20, stabilityCalculationMethod=1, initialVector=np.array([0,0,1]), initialMagnetization=np.array([0,0,1])):
    """Quality over a 200 x 200 grid of amplitude scaling x offset.

    Z[j, i] belongs to offset j and amplitude scaling i (as X, Y from np.meshgrid).
    Methods: 1 = PP/SS (final magnetization R m0 vs. target = on-resonance R m0, snapped to +-z
    if nearly along z), 2 = UR (total rotation vs. closest canonical rotation), 3 = polar angle
    (flip angle from +z) of the final magnetization R m0 in degrees.
    m0 = ``initialMagnetization`` (normalized); ``initialVector`` (noise axis) is not used here.
    """
    a = time.time()
    PS = np.asarray(PS, dtype=np.float64)
    initialVector = np.asarray(initialVector, dtype=np.float64)
    m0 = np.asarray(initialMagnetization, dtype=np.float64)
    m0 = m0 / np.linalg.norm(m0)

    offsetRange = 1e3 * np.linspace(-offsetRange_kHz/2, offsetRange_kHz/2, 200)

    scalingRange = np.linspace(1-scalingRange_percent/200,  1+scalingRange_percent/200, 200)
    scalingRange *= Umax

    X, Y = np.meshgrid(scalingRange, offsetRange)

    medianVec = rotate(propagate_batch(PS, T, Umax, 0.0), m0)
    Rtot_on_res = inverse_total_rotation(PS, T, Umax, 0.0, 1)

    if np.linalg.norm(medianVec[2]) > 5 * (np.linalg.norm(medianVec[1]) + np.linalg.norm(medianVec[0])):
        if medianVec[2] > 0:
            shallVec = np.array([0, 0, 1])
        elif medianVec[2] < 0:
            shallVec = np.array([0, 0, -1])
    else:
        shallVec = medianVec

    axis, angle, err_deg, R_target = closest_canonical_rotation(Rtot_on_res, return_matrix=True)

    Z = np.ones_like(X)
    if stabilityCalculationMethod == 1:
        Z = calculate_pulse_sequence_quality_batch(PS, T, l, X, Y, initialVector, shallVec, 1, initialMagnetization=m0)
    elif stabilityCalculationMethod == 2:
        Z = calculate_pulse_sequence_quality_batch(PS, T, l, X, Y, initialVector, shallVec, 2, target_Rotation=R_target)
    elif stabilityCalculationMethod == 3:
        M_end = rotate(propagate_batch(PS, T, X, Y), m0)
        Z = np.arctan2(np.sqrt(M_end[..., 0]**2 + M_end[..., 1]**2), M_end[..., 2]) * 180 / np.pi

    quality = float(np.mean(Z))

    b = time.time()
    print("Time: " + str(b-a))
    return X, Y, Z, quality, axis, angle

def _rotmat_from_axis_angle(axis, angle_deg):
    """Rodrigues' formula."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    th = np.deg2rad(angle_deg)
    c, s, C = np.cos(th), np.sin(th), 1 - np.cos(th)
    return np.array([
        [c + x*x*C,     x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s,   c + y*y*C,   y*z*C - x*s],
        [z*x*C - y*s,   z*y*C + x*s, c + z*z*C]
    ])

def _project_to_SO3(R):
    """Nearest proper rotation via SVD (handles small numerical drift)."""
    U, _, Vt = np.linalg.svd(R)
    Rproj = U @ Vt
    if np.linalg.det(Rproj) < 0:
        U[:, -1] *= -1
        Rproj = U @ Vt
    return Rproj

def _geodesic_angle(R):
    """Return rotation angle of R (in radians)."""
    c = (np.trace(R) - 1.0) / 2.0
    c = np.clip(c, -1.0, 1.0)
    return np.arccos(c)

def closest_canonical_rotation(R, project=True, return_matrix=False):
    """
    Find the closest rotation to R from:
      axes = {+/-x, +/-y, +/-z}, angles = {0°, 90°, 180°}.
    Returns: (axis_name, angle_deg, error_deg[, R_best])

    Notes:
      - For 0°, axis is irrelevant; returns axis_name='any'.
      - A 270° about +x, for example, maps to 90° about -x.
    """
    _AXES = {
    "+x": np.array([ 1., 0., 0.]),
    "-x": np.array([-1., 0., 0.]),
    "+y": np.array([ 0., 1., 0.]),
    "-y": np.array([ 0.,-1., 0.]),
    "+z": np.array([ 0., 0., 1.]),
    "-z": np.array([ 0., 0.,-1.]),
    }
    R = np.asarray(R, dtype=float)
    if R.shape != (3, 3):
        raise ValueError("R must be a 3x3 matrix.")
    if project:
        R = _project_to_SO3(R)

    candidates = [("any", None, 0, np.eye(3))]
    for name, axis in _AXES.items():
        for ang in (90, 180):
            candidates.append((name, axis, ang, _rotmat_from_axis_angle(axis, ang)))

    best = ("any", 0, np.inf, np.eye(3))
    for name, axis, ang, Rt in candidates:
        err = _geodesic_angle(Rt.T @ R)  # radians
        if err < best[2]:
            best = (name, ang, err, Rt)

    axis_name, angle_deg, err_rad, R_best = best
    if return_matrix:
        return axis_name, angle_deg, np.rad2deg(err_rad), R_best
    else:
        return axis_name, angle_deg, np.rad2deg(err_rad)

# Example:
# axis, angle, err_deg = closest_canonical_rotation(Rtot)
# print(axis, angle, f"{err_deg:.3f}°")
