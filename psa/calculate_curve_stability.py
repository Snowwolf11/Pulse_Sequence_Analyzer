import vectors
import totalRotMatrix

from psa.createVectors_Matrix import *
from psa.calculate_Rtot import *
from psa.calculate_pulse_sequence_quality import*

import time

import numpy as np
#46.5 vs 33
def calculate_curve_stability(PS, T, l, Umax, scalingRange_percent=20, offsetRange_kHz=20, stabilityCalculationMethod=1, initialVector=np.array([0,0,1]), language = "Rust"):
    a = time.time()
    
    offsetRange = 1e3 * np.linspace(-offsetRange_kHz/2, offsetRange_kHz/2, 200)
    
    scalingRange = np.linspace(1-scalingRange_percent/200,  1+scalingRange_percent/200, 200)
    scalingRange *= Umax
    
    #print("off: "+str(offsetRange))
    #print("amp: "+str(scalingRange))
    X, Y = np.meshgrid(scalingRange, offsetRange)
    Z = np.ones_like(X)
    counts = 0
    counto = 0
    countZ = 0
    
    if language == "Rust":
        VM =  vectors.createVectors_Matrix(PS.astype(np.float64),np.float64(T),np.float64(l),np.float64(Umax),np.float64(0),1,initialVector.astype(np.float64))
        Rtot_on_res = totalRotMatrix.create_Rtot(PS.astype(np.float64), np.float64(T), np.float64(Umax), np.float64(0.0), 1)
        medianVec = VM[-1,:]
    elif language == "Python":
        VM =  createVectors_Matrix_python(PS.astype(np.float64),np.float64(T),np.float64(l),np.float64(Umax),np.float64(0),1,initialVector.astype(np.float64))
        Rtot_on_res = create_Rtot_python(PS.astype(np.float64), np.float64(T), np.float64(Umax), np.float64(0.0), 1)
        medianVec = VM[-1]
    
    if np.linalg.norm(medianVec[2]) > 5 * (np.linalg.norm(medianVec[1]) + np.linalg.norm(medianVec[0])):
        if medianVec[2] > 0:
            shallVec = np.array([0, 0, 1])
        elif medianVec[2] < 0:
            shallVec = np.array([0, 0, -1])
    else:
        shallVec = medianVec
    
    axis, angle, err_deg, R_target = closest_canonical_rotation(Rtot_on_res, return_matrix=True)

    if stabilityCalculationMethod == 1: #app.SSButton.Value == 1:	#################################
        for (i,n) in np.ndenumerate(scalingRange):
            counts += 1
            for (j,n2) in np.ndenumerate(offsetRange):
                counto += 1
                if counto > len(offsetRange):
                    counto = 1
                Z[j,i] = calculate_pulse_sequence_quality(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(n), np.float64(n2), 1, initialVector.astype(np.float64), shallVec.astype(np.float64), 1, language=language)
                countZ += Z[j, i]
    
    if stabilityCalculationMethod == 2: #app.URButton.Value == 1:	###############################
        for (i,n) in np.ndenumerate(scalingRange):
            counts += 1
            for (j,n2) in np.ndenumerate(offsetRange):
                counto += 1
                if counto > len(offsetRange):
                    counto = 1
                Z[j,i] = calculate_pulse_sequence_quality(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(n), np.float64(n2), 1, initialVector.astype(np.float64), shallVec.astype(np.float64), 2, target_Rotation = R_target.astype(np.float64), language=language)
                countZ += Z[j, i]
                
                #(Rtot*ex
                
    if stabilityCalculationMethod == 3: #app.AngleButton.Value == 1:	###################################
        for (i,n) in np.ndenumerate(scalingRange):
            counts += 1
            for (j,n2) in np.ndenumerate(offsetRange):
                counto += 1
                if counto > len(offsetRange):
                    counto = 1
                if language == "Rust":
                    VM = vectors.createVectors_Matrix(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(n), np.float64(n2), 1, initialVector.astype(np.float64))
                    Z[j, i] = np.arctan2(np.sqrt(VM[-1,0]**2 + VM[-1,1]**2), VM[-1,2]) * 180 / np.pi
                elif language == "Python":
                    VM = createVectors_Matrix_python(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(n), np.float64(n2), 1, initialVector.astype(np.float64))
                    Z[j, i] = np.arctan2(np.sqrt(VM[-1][0]**2 + VM[-1][1]**2), VM[-1][2]) * 180 / np.pi
                countZ += Z[j, i]
    
    quality = countZ / (counto * counts)

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
