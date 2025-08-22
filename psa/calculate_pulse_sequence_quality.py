import vectors
import totalRotMatrix

from psa.createCoordinates_Matrix import *

import numpy as np

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

def calculate_pulse_sequence_quality(PS, T, l, Umax, off, inpoFact=1 , initialVector = np.array([0,0,1]), shallVec = np.array([0,0,1]), calcType = 1, target_Rotation = np.eye(3)):
    
    if calcType == 1:
       VM = vectors.createVectors_Matrix(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(Umax), np.float64(off), 1, initialVector.astype(np.float64))
       Q = np.dot(VM[-1,:], shallVec)
       
    elif calcType == 2:
       Rtot = totalRotMatrix.create_Rtot(PS.astype(np.float64), np.float64(T), np.float64(Umax), np.float64(off), 1)
       Q = 1/3 * (np.dot(Rtot[:,0], [1, 0, 0]) + np.dot(Rtot[:,1], [0, 1, 0]) + np.dot(Rtot[:,2], [0, 0, 1]))  # for 0° & 360° Pulse
       #Q = 1/3 * (np.dot(Rtot[:,0], [1, 0, 0]) + np.dot(Rtot[:,1], [0, 1, 0]) + np.dot(Rtot[:,2], [0, 0, 1]))  # for 0° & 360° Pulse
       #R_target = rotation_matrix([1,0,0], 180)
       Q = np.trace(target_Rotation.T @ Rtot) / 3.0
    
    elif calcType == 3:  #geschlossenheit
       CM = createCoordinates_Matrix(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(Umax), np.float64(off), 1, initialVector.astype(np.float64))
       Q = 2*np.sqrt(np.dot(CM[-1,:], CM[-1,:]))#/(l*(PS.shape[0]+1))

    return Q
