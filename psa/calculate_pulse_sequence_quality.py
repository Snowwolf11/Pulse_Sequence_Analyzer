from vectors import createVectors_Matrix
import totalRotMatrix

import numpy as np

def calculate_pulse_sequence_quality(PS, T, l, Umax, off, inpoFact=1 , initialVector = np.array([0,0,1]), shallVec = np.array([0,0,1]), calcType = 1):
    
    if calcType == 1:
       VM = createVectors_Matrix(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(Umax), np.float64(off), 1, initialVector.astype(np.float64))
       #print("shallvec: "+str(shallVec))
       #print("VM end: "+str(VM[-1,:]))
       #print("Q: "+str(np.dot(VM[-1,:], shallVec)))
       Q = np.dot(VM[-1,:], shallVec)
       
    elif calcType == 2:
       Rtot = totalRotMatrix.create_Rtot(PS.astype(np.float64), np.float64(T), np.float64(Umax), np.float64(off), 1)
       Q = 1/3 * (np.dot(Rtot[:,0], [1, 0, 0]) + np.dot(Rtot[:,1], [0, 1, 0]) + np.dot(Rtot[:,2], [0, 0, 1]))  # for 0° & 360° Pulse

    return Q
