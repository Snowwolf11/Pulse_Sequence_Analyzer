import vectors

from psa.createVectors_Matrix import *
from psa.calculate_Rtot import *
from psa.calculate_pulse_sequence_quality import*

import time

import numpy as np
#46.5 vs 33
def calculate_curve_stability(PS, T, l, Umax, scalingRange_percent=20, offsetRange_kHz=20, stabilityCalculationMethod=1, initialVector=np.array([0,0,1])):
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
    
    VM =  vectors.createVectors_Matrix(PS.astype(np.float64),np.float64(T),np.float64(l),np.float64(Umax),np.float64(0),1,initialVector.astype(np.float64))
    medianVec = VM[-1,:]
    
    if np.linalg.norm(medianVec[2]) > 5 * (np.linalg.norm(medianVec[1]) + np.linalg.norm(medianVec[0])):
        if medianVec[2] > 0:
            shallVec = np.array([0, 0, 1])
        elif medianVec[2] < 0:
            shallVec = np.array([0, 0, -1])
    else:
        shallVec = medianVec
    
    if stabilityCalculationMethod == 1: #app.SSButton.Value == 1:	#################################
        for (i,n) in np.ndenumerate(scalingRange):
            counts += 1
            for (j,n2) in np.ndenumerate(offsetRange):
                counto += 1
                if counto > len(offsetRange):
                    counto = 1
                Z[j,i] = calculate_pulse_sequence_quality(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(n), np.float64(n2), 1, initialVector.astype(np.float64), shallVec.astype(np.float64), 1)
                countZ += Z[j, i]
    
    if stabilityCalculationMethod == 2: #app.URButton.Value == 1:	###############################
        for (i,n) in np.ndenumerate(scalingRange):
            counts += 1
            for (j,n2) in np.ndenumerate(offsetRange):
                counto += 1
                if counto > len(offsetRange):
                    counto = 1
                Z[j,i] = calculate_pulse_sequence_quality(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(n), np.float64(n2), 1, initialVector.astype(np.float64), shallVec.astype(np.float64), 2)
                countZ += Z[j, i]
                
                #(Rtot*ex
                
    if stabilityCalculationMethod == 3: #app.AngleButton.Value == 1:	###################################
        for (i,n) in np.ndenumerate(scalingRange):
            counts += 1
            for (j,n2) in np.ndenumerate(offsetRange):
                counto += 1
                if counto > len(offsetRange):
                    counto = 1
                VM = vectors.createVectors_Matrix(PS.astype(np.float64), np.float64(T), np.float64(l), np.float64(n), np.float64(n2), 1, initialVector.astype(np.float64))
                Z[j, i] = np.arctan2(np.sqrt(VM[-1,0]**2 + VM[-1,1]**2), VM[-1,2]) * 180 / np.pi
                countZ += Z[j, i]
    
    quality = countZ / (counto * counts)

    b = time.time()
    print("Time: " + str(b-a))
    return X, Y, Z, quality
