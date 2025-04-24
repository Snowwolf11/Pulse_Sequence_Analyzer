from psa.calculate_pulse_sequence_quality import *
from psa.getPulseSequence import *

import vectors
import os
import numpy as np

def calculate_pulse_sequence_quality_images(GUI, dirname, Range, T, Umax, initialVector, calcType, changingVariable, Resolution):
    files = os.listdir(dirname)
    sorted_files = sorted(files)
    PSnames = [f for f in sorted_files if 'bruker' in f]
    PSnumber = len(PSnames)
    offsetRes = round(Range * Resolution/1000)
    offset = Range / offsetRes
    PS_1 = np.array(getPulseSequence(os.path.join(dirname, PSnames[-1])))
    VM_1 = vectors.createVectors_Matrix(PS_1.astype(np.float64), np.float64(T), np.float64(1.0), np.float64(Umax), np.float64(0.0), 1, initialVector.astype(np.float64))
    shallVec = round_vector(VM_1[-1,:], 45)#np.array([0, 0, -1])  # Define shallVec as a numpy array

    print("Original Vector:", VM_1[-1,:])
    print("Rounded Cartesian Vector:", shallVec)

    if changingVariable == 1:
        PI = np.ones((offsetRes, PSnumber))
        for n1 in range(PSnumber):
            PS = np.array(getPulseSequence(os.path.join(dirname, PSnames[n1])))
            GUI.text_area_set(text_area = GUI.info_error_text, text_str = f"{n1 + 1} von {PSnumber}", reset_bool = 0)
            for n2 in range(offsetRes):
                Q = calculate_pulse_sequence_quality(PS, T, 1, Umax + (-Range / 2 + (n2 - 1) * offset), 0, 1, np.array(initialVector), shallVec, calcType)
                #print("Q: "+str(Q))
                if abs(1 - Q) < 10**(-15):
                    Q = 1-10**(-15)
                PI[n2, n1] = np.log10(abs(1 - Q))

    if changingVariable == 2:
        PI = np.ones((offsetRes, PSnumber))
        for n1 in range(PSnumber):
            PS = np.array(getPulseSequence(os.path.join(dirname, PSnames[n1])))
            GUI.text_area_set(text_area = GUI.info_error_text, text_str = f"{n1 + 1} von {PSnumber}", reset_bool = 0)
            for n2 in range(offsetRes):
                Q = calculate_pulse_sequence_quality(PS, T, 1, Umax, -Range / 2 + (n2 - 1) * offset, 1, np.array(initialVector), shallVec, calcType)
                if abs(1 - Q) < 10**(-15):
                    Q = 1-10**(-15)
                #PI[n2, n1] = np.log10(abs(1-Q))
                PI[n2, n1] = np.log10(abs(1-Q))
    #PI = mean_deviation(PI, round(Resolution*4))
    #print(Resolution)
    return PI

def round_vector(v, nearest_degree):
    # Step 1: Convert to polar coordinates
    r = np.linalg.norm(v)
    theta = np.arctan2(v[1], v[0])
    phi = np.arccos(v[2] / r)

    # Step 2: Normalize the vector
    v_normalized = v / r

    # Step 3: Round angles to the nearest 45 degrees
    theta_deg = np.degrees(theta)
    phi_deg = np.degrees(phi)

    theta_rounded = np.round(theta_deg / nearest_degree) * nearest_degree
    phi_rounded = np.round(phi_deg / nearest_degree) * nearest_degree

    theta_rounded_rad = np.radians(theta_rounded)
    phi_rounded_rad = np.radians(phi_rounded)

    # Step 4: Convert back to Cartesian coordinates
    x_rounded = np.sin(phi_rounded_rad) * np.cos(theta_rounded_rad)
    y_rounded = np.sin(phi_rounded_rad) * np.sin(theta_rounded_rad)
    z_rounded = np.cos(phi_rounded_rad)

    v_rounded = np.array([x_rounded, y_rounded, z_rounded])

    return v_rounded

def mean_deviation(Q, k):
    n, m = Q.shape
    Q_av = np.zeros_like(Q, dtype=float)
    
    for i in range(n):
        for j in range(m):
            # Define the range considering the boundaries
            start = max(0, i - k)
            end = min(n, i + k + 1)
            
            # Extract the relevant sub-array
            sub_array = Q[start:end, j]
            
            # Calculate the mean deviation from Q[i, j]
            deviations = np.abs(sub_array - Q[i, j])
            mean_deviation = np.mean(deviations)
            
            # Assign the calculated mean deviation to Q_av
            Q_av[i, j] = np.log10(mean_deviation)
    
    return Q_av