import os

import numpy as np

from psa.calculate_pulse_sequence_quality import calculate_pulse_sequence_quality_batch
from psa.calculate_curve_stability import closest_canonical_rotation
from psa.getPulseSequence import getPulseSequence
from psa.rotations import propagate_batch, rotate, inverse_total_rotation


def calculate_pulse_sequence_quality_images(
    dirname,
    Range,
    T,
    Umax,
    initialVector,
    calcType,
    changingVariable,
    Resolution,
    progress_callback=None,
    initialMagnetization=np.array([0, 0, 1]),
):
    """log10(|1 - quality|) for every .bruker file in ``dirname`` (columns) over an amplitude
    (changingVariable 1) or offset (changingVariable 2) scan of width ``Range`` Hz (rows).

    calcType 1 (SS/PP): target = on-resonance final magnetization R m0 of the LAST file, rounded to
    45 degree steps; m0 = ``initialMagnetization``. calcType 2 (UR), 3 (closedness, noise axis
    ``initialVector``)."""
    files = os.listdir(dirname)
    sorted_files = sorted(files)
    PSnames = [f for f in sorted_files if 'bruker' in f]
    PSnumber = len(PSnames)
    offsetRes = round(Range * Resolution/1000)
    variation = np.linspace(-Range / 2, Range / 2, offsetRes)
    initialVector = np.asarray(initialVector, dtype=np.float64)
    m0 = np.asarray(initialMagnetization, dtype=np.float64)
    m0 = m0 / np.linalg.norm(m0)

    PS_1 = np.array(getPulseSequence(os.path.join(dirname, PSnames[-1])))
    M_1 = rotate(propagate_batch(PS_1, T, Umax, 0.0), m0)
    Rtot_on_res = inverse_total_rotation(PS_1, T, Umax, 0.0, 1)
    shallVec = round_vector(M_1, 45)
    print("Original Vector:", M_1)

    axis, angle, err_deg, R_target = closest_canonical_rotation(Rtot_on_res, return_matrix=True)

    print("Rounded Cartesian Vector:", shallVec)

    PI = np.ones((offsetRes, PSnumber))
    if changingVariable in (1, 2):
        for n1 in range(PSnumber):
            PS = np.array(getPulseSequence(os.path.join(dirname, PSnames[n1])))
            if progress_callback is not None:
                progress_callback(f"{n1 + 1} von {PSnumber}")
            if changingVariable == 1:
                Q = calculate_pulse_sequence_quality_batch(PS, T, 1, Umax + variation, 0.0, initialVector, shallVec, calcType, initialMagnetization=m0)
            else:
                Q = calculate_pulse_sequence_quality_batch(PS, T, 1, Umax, variation, initialVector, shallVec, calcType, target_Rotation=R_target, initialMagnetization=m0)
            Q = np.where(np.abs(1 - Q) < 10**(-15), 1 - 10**(-15), Q)
            PI[:, n1] = np.log10(np.abs(1 - Q))
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