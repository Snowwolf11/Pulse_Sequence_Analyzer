import numpy as np

def createCoordinates_HelixV2(PS,T,l,maximumAmplitude, offset, inpoFact, initialVector):
     # calculates the space curve using the rotation matrix method
    PS[PS == 0] = 1e-16
    N = len(PS)
    resolution = 2*inpoFact
    offset = 2 * np.pi * offset
    Rabi_frequencies = 2 * np.pi * PS[:,0] * maximumAmplitude/100
    phases = PS[:,1] *  np.pi / 180

    starting_tangent = initialVector
    starting_binormal = np.array([-np.cos(phases[1]), -np.sin(phases[1]), 0])
    starting_point = np.array([0, 0, 0])

    helix_array = []

    # Create helices
    for k in range(N):
        Omega_k = Rabi_frequencies[k]
        phi_k = phases[k]
        B_eff = np.array([Omega_k * np.cos(phi_k), Omega_k * np.sin(phi_k), offset])
        Omega_eff = np.linalg.norm(B_eff)
    
        radius = Omega_k / Omega_eff**2
        pitch = 2 * np.pi * offset / Omega_eff**2
        arc_length = T
    
        helix, ending_tangent, ending_binormal = create_helix(
            starting_point, starting_tangent, starting_binormal, radius, pitch, arc_length, -1, resolution
        )
    
        starting_point = helix[-1]
        starting_tangent = ending_tangent
        starting_binormal = ending_binormal
    
        if k < N - 1:
            theta = (phases[k + 1] - phases[k])
            K = skew_symmetric_matrix(starting_tangent)
            I = np.eye(3)
            R = I + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
            starting_binormal = R @ starting_binormal
    
        helix_array.append(helix.T)

    # Combine to one curve
    return np.hstack(helix_array).T*l

def skew_symmetric_matrix(v):
    """
    Compute the skew-symmetric matrix for a vector v.
    """
    K = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    return K

def create_helix(starting_point, starting_tangent, starting_binormal, radius, pitch, arc_length, chirality, resolution):
   # Create a helix given the parameters.

    t = np.linspace(0, arc_length / np.sqrt((radius * 2 * np.pi)**2 + pitch**2), resolution)
    standard_helix = np.array([
        radius * np.cos(2 * np.pi * t * chirality),
        radius * np.sin(2 * np.pi * t * chirality),
        pitch * t
    ])
    
    # Compute the starting tangent of the standard helix
    starting_tangent_standard = np.array([0, radius * 2 * np.pi * chirality, pitch])
    starting_tangent_standard /= np.linalg.norm(starting_tangent_standard)
    
    # Compute the starting normal of the standard helix
    starting_normal_standard = np.array([-radius * 4 * np.pi**2, 0, 0])
    starting_normal_standard /= np.linalg.norm(starting_normal_standard)
    
    # Compute the starting binormal of the standard helix
    starting_binormal_standard = np.cross(starting_tangent_standard, starting_normal_standard)
    starting_binormal_standard /= np.linalg.norm(starting_binormal_standard)
    
    # Compute the ending tangent of the standard helix
    ending_tangent_standard = np.array([
        -radius * 2 * np.pi * chirality * np.sin(2 * np.pi * t[-1] * chirality),
        radius * 2 * np.pi * chirality * np.cos(2 * np.pi * t[-1] * chirality),
        pitch
    ])
    ending_tangent_standard /= np.linalg.norm(ending_tangent_standard)
    
    # Compute the ending normal of the standard helix
    ending_normal_standard = np.array([
        -radius * 4 * np.pi**2 * np.cos(2 * np.pi * t[-1] * chirality),
        -radius * 4 * np.pi**2 * np.sin(2 * np.pi * t[-1] * chirality),
        0
    ])
    ending_normal_standard /= np.linalg.norm(ending_normal_standard)
    
    # Compute the ending binormal of the standard helix
    ending_binormal_standard = np.cross(ending_tangent_standard, ending_normal_standard)
    ending_binormal_standard /= np.linalg.norm(ending_binormal_standard)
    
    starting_normal = np.cross(starting_binormal, starting_tangent)
    
    # Compute the rotation matrix to align the standard helix with the given vectors
    R1 = np.column_stack((starting_normal_standard, starting_binormal_standard, starting_tangent_standard)).T
    R2 = np.column_stack((starting_normal, starting_binormal, starting_tangent))
    rotation_matrix = np.dot(R2,R1)
    
    # Apply the rotation to the standard helix
    rotated_helix = (rotation_matrix @ standard_helix).T
    
    # Translate the helix to the starting point
    helix = rotated_helix + starting_point - rotated_helix[0]
    
    # Compute the ending tangent and binormal for the rotated helix
    ending_tangent = rotation_matrix @ ending_tangent_standard
    ending_tangent /= np.linalg.norm(ending_tangent)
    ending_binormal = rotation_matrix @ ending_binormal_standard
    ending_binormal /= np.linalg.norm(ending_binormal)
    
    return helix, ending_tangent, ending_binormal
