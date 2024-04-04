#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

import numpy as np
cimport numpy as np

from cython.parallel import prange

def createVectors_Matrix_Cython(np.ndarray[double, ndim=2] PS, double T, double l, double Umax, double offset, int inpoFact, np.ndarray[double, ndim=1] initialVector):
    cdef int num_vectors = PS.shape[0] + 1
    cdef np.ndarray[double, ndim=2] m = np.empty((num_vectors, 3), dtype=np.float64)
    
    # Initialize the first vector
    m[0] = l * initialVector
    
    # Pre-calculate constants
    cdef double angle_factor = <double>(-2 * np.pi * T / inpoFact)
    cdef double offset_factor = <double>(Umax * inpoFact / 100)
    
    # Calculate rotation matrices
    cdef np.ndarray[double, ndim=3] rotation_matrices = np.zeros((PS.shape[0], 3, 3), dtype=np.float64)
    cdef int i
    cdef double Ux, Uy, Uz, norm_n, cosa, sina, mcosa
    cdef np.ndarray[double, ndim=1] n
    for i in range(PS.shape[0]):#, nogil=True):
        Ux = PS[i, 0] * offset_factor * np.cos(np.radians(PS[i, 1]))
        Uy = PS[i, 0] * offset_factor * np.sin(np.radians(PS[i, 1]))
        Uz = offset
        
        n = np.array([Ux, Uy, Uz])
        norm_n = np.linalg.norm(n)
        if norm_n != 0:
            n /= norm_n
            cosa = np.cos(angle_factor * norm_n)
            sina = np.sin(angle_factor * norm_n)
            mcosa = 1 - cosa
            rotation_matrices[i] = np.array([[n[0]**2 * mcosa + cosa, n[0] * n[1] * mcosa - n[2] * sina, n[0] * n[2] * mcosa + n[1] * sina],
                                              [n[0] * n[1] * mcosa + n[2] * sina, n[1]**2 * mcosa + cosa, n[1] * n[2] * mcosa - n[0] * sina],
                                              [n[2] * n[0] * mcosa - n[1] * sina, n[2] * n[1] * mcosa + n[0] * sina, n[2]**2 * mcosa + cosa]])
        else:
            rotation_matrices[i] = np.eye(3)

    # Apply rotation matrices
    cdef np.ndarray[double, ndim=1] vn = m[0]
    for i in range(PS.shape[0]):
        vn = np.dot(rotation_matrices[i], vn)
        m[i + 1] = vn

    return m
"""

import numpy as np
cimport numpy as np
from cython.parallel import prange

def createVectors_Matrix_Cython(np.ndarray[double, ndim=2] PS, double T, double l, double Umax, double offset, int inpoFact, np.ndarray[double, ndim=1] initialVector):
    cdef int num_vectors = PS.shape[0] + 1
    cdef np.ndarray[double, ndim=2] m = np.empty((num_vectors, 3), dtype=np.float64)
    
    # Initialize the first vector
    m[0] = l * initialVector
    
    # Pre-calculate constants
    cdef double angle_factor = <double>(-2 * np.pi * T / inpoFact)
    cdef double offset_factor = <double>(Umax * inpoFact / 100)
    
    # Calculate rotation matrices in parallel
    cdef np.ndarray[double, ndim=3] rotation_matrices = np.zeros((PS.shape[0], 3, 3), dtype=np.float64)
    cdef double Ux, Uy, Uz, norm_n, cosa, sina, mcosa
    cdef np.ndarray[double, ndim=1] n
    for i in prange(PS.shape[0], nogil=True):
        Ux = PS[i, 0] * offset_factor * np.cos(np.radians(PS[i, 1]))
        Uy = PS[i, 0] * offset_factor * np.sin(np.radians(PS[i, 1]))
        Uz = offset
        
        n = np.array([Ux, Uy, Uz])
        norm_n = np.linalg.norm(n)
        if norm_n != 0:
            n /= norm_n
            cosa = np.cos(angle_factor * norm_n)
            sina = np.sin(angle_factor * norm_n)
            mcosa = 1 - cosa
            rotation_matrices[i] = np.array([[n[0]**2 * mcosa + cosa, n[0] * n[1] * mcosa - n[2] * sina, n[0] * n[2] * mcosa + n[1] * sina],
                                              [n[0] * n[1] * mcosa + n[2] * sina, n[1]**2 * mcosa + cosa, n[1] * n[2] * mcosa - n[0] * sina],
                                              [n[2] * n[0] * mcosa - n[1] * sina, n[2] * n[1] * mcosa + n[0] * sina, n[2]**2 * mcosa + cosa]])
        else:
            rotation_matrices[i] = np.eye(3)

    # Apply rotation matrices
    cdef np.ndarray[double, ndim=1] vn = m[0]
    for i in prange(PS.shape[0], nogil=True):
        vn = np.dot(rotation_matrices[i], vn)
        m[i + 1] = vn

    return m
"""
