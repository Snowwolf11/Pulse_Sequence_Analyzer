#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
#import cProfile	#profiling
#import pstats	#profiling
import time

import numpy as np

#Etwa 6 mal langsamer als matlab

def _createVectors_Matrix(PS, T, l, Umax, offset, inpoFact, initialVector):
    #time_start = time.time()
    # Pre-allocate memory for the output matrix
    num_vectors = len(PS) + 1
    m = np.empty((num_vectors, 3))

    # Initialize the first vector
    m[0] = l * initialVector

    # Pre-calculate constants
    angle_factor = -2 * np.pi * T / inpoFact
    Umax_factor = Umax * inpoFact / 100

    vn = m[0]

    for i, v in enumerate(PS):
        Ux = v[0] * Umax_factor * np.cos(np.radians(v[1]))
        Uy = v[0] * Umax_factor * np.sin(np.radians(v[1]))
        Uz = offset

        n = np.array([Ux, Uy, Uz])

        norm_n = np.linalg.norm(n)
        if norm_n != 0:
            n /= norm_n
            cosa = np.cos(angle_factor * norm_n)
            sina = np.sin(angle_factor * norm_n)
            mcosa = 1 - cosa
            Rn = np.array([[n[0]**2 * mcosa + cosa, n[0] * n[1] * mcosa - n[2] * sina, n[0] * n[2] * mcosa + n[1] * sina],
                           [n[0] * n[1] * mcosa + n[2] * sina, n[1]**2 * mcosa + cosa, n[1] * n[2] * mcosa - n[0] * sina],
                           [n[2] * n[0] * mcosa - n[1] * sina, n[2] * n[1] * mcosa + n[0] * sina, n[2]**2 * mcosa + cosa]])
        else:
            Rn = np.eye(3)

        vn = np.dot(Rn, vn)
        m[i + 1] = vn
        
    #time_end = time.time()
    #print("Ellapsed time: " + str(time_end-time_start))
    return m

"""
#import cProfile	#profiling
#import pstats	#profiling
import time

import numpy as np

#Etwa 6 mal langsamer als matlab

def createVectors_Matrix(PS, T, l, Umax, offset, inpoFact, initialVector):
    time_start = time.time()
    # Pre-allocate memory for the output matrix
    num_vectors = len(PS) + 1
    m = np.empty((num_vectors, 3))

    # Initialize the first vector
    m[0] = l * initialVector

    # Pre-calculate constants
    angle_factor = -2 * np.pi * T / inpoFact
    Umax_factor = Umax * inpoFact / 100

    #vn = m[0]
	
    rotation_matrices = np.zeros((len(PS), 3, 3))
    for i, v in enumerate(PS):
        Ux = v[0] * Umax_factor * np.cos(np.radians(v[1]))
        Uy = v[0] * Umax_factor * np.sin(np.radians(v[1]))
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

        #vn = np.dot(Rn, vn)
        #m[i + 1] = vn
    
    vn = m[0]
    for i, Rn in enumerate(rotation_matrices):
        vn = np.dot(Rn, vn)
        m[i + 1] = vn    
    time_end = time.time()
    print("Ellapsed time: " + str(time_end-time_start))
    return m

"""
