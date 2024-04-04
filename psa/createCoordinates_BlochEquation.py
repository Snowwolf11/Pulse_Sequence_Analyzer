#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 18:44:02 2024

@author: leon
"""
from psa.solve_bloch_equation import *

import numpy as np

def createCoordinates_BlochEquation(PS, T, l, Umax, offset, inpoFact, initialVector):
    CM = np.ones((PS.shape[0] + 1, 3))
    M0 = initialVector
    mx0, my0, mz0 = M0
    cn = np.array([0, 0, 0])
    CM[0] = cn

    for ind in range(1, PS.shape[0] + 1):
        x0, y0, z0 = cn
        Mx0 = -mx0
        My0 = my0
        Mz0 = mz0
        phi = PS[ind - 1, 1] * (np.pi / 180)
        Om = 2 * np.pi * Umax * inpoFact / 100 * PS[ind - 1, 0]
        x0, y0, z0, mx0, my0, mz0 = solve_bloch_equation(x0, y0, z0, Mx0, My0, Mz0, phi, Om, 2 * np.pi * offset, T / inpoFact)
        cn = np.array([x0, y0, z0])
        direction = l * (cn - CM[ind - 1]) / np.linalg.norm(cn - CM[ind - 1])
        cn = CM[ind - 1] + direction
        CM[ind] = cn

    return CM