#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 18:28:37 2024

@author: leon
"""

import numpy as np
from scipy.integrate import odeint

def solve_bloch_equation(x0 = 0.0, y0 = 0.0, z0 = 0.0, Mx0 = 0.0, My0 = 0.0, Mz0 =1.0, phi = 0.0, Om = 2.0, off = 1.0, T = 10.0):
    
    w1 = 0.0
    wL = w1 - off

    M0 = [Mx0, My0, Mz0]
    t_span = [0, T]

    def diffM(M, t):
        return [wL * M[1] + Om * np.sin(w1 * t + phi) * M[2],
                -wL * M[0] + Om * np.cos(w1 * t + phi) * M[2],
                -Om * np.sin(w1 * t + phi) * M[0] - Om * np.cos(w1 * t + phi) * M[1]]

    t = np.linspace(0, T, 100)
    y = odeint(diffM, M0, t)

    y = np.column_stack((-y[:, 0], y[:, 1], y[:, 2]))

    AHT = np.zeros((y.shape[0] + 1, 3))
    cn = np.array([0.0, 0.0, 0.0])

    AHT[0] = cn
    for n in range(1, y.shape[0] + 1):
        cn += y[n - 1]
        AHT[n] = cn

    a1, a2, a3 = AHT[0]
    b1, b2, b3 = AHT[-1]

    xT = b1 + x0 - a1
    yT = b2 + y0 - a2
    zT = b3 + z0 - a3
    MxT, MyT, MzT = y[-1]

    return xT, yT, zT, MxT, MyT, MzT

[xT, yT, zT, MxT, MyT, MzT] = solve_bloch_equation()