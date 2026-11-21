#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
import numpy as np

from psa.rotations import tangent_vectors


def trapezoid_curve(VM):
    """Error-curve points from the M + 1 tangent step vectors t_0 ... t_M (already scaled to the
    step length): r_0 = 0, r_k = sum_{i<k} (t_i + t_{i+1}) / 2. Returns M + 1 points for M steps."""
    VM = np.asarray(VM, dtype=np.float64)
    CM = np.zeros((len(VM), 3))
    np.cumsum((VM[:-1] + VM[1:]) / 2, axis=0, out=CM[1:])
    return CM


def createCoordinates_Matrix(PS, T, l, Umax, offset, inpoFact, initialVector):
    """Error-curve points (rotation-matrix method, trapezoid rule), starting at the origin."""
    return trapezoid_curve(tangent_vectors(PS, T, l, Umax, offset, int(inpoFact), initialVector))
