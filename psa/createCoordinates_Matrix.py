#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
import numpy as np

from psa.rotations import tangent_vectors


def createCoordinates_Matrix(PS, T, l, Umax, offset, inpoFact, initialVector):
    """Error-curve points: cumulative sum of the tangent vectors, starting at the origin."""
    VM = tangent_vectors(PS, T, l, Umax, offset, int(inpoFact), initialVector)
    CM = np.zeros((len(VM) + 1, 3))
    np.cumsum(VM, axis=0, out=CM[1:])
    return CM
