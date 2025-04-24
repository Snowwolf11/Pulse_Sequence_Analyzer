#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:11:41 2024

@author: leon
"""
from psa.old.compute_helix_leo_van_damme import *

import numpy as np

def createCoordinates_HelixVanDamme(PS,T,l,maximumAmplitude, offset, inpoFact, initialVector):
           
   maximumAmplitude *= 2 * np.pi
   offset *= 2 * np.pi  # wtfff????
   
   # Replace zeros with small values to avoid division by zero
   PS[PS == 0] = 1e-16
   
   t = np.linspace(0, T, inpoFact + 1)
   X0 = np.array([[0],[0], [0]])
   m = np.zeros((inpoFact * PS.shape[0], 3))
   v0 = initialVector.reshape(3,1)
   vn = v0
   
   for i in range(PS.shape[0]):
       Ux = (PS[i, 0] / 100) * maximumAmplitude * np.cos(PS[i, 1] * (np.pi / 180))
       Uy = (PS[i, 0] / 100) * maximumAmplitude * np.sin(PS[i, 1] * (np.pi / 180))
       Uz = offset
       w = np.sqrt(Ux**2 + Uy**2 + Uz**2)
       n = np.transpose(np.array([[Ux, Uy, Uz]]) / w)
       X, v = compute_helix_leo_van_damme(w, vn, X0, n, t[0], t)
       
       X0 = X[:, -1].reshape(3,1)
       vn = v[:, -1].reshape(3,1)
       
       m[i * inpoFact : (i + 1) * inpoFact, :] = X[:, :-1].T      #sketch!!!
       t += T
   
   return l * m

    
   
       

