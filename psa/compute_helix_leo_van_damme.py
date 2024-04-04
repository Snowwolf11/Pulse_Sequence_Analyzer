#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:20:56 2024

@author: leon
"""

import numpy as np

def compute_helix_leo_van_damme(w, vk, Xk, nk, tk, t):
    # w: effective field
    # vk: starting Tangent
    # Xk: starting Point
    # nk: rotation axis (given by direction of effective field)
    # tk: initial time
    # t: time of interest
    
    cos_term = np.cos(w * (t - tk))
    sin_term = np.sin(w * (t - tk))
    
    v = (cos_term * vk +
         sin_term * np.cross(vk[:,0], nk[:,0]).reshape(3,1) +
         (1 - cos_term) * np.dot(nk[:,0], vk[:,0]) * nk)
    
    X = (sin_term / w * vk +
         (1 - cos_term) / w * np.cross(vk[:,0], nk[:,0]).reshape(3,1) +
         (t - tk - sin_term / w) * np.dot(nk[:,0], vk[:,0]) * nk +
         Xk)
    
    return X, v

#X is Helix with:
#       r = Om_s/w^2;
#       h = 2*pi*r*tan(m)
#       Axis_dir : nk (normalized)
#       Starting_tangent (for tk=0) : vk  (normalized)
#       vk_s := cross(n,cross(v,n)/norm(cross(v,n))) (normalized)
#       Om_s := w*norm(cross(v,n))      (field orthonormal to v)
