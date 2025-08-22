#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
import vectors

from psa.createVectors_Matrix import *
import numpy as np
import time

def createCoordinates_Matrix(PS,T,l,Umax,offset,inpoFact,initialVector):
    #all inputs eather floats or ndarrays
    #createCoordinateMatrix create a Matrix where each wor accords to a set of
    #3D Coordinates
    #   Detailed explanation goes here
    
    #a = time.time()
    VM = vectors.createVectors_Matrix(PS.astype(np.float64),np.float64(T),np.float64(l),np.float64(Umax),np.float64(offset),int(inpoFact),initialVector.astype(np.float64))
    #b = time.time()
    #print("Time: " + str(b-a))
    h = np.ones((np.shape(VM)[0]+1,3))
    cn=np.array([0,0,0])      #first Coordinates
    h[0,:]=cn
    for n in range(np.shape(VM)[0]):
        cn=cn+VM[n,:]      #Coordinates = Coordinates of point befor + the fitting Vector from VM
        h[n+1,:]=cn
    CM=h
    return CM
    