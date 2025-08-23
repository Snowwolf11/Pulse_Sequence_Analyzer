#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
from psa.createCoordinates_Matrix import *
from psa.old.createCoordinates_Helix import *
from psa.getPulseSequence import *
from psa.old.createCoordinates_HelixVanDamme import *
from psa.old.createCoordinates_BlochEquation import *
from psa.createCoordiantes_HelixV2 import *

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

import numpy as np


def  createCurve(PulseSequence,T=0.0000005,l=1,maximumAmplitude=10000,offset=8900, inpoFact=8,xExpand=0, calculationMethod=1, initialVector=np.array([0,0,1]), language = "Rust"):
            #print("offset: ", offset)
            #print("maxAmpl: ", maximumAmplitude)
            #createCurve creates a Curve from a raw pulse Sequence
            #   PulseSequence: Puls Sequence in form of a cell, a matrix or a char
            #                 (complete name of a .bruker file)
            #   T: Time of a single Puls (directly proportional to rotaion angle per
            #      puls) in seconds
            #   l: vector length per step
            #   maximumAmplitude: == 100# Puls intensity (==maximum Pulse power output)
            #                     in Hz
            #   offset: (y component of rotation axis) in Hz
            #   inpoFact: Interpolation Factor for smoothing the curve by showing
            #             intermediate rotation results
            #   xEpand:   :expands the curve in x Direction by adding a linearly
            #             increasing value
            #   path: logical (0:false 1:true) if path of vector-end should be plotted
            #   curvandtors: logical (0:false 1:true) if Curvature and Torsion should be plotted
            
            #   sample function call:
            #   createCurve('/Users/leon/Documents/MATLAB/Pulssequenzen/BIBOP_sorted_20kHz_noB1_rf10kHz/pulse0015.bruker',5*10^-7,1,10^4,0,10,0,1,1)
            #
            
           # profile on
           
            #reshape input
            #initialVector.resize((1,3))
            
            l=l*T/inpoFact     #That the total curve length stays the same if inpofact≠1
            PS = getPulseSequence(PulseSequence)
            PS_ori = PS
            PS=np.kron(PS, np.ones((inpoFact,1)));   #Each Element of the PS is exchanged by an itself * ones(inpoFact,1) to interpolate the Curve for better handability (it gets a better resolution)
            PS[:,0]=PS[:,0]/inpoFact;       #The Intensities of the PS get reduced by the factor inpoFact, to ensure the Curve stays the same (exept better resolution) for inpoFact≠0
            initialVector = initialVector/np.linalg.norm(initialVector)
            #VM=createVectorMatrix_app(app,PS,T,l,maximumAmplitude, offset, inpoFact); #A Matrix which contains all the vectors gets created
            if calculationMethod == 1:
              CM=createCoordinates_Matrix(PS,T,l,maximumAmplitude, offset, inpoFact, initialVector, language=language) #A Matrix which contains the Coordinate of Points of the Curve is created
            elif(calculationMethod == 2):
              CM = createCoordinates_HelixV2(PS_ori,T,l*inpoFact/T,maximumAmplitude, offset, inpoFact, initialVector)
            else:
              raise Exception("calculation Method must be either 1 (Rotation Matrix), 2 (Helix)")
            
            if not np.isreal(np.sum(CM)):
                raise ValueError("CM contains imaginary items!")
            if np.isnan(np.sum(CM)):
               raise ValueError("CM contains nan!")
            
            VM = calculate_VM(CM)
            
            #[l_Curv, l_Tors, CurvMat, TorsVec, CurvInt, TorsInt, TorsInt_abs] = calculateCurvatureAndTorsion(CM); #calculates the Curvature and The Torsion of the Curve
            
            #plot3DCurve(CM[:,0]+np.linspace(0,xExpand,np.shape(CM)[0]).transpose(), CM[:,1], CM[:,2], Title='AHT-Kurve')
            
            return CM, VM, PS

def calculate_VM(CM):
    size_CM = CM.shape[0]
    VM = np.ones((size_CM - 1, 3))
    for n in range(size_CM - 1):
        VM[n] = CM[n + 1] - CM[n]
    return VM
    
"""
def plot3DCurve(CurveX, CurveY, CurveZ, Title="", type="new"):
    #type: new or old (for new or old figure)
    if type == "new":
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        plot = ax.plot3D(CurveX, CurveY, CurveZ, 'blue')
        ax.set_title(Title)
        #ax.set_aspect('equal')
        ax.set_xlabel('x', labelpad=20)
        ax.set_ylabel('y', labelpad=20)
        ax.set_zlabel('z', labelpad=20)
        ax.set_box_aspect((abs(max(CurveX)-min(CurveX)), abs(max(CurveY)-min(CurveY)), abs(max(CurveZ)-min(CurveZ))))
    elif type == "old":
        ax = plt.gca()
        ax.plot3D(CurveX, CurveY, CurveZ)
    else:
        raise Exception("plot type must be new or old!")
    plt.show()
 """       
            
    
#createCurve("/Users/leon/Desktop/Physik/Glaser/Analyse_und_Visualisierung_von_robusten_Kontrollpulsen/Pulssequenzen/UR_Pulse/UR36020kHz_30B1_rf10kHz/pulse1000.bruker")
