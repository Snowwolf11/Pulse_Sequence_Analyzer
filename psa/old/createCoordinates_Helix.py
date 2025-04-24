#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""

from psa.old.createHelix import *
import numpy as np

def createCoordinates_Helix(PS,T,l,Umax,offset,inpoFact,initialVector):
           
    #all inputs eather floats or ndarrays
    #
    #TODO
    #checked all parameters --> calculaction correct --> mistake
    #probably in createHelix_app or numerical
    
    ###Zum debuggen###
    Distortion = np.ones((np.shape(PS)[0]+1))
    ##################
    CM = np.ones((np.shape(PS)[0]+1,3))
    startingPoint = np.array([0,0,0])
    startingTangent = initialVector
    CM[0,:] = startingPoint;
    resolution = 2                           #TODO: einstellbar machen
    testnum = 0                               #-0.08;# bei 3000 -0.22;#bei 800 #0.6345 bei 190, 0.5113 bei 3000
    #testnum2 = 0
    
    endingPoint = initialVector              #nur zum Vordefinieren, Wert hat noch keine Bedeutung
    endingTangent = initialVector            #nur zum Vordefinieren, Wert hat noch keine Bedeutung
            
    for ind, val in enumerate(PS):
        # calculate Parameters
        Om = Umax*inpoFact/100*val[0]           #B-field strength of pulse
        if Om == 0:
            Om=10e-15
        Om_s = np.sqrt(offset**2 + Om**2)              #effective B-field
        phi = val[1]*(np.pi/180)
        B = np.array([Om*np.cos(phi), Om*np.sin(phi), offset]) 
        n = B/Om_s

        #calc of radius
        off_eff = np.dot(startingTangent,B)
        O_eff = np.linalg.norm(np.cross(startingTangent,B))
        radius = O_eff/(Om_s**2) #=sin(atan(O_s/off_s))/w              
        #d = 2*R*np.cross(startingTangent,n)/np.linalg.norm(np.cross(startingTangent,n))    #diameter vector


        #radius = np.sin(np.arctan(Om/abs(offset)))/Om_s                      #radius doesnt have an influence on curve because its proportional to height
        if not np.isreal(radius):
            raise ValueError("radius n = %s has imaginary value!" %(ind+1))
              
        if ind > 1:
            startingPoint = endingPoint
            #startingPoint = penultimatePoint
            startingTangent = endingTangent
        if np.linalg.norm(B) != 0:
            if np.arccos(np.dot(B/np.linalg.norm(B), startingTangent/np.linalg.norm(startingTangent))) <= np.pi/2:
                chirality = 1
            elif np.arccos(np.dot(B/np.linalg.norm(B), startingTangent/np.linalg.norm(startingTangent))) > np.pi/2:
                chirality = -1
        elif np.linalg.norm(B) == 0:
            chirality = 1
        print(chirality)
        axisDir = chirality*B
        length = 2*np.pi*Om_s*radius*np.sqrt(1+(off_eff/O_eff)**2)*T/inpoFact
        #length = length*(resolution+testnum)/resolution
        #print(length)
        #CM(ind+1,:) = endingPoint;
        helix, endingPoint, endingTangent, endingNormal_new, startingNormal_new = createHelix(startingPoint, startingTangent, axisDir, radius, length, -1*chirality, resolution, 0)     #-chirality weil wir die rotation des frames, nicht des Vektors betrachten
        
        ###Zum debuggen###
        Distortion[ind+1] = l/np.linalg.norm(endingPoint-startingPoint);
        ##################
        endingPoint = startingPoint - l*(startingPoint-endingPoint)/np.linalg.norm(endingPoint-startingPoint)     #TODO: sketchy
        endingPoint = endingPoint.transpose()
        CM[ind+1,:] = endingPoint
        if not np.isreal(np.sum(endingPoint)):
            raise ValueError("endingPoint n = %s has imaginary value!" %(ind+1))
    return CM
       

