#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
import numpy as np

def createHelix(startingPoint, startingTangent, axisDir, radius, length, chirality, resolution, type):
            
            #length = 1
            
            # creates a numerical Helix (vor AHT-Curve -> Chirality inverted!!!) function with the parameters:
            #
            #       Input: (all inputs either ndarray or float)
            #       - radius: Radius of the Helix
            #       - startingPoint: starting point of the Helix
            #       - startingTangent: Tangent of the Helix at the startingPoint
            #       - length: arc length of the helix
            #       - resolution: Determines how many Coordinates the Helix consists
            #               off (resolution of helix)
            #       - chirality :+1 right-handed /-1 left-handed
            #       - axisDir: Vector parallel to the axis of the Helix
            #       - type: if(0): Normal, Tangent, ... isnt plotted
            #
            #       Output:
            #
            #       - Helix: Matrix, which contains the coordinates of the numerical Helix
            #               function as entries
            #       - endingPoint: Point at which the Helix Ends
            #       - endingNormal: normal Vector of the Helix at the ending Point
            #       - endingTangent: Tangent of tghe Helix at the endingPoint
            #       - startingNormal: Normal at beginning of Helix (Output for tests)
            #
            #       Formulas:
            #
            #       k=h/(2πr) (slope of the helix)
            #       s=2πr*sqrt(1+k^2)*t
            #       right handed: (chirality>1)
            #               a=cross(tp,n)   a:axis; tp:tangent projection; n:normal;
            #       left handed: (chirality<1)
            #               a=cross(n,tp)   a:axis; t:tangent projection; n:normal;
            #       ...(Wikipedia)
            #
            # sample function call:
            #       createHelix_V3([0;0;0],[0;1;1], [0;1;-1], 1, 10, 1, 1000, 1)
            
            # transform inpu into column-vectors
            #startingPoint = np.array([startingPoint[1]; startingPoint[2]; startingPoint[3]])
            #startingTangent = np.array([startingTangent[1]; startingTangent[2]; startingTangent[3]]);
            #axisDir = [axisDir(1); axisDir(2); axisDir(3)];
            
            # normalize Input
            startingTangent = startingTangent/np.linalg.norm(startingTangent)
            axisDir = axisDir/np.linalg.norm(axisDir)
            #calculate additional Parameters
            if chirality >= 0:
              #"right-handed"
              startingNormal = np.cross(startingTangent,axisDir)/np.linalg.norm(np.cross(startingTangent,axisDir))
              startingTangent_s = np.cross(axisDir,startingNormal)/np.linalg.norm(np.cross(axisDir, startingNormal))
            elif chirality < 0:
              #"left-handed"
              startingNormal = np.cross(axisDir, startingTangent)/np.linalg.norm(np.cross(axisDir, startingTangent))
              startingTangent_s = np.cross(startingNormal, axisDir)/np.linalg.norm(np.cross(startingNormal, axisDir))
            #slope = abs(sqrt(1-(dot(startingTangent, startingTangent_s))^2)/dot(startingTangent, startingTangent_s)); 
            #right???
            ###############################
            innerpr = np.dot(startingTangent, startingTangent_s)
            if abs(1-innerpr)<10e-15:
              innerpr = 1
            if abs(-1-innerpr)<10e-15:
              innerpr = -1
          
            """
            #{
            if(abs(innerpr)>1)
              if(abs(innerpr)<(1+10^-9))
                if(innerpr<0)
                  innerpr = -1;
                elseif(innerpr>1)
                  innerpr = 1;
                end
              else
                innerpr ##ok<NOPRT> 
                error('abs(dot(startingTangent, startingTangent_s))>1!')
              end
            end
            #}
        
            """
            
            slope =  np.tan(np.arccos(innerpr)) #sqrt(1-innerpr^2)/innerpr PROBLEM
            if isinstance(slope, complex):
                raise Exception('slope has imaginary value!')
            ###############################
            height = 2*np.pi*radius*slope
            T = length/(2*np.pi*radius*np.sqrt(1+slope**2)) #tmax of the Helix (t at ending point)
            t = np.linspace(0,T,resolution)
            
            #calculate transormation matrix and norm helix in z-Direction
            S = np.array([startingNormal[:], startingTangent_s[:], axisDir[:]]).transpose()
            if (1-np.linalg.norm(np.linalg.det(S))) > (10e9):
              raise Exception('det(S) ≠ ±1')
            helix_z =  np.array([radius*np.cos(2*np.pi*t), radius*np.sin(2*np.pi*t), height*t])
            
            #calculate helix
            helix = helix_z
            helix_0 = np.dot(S,helix[:,0])
            for n in range(len(t)):
              helix[:,n] = np.dot(S,helix_z[:, n]) + startingPoint - helix_0
            #[helix[:,n] = S*helix_z[:, n] + startingPoint - helix_0  for n,v in enumerate t]
            #generate Output
            endingPoint = helix[:, -1]
            endingTangent_num = helix[:,-1] - helix[:,-2]
            endingTangent_num = endingTangent_num/np.linalg.norm(endingTangent_num)
            T_end = length/(2*np.pi*radius*np.sqrt(1+(height/(2*np.pi*radius))**2))
            endingTangent = np.dot(S,np.array([-2*np.pi*radius*np.sin(2*np.pi*T_end), 2*np.pi*radius*np.cos(2*np.pi*T_end), height]).transpose())      #TODO
            endingTangent = endingTangent/np.linalg.norm(endingTangent)
            if chirality >= 0:
              endingNormal = np.cross(endingTangent,axisDir)/np.linalg.norm(np.cross(endingTangent,axisDir))
            elif(chirality < 0):
              endingNormal = np.cross(axisDir, endingTangent)/np.linalg.norm(np.cross(axisDir, endingTangent))
            #penultimatePoint = helix(:, end-1)
            return helix, endingPoint, endingTangent, endingNormal, startingNormal
        
        