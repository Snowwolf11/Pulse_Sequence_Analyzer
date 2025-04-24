#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
from psa.createCurve import *
import numpy as np
from matplotlib import pyplot as plt

#PS = "/Users/leon/Desktop/Physik/Glaser/Analyse_und_Visualisierung_von_robusten_Kontrollpulsen/Pulssequenzen/UR_Pulse/UR36020kHz_30B1_rf10kHz/pulse1400.bruker"
PS = "/Users/leon/Desktop/Physik/Glaser/Analyse_und_Visualisierung_von_robusten_Kontrollpulsen/Pulssequenzen/linearer_Phasengang/0-1440/pulse2000"
T=0.0000005
l=1
maximumAmplitude=10000
offset=0
inpoFact=1
xExpand=0
initialVector=np.array([0,0,1])

CM_Matrix = createCurve(PS,T,l,maximumAmplitude,offset, inpoFact,xExpand, 1, initialVector)
CM_Helix = createCurve(PS,T,l,maximumAmplitude,offset, inpoFact,xExpand, 2, initialVector)
CM_diff = CM_Matrix[0:-1]-CM_Helix
CM_Matrix = createCurve(PS,T,l,maximumAmplitude,offset, inpoFact,xExpand, 1, initialVector)
plot3DCurve(CM_Helix[:,0], CM_Helix[:,1], CM_Helix[:,2], type="old")
plot3DCurve(CM_diff[:,0], CM_diff[:,1], CM_diff[:,2], type="old")

CM_diffnorm = [np.linalg.norm(val) for val in CM_diff]

plt.figure()
plt.plot(np.linspace(0, np.shape(CM_diffnorm)[0]-1, np.shape(CM_diffnorm)[0]), CM_diffnorm)
plt.title("diffNorm")
plt.show()
