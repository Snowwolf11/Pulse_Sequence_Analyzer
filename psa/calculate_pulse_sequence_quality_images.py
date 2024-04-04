from psa.calculate_pulse_sequence_quality import *
from psa.getPulseSequence import *

import os
import numpy as np

def calculate_pulse_sequence_quality_images(GUI, dirname, Range, T, Umax, initialVector, calcType, changingVariable, Resolution):
    files = os.listdir(dirname)
    sorted_files = sorted(files)
    PSnames = [f for f in sorted_files if 'bruker' in f]
    PSnumber = len(PSnames)
    offsetRes = round(Range * Resolution/1000)
    offset = Range / offsetRes
    shallVec = np.array([0, 0, -1])  # Define shallVec as a numpy array

    print("step1")
    if changingVariable == 1:
        PI = np.ones((offsetRes, PSnumber))
        for n1 in range(PSnumber):
            PS = np.array(getPulseSequence(os.path.join(dirname, PSnames[n1])))
            GUI.text_area_set(text_area = GUI.info_error_text, text_str = f"{n1 + 1} von {PSnumber}", reset_bool = 0)
            for n2 in range(offsetRes):
                Q = calculate_pulse_sequence_quality(PS, T, 1, Umax + (-Range / 2 + (n2 - 1) * offset), 0, 1, np.array(initialVector), shallVec, calcType)
                #print("Q: "+str(Q))
                if abs(1 - Q) < 10**(-15):
                    Q = 1-10**(-15)
                PI[n2, n1] = np.log10(abs(1 - Q))

    if changingVariable == 2:
        PI = np.ones((offsetRes, PSnumber))
        for n1 in range(PSnumber):
            PS = np.array(getPulseSequence(os.path.join(dirname, PSnames[n1])))
            GUI.text_area_set(text_area = GUI.info_error_text, text_str = f"{n1 + 1} von {PSnumber}", reset_bool = 0)
            for n2 in range(offsetRes):
                Q = calculate_pulse_sequence_quality(PS, T, 1, Umax, -Range / 2 + (n2 - 1) * offset, 1, np.array(initialVector), shallVec, calcType)
                if abs(1 - Q) < 10**(-15):
                    Q = 1-10**(-15)
                PI[n2, n1] = np.log10(abs(1 - Q))

    return PI
