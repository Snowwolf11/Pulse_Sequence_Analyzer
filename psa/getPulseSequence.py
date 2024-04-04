#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leon
"""
#TODO: line 26
import re
import numpy as np

def getPulseSequence(PulseSequence):
    PS = Import_Bruker(PulseSequence)
    return PS

def Import_Bruker(filename):
    if isinstance(filename, str):
        with open(filename, "r") as f:
            text = f.read()
        text = re.sub(chr(13), chr(10), text)       
        text = re.sub(chr(10)+'{2}', chr(10), text)
        # Remove comments:
        #   char(10),
        #   then #,
        #   then any number of characters except char(10),
        #   then char(10)
        text = re.sub('\x0A#[^\x0A]*(?=\x0A)','', chr(10)+text+chr(10))
        PSlines = text.splitlines()
        while '' in PSlines:
            PSlines.remove('') 
        PS = np.empty((len(PSlines),2))
        for ind,line in enumerate(PSlines):
            PS[ind,:] = np.fromstring(PSlines[ind], dtype=float, sep=",")            #TODO: not gonna work??
        return np.array(PS)
    elif isinstance(filename, np.ndarray):
        PS = filename
        return PS
    else:
        raise Exception('unknown argument type for filename.')
        
