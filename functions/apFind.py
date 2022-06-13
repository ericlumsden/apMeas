import numpy as np
import apISI
from subFuncs import *
'''
The function 'apFind' is used to find all action potentials in a given trace and run the subsequent functions for analyzing each action potential
Action potentials are identified based on a threshold crossing which has a default value of -20.0 mV but can be changed
'''
def apFind(trace, threshold=-20.0):
    # Initiate empty arrays for all of the measurements
    apPeaks = []
    apRise = []
    apDecay = []
    apHalfWidth = []
    apAHPmin = []
    apAHPlen = []

    start = 0 # Start of threshold crossing for each AP
    for idx, val in enumerate(trace, buffer):
        if (val < threshold) and (start == 0):
            continue
        elif (val > threshold) and (start == 0):
            start = idx
        elif (val < threshold) and (start != 0):
            peakIDX = np.argmax[start:idx] + start
            subTrace = trace[int(peakIDX-buffer):int(peakIDX+buffer)]
            apPeaks.append(peakIDX)
            apRise.append(rise(subTrace, peakIDX))
            apDecay.append(decay(subTrace, peakIDX))
            apHalfWidth.append(halfWidth(subTrace, peakIDX))
            apAHPmin.append(ahpMin(subTrace, peakIDX))
            apAHPlen.append(ahpLen(subTrace, peakIDX))
            start = 0 # Need to reset for next AP

        apInterSpikeIntervals = apISI(apPeaks)
        
    return apPeaks, apInterSpikeIntervals, apRise, apDecay, apHalfWidth, apAHPmin, apAHPlen