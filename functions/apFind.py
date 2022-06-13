import numpy as np
'''
The function 'apFind' is used to find all action potentials in a given trace
Action potentials are identified based on a threshold crossing which has a default value of -20.0 mV but can be changed
This function returns an array of indices for action potential peaks
'''
def apFind(trace, threshold=-20.0):
    apPeaks = [] # Initiate an empty array for the action potential peaks
    start = 0 # Start of threshold crossing for each AP
    for idx, val in enumerate(trace):
        if (val < threshold) and (start == 0):
            continue
        elif (val > threshold) and (start == 0):
            start = idx
        elif (val < threshold) and (start != 0):
            apPeaks.append(np.argmax[start:idx])
            start = 0 # Need to reset for next AP
        
    return apPeaks