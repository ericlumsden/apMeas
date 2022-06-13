import numpy as np
'''
This function takes an array of spike indices (apSpikes) and returns an array of the inter-spike intervals (in s)
Sampling frequency (freq) is set to a value of 10kHz (my sampling frequency for current clamp recordings) but can be changed
'''
def apISI(apSpikes, freq=10000):
    isi_array = []
    for idx, val in enumerate(apSpikes):
        if (idx == 0) or (idx == len(apSpikes)-1): # Skip first and last values, we're looking at differences
            continue
        else:
            isi_array.append((val - apSpikes[idx-1]) / freq) # ISI is difference between the current value and previous value divided by the sampling frequency
        
    return isi_array