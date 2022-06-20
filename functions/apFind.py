import numpy as np
import apISI
from subFuncs import *
'''

The function 'apFind' is used to find all action potentials in a given trace and run the subsequent functions for analyzing each action potential

Action potentials are identified based on a threshold crossing which has a default value of -20.0 mV but can be changed

For other AP properties, a 'sub_trace' is utilized, which takes the index of the AP peak and +/- a 'buffer' amount of time in the trace... 
    ...the buffer is pre-set to 10ms, found by multiplying 0.01 by the recording frequency (10,000 Hz for my IC recordings)...
    ...that total stretch is made into an array and then passed into subsequent functions

'''
def apFind(trace, threshold=-20.0, buff=0.01, freq=10000)):
    # Initiate empty arrays for all of the measurements
    apPeaks = []
    apRise = []
    apDecay = []
    apHalfWidth = []
    apAHPmin = []
    apAHPlen = []
    apInterSpikeIntervals = []

    # Set buffer
    buffer = int(buff * freq)

    start = 0 # Start of threshold crossing for each AP
    for idx, val in enumerate(trace):
        ap_count = len(apPeaks)
        if (val > threshold) and (start == 0):
            start = idx
        elif (val < threshold) and (start != 0):
            # Find the peak of the AP based on the start/end of threshold crossing...
            peak_idx = np.argmax[start:idx] + start
            apPeaks.append(peak_idx)

            # Find inter-spike interval by taking difference in the indices of this AP and the previous one, then dividing by sampling frequency
            if (ap_count > 1):
                apInterSpikeIntervals.append((peak_idx - apPeaks[ap_count]) / freq)

            # Get the subtrace array based on the peak and the buffer
            sub_trace = trace[int(peak_idx - buffer):int(peak_idx + buffer)]

            # Use the max of the derivative to find the AP takeoff point (this will be used in sub functions)
            ap_takeoff = apTakeoff(sub_trace, freq) + start # have to include `+ start` as the function just finds the max of the derivative of the array passed
            # Then use the takeoff point to find where the AP returns to baseline
            ap_return = apReturn(sub_trace, (ap_takeoff - start), freq) + start

            apRise.append((peak_idx - ap_takeoff) / freq)
            apDecay.append((ap_return - peak_idx) / freq)
            apHalfWidth.append(halfWidth(sub_trace, peak_idx))
            apAHPmin.append(ahpMin(sub_trace, peak_idx))
            apAHPlen.append(ahpLen(sub_trace, peak_idx))

            start = 0 # Need to reset for next AP
        else:
            continue

    return apPeaks, apInterSpikeIntervals, apRise, apDecay, apHalfWidth, apAHPmin, apAHPlen