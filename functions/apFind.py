import numpy as np
from functions.subFuncs import *
import matplotlib.pyplot as plt
'''

The function 'apFind' is used to find all action potentials in a given trace and run the subsequent functions for analyzing each action potential

Action potentials are identified based on a threshold crossing which has a default value of -20.0 mV but can be changed

For other AP properties, a 'sub_trace' is utilized, which takes the index of the AP peak and +/- a 'buffer' amount of time in the trace... 
    ...the buffer is pre-set to 20ms, found by multiplying 0.02 by the recording frequency (10,000 Hz for my IC recordings)...
    ...that total stretch is made into an array and then passed into subsequent functions

'''
def apFind(trace, threshold=-20.0, buff=0.02, freq=10000):
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
        if (val > threshold) and (start == 0):
            start = idx
        elif (val < threshold) and (start != 0):
            # Find the peak of the AP based on the start/end of threshold crossing...
            peak_idx = np.argmax(trace[start:idx])

            # Find inter-spike interval by taking difference in the indices of this AP and the previous one, then dividing by sampling frequency
            if (len(apPeaks) > 1):
                apInterSpikeIntervals.append(((peak_idx + star) - apPeaks[-1]) / freq)

            apPeaks.append(peak_idx + start)
            # Get the subtrace array based on the peak and the buffer
            sub_trace = trace[int(peak_idx - buffer):int(peak_idx + buffer)]
            time = [x/freq for x in range(len(ap_array))]
            ahp_min = np.argmin(sub_trace)

            # Use the max of the derivative to find the AP takeoff point (this will be used in sub functions)
            ap_takeoff = apTakeoff(sub_trace, time) 
            half_height = (sub_trace[ap_takeoff] - sub_trace[peak_idx]) / 2
            # Then use the takeoff point to find where the AP returns to baseline
            ap_return = apReturn(sub_trace, peak_idx, ap_takeoff, ahp_min) 

            # Rise and decay are easily calculated by taking difference between peak and takeoff or return point, then dividing by sampling frequency
            apRise.append((peak_idx - ap_takeoff) / freq)
            apDecay.append((ap_return - peak_idx) / freq)

            # Half width, and after-hyperpolarization len/min are found using sub functions
            apHalfWidth.append(findHW(sub_trace, time, half_height, peak_idx, ap_takeoff, ap_return, freq))

            # Once all of the necessary components are discovered, the AHP measurements are very easy to collect, no function necessary
            apAHPmin.append(sub_trace[ap_return] - sub_trace[ahp_min])
            apAHPlen.append((ap_retun - ahp_min) / freq)

            start = 0 # Need to reset for next AP
        else:
            continue

    return apPeaks, apInterSpikeIntervals, apRise, apDecay, apHalfWidth, apAHPmin, apAHPlen

def runAPfind(trace, working_directory, sparrow_num, cell_num, exposure=False, freq=10000):
    # If there is an exposure we have to note that in the naming of the figure
    if exposure != False:
        cell_name = f'{cell_num}_{exposure}'
    else:
        cell_name = cell_num

    # Run the apFind function (above) with the trace as designated in apMeas
    apPeaks, apInterSpikeIntervals, apRise, apDecay, apHalfWidth, apAHPmin, apAHPlen = apFind(trace)
    # Make a dictionary of these values, to be returned to apMeas for saving as a json file
    ap_dict = {
        'apPeaks': [int(x) for x in apPeaks],
        'apISIs': [float(x) for x in apInterSpikeIntervals],
        'apRise': [float(x) for x in apRise],
        'apDecay': [float(x) for x in apDecay],
        'apHW': [float(x) for x in apHalfWidth],
        'apAHPmin': [float(x) for x in apAHPmin],
        'apAHPlen': [float(x) for x in apAHPlen]
    }

    # We will use the AP peak indices as our time points for the measured values
    apPeaks_tp = [x/freq for x in apPeaks]
    # Create a figure with 7 rows for the different measured values as well as a general trace plot
    fig, axs = plt.subplots(7,1)
    axs[0].scatter(apPeaks_tp[1:], apInterSpikeIntervals)
    axs[1].scatter(apPeaks_tp, apRise)
    axs[2].scatter(apPeaks_tp, apDecay)
    axs[3].scatter(apPeaks_tp, apHalfWidth)
    axs[4].scatter(apPeaks_tp, apAHPmin)
    axs[5].scatter(apPeaks_tp, apAHPlen)
    axs[6].plot([x/10000 for x in range(len(trace))], trace)

    # Save this figure in the sparrow # directory with the cell name as designated above
    plt.savefig(f'{working_directory}/{sparrow_num}/{cell_name}')

    # Return the dictionary to be saved as a json file in apMeas
    return ap_dict