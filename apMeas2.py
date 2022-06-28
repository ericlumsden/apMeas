import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf
import os
'''
New order of operations:
1. Find peak
2. Find takeoff point
3. Find first spot after peak that is same as takeoff
4. Find next spot after that that is the same as or greater than takeoff
5. Find half-width
6. Define window as -0.025 to that return to takeoff after AP
7. Find AHP min and length
8. Rise and decay? Last priorities
'''

def runAPfind(trace, threshold=-25.0, freq=10000, buffer=0.025):
    buff = int(buffer * freq)
    peakIndices = []
    peakAmps = []
    start = 0
    for idx, x in enumerate(trace):
        if (start == 0) and (x > threshold):
            start = idx
        elif (start != 0) and (x <= threshold):
            peak_idx = np.argmax(trace[int(start):idx]) + start
            peak_amp = trace[int(peak_idx)]

            peakIndices.append(peak_idx)
            peakAmps.append(peak_amp)

            sub_trace = trace[int(peak_idx - buff):int(peak_idx + buff)]
            time = [x/freq for x in range(len(sub_trace))]
            plt.figure(1)
            plt.plot(time, sub_trace)
            plt.show()
            plt.clf()

            start = 0
        else:
            continue

    return ap_dict
'''
In the main apMeas script I will load a csv file with all of the info for the traces to be analyzed

Qs that remain:
- Will I save each file as its own dataclass?
- How do I want to deal with spont vs induced APs?
'''

# First thing is to get the working directory where the csv and abf files are located...
# ...user input at command line or in-file edit?
# For now I will be testing using file saved in the apMeas directory
working_directory = '/media/eric/HDD/uw/VGCC'
#working_directory = os.getcwd()
freq = 10000 # If we need to change the sampling frequency do so here, and set it in the called functions
'''
Steps:
1. Load in the csv file
2. Iterate over rows
    NOTE: In the future I would like to skip over files already processed but for now let's continue
    2a. If the sparrow number does not have a folder, create on and continue
        2ai. If the trace is 'gap-free' analyze from start of sweep to Cd exposure start and then Cd exposure to the end
        2aii. Save the output as both a dataclass(?) and a graph into the sparrow numbers folder as cell 1, also save a copy to the main 'data' folder for easy analysis
        2aii. If the trace is 'cc-steps' analyze each sweep individually, saving each sweep as a dataclass in a subfolder, again putting a copy into the 'data' folder
3. Once all of this is done I will write a new script for doing summary data and figures
'''
csv_file_name = 'cadmiumBlock'
#csv_file_name = 'test'
df = pd.read_csv(f'{working_directory}/{csv_file_name}.csv')
df_len = len(df.index)

for index, row in df.iterrows():
    try:
        os.mkdir(f"{working_directory}/{row['sparrow_num']}/")
    except OSError:
        pass
    abf = pyabf.ABF(f"{working_directory}/{row['trace_num']}.abf")
    ap_dict = runAPfind(abf.sweepY)