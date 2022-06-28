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
            peak_i = np.argmax(sub_trace)
            time = [x/freq for x in range(len(sub_trace))]

            takeoff = np.argmax(np.diff(np.diff(np.diff(sub_trace[:peak_i]) / np.diff(time[:peak_i])) / np.diff(time[:int(peak_i-1)])) / np.diff(time[:int(peak_i - 2)]))

            half_height = ((sub_trace[peak_i] - sub_trace[takeoff]) / 2) + sub_trace[takeoff]
            window_min = np.argmin(sub_trace[peak_i:]) + peak_i
            post_peak_to = (np.interp(sub_trace[int(takeoff)], sub_trace[int(peak_i):int(window_min)][::-1], time[int(peak_i):int(window_min)][::-1])) * freq

            window_end = 0
            for idx, val in enumerate(trace[int(post_peak_to + peak_idx)+1:]):
                if val > sub_trace[int(takeoff)]:
                    window_end = idx + (post_peak_to) + 1
                    break
                else:
                    continue

            print(window_end, peak_idx)
            sub_trace = trace[int(peak_idx - buff):int(window_end + peak_idx)]
            time = [x/freq for x in range(len(sub_trace))]
            up = (np.interp(half_height, sub_trace[int(takeoff):int(peak_i)], time[int(takeoff):int(peak_i)]))
            down = (np.interp(half_height, sub_trace[int(peak_i):int(post_peak_to)+1][::-1], time[int(peak_i):int(post_peak_to)+1][::-1]))
            print(up, down)

            half_width = (down - up) / freq
            ahp_min = np.argmin(sub_trace[peak_i:]) + peak_i

            plt.figure(1)
            plt.plot(time, sub_trace)
            plt.scatter(time[peak_i], sub_trace[peak_i], marker='x', color='r')
            plt.scatter(up, half_height, marker='*', color='g')
            plt.scatter(down, half_height, marker='*', color='g')
            plt.scatter(time[ahp_min], sub_trace[ahp_min], marker='x', color='r')
            plt.scatter(time[takeoff], sub_trace[takeoff], marker='x', color='k')
            plt.scatter(post_peak_to/freq, sub_trace[takeoff], marker='x', color='k')
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