from functions.apFind import *
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf
import os

'''
In the main apMeas script I will load a csv file with all of the info for the traces to be analyzed

Qs that remain:
- Will I save each file as its own dataclass?
- How do I want to deal with spont vs induced APs?
'''

# First thing is to get the working directory where the csv and abf files are located...
# ...user input at command line or in-file edit?
# For now I will be testing using file saved in the apMeas directory
working_directory = os.getcwd()
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

df = pd.read_csv(f'{working_directory}/test.csv')
df_len = len(df.index)

for index, row in df.iterrows():
    os.mkdir(f"{working_directory}/{row['sparrow_num']}/")
    abf = pyabf.ABF(f"{working_directory}/{row['trace_num']}.abf")

    if (row['trace_type'] == "'gap-free'"):
        if row['cd_exp'] == True:
            for exposures in ['pre', 'post']:
                if exposures == 'pre':
                    trace = abf.sweepY[:int(row['cd_exp_start']*freq)]
                else:
                    trace = abf.sweepY[int(row['cd_exp_end']):]
                ap_dict = runAPfind(trace, working_directory, row['sparrow_num'], row['cell_num'], exposure=exposures)
                json.dump(ap_dict, open(f"{working_directory}/{row['sparrow_num']}/{row['cell_num']}_{exposures}.json", 'w'))
        else:
            trace = abf.sweepY
            ap_dict = runAPfind(trace, working_directory, row['sparrow_num'], row['cell_num'])
            json.dump(ap_dict, open(f"{working_directory}/{row['sparrow_num']}/{row['cell_num']}.json", 'w'))
    
    print(f'Finished row {index} of {df_len}')