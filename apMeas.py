from functions.apFind import *
import matplotlib.pyplot as plt
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

df = pd.read_csv(f'{working_directory}/test.csv')
abf = pyabf.ABF(f'{working_directory}/{df.trace_num.item()}.abf')

apPeaks, apInterSpikeIntervals, apRise, apDecay, apHalfWidth, apAHPmin, apAHPlen = apFind(abf.sweepY)
apPeaks_tp = [x/10000 for x in apPeaks]
print(len(apPeaks), len(apInterSpikeIntervals), len(apRise), len(apDecay), len(apHalfWidth), len(apAHPmin), len(apAHPlen))
fig, axs = plt.subplots(7,1)
axs[0].scatter(apPeaks_tp[1:-1], apInterSpikeIntervals)
axs[1].scatter(apPeaks_tp, apRise)
axs[2].scatter(apPeaks_tp, apDecay)
axs[3].scatter(apPeaks_tp, apHalfWidth)
axs[4].scatter(apPeaks_tp, apAHPmin)
axs[5].scatter(apPeaks_tp, apAHPlen)
axs[6].plot(abf.sweepX, abf.sweepY)

plt.savefig('./test.png')