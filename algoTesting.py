import pyabf
import numpy as np
import matplotlib.pyplot as plt

freq = 10000
abf = pyabf.ABF('./22323002.abf')
trace = abf.sweepY[int(127.9*freq):int(128.3*freq)]
time = [x/freq for x in range(len(trace))]
peak = np.argmax(trace)
ahp_peak = np.argmin(trace)

def findTakeoff(trace, time):
    # Returns the max value of the 3rd derivative of the trace
    return np.argmax(np.diff(np.diff(np.diff(trace) / np.diff(time)) / np.diff(time[:-1])) / np.diff(time[:-2]))

takeoff = findTakeoff(trace, time)
half_height = ((abs(trace[int(takeoff)]) + abs(trace[int(peak)])) / 2) + trace[int(takeoff)]

#3019, 3020
baseline_return = np.interp(trace[int(takeoff)], trace[int(peak):int(ahp_peak)][::-1], time[int(peak):int(ahp_peak)][::-1])
print(baseline_return)

def findHW(y, x, half_height, peak, takeoff, baseline_return):
    up = np.interp(half_height, y[int(takeoff):int(peak)], x[int(takeoff):int(peak)])
    down = np.interp(half_height, y[int(peak):int(baseline_return*freq)][::-1], x[int(peak):int(baseline_return*freq)][::-1])

    return [up, down, ((down - up) / freq)]

hw = findHW(trace, time, half_height, peak, takeoff, baseline_return)
print(hw)

plt.plot(time, trace)
plt.scatter(hw[0], half_height, marker='*', color='r')
plt.scatter(hw[1], half_height, marker='*', color='k')
plt.scatter(time[int(takeoff)], trace[int(takeoff)], marker='x', color='r')
plt.scatter(time[int(ahp_peak)], trace[int(ahp_peak)], marker='x', color='g')
plt.scatter(baseline_return, trace[int(takeoff)], marker='x', color='k')
plt.axhline(y=half_height, color='grey', linestyle='--')
plt.axhline(y=np.max(trace), color='grey', linestyle='--')
plt.axhline(y=trace[takeoff], color='grey', linestyle='--')
plt.xlim([0.295,0.315])
plt.savefig('algoTesting2.png')