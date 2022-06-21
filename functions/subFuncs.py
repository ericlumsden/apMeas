import numpy as np
'''
Rather than constantly iterating through each trace, I will find AP parameters when the peaks are discovered
Here I will define the subfunctions called in apFind at the time of AP peak discovery
'''

# To find the AP takeoff point, return the maximum value of the first derivative of the AP
def apTakeoff(ap_array, time):
    return np.argmax(np.diff(np.diff(np.diff(ap_array) / np.diff(time)) / np.diff(time[:-1])) / np.diff(time[:-1]))

# apReturn finds the index at which the AP returns to the same value as the takeoff point (for finding halfwidth and decay)
def apReturn(ap_array, peak, takeoff, ahp_min):
    return np.interp(ap_array[int(takeoff)], ap_array[int(peak):int(ahp_min)][::-1], time[int(peak):int(ahp_min)][::-1])

# I'm using baseline_return here in findHW because this function was initially written as a test function... COME BACK AND CHANGE THIS ONCE IT'S WORKING
def findHW(ap_array, time, half_height, peak, takeoff, baseline_return, freq):
    up = np.interp(half_height, ap_array[int(takeoff):int(peak)], time[int(takeoff):int(peak)])
    down = np.interp(half_height, ap_array[int(peak):int(baseline_return*freq)][::-1], time[int(peak):int(baseline_return*freq)][::-1])

    return ((down - up) / freq)