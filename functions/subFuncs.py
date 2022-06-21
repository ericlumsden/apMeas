import numpy as np
'''
Rather than constantly iterating through each trace, I will find AP parameters when the peaks are discovered
Here I will define the subfunctions called in apFind at the time of AP peak discovery
'''

# To find the AP takeoff point, return the maximum value of the first derivative of the AP
def apTakeoff(ap_array, freq):
    dVdt = np.diff(ap_array) / np.diff([x/freq for x in range(len(ap_array))])
    return np.argmax(dVdt)-1 

# apReturn finds the index at which the AP returns to the same value as the takeoff point (for finding halfwidth and decay)
def apReturn(ap_array, takeoff, freq):
    target = ap_array[int(takeoff)]
    bl_return = 0
    
    for idx, x in enumerate(ap_array[int(takeoff):]):
        if x > target:
            continue
        else:
            bl_return = idx
            break
    
    return bl_return

def halfWidth(ap_array, peak, takeoff, freq):
    # Find half height based on the voltage halfway between takeoff and peak
    half_height = ((ap_array[int(peak)] - ap_array[int(takeoff)]) / 2) + ap_array[int(takeoff)]
    half_rise = 0
    half_decay = 0

    # Find the points at which rise and decay cross half_height...
    for idx, x in enumerate(ap_array[:int(peak)]):
        if x > half_height:
            half_rise = idx
            break
        else:
            continue
    
    for idx, x in enumerate(ap_array[int(peak):]):
        if x < half_height:
            half_decay = idx
            break
        else:
            continue
    # Return the difference in these two points (divided by sampling frequency, to return in s)
    return ((half_decay - half_rise) / freq)

# All of the AHP measurements have been collapsed into this one function, as they primarily 
def ahpMeas(ap_array, ap_return, freq):
    ahp_min = np.min(ap_array)
    ahp_min_idx = np.argmin(ap_array)

    target = ap_array[int(ap_return)]
    ahp_end = 0
    for idx, x in enumerate(ap_array[int(ap_return):]):
        if x < target:
            continue
        elif x > target:
            ahp_end = idx
            break

    ahp_len = (ahp_end - ap_return) / freq

    return ahp_min, ahp_min_idx, ahp_len