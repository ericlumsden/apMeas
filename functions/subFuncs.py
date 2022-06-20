import numpy as np
'''
Rather than constantly iterating through each trace, I will find AP parameters when the peaks are discovered
Here I will define the subfunctions called in apFind at the time of AP peak discovery
'''

# To find the AP takeoff point, return the maximum value of the first derivative of the AP
def apTakeoff(ap_array, freq):
    dVdt = np.diff(ap_array) / np.diff(len(ap_array)/freq)
    return np.argmax(dVdt) 

# apReturn finds the index at which the AP returns to the same value as the takeoff point (for finding halfwidth and decay)
def apReturn(ap_array, takeoff, freq):
    target = ap_array[int(takeoff)]
    
    for idx, x in enumerate(ap_array[int(takeoff):]):
        if x > target:
            continue
        else:
            bl_return = idx
            break
    
    return bl_return

def halfWidth():
    return hw

def ahpMin():
    return ahpM

def ahpLen():
    return ahpL