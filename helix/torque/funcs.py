""" general functions for use in the torque module"""

import math

def to_one_eighty(x):
    x %= 360
    if x <= 180: return x
    return -360 + x
