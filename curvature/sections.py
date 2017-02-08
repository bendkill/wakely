"""
Script for processing steps from a given Nx2 csv file. Relevant functions:
parse: returns an array with the csv data.
map_curvatures: returns the list of points at which curvature spikes
    note that these points must be separated by at least 1 degree. This needs fixing.
plot: plot the csv data on a semilog scale with the sections added as vertical lines.

"""

import numpy as np
import math
import matplotlib.pyplot as plt

norm     = np.linalg.norm
subtract = np.subtract
sin      = math.sin
acos     = math.acos
log      = math.log

def parse(file_name):
    """for parsing a (2xN) .csv file into an (Nx2) numpy array."""
    lines = open(file_name).read().split('\n')
    temps = [float(point) for point in lines[0].split(',')]
    volts = [float(point) for point in lines[1].split(',')]
    return np.array( zip(temps,volts) )

def curvature(arr):
    """given a 3x2 array, [[x1,x2], [y1,y2], [z1,z2]] finds the curvature
    at the point y."""
    X  = norm(np.subtract(arr[1], arr[2]))
    Y  = norm(np.subtract(arr[0], arr[2]))
    Z  = norm(np.subtract(arr[0], arr[1]))
    th = (X**2 - Y**2 + Z**2) / (2*Z*X)
    return 2 * sin(acos(th)) / Y

def map_logs(arr):
    """maps log over the second elements in each point of arr,
    and Nx2 array"""
    return np.array([[point[0], log(point[1])] for point in arr])

def map_curvature(arr):
    """given an Nx2 array, map curvature over elements 1 to N-1, return
    curvs, an Nx2 array of xs and curvature, with 0 at the bounds."""
    curvs       = np.zeros_like(arr)
    curvs[0,0]  = arr[0,0]
    curvs[-1,0] = arr[-1,0]
    for i in range(1, len(arr)-1):
        curvs[i] = [arr[i,0], curvature(arr[i-1:i+2])]
    return curvs

def is_peak(arr):
    """A point is a peak if its y value is greater than the points to
    either side.
    So given a 5x2 np array of points, return true if the middle point
    is a peak and False otherwise."""
    if arr[0,1] < arr[1,1] < arr[2,1]:
        return True
    else:
        return False

def is_min(arr):
    if arr[0,1] > arr[1,1] > arr[2,1]:
        return True
    else:
        return False

def get_temp_with_max(temps_curvs):
    """takes an nx2 list with all the temperatures within 1 deg of each other.
    returns the temperature at which curvature is highest in that group"""
    curv = 0.
    for i in range(len(temps_curvs)):
        if temps_curvs[i][1] > curv:
            curv = temps_curvs[i][1]
            peak = temps_curvs[i][0]
    return peak
    
def filter_peaks(temps_curvs):
    """index through each point in curvs_temps. If two temps are within a degree
    acquire the group of temps thereafter, pick the max from that group."""
    peaks = []
    i = 0
    while i < len(temps_curvs)-2:
        dx = int(temps_curvs[i+1][0]) - int(temps_curvs[i][0])
        if dx > 1.:
            i += 1
            continue
        maybe_peaks = [[temps_curvs[i][0], temps_curvs[i][1]]]
        while dx <= 1. and i < len(temps_curvs)-2:
            maybe_peaks.append([temps_curvs[i+1][0], temps_curvs[i+1][1]])
            i += 1
            dx = int(temps_curvs[i+1][0]) - int(temps_curvs[i][0])
        peaks.append( get_temp_with_max(maybe_peaks) )
    return peaks

def lookup(temp,arr):
    """looks for temp in np arr (Nx2). If exists, return the second
    value at that point in arr. If doesn't exist, raise an error."""
    for i in range(len(arr)):
        if arr[i,0] == temp:
            return arr[i,1]
    print "error: temp not found"
    quit()

def filter_peaks_stat(peaks,arr):
    i = 0
    while i < len(peaks):
        if lookup(peaks[i],arr) < 10:
            del peaks[i]
            continue
        i += 1
    return peaks

def get_peaks(arr):
    """Given an Nx2 array (of temps and curvs), figure out if
    each value is a peak.
    If it is a peak, add its temp to a list of temps at which a
    curvature peak occurs.
    Filter that to just a list of temps with max or min peaks."""
    curvs_temps = []
    for i in range(1,len(arr)-1):
        if is_peak(arr[i-1:i+2]) or is_min(arr[i-1:i+2]):
            curvs_temps.append([arr[i,0], arr[i,1]])
    return filter_peaks(curvs_temps)

def peaks(fileName):
    """ takes a csv filename and returns an np array of temps for which changes
    occur"""
    data = parse(fileName)
    return filter_peaks_stat(get_peaks(map_curvature(data)), data)

def plot(fileName):
    """Take a 2xN csv file, parse it and plot the data.
    Use map_curvature and peaks to find the sections, then
    overlay those onto the plot."""
    data = parse(fileName)
    plt.semilogy([point[0] for point in data],
                 [point[1] for point in data])
    peaks_ = peaks(fileName)
    for peak in peaks_:
        plt.axvline(peak)
