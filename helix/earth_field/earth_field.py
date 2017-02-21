"""
Given the flight data file, find the Magnitude and direction of the
Earth field at each point in the data.
"""

import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np
import string
from time import time

flight_data_help = np.fromfile("flight_data.txt", dtype=float, sep=" ")
flight_data = np.zeros((len(flight_data_help) / 5, 5))

for i in range(0, len(flight_data)):
    flight_data[i,0] = flight_data_help[5*i+0]
    flight_data[i,1] = flight_data_help[5*i+1]
    flight_data[i,2] = flight_data_help[5*i+2]
    flight_data[i,3] = flight_data_help[5*i+3]
    flight_data[i,4] = flight_data_help[5*i+4]

"""
Need to convert from mission days to decimal calendar year, and then
find a way to feed the latitude, longitude, altitude, and decimal year
into the wmm_point.exe file and read the output. Then we somehow have
to parse this crap (below) into just the magnitude of the field F, as
well as its "Decl and Incl" ideally in an Nx3 array.


 Welcome to the World Magnetic Model (WMM) 2015 C-Program

              --- Model Release Date: 15 Dec 2014 ---
            --- Software Release Date: 21 Nov 2014 ---


 This program estimates the strength and direction of
 Earth's main Magnetic field for a given point/area.
 Enter h for help and contact information or c to continue.
 >

Please enter latitude
North Latitude positive, For example:
30, 30, 30 (D,M,S) or 30.508 (Decimal Degrees) (both are north)

Please enter longitude
East longitude positive, West negative.  For example:
-100.5 or -100, 30, 0 for 100.5 degrees west

Please enter height above mean sea level (in kilometers):
[For height above WGS-84 Ellipsoid prefix E, for example (E20.1)]

Please enter the decimal year or calendar date
 (YYYY.yyy, MM DD YYYY or MM/DD/YYYY):
2016.000

 Results For

Latitude    0.00N
Longitude    0.00E
Altitude:    0.00 Kilometers above mean sea level
Date:        2016.0

        Main Field            Secular Change
F    =      31901.4 +/- 152.0 nT         Fdot =  37.0    nT/yr
H    =      27665.5 +/- 133.0 nT         Hdot =  -2.0    nT/yr
X    =      27545.8 +/- 138.0 nT         Xdot =   3.6    nT/yr
Y    =      -2571.0 +/-  89.0 nT         Ydot =  59.9    nT/yr
Z    =     -15884.5 +/- 165.0 nT         Zdot = -77.8    nT/yr
Decl    =      -5 Deg -20 Min  (WEST) +/- 19 Min Ddot = 7.4    Min/yr
Incl    =     -29 Deg -52 Min  (UP)   +/- 13 Min Idot = -7.4    Min/yr

"""

# flight_data goes like:
# mission day, latitude, longitude, altitude, something

# field_data goes like: dec_year, latitude, longitude, altitude, decl, incl, magnitude
field_data = np.zeros((len(flight_data), 7))

def process_outputs(handle):
    # takes the output of the wmm_point.exe
    # returns len 3 np array of the desired info
    data = np.zeros(3)

    handle.seek(0)
    lines = handle.read().split('\n')
    for i in range(len(lines)):
        words = string.split(lines[-i])
        if len(words) < 4: continue
        if words[0] == "Incl":
            degs = float(words[2])
            mins = float(words[4])
            data[1] = degs + mins / 60.
        if words[0] == "Decl":
            degs = float(words[2])
            mins = float(words[4])
            data[0] = degs + mins / 60.
        if words[0] == "F":
            nTs = float(words[2])
            data[2] = nTs / 10**9
            break

    return data

data_string = ""

for i in range(len(flight_data)):

    # modify to depend on a start date + mission days
    dec_year = round(2018. + 359./365. + flight_data[i,0] / 365., 3)
    # print "i: %i, dec_year: %f" % (i, dec_year)
    field_data[i,0] = flight_data[i,0]
    field_data[i,1:4] = flight_data[i,1:4]

    inputs_string = "c\n" +\
                    str(flight_data[i,1]) + "\n" +\
                    str(flight_data[i,2]) + "\n" +\
                    str(flight_data[i,3]) + "\n" +\
                    str(dec_year) + "\n"

    inputs_file = open("wmm_inputs.txt", 'w+')
    inputs_file.write(inputs_string)
    inputs_file.seek(0)

    outputs_file = open("wmm_outputs.txt", 'w+')
    sp.call("./wmm_point.exe", stdin=inputs_file, stdout=outputs_file)
    field_data[i,4:7] = process_outputs(outputs_file)

outputs_file.close()
inputs_file.close()
    
# Save the data to csv and npy
np.savetxt("field_data.csv", field_data, delimiter=", ",
           fmt=["%.3f", "%.4f", "%.3f", "%.4f", "%.3f", "%.3f", "%.5e"],
           header="mission_day, lat, long, alt, B decl, B incl, B mag")
np.save("field_data.npy", field_data)

sp.call(["say", "I am done"])
