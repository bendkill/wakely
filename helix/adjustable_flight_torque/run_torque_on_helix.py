"""
Assumes that HELIX points at the sun at all times. Finds the vertical torque on
Helix from this.

theta is the LOCAL inclination. If the payload is level, this is
the inclination from the payload's vertical.

phi is the difference between the sun's declination, sun_decl, and the
Earth field declination, B_decl.

These are then fed into the function:

torque = -2 * B_mag * pi * R^2 * I * sin(theta) * cos(phi)

which yields the (local) vertical component of the torque on the
payload.

For this fixed system, have default values for the latitude and altitude of the payload,
keeping these fixed. --lat and --alt can be used to adjust these. -f
or --filename can be used to change the output file names (enter with no
extensions). The longitude is altered at the usual rate of the CREST
flight.

defaults:
lat = -80 deg
alt = 36000 (m)


"""

import string
import os
import subprocess as sp
import math
import numpy as np
import argparse
import ephem
from datetime import datetime as dt
from toYearFraction import toYearFraction

parser = argparse.ArgumentParser(description='Output angle to sun, \
torque on Helix.')
parser.add_argument('-l', '--lat', nargs=1, default=[-80], type=float,
                    help="set the constant latitude, in degrees")
parser.add_argument('-a', '--alt', nargs=1, default=[3600], type=float,
                    help="set the constant altitude, in meters")
parser.add_argument('-f', '--filename', nargs=1,
                    default=["fix_lat_alt_torque_data"],
                    help="set the output filename, no extension")
parser.add_argument('--crest', action='store_true', 
                    help="follow the original crest flight")
parser.add_argument('-s', '--startdate', nargs=1, default=["2018/12/30"], # 2018. + 364./365.
                    help="set the startdate for projections, yyyy/mm/dd, between 2015 and 2020")
                    
args = parser.parse_args()
global lat
global alt
global filename
global crest
global start_date
global crest_flight_data

lat = args.lat[0]
alt = args.alt[0]
filename = args.filename[0]
crest = args.crest
start_date = args.startdate[0]
crest_flight_data = np.load("crest_flight_data.npy")

def set_lat(new_lat):
    global lat
    lat = new_lat

def set_alt(new_alt):
    global alt
    alt = new_alt

################################################################################
# main torque function

def torque_on_helix(B_e, theta, phi):
    # Angles should be in radians
    R = 0.2286 # meters
    I = 1.437e6 # Amperes
    return 2 * B_e * math.pi * R**2 * I *\
        math.sqrt(math.cos(theta)**2 +
                  (math.sin(theta) * math.cos(phi))**2)

################################################################################
# functions for creating various data sets

def process_wmm_outputs(handle):
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

def make_field_data():
    if crest:
        return np.load("field_data_sets/crest.npy")[:,4:7]

    out_file_name = "field_data_sets/lat_%f_alt_%f_field_data.npy" % (lat, alt)
    if os.path.isfile(out_file_name):
        return np.load(out_file_name)[:,5:8]

    # mission day, lat, long, alt, (pressure?),
    # B local decl (deg), B local incl (deg), B mag (T)
    field_data = np.zeros((len(crest_flight_data), len(crest_flight_data[0]) + 3))
    field_data[:,0:5] = crest_flight_data

    print("Making field data; lat: %f, alt: %f..." % (lat, alt))

    start_year = toYearFraction(start_date) # default start_date = '2018/12/30'
    
    for i in range(len(crest_flight_data)):
        # set values according to args
        lon = crest_flight_data[i,2]
        
        # get the field data for the point
        dec_year = round(start_year + crest_flight_data[i,0] / 365., 3)
        inputs_string = "c\n" +\
                        str(lat) + "\n" +\
                        str(lon) + "\n" +\
                        str(alt) + "\n" +\
                        str(dec_year) + "\n"

        inputs_file = open("wmm_inputs.txt", 'w+')
        inputs_file.write(inputs_string)
        inputs_file.seek(0)

        outputs_file = open("wmm_outputs.txt", 'w+')
        sp.call("./wmm_point.exe", stdin=inputs_file, stdout=outputs_file)
        field_data[i,5:8] = process_wmm_outputs(outputs_file)

    outputs_file.close()
    inputs_file.close()

    np.save(out_file_name, field_data)
    
    return field_data[:,5:8]

def make_sun_data():
    if crest:
        return np.load('sun_data_sets/crest.npy')[:,5:7]

    out_file_name = "sun_data_sets/lat_%f_alt_%f_sun_data.npy" % (lat, alt)
    if os.path.isfile(out_file_name):
        return np.load(out_file_name)[:,5:7]
    
    # mission day, lat, long, alt, (pressure?), alt, azimuth
    sun_data = np.zeros((len(crest_flight_data), len(crest_flight_data[0]) + 2))
    
    start_day = ephem.date(start_date)
    sun = ephem.Sun()
    helix = ephem.Observer()
    helix.lat = radians(lat)
    helix.elevation = alt

    print("Making sun data...")

    for i in range(len(flight_data)):
        helix.date = start_date + flight_data[i,0]

        helix.lon = radians(flight_data[idx,2])
        # Can also specify temperature or pressure, if that's what the last
        # column specifies
        sun.compute(helix)
    
        sun_data[i,5] = degrees(sun.alt)
        sun_data[i,6] = degrees(sun.az)

    np.save(out_file_name, sun_data)

    return sun_data[:,5:7]

def make_torque_data():
    

################################################################################
# plot functions


################################################################################
# define and execute main, if this is being run as a script
# (might want to use plot, make_sun_data, and such in the python interpreter)

def main():
    # crest_flight_data has:
    # mission day, lat, long, alt, (pressure?)

    # mission day, lat, long, alt, (pressure?),
    # B local decl (deg), B local incl (deg), B mag (T)
    torque_data = np.zeros((len(crest_flight_data),
                            len(crest_flight_data) + 3 + 2))

    torque_data[:,0:5] = crest_flight_data
    torque_data[:,5:8] = make_field_data(crest_flight_data)
    torque_data[:,8:10] = make_sun_data(crest_flight_data)
    

if __name__ == "__main__": main()

# for i in range(len(sun_data)):
#     B_e = earth_field_data[i,6]
#     theta = math.radians(earth_field_data[i,5])
#     phi = math.radians(earth_field_data[i,4] - sun_data[i,5])
#     torque_data[i,9] = math.degrees(phi)
#     torque_data[i,10] = torque_on_helix(B_e, theta, phi)

# np.save("torque_data_fixed_sun.npy", torque_data)
# np.savetxt("torque_data_fixed_sun.csv", torque_data, delimiter=", ",
#            fmt=["%.3f", "%.4f", "%.3f", "%.4f",
#                 "%.3f", "%.3f", "%.5e", "%.5e",
#                 "%.3f", "%.3f", "%.3f"],
#            header="mission_day, lat, long, alt, B decl, B incl, " +\
#            "B mag (T), sun horizon alt, sun azimuth (E from N), phi " +\
#            "(sun to field, deg), torque (Nm)")
