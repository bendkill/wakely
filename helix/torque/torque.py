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
import matplotlib.pyplot as plt
import numpy as np
import argparse
import ephem
from datetime import datetime as dt
from toYearFraction import toYearFraction

parser = argparse.ArgumentParser(description='Output angle to sun, \
torque on Helix.')
parser.add_argument('-l', '--lat', nargs=1, default=[-80], type=float,
                    help="set the constant latitude, in degrees, default = -80")
parser.add_argument('-a', '--alt', nargs=1, default=[36000], type=float,
                    help="set the constant altitude, in meters, default = 36000")
parser.add_argument('-s', '--startdate', nargs=1, default=["2018/12/30"], # 2018. + 364./365.
                    help="set the startdate for projections, yyyy/mm/dd, between 2015 and 2020, do not change if following crest")
parser.add_argument('--crest', action='store_true', 
                    help="follow the original crest flight")
parser.add_argument('-p', '--plots', action='store_true',
                    help="generate plots for torque and sun data")
parser.add_argument('-o', '--overwrite', action='store_true',
                    help="overwrite sun and torque (but not field) data")
parser.add_argument('-f', '--field', action='store_true',
                    help="overwrite existing field data (time intensive)")

args = parser.parse_args()
global lat, alt, start_date, crest, plots
global overwrite, overwrite_field
global crest_flight_data, file_header

lat = args.lat[0]
alt = args.alt[0]
start_date = args.startdate[0]
crest = args.crest
plots = args.plots
overwrite = args.overwrite
overwrite_field = args.field
crest_flight_data = np.load("crest_flight_data.npy")

if crest: file_header = "crest"
else: file_header = "lat_%i_alt_%i" % (int(lat), int(alt))
                                       
def set_lat(new_lat):
    global lat
    lat = new_lat

def set_alt(new_alt):
    global alt
    alt = new_alt

def set_crest():
    global crest
    crest = not crest

def set_plots():
    global plots
    plots = not plots
    
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
            data[0] = (degs + mins / 60.) % 360
        if words[0] == "F":
            nTs = float(words[2])
            data[2] = nTs / 10**9
            break

    return data

def make_field_data():
    out_file_name = "field_data_sets/%s_field_data.npy" % file_header
    if os.path.isfile(out_file_name) and not overwrite_field:
        return np.load(out_file_name)[:,5:8]

    # mission day, lat, long, alt, (pressure?),
    # B local decl (deg), B local incl (deg), B mag (T)
    field_data = np.zeros((len(crest_flight_data), len(crest_flight_data[0]) + 3))
    field_data[:,0:5] = crest_flight_data[:,0:5]
    
    if not crest:                                   
        field_data[:,1] = lat
        field_data[:,3] = alt

    if crest: print("Making field data (crest)...")
    else: print("Making field data; lat: %f, alt: %f..." % (lat, alt))

    start_year = toYearFraction(start_date) # default start_date = '2018/12/30'
    
    for i in range(len(crest_flight_data)):
        # set values according to args
        local_lat = field_data[i,1]
        local_lon = field_data[i,2]
        local_alt = field_data[i,3]
        
        # get the field data for the point
        dec_year = round(start_year + field_data[i,0] / 365., 3)
        inputs_string = "c\n" +\
                        str(local_lat) + "\n" +\
                        str(local_lon) + "\n" +\
                        str(local_alt) + "\n" +\
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
    out_file_name = "sun_data_sets/%s_sun_data.npy" % file_header
    if os.path.isfile(out_file_name) and not overwrite:
        return np.load(out_file_name)[:,5:7]
    
    # mission day, lat, long, alt, (pressure?), alt, azimuth
    sun_data = np.zeros((len(crest_flight_data), len(crest_flight_data[0]) + 2))
    sun_data[:,0:5] = crest_flight_data[:,0:5]
    sun_data[:,1] = lat
    sun_data[:,3] = alt
    
    start_day = ephem.date(start_date)
    sun = ephem.Sun()
    helix = ephem.Observer()
    helix.lat = math.radians(lat)
    helix.elevation = alt

    print("Making sun data...")

    for i in range(len(crest_flight_data)):
        helix.date = start_day + crest_flight_data[i,0]

        helix.lon = math.radians(crest_flight_data[i,2])
        # Can also specify temperature or pressure, if that's what the last
        # column specifies
        sun.compute(helix)
    
        sun_data[i,5] = math.degrees(sun.alt) % 360
        sun_data[i,6] = math.degrees(sun.az) % 360

    np.save(out_file_name, sun_data)

    return sun_data[:,5:7]

def make_torque_data():
    out_file_name = "torque_data_sets/%s_torque_data.npy" % file_header
    if os.path.isfile(out_file_name) and not overwrite:
        return np.load(out_file_name)[:,10:12]

    print("Making torque data...")
    
    for i in range(len(torque_data)):
        B_e = torque_data[i,7]
        theta = math.radians(torque_data[i,6])
        phi = math.radians(torque_data[i,5] - torque_data[i,9]) # B_decl - sun_az
        torque_data[i,10] = math.degrees(phi) % 360
        torque_data[i,11] = torque_on_helix(B_e, theta, phi)

    np.save(out_file_name, torque_data)
    return torque_data[:,10:12]

################################################################################
# plot functions

def my_plot(xs, ys, fmt='r,', xlabel='', ylabel='',
            title='', filename='plots/plot.png'):
    fig = plt.figure(dpi=100)
    ax = fig.add_subplot(1,1,1)
    ax.plot(xs, ys, fmt)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.savefig(filename)
    plt.close()

def make_plots():

    # make plot headers
    if not os.path.isdir("plots/" + file_header):
        os.mkdir("plots/" + file_header)

    if crest: title_header = ", CREST flight path"
    else: title_header = ", Latitude: %i, Altitude: %i" % (int(lat), int(alt))

    print("Making plots...")
    
    # days, phis
    my_plot(torque_data[:,0], torque_data[:,10],
        xlabel="Mission Day",
        ylabel="Angle to Sun (deg)",
        title="Angle to Sun, HELIX Fixed to Earth Field" + title_header,
        filename="plots/%s/angle_to_sun.png" % file_header)

    # days, torques
    my_plot(torque_data[:,0], torque_data[:,11],
        xlabel="Mission Day",
        ylabel="Torque on HELIX (Nm)",
            title="Torque, given HELIX Fixed to Sun" + title_header,
        filename="plots/%s/fixed_to_sun_torques.png" % file_header)

    # days, alts
    my_plot(torque_data[:,0], torque_data[:,8],
        xlabel="Mission Day",
        ylabel="Altitude (deg)",
        title="Altitude Above Horizon of Sun" + title_header,
        filename="plots/%s/sun_alts.png" % file_header)

    # days, azims
    my_plot(torque_data[:,0], torque_data[:,9],
        xlabel="Mission Day",
        ylabel="Azimuth (deg)",
        title="Sun Azimuth (degrees East from North)" + title_header,
        filename="plots/%s/sun_azims.png" % file_header)

    # azims, alts
    my_plot(torque_data[:,9], torque_data[:,8],
        xlabel="Azimuth (deg)",
        ylabel="Altitude (deg)",
        title="Time Independent Path of Sun" + title_header,
        filename="plots/%s/sun_alts_vs_azims.png" % file_header)

    # days, lat
    my_plot(torque_data[:,0], torque_data[:,1],
        xlabel="Mission Day",
        ylabel="Latitude (deg)",
        title="Latitude of Payload" + title_header,
        filename="plots/%s/payload_lats.png" % file_header)

    # days, long
    my_plot(torque_data[:,0], torque_data[:,2],
        xlabel="Mission Day",
        ylabel="Longitude (deg)",
        title="Longitude of Payload" + title_header,
        filename="plots/%s/payload_longs.png" % file_header)

    # days, alt
    my_plot(torque_data[:,0], torque_data[:,3],
        xlabel="Mission Day",
        ylabel="Altitude (m)",
        title="Altitude of Payload" + title_header,
        filename="plots/%s/payload_alts.png" % file_header)

    # days, B local decl
    my_plot(torque_data[:,0], torque_data[:,5],
        xlabel="Mission Day",
        ylabel="B Declination (deg E of N)",
        title="Local Declination of Earth's Magnetic Field" + title_header,
        filename="plots/%s/decl_vs_days.png" % file_header)

    # long, B local decl
    my_plot(torque_data[:,2], torque_data[:,5],
        xlabel="Longitude (deg)",
        ylabel="B Declination (deg E of N)",
        title="Declination of B Field vs Longitude" + title_header,
        filename="plots/%s/decl_vs_longs.png" % file_header)
    
    # days, B mag
    my_plot(torque_data[:,0], torque_data[:,7],
        xlabel="Mission Day",
        ylabel="B Magnitude (T)",
        title="Magnitude of Earth's Magnetic Field" + title_header,
        filename="plots/%s/earth_field_mag.png" % file_header)

    
    
################################################################################
# define and execute main, if this is being run as a script
# (might want to use plot, make_sun_data, and such in the python interpreter)

def main():
    # crest_flight_data has:
    # mission day, lat, long, alt, (pressure?)

    # mission day, lat, long, alt, (pressure?),
    # B local decl (deg), B local incl (deg), B mag (T)
    # sun alt (deg), sun azimuth (deg)
    # phi (helix to sun), torque on helix (Nm)
    global torque_data
    torque_data = np.zeros((len(crest_flight_data),
                            len(crest_flight_data[0]) + 3 + 2 + 2))

    torque_data[:,0:5] = crest_flight_data
    if not crest:
        torque_data[:,1] = lat
        torque_data[:,3] = alt
    torque_data[:,5:8] = make_field_data()
    torque_data[:,8:10] = make_sun_data()
    torque_data[:,10:12] = make_torque_data()

    if plots: make_plots()
    
if __name__ == "__main__": main()
