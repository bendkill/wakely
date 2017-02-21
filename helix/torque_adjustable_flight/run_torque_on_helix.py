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

import math
import numpy as np
import argparse
import ephem

parser = argparse.ArgumentParser(description='Output angle to sun, \
torque on Helix.')
parser.add_argument('-l', '--lat', nargs=1, default=[-80], type=float,
                    help="set the constant latitude")
parser.add_argument('-a', '--alt', nargs=1, default=[3600], type=float,
                    help="set the constant altitude")
parser.add_argument('-f', '--filename', nargs=1,
                    default=["fix_lat_alt_torque_data"],
                    help="set the output filename, no extension")

args = parser.parse_args()
lat = args.lat[0]
alt = args.alt[0]
filename = args.filename[0]

def torque_on_helix(B_e, theta, phi):
    # Angles should be in radians
    R = 0.2286 # meters
    I = 1.437e6 # Amperes
    return 2 * B_e * math.pi * R**2 * I *\
        math.sqrt(math.cos(theta)**2 +
                  (math.sin(theta) * math.cos(phi))**2)

# # mission day, lat, long, alt, (pressure?), horizon alts (deg), azimuth (deg)
# sun_data = np.load("sun_data_alts.npy")

# mission_day, lat, long, alt, B decl, B incl, B mag
earth_field_data = np.load("field_data.npy")

# # mission_day, lat, long, alt, B decl, B incl, B mag, sun horizon alt
# # (deg), sun azimuth (deg), phi, torque (Nm)
# torque_data = np.zeros((len(sun_data), 11))
# torque_data[:,:7] = earth_field_data[:,:7]
# torque_data[:,7:9] = sun_data[:,5:7]

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
