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

"""

import math
import numpy as np

def torque_on_helix(B_e, theta, phi):
    # Angles should be in radians
    R = 0.2286 # meters
    I = 1.437e6 # Amperes
    return 2 * B_e * math.pi * R**2 * I *\
        math.sqrt(math.cos(theta)**2 +
                  (math.sin(theta) * math.cos(phi))**2)

# mission day, lat, long, alt, (pressure?), horizon alts (deg), azimuth (deg)
sun_data = np.load("sun_data_alts.npy")

# mission_day, lat, long, alt, B decl, B incl, B mag
earth_field_data = np.load("field_data.npy")

# mission_day, lat, long, alt, B decl, B incl, B mag, sun horizon alt
# (deg), sun azimuth (deg), torque (Nm)
torque_data = np.zeros((len(sun_data), 10))
torque_data[:,:7] = earth_field_data[:,:7]
torque_data[:,7:9] = sun_data[:,5:7]

for i in range(len(sun_data)):
    B_e = earth_field_data[i,6]
    theta = math.radians(earth_field_data[i,5])
    phi = math.radians(earth_field_data[i,4] - sun_data[i,5])
    torque_data[i,9] = torque_on_helix(B_e, theta, phi)

np.save("torque_data.npy", torque_data)
