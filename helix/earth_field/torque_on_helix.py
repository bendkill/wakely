# Benjamin Killeen
# January 26, 2017

import numpy as np
import math
import matplotlib.pyplot as plt

def torque_on_helix(B_e, theta, phi):
    # Angles should be in radians
    R = 0.2286 # meters
    I = 1.437e6 # Amperes
    return 2 * B_e * math.pi * R**2 * I *\
        math.sqrt(math.cos(theta)**2 +
                  (math.sin(theta) * math.cos(phi))**2)

# dec_year, lat, long, alt, B decl, B incl, B mag
field_data = np.load("field_data.npy")

# dec_year, lat, long, alt, B decl, B incl, B mag, max torque
# max torque on Helix, torque assuming decl = 90
torque_data = np.zeros((len(field_data), len(field_data[0]) + 1))

for i in range(len(field_data)):
    torque_data[i,0:7] = field_data[i]
    
    # B incl - (90 - lat)
    B_incl_local = math.radians(field_data[i,5] - (90 - field_data[i,1])) 
    B_decl = math.radians(90)

    # field_data[i,6] = B mag
    torque_data[i,7] = torque_on_helix(field_data[i,6], B_incl_local, B_decl)

# Save the data to csv and npy
np.savetxt("torque_data.csv", torque_data, delimiter=", ",
           fmt=["%.3f", "%.4f", "%.3f", "%.4f",
                "%.3f", "%.3f", "%.5e", "%.5e"],
           header="mission_day, lat, long, alt, B decl, B incl, " +\
           "B mag (T), max torque (B_decl = 90)")
np.save("torque_data.npy", torque_data)

# sp.call(["say", "I am done"])

