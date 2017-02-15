import numpy as np
import ephem
from math import degrees
from math import radians

# Expects an Nx5 table with mission day, lat, long, alt, (pressure?)
flight_data = np.load("flight_data.npy")

# mission day, lat, long, alt (km), (pressure?), topocentric ra, dec
sun_data = np.zeros((len(flight_data), len(flight_data[0]) + 2))

start_date = ephem.date("2019/12/30")
sun = ephem.Sun()
helix = ephem.Observer()

for i in range(len(flight_data)):
    idx = i                    # change to be stationary or i
    helix.date = start_date + flight_data[i,0]
    # print(flight_data[idx,1])
    helix.lat = radians(flight_data[idx,1])
    helix.lon = radians(flight_data[idx,2])
    # convert to meters:
    helix.elevation = flight_data[idx,3] * 1000
    # Can also specify temperature or pressure, if that's what the last
    # column specifies
    sun.compute(helix)

    sun_data[i,0:5] = flight_data[i]
    sun_data[i,5] = degrees(sun.alt)
    sun_data[i,6] = degrees(sun.az)

np.save("sun_data_alts.npy", sun_data)
np.savetxt("sun_data_alts.csv", sun_data, delimiter=", ",
           header="mission day, lat, long, alt, (pressure?), horizon alts, azimuth",
           fmt=["%.3f", "%.4f", "%.3f", "%.4f", "%.3f", "%.5f", "%.5f"]
           )

