import numpy as np
import ephem
from sys import argv

script, filename = argv

# Expects an Nx5 table with mission day, lat, long, alt, (pressure?)
flight_data = np.load(filename)
# mission day, lat, long, alt (km), (pressure?), horizon alt, azimuth
sun_data = np.zeros((len(flight_data), len(flight_data[0]) + 2))

start_date = ephem.date("2019/12/30")
sun = ephem.Sun()
helix = ephem.Observer()

for i in range(len(flight_data)):
    helix.date = start_date + flight_data[i,0]
    helix.lat = flight_data[i,1]
    helix.lon = flight_data[i,2]
    helix.elevation = flight_data[i,3] * 1000 # convert to meters
    # Can also specify temperature or pressure, if that's what the last
    # column specifies
    sun.compute(helix)

    sun_data[i,0:5] = flight_data[i]
    sun_data[i,5] = sun.alt
    sun_data[i,6] = sun.az

print(sun_data)
np.save("sun_data.npy", sun_data)
np.savetxt("sun_data.csv", sun_data, delimiter=", ",
           header="mission day, lat, long, alt, (pressure?), horizon alt, azimuth",
           fmt=["%.3f", "%.4f", "%.3f", "%.4f", "%.3f", "%.5f", "%.5f"]
           )


