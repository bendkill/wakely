import numpy as np
import ephem
from sys import argv

script, filename = argv

# Expects an Nx5 table with mission day, lat, long, alt
flight_data = np.load(filename)
# mission day, lat, long, alt, 
sun_data = np.zeros((len(flight_data), 4))

start_date = ephem.date("2019/12/30")
                     
for 
