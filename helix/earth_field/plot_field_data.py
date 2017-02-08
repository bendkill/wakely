import matplotlib.pyplot as plt
import numpy as np

# dec_year, lat, long, alt, B decl, B incl, B mag, max torque
torque_data = np.load("torque_data.npy")

days = np.zeros(len(torque_data))
decls = np.zeros(len(torque_data))
torques = np.zeros(len(torque_data))
for i in range(len(torque_data)):
    days[i] = torque_data[i,0]
    decls[i] = torque_data[i,4]
    torques[i] = torque_data[i,7]

fig = plt.figure(dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.plot(days, torques, 'r,')
ax.set_xlabel("Projected Mission Day (Days)")
ax.set_ylabel("Max Torque on Helix (Nm)")
ax.set_title("Maximum Torque on Helix over Mission Duration")
plt.savefig("plots/max_torque.png")
plt.close()

fig = plt.figure(dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.plot(days, decls, 'b,')
ax.set_xlabel("Projected Mission Day (Days)")
ax.set_ylabel("Declination of Field (Degrees)")
ax.set_title("Declination of Field over Mission Duration")
plt.savefig("plots/aligned_declination.png")
plt.close()
