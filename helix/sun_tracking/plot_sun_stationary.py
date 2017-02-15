import numpy as np
import matplotlib.pyplot as plt

# mission day, lat, long, alt (km), (pressure?), horizon alt, azimuth
sun_data = np.load("sun_data_alts.npy")
days = sun_data[:,0]
alts = sun_data[:,5]
azims = sun_data[:,6]
helix_alts = sun_data[:,3]

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

my_plot(days, alts,
        xlabel="Projected Mission Day",
        ylabel="Altitude (deg)",
        title="Altitude Above Horizon of Sun for a Fixed Point",
        filename="plots/stationary_sun_alts.png")

my_plot(days, azims,
        xlabel="Projected Mission Day",
        ylabel="Azimuth (deg)",
        title="Sun Azimuth (degrees East from North) for a Fixed Point",
        filename="plots/stationary_sun_azims.png")

my_plot(azims, alts,
        xlabel="Azimuth (deg)",
        ylabel="Altitude (deg)",
        title="Time Independent Path of Sun for a Fixed Point",
        filename="plots/stationary_sun_alts_vs_azims.png")

