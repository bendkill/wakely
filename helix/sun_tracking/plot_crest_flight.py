import numpy as np
import matplotlib.pyplot as plt

# mission day, lat, long, alt (km), (pressure?), horizon alt, azimuth
sun_data = np.load("sun_data_alts.npy")
days = sun_data[:,0]
alts = sun_data[:,5]
azims = sun_data[:,6]

crest_lats = sun_data[:,1]
crest_longs = sun_data[:,2]
crest_alts = sun_data[:,3]      # in km

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

def double_plot(xs, ys1, ys2, fmt='r,', xlabel='', ylabel1='', ylabel2='',
                title='', filename='plots/plot.png'):
    fig = plt.figure(dpi=100)
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(xs, ys1, fmt)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1)
    ax1.set_title(title)

    ax2 = ax1.twinx()
    ax2.plot(xs, ys2, 'b,')
    ax2.set_ylabel(ylabel2)
    plt.savefig(filename)
    plt.close()

double_plot(days, alts, crest_alts,
            xlabel="Mission Day",
            ylabel1="Altitude of Sun over Horizon (deg)",
            ylabel2="Altitude of CREST (km)",
            title="Double Plot",
            filename="plots/double_plot.png")
    
my_plot(days, alts,
        xlabel="Mission Day",
        ylabel="Altitude (deg)",
        title="Altitude Above Horizon of Sun For Projected Flight Path",
        filename="plots/crest_sun_alts_crest.png")

my_plot(days, azims,
        xlabel="Mission Day",
        ylabel="Azimuth (deg)",
        title="Sun Azimuth (degres East from North) over CREST Flight Path",
        filename="plots/crest_sun_azims.png")

my_plot(azims, alts,
        xlabel="Azimuth (deg)",
        ylabel="Altitude (deg)",
        title="Time Independent Path of Sun over Projected Flight Data",
        filename="plots/crest_sun_alts_vs_azims.png")

my_plot(days, crest_lats,
        xlabel="Mission Day",
        ylabel="CREST Latitude (deg)",
        title="Latitude of CREST over Flight",
        filename="plots/crest_payload_lats.png")

my_plot(days, crest_longs,
        xlabel="Mission Day",
        ylabel="CREST Longitude (deg)",
        title="Longitudes of CREST over Flight",
        filename="plots/crest_payload_longs.png")

my_plot(days, crest_alts,
        xlabel="Projected Mission Day",
        ylabel="CREST Altitude (km)",
        title="Altitude of Crest over Flight",
        filename="plots/crest_payload_alts.png")
