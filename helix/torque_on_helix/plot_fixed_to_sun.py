import numpy as np
import matplotlib.pyplot as plt

# mission_day, lat, long, alt, B decl, B incl, B mag, sun horizon alt
# (deg), sun azimuth (deg), torque (Nm)
torque_data = np.load("torque_data.npy")
days = torque_data[:,0]
torques = torque_data[:,9]
sun_azims = torque_data[:,8]

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
 
my_plot(days, torques,
        xlabel="Projected Mission Day",
        ylabel="Torque on HELIX (Nm)",
        title="Torque on HELIX over Projected Flight, Fixed to Sun",
        filename="plots/fixed_to_sun_torques.png")

double_plot(days, torques, sun_azims,
            xlabel="Projected Mission Day",
            ylabel1="Torque on HELIX (Nm)",
            ylabel2="Azimuth of Sun (deg)",
            title="Comparison of Torque on Helix with Azimuth of Sun",
            filename="plots/torque_sun_comparison.png")
