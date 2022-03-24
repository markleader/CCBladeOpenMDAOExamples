import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# File for use_twist=False, use_RPM=False:
dvs_sf_0 = pd.read_csv('dvs_using_chord_sf_0.0.csv')
dvs_sf_1 = pd.read_csv('dvs_using_chord_sf_1.0.csv')

forces_sf_0 = pd.read_csv('aero_forces_using_chord_sf_0.0.csv')
forces_sf_1 = pd.read_csv('aero_forces_using_chord_sf_1.0.csv')

stress_sf_0 = pd.read_csv('stress_using_chord_sf_0.0.csv')
stress_sf_1 = pd.read_csv('stress_using_chord_sf_1.0.csv')

disp_sf_0 = pd.read_csv('disp_using_chord_sf_0.0.csv')
disp_sf_1 = pd.read_csv('disp_using_chord_sf_1.0.csv')

# File for use_twist=False, use_RPM=True:
dvs_omega_sf_0 = pd.read_csv('dvs_using_omega_sf_0.0.csv')
dvs_omega_sf_1 = pd.read_csv('dvs_using_omega_sf_1.0.csv')
dvs_omega_sf_2 = pd.read_csv('dvs_using_omega_sf_2.0.csv')
dvs_omega_sf_3 = pd.read_csv('dvs_using_omega_sf_3.0.csv')

forces_omega_sf_0 = pd.read_csv('aero_forces_using_omega_sf_0.0.csv')
forces_omega_sf_1 = pd.read_csv('aero_forces_using_omega_sf_1.0.csv')
forces_omega_sf_2 = pd.read_csv('aero_forces_using_omega_sf_2.0.csv')
forces_omega_sf_3 = pd.read_csv('aero_forces_using_omega_sf_3.0.csv')

stress_omega_sf_0 = pd.read_csv('stress_using_omega_sf_0.0.csv')
stress_omega_sf_1 = pd.read_csv('stress_using_omega_sf_1.0.csv')
stress_omega_sf_2 = pd.read_csv('stress_using_omega_sf_2.0.csv')
stress_omega_sf_3 = pd.read_csv('stress_using_omega_sf_3.0.csv')

disp_omega_sf_0 = pd.read_csv('disp_using_omega_sf_0.0.csv')
disp_omega_sf_1 = pd.read_csv('disp_using_omega_sf_1.0.csv')
disp_omega_sf_2 = pd.read_csv('disp_using_omega_sf_2.0.csv')
disp_omega_sf_3 = pd.read_csv('disp_using_omega_sf_3.0.csv')

# File for use_twist=True, use_RPM=False:
dvs_theta_sf_0 = pd.read_csv('dvs_using_theta_sf_0.0.csv')
dvs_theta_sf_1 = pd.read_csv('dvs_using_theta_sf_1.0.csv')
dvs_theta_sf_2 = pd.read_csv('dvs_using_theta_sf_2.0.csv')
dvs_theta_sf_3 = pd.read_csv('dvs_using_theta_sf_3.0.csv')
dvs_theta_sf_4 = pd.read_csv('dvs_using_theta_sf_4.0.csv')

forces_theta_sf_0 = pd.read_csv('aero_forces_using_theta_sf_0.0.csv')
forces_theta_sf_1 = pd.read_csv('aero_forces_using_theta_sf_1.0.csv')
forces_theta_sf_2 = pd.read_csv('aero_forces_using_theta_sf_2.0.csv')
forces_theta_sf_3 = pd.read_csv('aero_forces_using_theta_sf_3.0.csv')
forces_theta_sf_4 = pd.read_csv('aero_forces_using_theta_sf_4.0.csv')

stress_theta_sf_0 = pd.read_csv('stress_using_theta_sf_0.0.csv')
stress_theta_sf_1 = pd.read_csv('stress_using_theta_sf_1.0.csv')
stress_theta_sf_2 = pd.read_csv('stress_using_theta_sf_2.0.csv')
stress_theta_sf_3 = pd.read_csv('stress_using_theta_sf_3.0.csv')
stress_theta_sf_4 = pd.read_csv('stress_using_theta_sf_4.0.csv')

disp_theta_sf_0 = pd.read_csv('disp_using_theta_sf_0.0.csv')
disp_theta_sf_1 = pd.read_csv('disp_using_theta_sf_1.0.csv')
disp_theta_sf_2 = pd.read_csv('disp_using_theta_sf_2.0.csv')
disp_theta_sf_3 = pd.read_csv('disp_using_theta_sf_3.0.csv')
disp_theta_sf_4 = pd.read_csv('disp_using_theta_sf_4.0.csv')

# File for use_twist=True, use_RPM=True:
dvs_theta_omega_sf_0 = pd.read_csv('dvs_using_theta_omega_sf_0.0.csv')
dvs_theta_omega_sf_1 = pd.read_csv('dvs_using_theta_omega_sf_1.0.csv')
dvs_theta_omega_sf_2 = pd.read_csv('dvs_using_theta_omega_sf_2.0.csv')
dvs_theta_omega_sf_3 = pd.read_csv('dvs_using_theta_omega_sf_3.0.csv')
dvs_theta_omega_sf_4 = pd.read_csv('dvs_using_theta_omega_sf_4.0.csv')

forces_theta_omega_sf_0 = pd.read_csv('aero_forces_using_theta_omega_sf_0.0.csv')
forces_theta_omega_sf_1 = pd.read_csv('aero_forces_using_theta_omega_sf_1.0.csv')
forces_theta_omega_sf_2 = pd.read_csv('aero_forces_using_theta_omega_sf_2.0.csv')
forces_theta_omega_sf_3 = pd.read_csv('aero_forces_using_theta_omega_sf_3.0.csv')
forces_theta_omega_sf_4 = pd.read_csv('aero_forces_using_theta_omega_sf_4.0.csv')

stress_theta_omega_sf_0 = pd.read_csv('stress_using_theta_omega_sf_0.0.csv')
stress_theta_omega_sf_1 = pd.read_csv('stress_using_theta_omega_sf_1.0.csv')
stress_theta_omega_sf_2 = pd.read_csv('stress_using_theta_omega_sf_2.0.csv')
stress_theta_omega_sf_3 = pd.read_csv('stress_using_theta_omega_sf_3.0.csv')
stress_theta_omega_sf_4 = pd.read_csv('stress_using_theta_omega_sf_4.0.csv')

disp_theta_omega_sf_0 = pd.read_csv('disp_using_theta_omega_sf_0.0.csv')
disp_theta_omega_sf_1 = pd.read_csv('disp_using_theta_omega_sf_1.0.csv')
disp_theta_omega_sf_2 = pd.read_csv('disp_using_theta_omega_sf_2.0.csv')
disp_theta_omega_sf_3 = pd.read_csv('disp_using_theta_omega_sf_3.0.csv')
disp_theta_omega_sf_4 = pd.read_csv('disp_using_theta_omega_sf_4.0.csv')


# Set the x1 variable and set the colors that we will use
span = 12.0
x0 = 0.2*span
nelems = len(dvs_sf_1)
xe = np.linspace(x0, span, nelems)

n = 8
colors = ["#ffebee", "#ffcdd2", "#ef9a9a", "#e57373", "#ef5350", "#f44336", "#e53935", "#d32f2f", "#c62828", "#b71c1c"]  #pl.cm.Blues(np.linspace(0, 1, n))
colors_twist = ["#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#ef6c00", "#e65100"]
colors_rpm = ["#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#283593", "#1a237e"]
colors_twist_rpm = ["#e0f2f1", "#b2dfdb", "#80cbc4", "#4db6ac", "#26a69a", "#009688", "#00897b", "#00796b", "#00695c", "#004d40"]



"""
Plot design variables
"""



# Plot dvs for use_twist=False, use_RPM=False:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, dvs_sf_0.chord, color=colors[1], label="sf = 0")
ax0.plot(xe, dvs_sf_1.chord, color=colors[3], label="sf = 1")

ax1.plot(xe, dvs_sf_0.theta, color=colors[1])
ax1.plot(xe, dvs_sf_1.theta, color=colors[3])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("No twist, no RPM")
ax0.set_ylabel('Chord (in.)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Twist (deg.)')

ax0.legend()

plt.savefig("dvs_sf_compare.png", transparent=False, dpi=300)

# Plot dvs for use_twist=False, use_RPM=True:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, dvs_omega_sf_0.chord, color=colors_rpm[1], label="sf = 0")
ax0.plot(xe, dvs_omega_sf_1.chord, color=colors_rpm[3], label="sf = 1")
ax0.plot(xe, dvs_omega_sf_2.chord, color=colors_rpm[5], label="sf = 2")
ax0.plot(xe, dvs_omega_sf_3.chord, color=colors_rpm[7], label="sf = 3")

ax1.plot(xe, dvs_omega_sf_0.theta, color=colors_rpm[1])
ax1.plot(xe, dvs_omega_sf_1.theta, color=colors_rpm[3])
ax1.plot(xe, dvs_omega_sf_2.theta, color=colors_rpm[5])
ax1.plot(xe, dvs_omega_sf_3.theta, color=colors_rpm[7])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("No twist, RPM")
ax0.set_ylabel('Chord (in.)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Twist (deg.)')

ax0.legend()

plt.savefig("dvs_omega_sf_compare.png", transparent=False, dpi=300)

# Plot dvs for use_twist=True, use_RPM=False:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, dvs_theta_sf_0.chord, color=colors_twist[1], label="sf = 0")
ax0.plot(xe, dvs_theta_sf_1.chord, color=colors_twist[3], label="sf = 1")
ax0.plot(xe, dvs_theta_sf_2.chord, color=colors_twist[5], label="sf = 2")
ax0.plot(xe, dvs_theta_sf_3.chord, color=colors_twist[7], label="sf = 3")
ax0.plot(xe, dvs_theta_sf_4.chord, color=colors_twist[9], label="sf = 4")

ax1.plot(xe, dvs_theta_sf_0.theta, color=colors_twist[1])
ax1.plot(xe, dvs_theta_sf_1.theta, color=colors_twist[3])
ax1.plot(xe, dvs_theta_sf_2.theta, color=colors_twist[5])
ax1.plot(xe, dvs_theta_sf_3.theta, color=colors_twist[7])
ax1.plot(xe, dvs_theta_sf_4.theta, color=colors_twist[9])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("Twist, no RPM")
ax0.set_ylabel('Chord (in.)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Twist (deg.)')

ax0.legend()

plt.savefig("dvs_theta_sf_compare.png", transparent=False, dpi=300)

# Plot dvs for use_twist=True, use_RPM=True:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, dvs_theta_omega_sf_0.chord, color=colors_twist_rpm[1], label="sf = 0")
ax0.plot(xe, dvs_theta_omega_sf_1.chord, color=colors_twist_rpm[3], label="sf = 1")
ax0.plot(xe, dvs_theta_omega_sf_2.chord, color=colors_twist_rpm[5], label="sf = 2")
ax0.plot(xe, dvs_theta_omega_sf_3.chord, color=colors_twist_rpm[7], label="sf = 3")
ax0.plot(xe, dvs_theta_omega_sf_4.chord, color=colors_twist_rpm[9], label="sf = 4")

ax1.plot(xe, dvs_theta_omega_sf_0.theta, color=colors_twist_rpm[1])
ax1.plot(xe, dvs_theta_omega_sf_1.theta, color=colors_twist_rpm[3])
ax1.plot(xe, dvs_theta_omega_sf_2.theta, color=colors_twist_rpm[5])
ax1.plot(xe, dvs_theta_omega_sf_3.theta, color=colors_twist_rpm[7])
ax1.plot(xe, dvs_theta_omega_sf_4.theta, color=colors_twist_rpm[9])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("Twist, RPM")
ax0.set_ylabel('Chord (in.)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Twist (deg.)')

ax0.legend()

plt.savefig("dvs_theta_omega_sf_compare.png", transparent=False, dpi=300)



"""
Plot aerodynamic forces
"""



# Plot the aerodynamic forces for use_twist=False, use_RPM=False:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, forces_sf_0.Np, color=colors[1], label='sf = 0')
ax0.plot(xe, forces_sf_1.Np, color=colors[3], label='sf = 1')

ax1.plot(xe, forces_sf_0.Tp, color=colors[1], label='sf = 0')
ax1.plot(xe, forces_sf_1.Tp, color=colors[3], label='sf = 1')

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("No twist, no RPM")
ax0.set_ylabel('Normal force (N/m)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Transverse force (N/m)')

ax0.legend()

plt.savefig("aero_forces_sf_compare.png", transparent=False, dpi=300)

# Plot the aerodynamic forces for use_twist=False, use_RPM=True:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, forces_omega_sf_0.Np, color=colors_rpm[1], label='sf = 0')
ax0.plot(xe, forces_omega_sf_1.Np, color=colors_rpm[3], label='sf = 1')
ax0.plot(xe, forces_omega_sf_2.Np, color=colors_rpm[5], label='sf = 2')
ax0.plot(xe, forces_omega_sf_3.Np, color=colors_rpm[7], label='sf = 3')

ax1.plot(xe, forces_omega_sf_0.Tp, color=colors_rpm[1])
ax1.plot(xe, forces_omega_sf_1.Tp, color=colors_rpm[3])
ax1.plot(xe, forces_omega_sf_2.Tp, color=colors_rpm[5])
ax1.plot(xe, forces_omega_sf_3.Tp, color=colors_rpm[7])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("No twist, RPM")
ax0.set_ylabel('Normal force (N/m)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Transverse force (N/m)')

ax0.legend()

plt.savefig("aero_forces_omega_sf_compare.png", transparent=False, dpi=300)

# Plot the aerodynamic forces for use_twist=True, use_RPM=False:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, forces_theta_sf_0.Np, color=colors_twist[1], label='sf = 0')
ax0.plot(xe, forces_theta_sf_1.Np, color=colors_twist[3], label='sf = 1')
ax0.plot(xe, forces_theta_sf_2.Np, color=colors_twist[5], label='sf = 2')
ax0.plot(xe, forces_theta_sf_3.Np, color=colors_twist[7], label='sf = 3')
ax0.plot(xe, forces_theta_sf_4.Np, color=colors_twist[9], label='sf = 4')

ax1.plot(xe, forces_theta_sf_0.Tp, color=colors_twist[1])
ax1.plot(xe, forces_theta_sf_1.Tp, color=colors_twist[3])
ax1.plot(xe, forces_theta_sf_2.Tp, color=colors_twist[5])
ax1.plot(xe, forces_theta_sf_3.Tp, color=colors_twist[7])
ax1.plot(xe, forces_theta_sf_4.Tp, color=colors_twist[9])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("Twist, no RPM")
ax0.set_ylabel('Normal force (N/m)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Transverse force (N/m)')

ax0.legend()

plt.savefig("aero_forces_theta_sf_compare.png", transparent=False, dpi=300)

# Plot the aerodynamic forces for use_twist=True, use_RPM=True:
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, forces_theta_omega_sf_0.Np, color=colors_twist_rpm[1], label='sf = 0')
ax0.plot(xe, forces_theta_omega_sf_1.Np, color=colors_twist_rpm[3], label='sf = 1')
ax0.plot(xe, forces_theta_omega_sf_2.Np, color=colors_twist_rpm[5], label='sf = 2')
ax0.plot(xe, forces_theta_omega_sf_3.Np, color=colors_twist_rpm[7], label='sf = 3')
ax0.plot(xe, forces_theta_omega_sf_4.Np, color=colors_twist_rpm[9], label='sf = 4')

ax1.plot(xe, forces_theta_omega_sf_0.Tp, color=colors_twist_rpm[1])
ax1.plot(xe, forces_theta_omega_sf_1.Tp, color=colors_twist_rpm[3])
ax1.plot(xe, forces_theta_omega_sf_2.Tp, color=colors_twist_rpm[5])
ax1.plot(xe, forces_theta_omega_sf_3.Tp, color=colors_twist_rpm[7])
ax1.plot(xe, forces_theta_omega_sf_4.Tp, color=colors_twist_rpm[9])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax0.grid(True)
ax1.grid(True)

for label in ax0.get_xticklabels():
    label.set_visible(False)

ax0.set_title("Twist, RPM")
ax0.set_ylabel('Normal force (N/m)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Transverse force (N/m)')

ax0.legend()

plt.savefig("aero_forces_theta_omega_sf_compare.png", transparent=False, dpi=300)



"""
Plot stress distribution
"""



# Plot the stress distribution for use_twist=False, use_RPM=False:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, stress_sf_0.sigma, color=colors[1], label='sf = 0')
ax.plot(xe, stress_sf_1.sigma, color=colors[3], label='sf = 1')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("No twist, no RPM")
ax.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("stress_sf_compare.png", transparent=False, dpi=300)

# Plot the stress distribution for use_twist=False, use_RPM=True:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, stress_omega_sf_0.sigma, color=colors_rpm[1], label='sf = 0')
ax.plot(xe, stress_omega_sf_1.sigma, color=colors_rpm[3], label='sf = 1')
ax.plot(xe, stress_omega_sf_2.sigma, color=colors_rpm[5], label='sf = 2')
ax.plot(xe, stress_omega_sf_3.sigma, color=colors_rpm[7], label='sf = 3')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("No twist, RPM")
ax.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("stress_omega_sf_compare.png", transparent=False, dpi=300)

# Plot the stress distribution for use_twist=True, use_RPM=False:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, stress_theta_sf_0.sigma, color=colors_twist[1], label='sf = 0')
ax.plot(xe, stress_theta_sf_1.sigma, color=colors_twist[3], label='sf = 1')
ax.plot(xe, stress_theta_sf_2.sigma, color=colors_twist[5], label='sf = 2')
ax.plot(xe, stress_theta_sf_3.sigma, color=colors_twist[7], label='sf = 3')
ax.plot(xe, stress_theta_sf_4.sigma, color=colors_twist[9], label='sf = 4')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("Twist, no RPM")
ax.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("stress_theta_sf_compare.png", transparent=False, dpi=300)

# Plot the stress distribution for use_twist=True, use_RPM=True:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, stress_theta_omega_sf_0.sigma, color=colors_twist_rpm[1], label='sf = 0')
ax.plot(xe, stress_theta_omega_sf_1.sigma, color=colors_twist_rpm[3], label='sf = 1')
ax.plot(xe, stress_theta_omega_sf_2.sigma, color=colors_twist_rpm[5], label='sf = 2')
ax.plot(xe, stress_theta_omega_sf_3.sigma, color=colors_twist_rpm[7], label='sf = 3')
ax.plot(xe, stress_theta_omega_sf_4.sigma, color=colors_twist_rpm[9], label='sf = 4')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("Twist, RPM")
ax.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("stress_theta_omega_sf_compare.png", transparent=False, dpi=300)

"""
Plot displacement distribution
"""



# Plot the displacement distribution for use_twist=False, use_RPM=False:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, disp_sf_0.u, color=colors[1], label='sf = 0')
ax.plot(xe, disp_sf_1.u, color=colors[3], label='sf = 1')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("No twist, no RPM")
ax.set_ylabel("Displacement (in)")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("displacement_sf_compare.png", transparent=False, dpi=300)

# Plot the displacement distribution for use_twist=False, use_RPM=True:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, disp_omega_sf_0.u, color=colors_rpm[1], label='sf = 0')
ax.plot(xe, disp_omega_sf_1.u, color=colors_rpm[3], label='sf = 1')
ax.plot(xe, disp_omega_sf_2.u, color=colors_rpm[5], label='sf = 2')
ax.plot(xe, disp_omega_sf_3.u, color=colors_rpm[7], label='sf = 3')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("No twist, RPM")
ax.set_ylabel("Displacement (in)")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("displacement_omega_sf_compare.png", transparent=False, dpi=300)

# Plot the displacement distribution for use_twist=True, use_RPM=False:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, disp_theta_sf_0.u, color=colors_twist[1], label='sf = 0')
ax.plot(xe, disp_theta_sf_1.u, color=colors_twist[3], label='sf = 1')
ax.plot(xe, disp_theta_sf_2.u, color=colors_twist[5], label='sf = 2')
ax.plot(xe, disp_theta_sf_3.u, color=colors_twist[7], label='sf = 3')
ax.plot(xe, disp_theta_sf_4.u, color=colors_twist[9], label='sf = 4')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("Twist, no RPM")
ax.set_ylabel("Displacement (in)")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("displacement_theta_sf_compare.png", transparent=False, dpi=300)

# Plot the displacement distribution for use_twist=True, use_RPM=True:
fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()

ax.plot(xe, disp_theta_omega_sf_0.u, color=colors_twist_rpm[1], label='sf = 0')
ax.plot(xe, disp_theta_omega_sf_1.u, color=colors_twist_rpm[3], label='sf = 1')
ax.plot(xe, disp_theta_omega_sf_2.u, color=colors_twist_rpm[5], label='sf = 2')
ax.plot(xe, disp_theta_omega_sf_3.u, color=colors_twist_rpm[7], label='sf = 3')
ax.plot(xe, disp_theta_omega_sf_4.u, color=colors_twist_rpm[9], label='sf = 4')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_title("Twist, RPM")
ax.set_ylabel("Displacement (in)")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("displacement_theta_omega_sf_compare.png", transparent=False, dpi=300)