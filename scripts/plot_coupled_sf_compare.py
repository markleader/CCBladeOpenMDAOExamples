import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

dvs_sf_0 = pd.read_csv('coupled_chord_theta_sf_0.0.csv')
dvs_sf_1 = pd.read_csv('coupled_chord_theta_sf_1.0.csv')
dvs_sf_2 = pd.read_csv('coupled_chord_theta_sf_2.0.csv')
dvs_sf_3 = pd.read_csv('coupled_chord_theta_sf_3.0.csv')
dvs_sf_4 = pd.read_csv('coupled_chord_theta_sf_4.0.csv')

forces_sf_0 = pd.read_csv('coupled_aero_forces_sf_0.0.csv')
forces_sf_1 = pd.read_csv('coupled_aero_forces_sf_1.0.csv')
forces_sf_2 = pd.read_csv('coupled_aero_forces_sf_2.0.csv')
forces_sf_3 = pd.read_csv('coupled_aero_forces_sf_3.0.csv')
forces_sf_4 = pd.read_csv('coupled_aero_forces_sf_4.0.csv')

span = 12.0
x0 = 0.2*span
nelems = len(dvs_sf_0)
xe = np.linspace(x0, span, nelems)

n = 8
colors = pl.cm.Blues(np.linspace(0, 1, n))
#colors = ["#364fc7", "#3b5bdb", "#4263eb", "#4c6ef5", "#5c7cfa"]

fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, dvs_sf_0.chord, color=colors[-1], label='sf = 0')
ax0.plot(xe, dvs_sf_1.chord, color=colors[-2], label='sf = 1')
ax0.plot(xe, dvs_sf_2.chord, color=colors[-3], label='sf = 2')
ax0.plot(xe, dvs_sf_3.chord, color=colors[-4], label='sf = 3')
ax0.plot(xe, dvs_sf_4.chord, color=colors[-5], label='sf = 4')

ax1.plot(xe, dvs_sf_0.theta, color=colors[-1], label='sf = 0')
ax1.plot(xe, dvs_sf_1.theta, color=colors[-2], label='sf = 1')
ax1.plot(xe, dvs_sf_2.theta, color=colors[-3], label='sf = 2')
ax1.plot(xe, dvs_sf_3.theta, color=colors[-4], label='sf = 3')
ax1.plot(xe, dvs_sf_4.theta, color=colors[-5], label='sf = 4')

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

ax0.set_ylabel('Chord (in.)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Twist (deg.)')

ax0.legend()

plt.savefig("coupled_dvs_sf_compare.pdf", transparent=False)#, dpi=300)

# Plot the aerodynamic forces
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe, forces_sf_0.Np, color=colors[-1], label='sf = 0')
ax0.plot(xe, forces_sf_1.Np, color=colors[-2], label='sf = 1')
ax0.plot(xe, forces_sf_2.Np, color=colors[-3], label='sf = 2')
ax0.plot(xe, forces_sf_3.Np, color=colors[-4], label='sf = 3')
ax0.plot(xe, forces_sf_4.Np, color=colors[-5], label='sf = 4')

ax1.plot(xe, forces_sf_0.Tp, color=colors[-1], label='sf = 0')
ax1.plot(xe, forces_sf_1.Tp, color=colors[-2], label='sf = 1')
ax1.plot(xe, forces_sf_2.Tp, color=colors[-3], label='sf = 2')
ax1.plot(xe, forces_sf_3.Tp, color=colors[-4], label='sf = 3')
ax1.plot(xe, forces_sf_4.Tp, color=colors[-5], label='sf = 4')

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

ax0.set_ylabel('Normal force (N)')
ax1.set_xlabel(r'$x_1$ (in)')
ax1.set_ylabel('Transverse force (N)')

ax0.legend()

plt.savefig("aero_forces_sf_compare.pdf", transparent=False)#, dpi=300)