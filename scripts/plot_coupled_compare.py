import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df_aero = pd.read_csv('coupled_chord_theta_sf_0.0.csv')
df_struct = pd.read_csv('chord_theta_SNOPT_w_splines_no_ks_sf_4.0.csv')
df_coupled = pd.read_csv('coupled_chord_theta_sf_4.0.csv')

span = 12.0
x0 = 0.2*span
xe1 = np.linspace(x0, span, len(df_aero))
xe2 = np.linspace(x0, span, len(df_struct))

colors = ["tab:blue", "tab:red", "tab:purple"]
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
ax0.plot(xe1, df_aero.chord, color=colors[0], label='Aerodynamic')
ax0.plot(xe2, df_struct.chord, color=colors[1], label='Structural')
ax0.plot(xe1, df_coupled.chord, color=colors[2], label='Coupled')

ax1.plot(xe1, df_aero.theta, color=colors[0], label='Aerodynamic')
ax1.plot(xe2, df_struct.theta, color=colors[1], label='Structural')
ax1.plot(xe1, df_coupled.theta, color=colors[2], label='Coupled')

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

plt.savefig("chord_theta_problem_compare.pdf", transparent=False)#, dpi=300)