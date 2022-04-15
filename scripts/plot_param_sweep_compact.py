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
colors = ["#ffebee", "#ffcdd2", "#ef9a9a", "#e57373", "#ef5350", "#f44336", "#e53935", "#d32f2f", "#c62828", "#b71c1c"]
colors_twist = ["#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#ef6c00", "#e65100"]
colors_rpm = ["#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#283593", "#1a237e"]
colors_twist_rpm = ["#e0f2f1", "#b2dfdb", "#80cbc4", "#4db6ac", "#26a69a", "#009688", "#00897b", "#00796b", "#00695c", "#004d40"]
legend_colors = ["#eceff1", "#cfd8dc", "#b0bec5", "#90a4ae", "#78909c", "#607d8b", "#546e7a", "#455a64", "#37474f", "#263238"]

ax_label_fs = 6 # Fontsize for axis labels
ax_fs = 6 # Fontsize for axis values

plt.rc("xtick", labelsize=ax_fs)
plt.rc("ytick", labelsize=ax_fs)

"""
Plot design variables
"""


chord_lower = 0.0
chord_upper = 3.0
theta_lower = 15.0
theta_upper = 75.0
xlower = 2.0
xupper = 12.0

xax_pos = 10.0
yax_pos = 1.75

ytick_vals_theta = [15.0, 30.0, 45.0, 60.0, 75.0]
ytick_vals_chord = [0.0, 1.0, 2.0, 3.0]
xtick_vals = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]

# Plot chord for use_twist=False, use_RPM=False:
ax1 = plt.subplot(4, 2, 1, xlim=(xlower, xupper), ylim=(chord_lower, chord_upper), xticks=xtick_vals, yticks=ytick_vals_chord)
ax1.plot(xe, dvs_sf_0.chord, color=colors[1], clip_on=False)
ax1.plot(xe, dvs_sf_1.chord, color=colors[3], clip_on=False)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.spines["left"].set_position(("data", yax_pos))
ax1.spines["bottom"].set_position(("data", xax_pos))
ax1.grid(True)
for label in ax1.get_xticklabels():
    label.set_visible(False)

#ax1.set_title("No twist, no RPM")
ax1.set_ylabel('Chord (in.)', fontsize=ax_label_fs)

# Plot chord for use_twist=False, use_RPM=True:
ax2 = plt.subplot(4, 2, 2, xlim=(xlower, xupper), ylim=(chord_lower, chord_upper), sharey=ax1, xticks=xtick_vals, yticks=ytick_vals_chord)
ax2.plot(xe, dvs_omega_sf_0.chord, color=colors_rpm[1], clip_on=False)
ax2.plot(xe, dvs_omega_sf_1.chord, color=colors_rpm[3], clip_on=False)
ax2.plot(xe, dvs_omega_sf_2.chord, color=colors_rpm[5], clip_on=False)
ax2.plot(xe, dvs_omega_sf_3.chord, color=colors_rpm[7], clip_on=False)

ax2.spines['left'].set_visible(False)
ax2.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.spines["left"].set_position(("data", yax_pos))
ax2.spines["bottom"].set_position(("data", xax_pos))
ax2.grid(True)
for label in ax2.get_xticklabels():
    label.set_visible(False)
for label in ax2.get_yticklabels():
    label.set_visible(False)

# Plot twist for use_twist=False, use_RPM=False:
ax3 = plt.subplot(4, 2, 3, xlim=(xlower, xupper), ylim=(theta_lower, theta_upper), sharex=ax1, xticks=xtick_vals, yticks=ytick_vals_theta)
ax3.plot(xe, dvs_sf_0.theta, color=colors[1], clip_on=False)
ax3.plot(xe, dvs_sf_1.theta, color=colors[3], clip_on=False)

ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax3.spines["left"].set_position(("data", yax_pos))
ax3.spines["bottom"].set_position(("data", xax_pos))
ax3.set_xticks(xtick_vals)
ax3.grid(True)
for label in ax3.get_xticklabels():
    label.set_visible(False)

ax3.set_ylabel('Twist (deg.)', fontsize=ax_label_fs)

# Plot twist for use_twist=False, use_RPM=True:
ax4 = plt.subplot(4, 2, 4, xlim=(xlower, xupper), ylim=(theta_lower, theta_upper), sharex=ax2, sharey=ax3, xticks=xtick_vals, yticks=ytick_vals_theta)
ax4.plot(xe, dvs_omega_sf_0.theta, color=colors_rpm[1], clip_on=False)
ax4.plot(xe, dvs_omega_sf_1.theta, color=colors_rpm[3], clip_on=False)
ax4.plot(xe, dvs_omega_sf_2.theta, color=colors_rpm[5], clip_on=False)
ax4.plot(xe, dvs_omega_sf_3.theta, color=colors_rpm[7], clip_on=False)

ax4.spines['left'].set_visible(False)
ax4.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax4.spines["left"].set_position(("data", yax_pos))
ax4.spines["bottom"].set_position(("data", xax_pos))
ax4.set_xticks(xtick_vals)
ax4.grid(True)
for label in ax4.get_xticklabels():
    label.set_visible(False)
for label in ax4.get_yticklabels():
    label.set_visible(False)

# Plot chord for use_twist=True, use_RPM=False:
ax5 = plt.subplot(4, 2, 5, xlim=(xlower, xupper), ylim=(chord_lower, chord_upper), sharex=ax1, xticks=xtick_vals, yticks=ytick_vals_chord)
ax5.plot(xe, dvs_theta_sf_0.chord, color=colors_twist[1], clip_on=False)
ax5.plot(xe, dvs_theta_sf_1.chord, color=colors_twist[3], clip_on=False)
ax5.plot(xe, dvs_theta_sf_2.chord, color=colors_twist[5], clip_on=False)
ax5.plot(xe, dvs_theta_sf_3.chord, color=colors_twist[7], clip_on=False)
ax5.plot(xe, dvs_theta_sf_4.chord, color=colors_twist[9], clip_on=False)

ax5.spines['right'].set_visible(False)
ax5.spines['top'].set_visible(False)
ax5.spines['bottom'].set_visible(False)
ax5.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax5.yaxis.set_ticks_position('left')
ax5.xaxis.set_ticks_position('bottom')
ax5.spines["left"].set_position(("data", yax_pos))
ax5.spines["bottom"].set_position(("data", xax_pos))
ax5.grid(True)
for label in ax5.get_xticklabels():
    label.set_visible(False)

ax5.set_ylabel('Chord (in.)', fontsize=ax_label_fs)

# Plot chord for use_twist=True, use_RPM=True:
ax6 = plt.subplot(4, 2, 6, xlim=(xlower, xupper), ylim=(chord_lower, chord_upper), sharex=ax2, sharey=ax5, xticks=xtick_vals, yticks=ytick_vals_chord)
ax6.plot(xe, dvs_theta_omega_sf_0.chord, color=colors_twist_rpm[1], clip_on=False)
ax6.plot(xe, dvs_theta_omega_sf_1.chord, color=colors_twist_rpm[3], clip_on=False)
ax6.plot(xe, dvs_theta_omega_sf_2.chord, color=colors_twist_rpm[5], clip_on=False)
ax6.plot(xe, dvs_theta_omega_sf_3.chord, color=colors_twist_rpm[7], clip_on=False)
ax6.plot(xe, dvs_theta_omega_sf_4.chord, color=colors_twist_rpm[9], clip_on=False)

ax6.spines['left'].set_visible(False)
ax6.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax6.spines['right'].set_visible(False)
ax6.spines['top'].set_visible(False)
ax6.spines['bottom'].set_visible(False)
ax6.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax6.yaxis.set_ticks_position('left')
ax6.xaxis.set_ticks_position('bottom')
ax6.spines["left"].set_position(("data", yax_pos))
ax6.spines["bottom"].set_position(("data", xax_pos))

ax6.grid(True)
for label in ax6.get_xticklabels():
    label.set_visible(False)
for label in ax6.get_yticklabels():
    label.set_visible(False)

# Plot twist for use_twist=True, use_RPM=False:
ax7 = plt.subplot(4, 2, 7, xlim=(xlower, xupper), ylim=(theta_lower, theta_upper), sharex=ax1, xticks=xtick_vals, yticks=ytick_vals_theta)
ax7.plot(xe, dvs_theta_sf_0.theta, color=colors_twist[1], clip_on=False)
ax7.plot(xe, dvs_theta_sf_1.theta, color=colors_twist[3], clip_on=False)
ax7.plot(xe, dvs_theta_sf_2.theta, color=colors_twist[5], clip_on=False)
ax7.plot(xe, dvs_theta_sf_3.theta, color=colors_twist[7], clip_on=False)
ax7.plot(xe, dvs_theta_sf_4.theta, color=colors_twist[9], clip_on=False)

ax7.spines['right'].set_visible(False)
ax7.spines['top'].set_visible(False)
ax7.spines["left"].set_position(("data", yax_pos))
ax7.spines["bottom"].set_position(("data", xax_pos))
ax7.set_xticks(xtick_vals)
ax7.grid(True)

ax7.set_xlabel('x (in.)', fontsize=ax_label_fs)
ax7.set_ylabel('Twist (deg.)', fontsize=ax_label_fs)

# Plot twist for use_twist=True, use_RPM=True:
ax8 = plt.subplot(4, 2, 8, xlim=(xlower, xupper), ylim=(theta_lower, theta_upper), sharex=ax1, sharey=ax7, xticks=xtick_vals, yticks=ytick_vals_theta)
ax8.plot(xe, dvs_theta_omega_sf_0.theta, color=colors_twist_rpm[1], clip_on=False)
ax8.plot(xe, dvs_theta_omega_sf_1.theta, color=colors_twist_rpm[3], clip_on=False)
ax8.plot(xe, dvs_theta_omega_sf_2.theta, color=colors_twist_rpm[5], clip_on=False)
ax8.plot(xe, dvs_theta_omega_sf_3.theta, color=colors_twist_rpm[7], clip_on=False)
ax8.plot(xe, dvs_theta_omega_sf_4.theta, color=colors_twist_rpm[9], clip_on=False)

ax8.spines['left'].set_visible(False)
ax8.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax8.spines['right'].set_visible(False)
ax8.spines['top'].set_visible(False)
ax8.spines["left"].set_position(("data", yax_pos))
ax8.spines["bottom"].set_position(("data", xax_pos))
ax8.set_xticks(xtick_vals)
ax8.grid(True)
for label in ax8.get_yticklabels():
    label.set_visible(False)

ax8.set_xlabel('x (in.)', fontsize=ax_label_fs)

plt.subplots_adjust(hspace=0.175, wspace=0.075)

l1 = plt.legend(
    title="Safety factor",
    bbox_to_anchor=(0.8, 4.95, 0.2, 0.05),
    loc="upper right",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=5,
    handles=[
        mpatches.Patch(color=legend_colors[1], label="0"),
        mpatches.Patch(color=legend_colors[3], label="1"),
        mpatches.Patch(color=legend_colors[5], label="2"),
        mpatches.Patch(color=legend_colors[7], label="3"),
        mpatches.Patch(color=legend_colors[9], label="4")
    ],
)

plt.gca().add_artist(l1)

l2 = plt.legend(
    title="Design variables",
    bbox_to_anchor=(-1.05, 4.95, 0.2, 0.05),
    loc="upper left",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=2,
    handles=[
        mpatches.Patch(color=colors[5], label="No twist, no RPM"),
        mpatches.Patch(color=colors_twist[5], label="Twist, no RPM"),
        mpatches.Patch(color=colors_rpm[5], label="No twist, RPM"),
        mpatches.Patch(color=colors_twist_rpm[5], label="Twist, RPM")
    ],
)

plt.savefig("dvs_compact_compare.png", transparent=False, dpi=300, bbox_inches='tight')


"""
Plot aerodynamic forces
"""


Np_lower = 0.0
Np_upper = 400.0
Tp_lower = 0.0
Tp_upper = 200.0

xax_pos = -10.0

ytick_vals_Np = [0.0, 100.0, 200.0, 300.0, 400.0]
ytick_vals_Tp = [0.0, 50.0, 100.0, 150.0, 200.0]

# Plot Np for use_twist=False, use_RPM=False:
ax1 = plt.subplot(4, 2, 1, xlim=(xlower, xupper), ylim=(Np_lower, Np_upper), xticks=xtick_vals, yticks=ytick_vals_Np)
ax1.plot(xe, forces_sf_0.Np, color=colors[1], clip_on=False)
ax1.plot(xe, forces_sf_1.Np, color=colors[3], clip_on=False)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.spines["left"].set_position(("data", yax_pos))
ax1.spines["bottom"].set_position(("data", xax_pos))
ax1.grid(True)
for label in ax1.get_xticklabels():
    label.set_visible(False)

ax1.set_ylabel('Normal force (N/m)', fontsize=ax_label_fs)

# Plot Np for use_twist=False, use_RPM=True:
ax2 = plt.subplot(4, 2, 2, xlim=(xlower, xupper), ylim=(Np_lower, Np_upper), sharey=ax1, xticks=xtick_vals, yticks=ytick_vals_Np)
ax2.plot(xe, forces_omega_sf_0.Np, color=colors_rpm[1], clip_on=False)
ax2.plot(xe, forces_omega_sf_1.Np, color=colors_rpm[3], clip_on=False)
ax2.plot(xe, forces_omega_sf_2.Np, color=colors_rpm[5], clip_on=False)
ax2.plot(xe, forces_omega_sf_3.Np, color=colors_rpm[7], clip_on=False)

ax2.spines['left'].set_visible(False)
ax2.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.spines["left"].set_position(("data", yax_pos))
ax2.spines["bottom"].set_position(("data", xax_pos))
ax2.grid(True)
for label in ax2.get_xticklabels():
    label.set_visible(False)
for label in ax2.get_yticklabels():
    label.set_visible(False)

# Plot Tp for use_twist=False, use_RPM=False:
ax3 = plt.subplot(4, 2, 3, xlim=(xlower, xupper), ylim=(Tp_lower, Tp_upper), sharex=ax1, xticks=xtick_vals, yticks=ytick_vals_Tp)
ax3.plot(xe, forces_sf_0.Tp, color=colors[1], clip_on=False)
ax3.plot(xe, forces_sf_1.Tp, color=colors[3], clip_on=False)

ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax3.spines["left"].set_position(("data", yax_pos))
ax3.spines["bottom"].set_position(("data", xax_pos))
ax3.set_xticks(xtick_vals)
ax3.grid(True)
for label in ax3.get_xticklabels():
    label.set_visible(False)

ax3.set_ylabel('Transverse force (N/m)', fontsize=ax_label_fs)

# Plot Tp for use_twist=False, use_RPM=True:
ax4 = plt.subplot(4, 2, 4, xlim=(xlower, xupper), ylim=(Tp_lower, Tp_upper), sharex=ax2, sharey=ax3, xticks=xtick_vals, yticks=ytick_vals_Tp)
ax4.plot(xe, forces_omega_sf_0.Tp, color=colors_rpm[1], clip_on=False)
ax4.plot(xe, forces_omega_sf_1.Tp, color=colors_rpm[3], clip_on=False)
ax4.plot(xe, forces_omega_sf_2.Tp, color=colors_rpm[5], clip_on=False)
ax4.plot(xe, forces_omega_sf_3.Tp, color=colors_rpm[7], clip_on=False)

ax4.spines['left'].set_visible(False)
ax4.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax4.spines["left"].set_position(("data", yax_pos))
ax4.spines["bottom"].set_position(("data", xax_pos))
ax4.set_xticks(xtick_vals)
ax4.grid(True)
for label in ax4.get_xticklabels():
    label.set_visible(False)
for label in ax4.get_yticklabels():
    label.set_visible(False)

# Plot Np for use_twist=True, use_RPM=False:
ax5 = plt.subplot(4, 2, 5, xlim=(xlower, xupper), ylim=(Np_lower, Np_upper), sharex=ax1, xticks=xtick_vals, yticks=ytick_vals_Np)
ax5.plot(xe, forces_theta_sf_0.Np, color=colors_twist[1], clip_on=False)
ax5.plot(xe, forces_theta_sf_1.Np, color=colors_twist[3], clip_on=False)
ax5.plot(xe, forces_theta_sf_2.Np, color=colors_twist[5], clip_on=False)
ax5.plot(xe, forces_theta_sf_3.Np, color=colors_twist[7], clip_on=False)
ax5.plot(xe, forces_theta_sf_4.Np, color=colors_twist[9], clip_on=False)

ax5.spines['right'].set_visible(False)
ax5.spines['top'].set_visible(False)
ax5.spines['bottom'].set_visible(False)
ax5.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax5.yaxis.set_ticks_position('left')
ax5.xaxis.set_ticks_position('bottom')
ax5.spines["left"].set_position(("data", yax_pos))
ax5.spines["bottom"].set_position(("data", xax_pos))
ax5.grid(True)
for label in ax5.get_xticklabels():
    label.set_visible(False)

#ax5.set_title("Twist, no RPM")
ax5.set_ylabel('Normal force (N/m)', fontsize=ax_label_fs)

# Plot Np for use_twist=True, use_RPM=True:
ax6 = plt.subplot(4, 2, 6, xlim=(xlower, xupper), ylim=(Np_lower, Np_upper), sharex=ax2, sharey=ax5, xticks=xtick_vals, yticks=ytick_vals_Np)
ax6.plot(xe, forces_theta_omega_sf_0.Np, color=colors_twist_rpm[1], clip_on=False)
ax6.plot(xe, forces_theta_omega_sf_1.Np, color=colors_twist_rpm[3], clip_on=False)
ax6.plot(xe, forces_theta_omega_sf_2.Np, color=colors_twist_rpm[5], clip_on=False)
ax6.plot(xe, forces_theta_omega_sf_3.Np, color=colors_twist_rpm[7], clip_on=False)
ax6.plot(xe, forces_theta_omega_sf_4.Np, color=colors_twist_rpm[9], clip_on=False)

ax6.spines['left'].set_visible(False)
ax6.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax6.spines['right'].set_visible(False)
ax6.spines['top'].set_visible(False)
ax6.spines['bottom'].set_visible(False)
ax6.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax6.yaxis.set_ticks_position('left')
ax6.xaxis.set_ticks_position('bottom')
ax6.spines["left"].set_position(("data", yax_pos))
ax6.spines["bottom"].set_position(("data", xax_pos))

ax6.grid(True)
for label in ax6.get_xticklabels():
    label.set_visible(False)
for label in ax6.get_yticklabels():
    label.set_visible(False)

# Plot Tp for use_twist=True, use_RPM=False:
ax7 = plt.subplot(4, 2, 7, xlim=(xlower, xupper), ylim=(Tp_lower, Tp_upper), sharex=ax1, xticks=xtick_vals, yticks=ytick_vals_Tp)
ax7.plot(xe, forces_theta_sf_0.Tp, color=colors_twist[1], clip_on=False)
ax7.plot(xe, forces_theta_sf_1.Tp, color=colors_twist[3], clip_on=False)
ax7.plot(xe, forces_theta_sf_2.Tp, color=colors_twist[5], clip_on=False)
ax7.plot(xe, forces_theta_sf_3.Tp, color=colors_twist[7], clip_on=False)
ax7.plot(xe, forces_theta_sf_4.Tp, color=colors_twist[9], clip_on=False)

ax7.spines['right'].set_visible(False)
ax7.spines['top'].set_visible(False)
ax7.spines["left"].set_position(("data", yax_pos))
ax7.spines["bottom"].set_position(("data", xax_pos))
ax7.set_xticks(xtick_vals)
ax7.grid(True)

ax7.set_xlabel('x (in.)', fontsize=ax_label_fs)
ax7.set_ylabel('Transverse force (N/m)', fontsize=ax_label_fs)

# Plot Tp for use_twist=True, use_RPM=True:
ax8 = plt.subplot(4, 2, 8, xlim=(xlower, xupper), ylim=(Tp_lower, Tp_upper), sharex=ax1, sharey=ax7, xticks=xtick_vals, yticks=ytick_vals_Tp)
ax8.plot(xe, forces_theta_omega_sf_0.Tp, color=colors_twist_rpm[1], clip_on=False)
ax8.plot(xe, forces_theta_omega_sf_1.Tp, color=colors_twist_rpm[3], clip_on=False)
ax8.plot(xe, forces_theta_omega_sf_2.Tp, color=colors_twist_rpm[5], clip_on=False)
ax8.plot(xe, forces_theta_omega_sf_3.Tp, color=colors_twist_rpm[7], clip_on=False)
ax8.plot(xe, forces_theta_omega_sf_4.Tp, color=colors_twist_rpm[9], clip_on=False)

ax8.spines['left'].set_visible(False)
ax8.tick_params(axis='y', direction='out', length=0.0, width=0.0)
ax8.spines['right'].set_visible(False)
ax8.spines['top'].set_visible(False)
ax8.spines["left"].set_position(("data", yax_pos))
ax8.spines["bottom"].set_position(("data", xax_pos))
ax8.set_xticks(xtick_vals)
ax8.grid(True)
for label in ax8.get_yticklabels():
    label.set_visible(False)

ax8.set_xlabel('x (in.)', fontsize=ax_label_fs)

plt.subplots_adjust(hspace=0.175, wspace=0.075)

l1 = plt.legend(
    title="Safety factor",
    bbox_to_anchor=(0.8, 4.95, 0.2, 0.05),
    loc="upper right",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=5,
    handles=[
        mpatches.Patch(color=legend_colors[1], label="0"),
        mpatches.Patch(color=legend_colors[3], label="1"),
        mpatches.Patch(color=legend_colors[5], label="2"),
        mpatches.Patch(color=legend_colors[7], label="3"),
        mpatches.Patch(color=legend_colors[9], label="4")
    ],
)

plt.gca().add_artist(l1)

l2 = plt.legend(
    title="Design variables",
    bbox_to_anchor=(-1.05, 4.95, 0.2, 0.05),
    loc="upper left",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=2,
    handles=[
        mpatches.Patch(color=colors[5], label="No twist, no RPM"),
        mpatches.Patch(color=colors_twist[5], label="Twist, no RPM"),
        mpatches.Patch(color=colors_rpm[5], label="No twist, RPM"),
        mpatches.Patch(color=colors_twist_rpm[5], label="Twist, RPM")
    ],
)

plt.savefig("aero_forces_compact_compare.png", transparent=False, dpi=300, bbox_inches='tight')



"""
Plot stress distribution
"""


sigma_lower = 0.1
sigma_upper = 1.0

xax_pos = -0.1
yax_pos = 1.85

ytick_vals = [0.0, 0.25, 0.33, 0.5, 1.0]
ytick_labels = ["0.0", "0.25", "", "0.50", "1.00"]

# Plot stress for use_twist=False, use_RPM=False:
ax1 = plt.subplot(4, 1, 1, xlim=(xlower, xupper), ylim=(sigma_lower, sigma_upper), xticks=xtick_vals, yticks=ytick_vals, yticklabels=ytick_labels)
ax1.plot(xe, stress_sf_0.sigma, color=colors[1])
ax1.plot(xe, stress_sf_1.sigma, color=colors[3])

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.spines["left"].set_position(("data", yax_pos))
ax1.spines["bottom"].set_position(("data", xax_pos))
ax1.grid(True)
for label in ax1.get_xticklabels():
    label.set_visible(False)

ax1.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$", fontsize=ax_label_fs)

# Plot stress for use_twist=False, use_RPM=True:
ax2 = plt.subplot(4, 1, 2, xlim=(xlower, xupper), ylim=(sigma_lower, sigma_upper), xticks=xtick_vals, yticks=ytick_vals, yticklabels=ytick_labels, sharex=ax1)
ax2.plot(xe, stress_omega_sf_0.sigma, color=colors_rpm[1])
ax2.plot(xe, stress_omega_sf_1.sigma, color=colors_rpm[3])
ax2.plot(xe, stress_omega_sf_2.sigma, color=colors_rpm[5])
ax2.plot(xe, stress_omega_sf_3.sigma, color=colors_rpm[7])

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.spines["left"].set_position(("data", yax_pos))
ax2.spines["bottom"].set_position(("data", xax_pos))
ax2.grid(True)
for label in ax2.get_xticklabels():
    label.set_visible(False)

ax2.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$", fontsize=ax_label_fs)

# Plot stress for use_twist=True, use_RPM=False:
ax3 = plt.subplot(4, 1, 3, xlim=(xlower, xupper), ylim=(sigma_lower, sigma_upper), xticks=xtick_vals, yticks=ytick_vals, yticklabels=ytick_labels, sharex=ax1)
ax3.plot(xe, stress_theta_sf_0.sigma, color=colors_twist[1])
ax3.plot(xe, stress_theta_sf_1.sigma, color=colors_twist[3])
ax3.plot(xe, stress_theta_sf_2.sigma, color=colors_twist[5])
ax3.plot(xe, stress_theta_sf_3.sigma, color=colors_twist[7])
ax3.plot(xe, stress_theta_sf_4.sigma, color=colors_twist[9])

ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')
ax3.spines["left"].set_position(("data", yax_pos))
ax3.spines["bottom"].set_position(("data", xax_pos))
ax3.grid(True)
for label in ax3.get_xticklabels():
    label.set_visible(False)

ax3.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$", fontsize=ax_label_fs)

# Plot stress for use_twist=True, use_RPM=True:
ax4 = plt.subplot(4, 1, 4, xlim=(xlower, xupper), ylim=(sigma_lower, sigma_upper), xticks=xtick_vals, yticks=ytick_vals, yticklabels=ytick_labels, sharex=ax1)
ax4.plot(xe, stress_theta_omega_sf_0.sigma, color=colors_twist_rpm[1])
ax4.plot(xe, stress_theta_omega_sf_1.sigma, color=colors_twist_rpm[3])
ax4.plot(xe, stress_theta_omega_sf_2.sigma, color=colors_twist_rpm[5])
ax4.plot(xe, stress_theta_omega_sf_3.sigma, color=colors_twist_rpm[7])
ax4.plot(xe, stress_theta_omega_sf_4.sigma, color=colors_twist_rpm[9])

ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.yaxis.set_ticks_position('left')
ax4.xaxis.set_ticks_position('bottom')
ax4.spines["left"].set_position(("data", yax_pos))
ax4.spines["bottom"].set_position(("data", xax_pos))
ax4.grid(True)

ax4.set_xlabel('x (in.)', fontsize=ax_label_fs)
ax4.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$", fontsize=ax_label_fs)

plt.subplots_adjust(hspace=0.175)

l1 = plt.legend(
    title="Safety factor",
    bbox_to_anchor=(0.8, 5.3, 0.2, 0.05),
    loc="upper right",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=5,
    handles=[
        mpatches.Patch(color=legend_colors[1], label="0"),
        mpatches.Patch(color=legend_colors[3], label="1"),
        mpatches.Patch(color=legend_colors[5], label="2"),
        mpatches.Patch(color=legend_colors[7], label="3"),
        mpatches.Patch(color=legend_colors[9], label="4")
    ],
)

plt.gca().add_artist(l1)

l2 = plt.legend(
    title="Design variables",
    bbox_to_anchor=(0.0, 5.3, 0.2, 0.05),
    loc="upper left",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=1,
    handles=[
        mpatches.Patch(color=colors[5], label="No twist, no RPM"),
        mpatches.Patch(color=colors_rpm[5], label="No twist, RPM"),
        mpatches.Patch(color=colors_twist[5], label="Twist, no RPM"),
        mpatches.Patch(color=colors_twist_rpm[5], label="Twist, RPM")
    ],
)

plt.savefig("stress_compact_compare.png", transparent=False, dpi=300, bbox_inches='tight')



"""
Plot displacement distribution
"""


disp_lower = 0.1
disp_upper = 2.0

xax_pos = -0.2

ytick_vals = [0.0, 0.5, 1.0, 1.5, 2.0]

# Plot displacement for use_twist=False, use_RPM=False:
ax1 = plt.subplot(4, 1, 1, xlim=(xlower, xupper), ylim=(disp_lower, disp_upper), xticks=xtick_vals, yticks=ytick_vals)
ax1.plot(xe, disp_sf_0.u, color=colors[1], clip_on=False)
ax1.plot(xe, disp_sf_1.u, color=colors[3], clip_on=False)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.spines["left"].set_position(("data", yax_pos))
ax1.spines["bottom"].set_position(("data", xax_pos))
ax1.grid(True)
for label in ax1.get_xticklabels():
    label.set_visible(False)

ax1.set_ylabel("Displacement (in.)", fontsize=ax_label_fs)

# Plot displacement for use_twist=False, use_RPM=True:
ax2 = plt.subplot(4, 1, 2, xlim=(xlower, xupper), ylim=(disp_lower, disp_upper), xticks=xtick_vals, yticks=ytick_vals, sharex=ax1)
ax2.plot(xe, disp_omega_sf_0.u, color=colors_rpm[1], clip_on=False)
ax2.plot(xe, disp_omega_sf_1.u, color=colors_rpm[3], clip_on=False)
ax2.plot(xe, disp_omega_sf_2.u, color=colors_rpm[5], clip_on=False)
ax2.plot(xe, disp_omega_sf_3.u, color=colors_rpm[7], clip_on=False)

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.spines["left"].set_position(("data", yax_pos))
ax2.spines["bottom"].set_position(("data", xax_pos))
ax2.grid(True)
for label in ax2.get_xticklabels():
    label.set_visible(False)

ax2.set_ylabel("Displacement (in.)", fontsize=ax_label_fs)

# Plot displacement for use_twist=True, use_RPM=False:
ax3 = plt.subplot(4, 1, 3, xlim=(xlower, xupper), ylim=(disp_lower, disp_upper), xticks=xtick_vals, yticks=ytick_vals, sharex=ax1)
ax3.plot(xe, disp_theta_sf_0.u, color=colors_twist[1], clip_on=False)
ax3.plot(xe, disp_theta_sf_1.u, color=colors_twist[3], clip_on=False)
ax3.plot(xe, disp_theta_sf_2.u, color=colors_twist[5], clip_on=False)
ax3.plot(xe, disp_theta_sf_3.u, color=colors_twist[7], clip_on=False)
ax3.plot(xe, disp_theta_sf_4.u, color=colors_twist[9], clip_on=False)

ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.tick_params(axis='x', direction='out', length=0.0, width=0.0)
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')
ax3.spines["left"].set_position(("data", yax_pos))
ax3.spines["bottom"].set_position(("data", xax_pos))
ax3.grid(True)
for label in ax3.get_xticklabels():
    label.set_visible(False)

ax3.set_ylabel("Displacement (in.)", fontsize=ax_label_fs)

# Plot displacement for use_twist=True, use_RPM=True:
ax4 = plt.subplot(4, 1, 4, xlim=(xlower, xupper), ylim=(disp_lower, disp_upper), xticks=xtick_vals, yticks=ytick_vals, sharex=ax1)
ax4.plot(xe, disp_theta_omega_sf_0.u, color=colors_twist_rpm[1], clip_on=False)
ax4.plot(xe, disp_theta_omega_sf_1.u, color=colors_twist_rpm[3], clip_on=False)
ax4.plot(xe, disp_theta_omega_sf_2.u, color=colors_twist_rpm[5], clip_on=False)
ax4.plot(xe, disp_theta_omega_sf_3.u, color=colors_twist_rpm[7], clip_on=False)
ax4.plot(xe, disp_theta_omega_sf_4.u, color=colors_twist_rpm[9], clip_on=False)

ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.yaxis.set_ticks_position('left')
ax4.xaxis.set_ticks_position('bottom')
ax4.spines["left"].set_position(("data", yax_pos))
ax4.spines["bottom"].set_position(("data", xax_pos))
ax4.grid(True)

ax4.set_xlabel('x (in.)', fontsize=ax_label_fs)
ax4.set_ylabel("Displacement (in.)", fontsize=ax_label_fs)

plt.subplots_adjust(hspace=0.175)

l1 = plt.legend(
    title="Safety factor",
    bbox_to_anchor=(0.8, 5.3, 0.2, 0.05),
    loc="upper right",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=5,
    handles=[
        mpatches.Patch(color=legend_colors[1], label="0"),
        mpatches.Patch(color=legend_colors[3], label="1"),
        mpatches.Patch(color=legend_colors[5], label="2"),
        mpatches.Patch(color=legend_colors[7], label="3"),
        mpatches.Patch(color=legend_colors[9], label="4")
    ],
)

plt.gca().add_artist(l1)

l2 = plt.legend(
    title="Design variables",
    bbox_to_anchor=(0.0, 5.3, 0.2, 0.05),
    loc="upper left",
    borderaxespad=0.0,
    fontsize="xx-small",
    title_fontsize="xx-small",
    frameon=False,
    ncol=1,
    handles=[
        mpatches.Patch(color=colors[5], label="No twist, no RPM"),
        mpatches.Patch(color=colors_rpm[5], label="No twist, RPM"),
        mpatches.Patch(color=colors_twist[5], label="Twist, no RPM"),
        mpatches.Patch(color=colors_twist_rpm[5], label="Twist, RPM")
    ],
)

plt.savefig("displacement_compact_compare.png", transparent=False, dpi=300, bbox_inches='tight')

