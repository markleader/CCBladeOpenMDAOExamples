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

#ax2.set_title("No twist, RPM")

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

#ax5.set_title("Twist, no RPM")
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

#ax6.set_title("Twist, RPM")

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

ax7.set_xlabel(r'$x_1$ (in)', fontsize=ax_label_fs)
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

ax8.set_xlabel(r'$x_1$ (in)', fontsize=ax_label_fs)

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