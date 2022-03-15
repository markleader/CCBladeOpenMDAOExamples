import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

dvs_sf_0 = pd.read_csv('coupled_stress_sf_0.0.csv')
dvs_sf_1 = pd.read_csv('coupled_stress_sf_1.0.csv')
dvs_sf_2 = pd.read_csv('coupled_stress_sf_2.0.csv')
dvs_sf_3 = pd.read_csv('coupled_stress_sf_3.0.csv')
dvs_sf_4 = pd.read_csv('coupled_stress_sf_4.0.csv')

span = 12.0
x0 = 0.2*span
nelems = len(dvs_sf_0)
xe = np.linspace(x0, span, nelems)

n = 8
colors = pl.cm.Blues(np.linspace(0, 1, n))
#colors = ["#364fc7", "#3b5bdb", "#4263eb", "#4c6ef5", "#5c7cfa"]

fig = plt.figure(figsize=(8,4), constrained_layout=True)
ax = plt.subplot()
ax.plot(xe, dvs_sf_0.sigma, color=colors[-1], label='sf = 0')
ax.plot(xe, dvs_sf_1.sigma, color=colors[-2], label='sf = 1')
ax.plot(xe, dvs_sf_2.sigma, color=colors[-3], label='sf = 2')
ax.plot(xe, dvs_sf_3.sigma, color=colors[-4], label='sf = 3')
ax.plot(xe, dvs_sf_4.sigma, color=colors[-5], label='sf = 4')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(True)

for label in ax.get_xticklabels():
    label.set_visible(False)

ax.set_ylabel(r"$\sigma_{vM}(x)/\sigma_y$")
ax.set_xlabel(r'$x_1$ (in)')

ax.legend()

plt.savefig("coupled_stress_sf_compare.pdf", transparent=False)#, dpi=300)