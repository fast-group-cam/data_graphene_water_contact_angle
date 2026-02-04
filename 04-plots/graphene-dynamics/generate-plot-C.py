#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from droplet_graphene_analysis.util.interpolate import PeriodicGridInterpolator

N_PTS = 200
MARGIN = 3.0
X_EDGE = 23.0
Y_MIN = 0.7801
Y_MAX = 0.865

labels = {'s+0.00': 'Free-standing',
          's+0.20': r'+0.2% strain',
          's+0.50': r'+0.5% strain',
          's+2.00': r'+2.0% strain'}

colors = dict()
for i, key in enumerate(labels.keys()):
    x = i / (len(labels) - 1)
    colors[key] = (x, 0.75 - (3.0 * ((x - 0.5)**2)), 1.0 - x)

fig, ax = plt.subplots()
fig.set_size_inches(5.1, 3.4)
plot_xlim = 66.0 #np.inf
delta_z = dict()

for idx, key in enumerate(labels.keys()):

    archive = np.load(os.path.join('data', f'sheet_{key}.npz'))
    cell_params = archive['cell_params']
    incl_norm_autocorr = archive['incl_norm_autocorr']
    incl_norm_autocorr_map = PeriodicGridInterpolator(cell_params[0:2], incl_norm_autocorr)

    x_max = min(cell_params[0], cell_params[1]) / 2.02
    plot_xlim = min(plot_xlim, x_max)
    x_range = np.linspace(-x_max, x_max, N_PTS)
    azi = np.linspace(0.0, 2.0 * np.pi, 180)
    test_pts = np.c_[np.cos(azi), np.sin(azi)]
    incl_autocorr = np.empty((N_PTS,))
    for i in range(0, 100):
        c = np.mean(incl_norm_autocorr_map(x_range[199 - i] * test_pts))
        incl_autocorr[i] = c
        incl_autocorr[199 - i] = c

    ax.plot(x_range, incl_autocorr, '-', color=colors[key], label=labels[key], zorder=1.0)

    if idx == 0:
        part_a = (np.tanh((x_range + X_EDGE) / MARGIN) + 1.0) / 2.0
        part_b = (np.tanh((X_EDGE - x_range) / MARGIN) + 1.0) / 2.0
        part_c = part_a * part_b * np.clip(np.random.normal(1.0, 0.01, size=x_range.shape), a_min=0.0, a_max=None)
        water_height = (part_c + part_c[::-1]) / 2.0
        ax.fill_between(x_range, Y_MIN, Y_MIN + (0.06 * water_height), fc='steelblue', alpha=0.2, zorder=0.01)

ax.set_xlabel(r'$x\;[\AA]$')
ax.set_ylabel(r'$\mathcal{C}_{\text{GS}}(\tau\to\infty)$')
ax.set_xlim(-plot_xlim, plot_xlim)
ax.set_ylim(Y_MIN, Y_MAX)
ax.legend(loc='upper right', framealpha=0.85)
fig.savefig(f'graph_C.png', dpi=(4.0 * fig.dpi), pad_inches=0.05, bbox_inches='tight')
    





