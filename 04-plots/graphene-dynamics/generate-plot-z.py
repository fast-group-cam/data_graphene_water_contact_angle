#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from droplet_graphene_analysis.util.interpolate import PeriodicGridInterpolator

N_PTS = 200
Z_GAP = 2.0
STRETCH_FACTOR = 5.1
DROPLET_OFFSET = 1.1

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
    mean_sheet = archive['mean_sheet']
    mean_heightmap = PeriodicGridInterpolator(cell_params[0:2], mean_sheet)

    x_max = min(cell_params[0], cell_params[1]) / 2.02
    plot_xlim = min(plot_xlim, x_max)
    x_range = np.linspace(-x_max, x_max, N_PTS)
    azi = np.linspace(0.0, 2.0 * np.pi, 180)
    test_pts = np.c_[np.cos(azi), np.sin(azi)]
    z_coords = np.empty((N_PTS,))
    for i in range(0, 100):
        z = np.mean(mean_heightmap(x_range[199 - i] * test_pts))
        z_coords[i] = z
        z_coords[199 - i] = z
    z_offset = z_coords[0] + (Z_GAP * idx)
    delta_z[key] = z_coords[0] - z_coords[N_PTS//2]

    ax.plot(x_range, z_coords - z_offset, '-', color=colors[key], label=labels[key], zorder=1.0)
    if idx == 0:
        interface_points = np.load(os.path.join('data', f'interface_slice_{key}.npy'))
        ax.add_patch(Polygon(np.c_[interface_points[:,0], interface_points[:,2] - z_offset - DROPLET_OFFSET], closed=True, fc='steelblue', ec='none', alpha=0.25))

ax.set_aspect(STRETCH_FACTOR)
ax.set_xlabel(r'$x\;[\AA]$')
ax.set_ylabel(r'$\langle z\rangle\;[\AA]$')
ax.set_xlim(-plot_xlim, plot_xlim)
ax.set_ylim(-(3.8 * Z_GAP), Z_GAP)
ax.set_yticks([-(3.0 * Z_GAP) + (2.0 * i) for i in range(4)], [str(2 * i) for i in range(4)])
#ax.set_yticks([])
for idx, key in enumerate(labels.keys()):
    ax.text(-0.99 * plot_xlim, 0.1 - (Z_GAP * idx), labels[key], color=colors[key], ha='left', va='bottom', fontsize=9.5)
    ax.text(0.99 * plot_xlim, 0.1 - (Z_GAP * idx), r'$\Delta\langle z\rangle \,=\, ' + f'{delta_z[key]:.2f}' + r'\,\AA$', color=colors[key], ha='right', va='bottom', fontsize=9.5)
fig.savefig(f'graph_z.png', dpi=(4.0 * fig.dpi), pad_inches=0.05, bbox_inches='tight')
    





