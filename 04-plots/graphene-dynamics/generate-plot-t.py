#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from droplet_graphene_analysis.util.interpolate import PeriodicGridInterpolator

N_PTS = 200

color = (0.0, 0.0, 1.0)

fig, ax = plt.subplots()
fig.set_size_inches(2.6, 1.7)
delta_z = dict()

archive = np.load(os.path.join('data', f'sheet_s+0.00.npz'))
cell_params = archive['cell_params']
incl_norm_autocorr = archive['incl_norm_autocorr']
incl_norm_autocorr_map = PeriodicGridInterpolator(cell_params[0:2], incl_norm_autocorr)

x_max = min(cell_params[0], cell_params[1]) / 2.02
x_range = np.linspace(-x_max, x_max, N_PTS)
azi = np.linspace(0.0, 2.0 * np.pi, 180)
test_pts = np.c_[np.cos(azi), np.sin(azi)]
incl_autocorr = np.empty((N_PTS,))
for i in range(0, 100):
    c = np.mean(incl_norm_autocorr_map(x_range[199 - i] * test_pts))
    incl_autocorr[i] = c
    incl_autocorr[199 - i] = c

C_edge = np.max(incl_autocorr)
C_dry = (incl_autocorr[0] + incl_autocorr[199]) / 2.0
k_edge = 8.213
k_dry = 11.115

def exp_decay(x, k, c):
    return ((1 - c) * np.exp(-k * x)) + c

tau_range = np.linspace(0.0, 1.0, 200)
ax.plot(tau_range, exp_decay(tau_range, k_edge, C_edge), '-', color=color)
ax.plot(tau_range, exp_decay(tau_range, k_dry, C_dry), '-', color=color)

ax.text(0.9, C_edge + 0.01, 'At droplet edge', color=color, ha='right', va='bottom')
ax.text(0.9, C_dry + 0.01, 'Far away', color=color, ha='right', va='bottom')

ax.set_xlim(0.0, 0.92)
ax.set_ylim(0.76, 1.02)
ax.set_xlabel(r'$\tau\;[ps]$')
ax.set_ylabel(r'$\mathcal{C}_{\text{GS}}(\tau)$')
fig.savefig(f'graph_t.png', dpi=(4.0 * fig.dpi), pad_inches=0.05, bbox_inches='tight')

    





