#! /usr/bin/env python

import os
import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

PATTERN = '../../02-strained/$$$$/run_prod/contact-angle/results.ini'
ZERO_FILE = '../../01-unstrained/mace/1000-molecules/run_prod/contact-angle/results.ini'
PATTERN_ANISO = '../../02-strained/$$$$/run_prod/contact-angle-aniso/results.ini'
PATTERN_ANISO_DIST = '../../02-strained/$$$$/run_prod/contact-angle-aniso/sphere_angles.npy'
STRAINS = [-2.0, -1.5, -1.0, -0.7, -0.5, -0.4, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.2, 0.5, 1.0, 1.5, 2.0]
STRAINS_RANDOM = STRAINS[9:]
STRAINS_COHERENT = STRAINS[:9]

ang_random = list()
ang_random_err = list()
ang_coherent_min = list()
ang_coherent_max = list()
ang_coherent_dist = list()

for strain in STRAINS_RANDOM:
    if strain == 0.0:
        filename = ZERO_FILE
    else:
        folder_name = f's{strain:+.2f}'
        filename = PATTERN.replace('$$$$', folder_name)
    if os.path.isfile(filename):
        config = configparser.ConfigParser()
        config.read(filename)
        if 'Block-Averaged Interface' in config:
            a = config.getfloat('Block-Averaged Interface', 'Contact angle [deg], mean of block means', fallback=None)
            b = config.getfloat('Block-Averaged Interface', 'Contact angle [deg], uncertainty', fallback=None)
            if None not in (a, b):
                ang_random.append(a)
                ang_random_err.append(b)
            else:
                raise RuntimeError(f'Error while reading {filename}: field unreadable')
        else:
            raise RuntimeError(f'Error while reading {filename}: missing [Block-Averaged Interface] section')
    else:
        raise RuntimeError(f'{filename} not found!')
    
for strain in STRAINS_COHERENT:
    folder_name = f's{strain:+.2f}'
    filename = PATTERN_ANISO.replace('$$$$', folder_name)
    if os.path.isfile(filename):
        config = configparser.ConfigParser()
        config.read(filename)
        if 'Time-Averaged Interface' in config:
            a = config.getfloat('Time-Averaged Interface', 'Dist. of contact angles, min [deg]', fallback=None)
            b = config.getfloat('Time-Averaged Interface', 'Dist. of contact angles, max [deg]', fallback=None)
            if None not in (a, b):
                ang_coherent_min.append(a)
                ang_coherent_max.append(b)
            else:
                raise RuntimeError(f'Error while reading {filename}: field unreadable')
        else:
            raise RuntimeError(f'Error while reading {filename}: missing [Time-Averaged Interface] section')
    else:
        raise RuntimeError(f'{filename} not found!')
    filename = PATTERN_ANISO_DIST.replace('$$$$', folder_name)
    if os.path.isfile(filename):
        ang_coherent_dist.append(np.load(filename))
    else:
        raise RuntimeError(f'{filename} not found!')
    
strains_random = np.array(STRAINS_RANDOM)
ang_random = np.array(ang_random)
ang_random_err = np.array(ang_random_err)

strains_coherent = np.array(STRAINS_COHERENT)
ang_coherent_min = np.array(ang_coherent_min)
ang_coherent_max = np.array(ang_coherent_max)

fig, ax = plt.subplots()
fig.set_size_inches(10.5, 4)

y_min = 66.3
y_max = 87.0
y_text = (0.9 * y_max) + (0.1 * y_min)
grey = (0.5, 0.5, 0.5)
transition_strain = (np.min(strains_random) + np.max(strains_coherent)) / 2.0
ax.plot((transition_strain, transition_strain), (y_min, y_max), '--', color=grey, zorder=0.01)
ax.annotate('Long-ranged\ncoherent\nrippling', (-1.45, y_text), (-0.73, y_text), xycoords='data', textcoords='data', arrowprops={'arrowstyle':'-|>', 'color': grey}, color=grey, ha='left', va='center', style='italic')
ax.annotate('Random\nrippling', (0.9, y_text), (-0.1, y_text), xycoords='data', textcoords='data', arrowprops={'arrowstyle':'-|>', 'color': grey}, color=grey, ha='left', va='center', style='italic')

guide_x = np.linspace(0.0, 2.23, 200)
guide_y = 74.3 + (2.625 * (guide_x**2))
guide_pts = np.array([guide_x, guide_y]).T.reshape(-1, 1, 2)
guide_segments = np.concatenate([guide_pts[:-1], guide_pts[1:]], axis=1)
lc = LineCollection(guide_segments, zorder=0.0)
lc.set_colors([('steelblue', min(i / 150, 0.7) if (i//5)%2 == 0 else 0.0) for i in range(199)])
ax.add_collection(lc)

vp = ax.violinplot(ang_coherent_dist, strains_coherent, widths=0.05, showmeans=True, showextrema=False)
for pc in vp['bodies']:
    pc.set_facecolor('salmon')
    pc.set_edgecolor('red')
    pc.set_alpha(0.5)
vp['cmeans'].set_edgecolor('red')
artist0 = ax.scatter(strains_coherent, ang_coherent_max, marker='^', fc='red', ec='red')
artist1 = ax.scatter(strains_coherent, ang_coherent_min, marker='v', fc='red', ec='red')
artist2 = ax.errorbar(strains_random, ang_random, yerr=ang_random_err, color='steelblue', fmt='o', mfc='none', mec='steelblue')

ax.set_xlabel(r'Applied strain $\epsilon$ [%]')
ax.set_ylabel(r'Microscopic contact angle $\theta\,[\degree]$')
ax.set_xlim(-2.23, 2.23)
ax.set_ylim(y_min, y_max)

ax.legend([artist0, artist1, artist2], ['Max. anisotropic contact angle', 'Min. anisotropic contact angle', 'Isotropic (mean) contact angle'], loc='lower right')

fig.savefig('strain-contact-angle.png', dpi=(3.2*fig.dpi), bbox_inches='tight', pad_inches=0.05)
