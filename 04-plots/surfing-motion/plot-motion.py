#! /usr/bin/env python

import time
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from droplet_graphene_analysis.util import elapsed_time
from droplet_graphene_analysis.util.droplet.contact_angle import find_spherical_cap_aniso
from droplet_graphene_analysis.util.graphene.sheet import regularized_heightmap
from droplet_graphene_analysis.util.interpolate import PeriodicGridInterpolator

TIME_INTERVAL = 80
OUTPUT_FILE = 'surfing-motion.png'

if __name__ == '__main__':
        
    fig, ax = plt.subplots(1, 4, figsize=(10, 3), layout='constrained')

    time_start = time.time()
    loaded = np.load('frames.npz')
    cell_params = loaded['cell_params']
    waters = loaded['waters']
    carbons = loaded['carbons']
    shifts = loaded['shifts']
    print(f'Read trajectory in {elapsed_time(time_start)}')

    cell_xy = cell_params[0:2]
    z_width = max(np.max(carbons[:,:,2]) - np.min(carbons[:,:,2]), 1.0)
    azi = np.linspace(0.0, 2.0 * np.pi, 200, endpoint=False)
    search_dirs = np.c_[np.cos(azi), np.sin(azi)]
    sheet_N = int(10.0 * np.sqrt(carbons.shape[1]))
    dx = cell_params[0] / sheet_N
    dy = cell_params[1] / sheet_N
    extent = (-0.5 * cell_params[0], 0.5 * cell_params[0], -0.5 * cell_params[1], 0.5 * cell_params[1])

    for f in range(4):

        time_start = time.time()
        graphene_sheet = regularized_heightmap(carbons[f], cell_xy, sheet_N)
        heightmap = PeriodicGridInterpolator(cell_xy, graphene_sheet)
        spherical_cap = find_spherical_cap_aniso(waters[f], cell_params, heightmap)
        sphere_r = spherical_cap['r']
        sphere_c = spherical_cap['c']

        CoM = np.mean(waters[f], axis=0)
        inters = np.empty((200, 3))
        for i, search_dir in enumerate(search_dirs):
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    floor = heightmap(sphere_c[0:2])[0]
                    for _ in range(15):
                        a = np.sqrt((sphere_r**2) - ((sphere_c[2] - floor)**2))
                        test_pt = (a * search_dir) + sphere_c[0:2]
                        floor = heightmap(test_pt)[0]
                    a = np.sqrt((sphere_r**2) - ((sphere_c[2] - floor)**2))
                    inters[i,0:2] = (a * search_dir) + sphere_c[0:2]
                    inters[i,2] = heightmap(inters[i,0:2])[0]
                except RuntimeWarning:
                    inters[i] = CoM

        CoM = CoM[0:2] + shifts[f,0:2]
        inters = inters[:,0:2] + shifts[f,0:2]
        carbon_xy = carbons[f,:,0:2] + shifts[f,0:2]
        carbon_z = carbons[f,:,2]

        CoM -= cell_xy * np.round(CoM / cell_xy)
        inters -= cell_xy * np.round(inters / cell_xy)
        carbon_xy -= cell_xy * np.round(carbon_xy / cell_xy)

        shift_quant_x = int(shifts[f,0] / dx)
        shift_quant_y = int(shifts[f,1] / dy)
        graphene_sheet = np.roll(graphene_sheet, shift_quant_x, axis=0)
        graphene_sheet = np.roll(graphene_sheet, shift_quant_y, axis=1)
        scattercolors = np.empty(graphene_sheet.shape + (3,))
        z_devs = (graphene_sheet - np.min(graphene_sheet)) / z_width
        z_devs += 0.5 * (1.0 - np.max(z_devs))
        scattercolors[:,:,0] = np.clip(z_devs, a_min=0.0, a_max=1.0)
        scattercolors[:,:,1] = np.clip((z_devs - np.square(z_devs)) / 0.3, a_min=0.0, a_max=1.0)
        scattercolors[:,:,2] = np.clip(1.0 - z_devs, a_min=0.0, a_max=1.0)

        ax[f].imshow(np.swapaxes(scattercolors, 0, 1), origin='lower', extent=extent)
        ax[f].fill(inters[:,0], inters[:,1], color=(0.275, 0.510, 0.706, 0.4))
        inters = np.concat((inters, np.atleast_2d(inters[0])))
        ax[f].plot(inters[:,0], inters[:,1], '-', color='black')
        ax[f].plot((CoM[0],), (CoM[1],), 'o', color='black')

        ax[f].set_title(f't = {(f * TIME_INTERVAL):.0f} ps')
        ax[f].set_xlabel(r'x [$\AA$]')
        if f == 0:
            ax[f].set_ylabel(r'y [$\AA$]')
        ax[f].set_aspect('equal')

        print(f'Plotted frame #{f+1} in {elapsed_time(time_start)}')

    x = np.linspace(0.0, 1.0, 256)
    y = np.clip((x - (x**2)) / 0.3, a_min=0.0, a_max=1.0)
    cdict = {'red': [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)],
             'green': [(x[i], y[i], y[i]) for i in range(256)],
             'blue': [(0.0, 1.0, 1.0), (1.0, 0.0, 0.0)]}
    cmap = LinearSegmentedColormap('custom_colormap', cdict)
    fig.suptitle('Motion of droplet CoM')
    fig.tight_layout()
    fig.colorbar(ScalarMappable(norm=Normalize(-z_width, z_width), cmap=cmap),
                 ax=ax, shrink=0.6, orientation='vertical', label=r'z [$\AA$]')
    fig.savefig(OUTPUT_FILE, dpi=(3 * fig.dpi), bbox_inches='tight', pad_inches=0.05)

