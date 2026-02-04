#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.cm import ScalarMappable
from droplet_graphene_analysis.util.droplet.coarse_grain import BULK_DENSITY

TRAJECTORY_FILES = '../../01-unstrained/mace/2000-molecules/run_prod/nvt_prod_$$$$.lammpstrj'
INDEX = '0::10'
N_PLOT_BINS = 201

#==================================================================================================

if not os.path.isfile('cache.npz'):

    import time
    from droplet_graphene_analysis.util import read_droplet_trajectory, elapsed_time
    from droplet_graphene_analysis.util.droplet.coarse_grain import coarse_grained_density, SLICING_CUTOFF, COARSE_GRAIN_LENGTH
    from droplet_graphene_analysis.util.droplet.contact_angle import find_spherical_cap
    from droplet_graphene_analysis.util.graphene.sheet import raw_heightmap
    from droplet_graphene_analysis.util.interpolate import PeriodicGridInterpolator

    trajectories = list()
    for i in range(50):
        filename = TRAJECTORY_FILES.replace('$$$$', str(i))
        if os.path.isfile(filename):
            trajectories.append(filename)
    if len(trajectories) == 0:
        raise RuntimeError('No files found!')
    
    time_start = time.time()
    cell_params, waters, carbons, _ = read_droplet_trajectory(trajectories, index=INDEX)
    N_frames = waters.shape[0]
    print(f'Read {N_frames} frames in {elapsed_time(time_start)}.')

    time_start = time.time()
    sheet_Nx = int(np.ceil(1.5 * cell_params[0]))
    sheet_Ny = int(np.ceil(1.5 * cell_params[1]))
    sheets = np.empty((N_frames, sheet_Nx, sheet_Ny), dtype=float)
    for f in range(N_frames):
        sheets[f] = raw_heightmap(carbons[f], cell_params[0:2], (sheet_Nx, sheet_Ny))
    mean_heightmap = PeriodicGridInterpolator(cell_params[0:2], np.mean(sheets, axis=0))
    plot_zmin = mean_heightmap(np.array([0.0, 0.0]))[0] - 5.0
    print(f'Calculated mean heightmap in {elapsed_time(time_start)}.')

    time_start = time.time()
    spherical_cap = find_spherical_cap(waters, cell_params, mean_heightmap)
    sphere_r = spherical_cap['r']
    sphere_z = spherical_cap['z']
    sphere_a = spherical_cap['a']
    contact_angle = spherical_cap['angle']
    plot_zmax = sphere_z + sphere_r + 8.0
    print(f'Calculated best-fit sphere in {elapsed_time(time_start)}.')

    time_start = time.time()
    plot_xmax = 1.5 * sphere_a
    carbon_line = np.empty((N_PLOT_BINS, 2), dtype=float)
    carbon_line[:,0] = np.linspace(-plot_xmax, plot_xmax, N_PLOT_BINS, endpoint=True)
    theta = np.linspace(0.0, 2.0 * np.pi, N_PLOT_BINS, endpoint=False)
    azi = np.c_[np.cos(theta), np.sin(theta)]
    for i in range(N_PLOT_BINS):
        carbon_line[i,1] = np.mean(mean_heightmap(carbon_line[i,0] * azi))
    print(f'Calculated mean carbon line in {elapsed_time(time_start)}.')

    time_start = time.time()
    slice_width = SLICING_CUTOFF * COARSE_GRAIN_LENGTH
    skip_interval = int(N_frames // 10)
    waters_skipped = waters[::skip_interval]
    N_frames_reduced = waters_skipped.shape[0]
    waters_skipped = waters_skipped.reshape(-1, 3)
    sliced = waters_skipped[waters_skipped[:,1] < slice_width]
    sliced = sliced[sliced[:,1] > -slice_width]
    x_space = np.linspace(-plot_xmax, plot_xmax, N_PLOT_BINS, endpoint=True)
    z_space = np.linspace(plot_zmin, plot_zmax, N_PLOT_BINS, endpoint=True)
    xx, zz = np.meshgrid(x_space, z_space)
    testpoints = np.column_stack((xx.ravel(), np.zeros(N_PLOT_BINS**2), zz.ravel()))
    densities = coarse_grained_density(testpoints, sliced) / N_frames_reduced
    densities = np.reshape(densities, (N_PLOT_BINS, N_PLOT_BINS))
    print(f'Calculated water densities in {elapsed_time(time_start)}.')

    np.savez_compressed('cache.npz', plot_xmax=plot_xmax, plot_zmin=plot_zmin, plot_zmax=plot_zmax,
                        sphere_r=sphere_r, sphere_z=sphere_z, sphere_a=sphere_a,
                        contact_angle=contact_angle, carbon_line=carbon_line, densities=densities)

#==================================================================================================

fig, ax = plt.subplots()
fig.set_size_inches(7.0, 3.5)

loaded = np.load('cache.npz')
plot_xmax = loaded['plot_xmax']
plot_zmin = loaded['plot_zmin']
plot_zmax = loaded['plot_zmax']
sphere_r = loaded['sphere_r']
sphere_z = loaded['sphere_z']
sphere_a = loaded['sphere_a']
contact_angle = loaded['contact_angle']
carbon_line = loaded['carbon_line']
densities = loaded['densities']

x_pad = plot_xmax / (N_PLOT_BINS - 1)
z_pad = 0.5 * (plot_zmax - plot_zmin) / (N_PLOT_BINS - 1)
colors = np.zeros((N_PLOT_BINS, N_PLOT_BINS, 4), dtype=float)
colors[:,:,0] = np.clip((densities / BULK_DENSITY) - 1.0, a_min=0.0, a_max=1.0)
colors[:,:,2] = np.clip(2.0 - (densities / BULK_DENSITY), a_min=0.0, a_max=1.0)
colors[:,:,3] = np.clip((densities / BULK_DENSITY), a_min=0.0, a_max=1.0)

ax.imshow(colors, origin='lower', extent=(-(plot_xmax + x_pad), plot_xmax + x_pad, plot_zmin - z_pad, plot_zmax + z_pad), zorder=0.0)
ax.plot(carbon_line[:,0], carbon_line[:,1], '-', color=(0.5, 0.5, 0.5), zorder=0.5)

max_phi = np.arcsin(np.clip(sphere_a / sphere_r, -1.0, 1.0))
phi = np.linspace(-max_phi, max_phi, N_PLOT_BINS)
ax.plot(sphere_r * np.sin(phi), sphere_z + (sphere_r * np.cos(phi)), '-', color=(1.0, 0.5, 0.0), zorder=1.0)

intersection_z = sphere_z + (sphere_r * np.cos(max_phi))
gradient = (np.interp(sphere_a + 0.1, carbon_line[:,0], carbon_line[:,1]) - np.interp(sphere_a - 0.1, carbon_line[:,0], carbon_line[:,1])) / 0.2
min_phi = -np.arctan(gradient)
phi = np.linspace(min_phi, max_phi, N_PLOT_BINS//2)
ax.plot((-sphere_a + (20.0 * np.cos(max_phi)), -sphere_a), (intersection_z + (20.0 * np.sin(max_phi)), intersection_z), '--', color='black', zorder=1.5)
ax.plot((sphere_a - (20.0 * np.cos(max_phi)), sphere_a), (intersection_z + (20.0 * np.sin(max_phi)), intersection_z), '--', color='black', zorder=1.5)
ax.plot(-sphere_a + (6.0 * np.cos(phi)), intersection_z + (6.0 * np.sin(phi)), '-', color='black', zorder=1.5)
ax.plot(sphere_a - (6.0 * np.cos(phi)), intersection_z + (6.0 * np.sin(phi)), '-', color='black', zorder=1.5)

ax.set_xlabel(r'$x\;[\AA]$')
ax.set_ylabel(r'$z\;[\AA]$')
ax.set_xlim(-plot_xmax, plot_xmax)
ax.set_ylim(plot_zmin, plot_zmax)
ax.set_aspect('equal')

ax.text(0.985, 0.975, (r'$\theta\;=\;' + f'{contact_angle:.1f}' + r'\degree$'), ha='right', va='top', transform=ax.transAxes)

cmap = mplc.LinearSegmentedColormap('water_density', {'red':   [(0.0, 1.0, 1.0), (0.6667, 0.0, 0.0), (1.0, 1.0, 1.0)],
                                                      'green': [(0.0, 1.0, 1.0), (0.6667, 0.0, 0.0), (1.0, 0.0, 0.0)],
                                                      'blue':  [(0.0, 1.0, 1.0), (0.6667, 1.0, 1.0), (1.0, 0.0, 0.0)]})
norm = mplc.Normalize(vmin=0.0, vmax=(1.5 * BULK_DENSITY))
fig.colorbar(ScalarMappable(norm, cmap), ax=ax, label=r'$\langle\rho(\mathbf{r})\rangle_{t}\;\;[\AA^{-3}$]', fraction=0.02, pad=0.04)

fig.savefig(f'best-fit-sphere.png', dpi=(3.0 * fig.dpi), pad_inches=0.05, bbox_inches='tight')
    





