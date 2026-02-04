#! /usr/bin/env python

import os
import time
import numpy as np
from droplet_graphene_analysis.util import read_droplet_trajectory, elapsed_time
from droplet_graphene_analysis.util.graphene.sheet import raw_heightmap
from droplet_graphene_analysis.util.graphene.grid import generate_grid

ROOT_PATTERN = '../../02-strained/$$$$/run_prod'
TRAJ_FILENAME = 'nvt_prod_$$$$.lammpstrj'
ZERO_ROOT = '../../01-unstrained/mace/1000-molecules/run_prod'
STRAINS = [-2.0, -1.5, -1.0, -0.7, -0.5, -0.4, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.2, 0.5, 1.0, 1.5, 2.0]

CALC_POINTS = 200 # Must be even!

if __name__ == '__main__':

    for strain in STRAINS:

        if strain == 0.0:
            folder_name = 's+0.00'
            root_dir = ZERO_ROOT
        else:
            folder_name = f's{strain:+.2f}'
            root_dir = ROOT_PATTERN.replace('$$$$', folder_name)
        
        if not os.path.isfile(os.path.join('cache', f'fourier_{folder_name}.npz')):

            time_start = time.time()
            traj_files = list()
            for i in range(100):
                possible_traj = os.path.join(root_dir, TRAJ_FILENAME.replace('$$$$', str(i)))
                if os.path.isfile(possible_traj):
                    traj_files.append(possible_traj)
            if len(traj_files) < 1:
                print(f'[{folder_name}] ERROR: could not find any trajectory files\n')
                continue
            
            cell_params, _, carbons, _ = read_droplet_trajectory(traj_files, index=':')
            print(f'[{folder_name}] read trajectory files in {elapsed_time(time_start)}')
            time_start = time.time()

            N_frames = carbons.shape[0]
            heightmaps = np.empty((N_frames, CALC_POINTS, CALC_POINTS))
            for f in range(N_frames):
                heightmaps[f] = raw_heightmap(carbons[f], cell_params[0:2], CALC_POINTS)
            grid_pts = generate_grid(CALC_POINTS, cell_params[0:2])
            print(f'[{folder_name}] generated heightmaps in {elapsed_time(time_start)}')
            time_start = time.time()

            phase = np.exp(2j * np.pi * grid_pts[:,:,0] / cell_params[0])
            cx = np.sum(heightmaps * phase[None,:,:], axis=(1,2)) / (CALC_POINTS**2)

            phase = np.exp(2j * np.pi * grid_pts[:,:,1] / cell_params[1])
            cy = np.sum(heightmaps * phase[None,:,:], axis=(1,2)) / (CALC_POINTS**2)

            np.savez_compressed(os.path.join('cache', f'fourier_{folder_name}.npz'), cx=cx, cy=cy)
            print(f'[{folder_name}] performed Fourier transforms in {elapsed_time(time_start)}')

        else:
            print(f'[{folder_name}] skipped as coefficients are already in cache')

