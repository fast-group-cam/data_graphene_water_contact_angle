#! /usr/bin/env python

import sys
import os
import time
import numpy as np
from droplet_graphene_analysis.util import read_droplet_trajectory, elapsed_time
from droplet_graphene_analysis.util.graphene.angle import inclination_norm_inf_autocor

FILENAMES = ['run/nvt_prod_0.lammpstrj', 'run/nvt_prod_1.lammpstrj']
CALC_POINTS = 120 # Must be even!
OUTPUT_FILE = 'inf_autoc.txt'

if __name__ == '__main__':

    for file in FILENAMES:
        if not os.path.isfile(file):
            print(f'[ERROR] {file} not found!')
            sys.exit()
    
    time_start = time.time()
    cell_params, _, carbons, _ = read_droplet_trajectory(FILENAMES, index=':')
    inf_autoc = inclination_norm_inf_autocor(carbons, cell_params[0:2], 30, calc_points=CALC_POINTS)
    N_x, N_y = inf_autoc.shape
    mean_val = np.mean(inf_autoc)
    std_val = np.std(inf_autoc)

    with open(OUTPUT_FILE, 'w') as stream:
        stream.write('Long-time inclination correlation:\n')
        stream.write(f'Mean = {mean_val}\n')
        stream.write(f'Std. Dev. = {std_val}\n')
        stream.write(f'Uncertainty = {std_val / np.sqrt((N_x * N_y) - 1)}\n')
    print(f'Program complete in {elapsed_time(time_start)}')

