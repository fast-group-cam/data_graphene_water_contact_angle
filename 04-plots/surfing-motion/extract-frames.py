#! /usr/bin/env python

import time
import numpy as np
from droplet_graphene_analysis.util import read_droplet_trajectory, elapsed_time

TRAJ_FILE = ['../../02-strained/s-1.00/run_prod/nvt_prod_5.lammpstrj',
             '../../02-strained/s-1.00/run_prod/nvt_prod_6.lammpstrj',
             '../../02-strained/s-1.00/run_prod/nvt_prod_7.lammpstrj',
             '../../02-strained/s-1.00/run_prod/nvt_prod_8.lammpstrj']

if __name__ == '__main__':

    time_start = time.time()
    cell_params, waters, carbons, _, shifts = read_droplet_trajectory(TRAJ_FILE, index=':', return_shift_trajectory=True)
    waters = waters[500:3700:800]
    carbons = carbons[500:3700:800]
    shifts = shifts[500:3700:800]
    print(f'Read trajectory in {elapsed_time(time_start)}')

    np.savez_compressed('frames.npz', cell_params=cell_params, waters=waters, carbons=carbons, shifts=shifts)

