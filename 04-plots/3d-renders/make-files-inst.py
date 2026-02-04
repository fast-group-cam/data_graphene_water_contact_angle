#! /usr/bin/env python

import time
import numpy as np
import ase.io
from ase import Atoms
from droplet_graphene_analysis.util import elapsed_time, center_coordinates

SOURCES = ['../../01-unstrained/mace/1000-molecules/run_prod/nvt_prod_final.lammps-data',
           '../../02-strained/s-2.00/run_prod/nvt_prod_final.lammps-data',
           '../../02-strained/s+2.00/run_prod/nvt_prod_final.lammps-data']

DESTINATIONS = ['free-inst', 'compressed-inst', 'stretched-inst']

def transfer(src, dest):

    system = ase.io.read(src)
    if np.array_equal(np.unique(system.numbers), [1, 2, 3]):
        system.numbers[system.numbers == 1] = 6
        system.numbers[system.numbers == 2] = 1
        system.numbers[system.numbers == 3] = 8
    cell_params = system.cell.cellpar()[0:3]

    _, carbons, _ = center_coordinates(system, cell_params)
    atoms = Atoms(['C',] * carbons.shape[0], carbons, cell=cell_params, pbc=True)
    atoms.positions += cell_params / 2.0
    ase.io.write(dest + '.xyz', atoms)

if __name__ == '__main__':
    for src, dest in zip(SOURCES, DESTINATIONS):
        print(f'Transferring from "{src}" to "{dest}"...', end='')
        time_start = time.time()
        transfer(src, dest)
        print(f'done in {elapsed_time(time_start)}.')
    print('All done.')
