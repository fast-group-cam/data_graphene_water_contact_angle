#! /usr/bin/env python

import os
import time
import numpy as np
import ase.io
from ase import Atoms
from droplet_graphene_analysis.util import elapsed_time, center_coordinates

SOURCES = {'01-droplet-dynamic-graphene.xyz':    ('../../01-unstrained/mace/4680-molecules/run_prod/nvt_prod_final.lammps-data', True, True),
           '02-droplet-flat-graphene.xyz':       ('../../01-unstrained/mace/fixed-4680-molecules/run_prod/nvt_prod_final.lammps-data', True, True),
           '03-droplet-FF-flat-graphene.xyz':    ('../../01-unstrained/spce/4680-molecules/run_prod/nvt_prod_final.lammps-data', True, True),
           '04-droplet-stretched-graphene.xyz':  ('../../02-strained/s+2.00/run_prod/nvt_prod_final.lammps-data', True, True),
           '05-droplet-compressed-graphene.xyz': ('../../02-strained/s-2.00/run_prod/nvt_prod_final.lammps-data', True, True),
           '06-dry-graphene.xyz':                ('../../03-others/graphene_dry/run/nvt_prod_final.lammps-data', True, True),
           '07-wet-graphene.xyz':                ('../../03-others/graphene_wet/run/nvt_prod_final.lammps-data', True, True),
           '08-bulk-water-nvt.xyz':              ('../../03-others/water/bulk-nvt/run/nvt_prod_final.lammps-data', False, False),
           '09-bulk-water-npt.xyz':              ('../../03-others/water/bulk-npt/run/npt_prod_final.lammps-data', False, False),
           '10-water-slab.xyz':                  ('../../03-others/water/slab-4680/run/nvt_prod_final_1.lammps-data', True, False)}

if __name__ == '__main__':

    for dest, (src, perform_shift, semi_shift) in SOURCES.items():

        if os.path.isfile(dest):
            print(f'"{dest}" already exists, skipped')
        else:
            if os.path.isfile(src):
                print(f'Transferring from "{src}" to "{dest}"...', end='')
                time_start = time.time()
                atoms = ase.io.read(src)
                if np.array_equal(np.unique(atoms.numbers), [1, 2, 3]):
                    atoms.numbers[atoms.numbers == 1] = 6
                    atoms.numbers[atoms.numbers == 2] = 1
                    atoms.numbers[atoms.numbers == 3] = 8
                elif np.array_equal(np.unique(atoms.numbers), [1, 2]):
                    atoms.numbers[atoms.numbers == 2] = 8
                cell_params = atoms.cell.cellpar()[0:3]
                if perform_shift:
                    oxygens, carbons, hydrogens = center_coordinates(atoms, cell_params)
                    shift = np.array((0.0, 0.0, -0.25 * cell_params[2])) if semi_shift else np.zeros(3)
                    oxygens += (0.5 * cell_params) + shift
                    carbons += (0.5 * cell_params) + shift
                    hydrogens += (0.5 * cell_params) + shift
                else:
                    oxygens = atoms.positions[atoms.numbers == 8]
                    carbons = atoms.positions[atoms.numbers == 6]
                    hydrogens = atoms.positions[atoms.numbers == 1]
                    oxygens = np.mod(oxygens, cell_params)
                    carbons = np.mod(carbons, cell_params)
                    hydrogens = np.mod(hydrogens, cell_params)
                N_oxy = oxygens.shape[0]
                N_car = carbons.shape[0]
                N_hyd = hydrogens.shape[0]
                atoms = Atoms(symbols=(['C',] * N_car + ['H',] * N_hyd + ['O',] * N_oxy),
                              positions=np.concat((carbons, hydrogens, oxygens), axis=0),
                              cell=cell_params, pbc=True)
                ase.io.write(dest, atoms)
                print(f'done in {elapsed_time(time_start)}.')
            else:
                print(f'"{src}" not found, unable to transfer to "{dest}"!')

    print('All done.\n')
