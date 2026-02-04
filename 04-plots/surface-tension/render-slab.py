#! /usr/bin/env python

import numpy as np
import ase.io
from ase import Atoms

INPUT_FILE = '../../03-others/water/slab-520/equi_final.lammps-data'
OUTPUT_FILE = 'render.xyz'

if __name__ == '__main__':

    atoms = ase.io.read(INPUT_FILE)
    need_to_reassign = np.array_equal(np.unique(atoms.numbers), (1, 2))
    hydrogens = atoms.positions[atoms.numbers == 1]
    oxygens = atoms.positions[atoms.numbers == (2 if need_to_reassign else 8)]
    cell_params = atoms.cell.cellpar()[0:3]

    hydrogens -= cell_params * np.round(hydrogens / cell_params)
    oxygens -= cell_params * np.round(oxygens / cell_params)

    hydrogens += 0.5 * cell_params
    oxygens += 0.5 * cell_params

    NH = hydrogens.shape[0]
    NO = oxygens.shape[0]

    dest = Atoms(['O',] * NO + ['H',] * NH, np.concat((oxygens, hydrogens), axis=0), cell=cell_params, pbc=True)
    ase.io.write(OUTPUT_FILE, dest)

