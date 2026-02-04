#! /usr/bin/env python

import shutil
import numpy as np
import ase
import ase.io.lammpsdata
from ase.build import make_supercell
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

system = ase.io.read('start.pdb')
water = system[[atom.index for atom in system if atom.symbol=='O' or atom.symbol=='H']]

height = 2.0 + np.max(water.positions[:,2]) - np.min(water.positions[:,2])

P = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
super_water = make_supercell(water, P, tol=1e-5, wrap=False)

MaxwellBoltzmannDistribution(super_water, temperature_K=300)
cell_length_x = super_water.cell.cellpar()[0]
cell_length_y = super_water.cell.cellpar()[1]
box_length = max(max(cell_length_x, cell_length_y), height)
super_water.set_cell(np.array([[box_length, 0, 0],
                               [0, box_length, 0],
                               [0, 0, box_length]]), scale_atoms=False)
super_water.positions = super_water.positions % np.array((box_length, box_length, box_length))

ase.io.lammpsdata.write_lammps_data('start.lammps-data', super_water, specorder=['H', 'O'],
                                    masses=True, velocities=True, units='metal', atom_style='full')
shutil.move('start.lammps-data', '../start.lammps-data')

