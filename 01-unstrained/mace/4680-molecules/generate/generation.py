#! /usr/bin/env python

import shutil
import numpy as np
import ase
import ase.io.lammpsdata
from ase.build import make_supercell
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

system = ase.io.read('start.pdb')
carbon = system[[atom.index for atom in system if atom.symbol=='C']]
water = system[[atom.index for atom in system if atom.symbol=='O' or atom.symbol=='H']]

P = np.array([[12, 0, 0], [0, 12, 0], [0, 0, 1]])
super_carbon = make_supercell(carbon, P, tol=1e-5)

P = np.array([[6, 0, 0], [0, 6, 0], [0, 0, 1]])
super_water = make_supercell(water, P, tol=1e-5, wrap=False)

cell_length_x = system.cell.cellpar()[0]
cell_length_y = system.cell.cellpar()[1]
water_height = np.max(water.positions[:,2]) - np.min(water.positions[:,2])
super_water.positions += np.array([3 * cell_length_x, 3 * cell_length_y, 0])

system = super_carbon + super_water
MaxwellBoltzmannDistribution(system, temperature_K=300)
cell_length_x = system.cell.cellpar()[0]
cell_length_y = system.cell.cellpar()[1]
system.set_cell(np.array([[cell_length_x, 0, 0],
                          [0, cell_length_y, 0],
                          [0, 0, np.sqrt(cell_length_x * cell_length_y)]]), scale_atoms=False)

ase.io.lammpsdata.write_lammps_data('start.lammps-data', system, specorder=['C', 'H', 'O'],
                                    masses=True, velocities=True, units='metal', atom_style='full')
shutil.move('start.lammps-data', '../start.lammps-data')

