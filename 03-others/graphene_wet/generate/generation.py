#! /usr/bin/env python

import shutil
import numpy as np
import ase
import ase.io.lammpsdata
from ase.build import make_supercell
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

system = ase.io.read('start.pdb')
carbon = system.positions[system.symbols=='C']

P = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 1]])
supersys = make_supercell(system, P, tol=1e-5, wrap=False)

cell_length_x = supersys.cell.cellpar()[0]
cell_length_y = supersys.cell.cellpar()[1]

supersys_carbon = supersys[[atom.index for atom in supersys if atom.symbol=='C']]
mean_z_coord = np.mean(supersys_carbon.positions[:,2])

MaxwellBoltzmannDistribution(supersys, temperature_K=300)
supersys.positions -= np.array([0.0, 0.0, mean_z_coord])
supersys.set_cell(np.array([[cell_length_x, 0, 0],
                            [0, cell_length_y, 0],
                            [0, 0, 100.0]]), scale_atoms=False)

ase.io.lammpsdata.write_lammps_data('start.lammps-data', supersys, specorder=['C', 'H', 'O'],
                                    masses=True, velocities=True, units='metal', atom_style='full')
shutil.move('start.lammps-data', '../start.lammps-data')

