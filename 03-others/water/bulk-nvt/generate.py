#! /usr/bin/env python

import numpy as np
import ase
import ase.io.lammpsdata
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

CELL_LENGTH = 12.5
BOND_LENGTH = 1.0
BOND_ANGLE = 109.47 * (np.pi / 180.0)
N_WATER = 64

#oxygens = np.random.random((N_WATER, 3)) * CELL_LENGTH

oxygens = np.empty((N_WATER, 3), dtype=float)
for i in range(4):
    for j in range(4):
        for k in range(4):
            oxygens[(16 * i) + (4 * j) + k] = np.array([i, j, k], dtype=float) * CELL_LENGTH / 4

theta = np.random.random(N_WATER) * np.pi
phi = np.random.random(N_WATER) * 2 * np.pi
hydrogen_1 = np.c_[np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)] * BOND_LENGTH + oxygens
hydrogen_2 = np.c_[np.sin(theta + BOND_ANGLE) * np.cos(phi), np.sin(theta + BOND_ANGLE) * np.sin(phi), np.cos(theta + BOND_ANGLE)] * BOND_LENGTH + oxygens

system = ase.Atoms(symbols=((['O',] * N_WATER) + (['H', 'H'] * N_WATER)),
                   positions=np.concat((oxygens, hydrogen_1, hydrogen_2), axis=0),
                   cell=[CELL_LENGTH, CELL_LENGTH, CELL_LENGTH],
                   pbc=True)

MaxwellBoltzmannDistribution(system, temperature_K=300)

ase.io.lammpsdata.write_lammps_data('start.lammps-data', system, specorder=['H', 'O'],
                                    masses=True, velocities=True, units='metal', atom_style='full')

