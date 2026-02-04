#! /usr/bin/env python

#==================================================================================================
# Python script for flattening the graphene sheet in a file. Use as:
#
#     python flatten-graphene.py <input_file> -o <output_file>
#
#==================================================================================================

import os
import argparse
import numpy as np
import ase.io
import ase.io.lammpsdata
from droplet_graphene_analysis.util.graphene import generate_sheet

if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(prog='delete-waters.py')
    parser.add_argument('input_file')
    parser.add_argument('-o', '--output', default=None)
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        raise RuntimeError(f'File "{args.input_file}" not found!')
    if args.output is None:
        raise RuntimeError('Output file must be specified!')

    # Read input and separate system by elements
    system = ase.io.read(args.input_file)
    hydrogen = system[system.symbols=='H']
    oxygen = system[system.symbols=='O']
    current_count = len(oxygen)

    # Flatten graphene
    cell = system.cell.cellpar()[0:2]
    carbon = generate_sheet(cell[0], cell[1], origin='corner')

    # Shift water to 1A above the graphene
    lowest_z = np.min(hydrogen.positions[:,2])
    lowest_z = min(lowest_z, np.min(oxygen.positions[:,2]))
    hydrogen.positions[:,2] -= (lowest_z - 1.0)
    oxygen.positions[:,2] -= (lowest_z - 1.0)

    # Write system to output
    system = carbon + hydrogen + oxygen
    ase.io.lammpsdata.write_lammps_data(args.output, system, specorder=['C', 'H', 'O'],
                                        masses=True, velocities=True, units='metal',
                                        atom_style='full')
