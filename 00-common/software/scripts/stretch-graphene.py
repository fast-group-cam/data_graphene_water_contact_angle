#! /usr/bin/env python

#==================================================================================================
# Python script for stretching the graphene sheet in a file. Use as:
#
#     python stretch-graphene.py <input_file> -s <strain> -o <output_file>
#
#==================================================================================================

import os
import argparse
import numpy as np
import ase.io
import ase.io.lammpsdata

if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(prog='stretch-graphene.py')
    parser.add_argument('input_file')
    parser.add_argument('-s', '--strain', default=0.0, type=float)
    parser.add_argument('-o', '--output', default=None)
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        raise RuntimeError(f'File "{args.input_file}" not found!')
    if args.output is None:
        raise RuntimeError('Output file must be specified!')

    # Read input
    system = ase.io.read(args.input_file)

    # Stretch graphene
    cell = system.cell.cellpar()[0:3]
    factor = 1.0 + args.strain
    system.positions[system.symbols=='C',0] *= factor
    system.positions[system.symbols=='C',1] *= factor
    system.set_cell(np.array([[factor * cell[0], 0, 0],
                              [0, factor * cell[1], 0],
                              [0, 0, cell[2]]]), scale_atoms=False)

    # Write system to output
    ase.io.lammpsdata.write_lammps_data(args.output, system, specorder=['C', 'H', 'O'],
                                        masses=True, velocities=True, units='metal',
                                        atom_style='full')
