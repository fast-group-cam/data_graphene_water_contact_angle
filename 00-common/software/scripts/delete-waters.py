#! /usr/bin/env python

#==================================================================================================
# Python script for deleting water molecules from a file. Use as:
#
#     python delete-waters.py <input_file> -o <output_file> -N <target_number>
#
#==================================================================================================

import os
import argparse
import numpy as np
import ase.io
import ase.io.lammpsdata

if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(prog='delete-waters.py')
    parser.add_argument('input_file')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-N', '--N_target', type=int, default=None)
    parser.add_argument('--N_remove', type=int, default=None)
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        raise RuntimeError(f'File "{args.input_file}" not found!')
    if args.output is None:
        raise RuntimeError('Output file must be specified!')

    # Read input and separate system by elements
    system = ase.io.read(args.input_file)
    carbon = system[system.symbols=='C']
    hydrogen = system[system.symbols=='H']
    oxygen = system[system.symbols=='O']
    current_count = len(oxygen)

    # Check how many molecules should be removed
    if args.N_target is not None:
        if args.N_remove is not None:
            raise RuntimeError('Only one of --N_target or --N_remove should be specified!')
        if args.N_target < 1:
            raise RuntimeError('Target number of water molecules must be at least 1!')
        if current_count < args.N_target:
            raise RuntimeError(f'Existing number of water molecules ({current_count}) already ' +
                               f'less than target number ({args.N_target})!')
        N_to_remove = current_count - args.N_target
    elif args.N_remove is not None:
        if args.N_remove < 0:
            raise RuntimeError('Number of water molecules to remove must be positive!')
        N_to_remove = args.N_remove
    else:
        raise RuntimeError('Either --N_target or --N_remove must be specified!')

    # Remove random water molecules
    for i in range(N_to_remove):

        # Select random oxygen atom and remove it
        idx_to_remove = np.random.randint(0, len(oxygen))
        remaining_idxs = list(range(0, idx_to_remove)) + list(range(idx_to_remove + 1, len(oxygen)))
        pos_to_remove = oxygen[idx_to_remove].position
        oxygen = oxygen[remaining_idxs]

        # Remove two closest hydrogen atoms to the removed oxygen atom
        for j in range(2):
            hydrogen_dists = np.sum(np.square(hydrogen.positions - pos_to_remove), axis=-1)
            idx_to_remove = np.argmin(hydrogen_dists)
            remaining_idxs = list(range(0, idx_to_remove)) + list(range(idx_to_remove + 1, len(hydrogen)))
            hydrogen = hydrogen[remaining_idxs]

    # Write system to output
    system = carbon + hydrogen + oxygen
    ase.io.lammpsdata.write_lammps_data(args.output, system, specorder=['C', 'H', 'O'],
                                        masses=True, velocities=True, units='metal',
                                        atom_style='full')
