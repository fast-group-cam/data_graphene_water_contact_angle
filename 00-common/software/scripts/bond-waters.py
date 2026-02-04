#! /usr/bin/env python

#==================================================================================================
# Python script for adding bonds to water molecules in a file. Use as:
#
#     python bond-waters.py <input_file> -o <output_file> -N <target_number>
#
#==================================================================================================

import os
import argparse
import numpy as np
import ase.io
import ase.io.lammpsdata

R_OH = 1.0
THETA = 109.47
COS_T = np.cos(THETA * np.pi / 180)
SIN_T = np.sin(THETA * np.pi / 180)

if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(prog='bond-waters.py')
    parser.add_argument('input_file')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('--fix', action='store_true')
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        raise RuntimeError(f'File "{args.input_file}" not found!')
    if args.output is None:
        raise RuntimeError('Output file must be specified!')

    # Read input and track which indices need to be bonded
    system = ase.io.read(args.input_file)
    N_atoms = len(system)
    N_oxy = np.sum(system.symbols == 'O')
    hydrogen_positions = system.positions[system.symbols=='H']
    bondsi = list()
    for i, atom in enumerate(system):
        if atom.symbol == 'O':

            # Find closest two hydrogen atoms
            hydrogen_dists = np.sum(np.square(hydrogen_positions - atom.position), axis=-1)
            target_positions = hydrogen_positions[np.argpartition(hydrogen_dists, 1)[0:2]]
            match = np.all(np.equal(system.positions, np.broadcast_to(target_positions[0], (N_atoms, 3))), axis=-1)
            j = np.argwhere(match)[0][0]
            match = np.all(np.equal(system.positions, np.broadcast_to(target_positions[1], (N_atoms, 3))), axis=-1)
            k = np.argwhere(match)[0][0]

            # If necessary, fix lengths & angles
            if args.fix:
                disp_0 = target_positions[0] - atom.position
                axis_0 = disp_0 / np.linalg.norm(disp_0)
                system.positions[j] = atom.position + (R_OH * axis_0)
                disp_1 = target_positions[1] - atom.position
                axis_1 = np.cross(axis_0, disp_1)
                axis_1 /= np.linalg.norm(axis_1)
                axis_2 = np.cross(axis_1, axis_0)
                system.positions[k] = atom.position + (R_OH * axis_0 * COS_T) + (R_OH * axis_2 * SIN_T)

            bondsi.append(f'{j}(1),{k}(1)')

        else:
            bondsi.append('_')

    # Add bonds and write system to output
    system.arrays['bonds'] = bondsi
    ase.io.lammpsdata.write_lammps_data(args.output + '_tmp', system, specorder=['C', 'H', 'O'],
                                        masses=True, velocities=True, units='metal',
                                        bonds=True, atom_style='full')
    
    # Add angles to output
    counter = 1
    with open(args.output + '_tmp', 'a') as dest:
        dest.write('\nAngles\n')
        for i, bond in enumerate(bondsi):
            if bond != '_':
                left, right = bond.split(',')
                j = int(left.removesuffix('(1)'))
                k = int(right.removesuffix('(1)'))
                dest.write(f'\n{counter} 1 {j+1} {i+1} {k+1}')
                counter += 1
    counter -= 1

    # Add angle count to output
    with open(args.output + '_tmp', 'r') as src:
        with open(args.output, 'w') as dest:
            for _ in range(7):
                dest.write(src.readline())
            dest.write(f'{counter} angles\n')
            dest.write('1 angle types\n')
            for line in src.readlines():
                dest.write(line)
    os.remove(args.output + '_tmp')
