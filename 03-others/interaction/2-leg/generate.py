import os
import argparse
import numpy as np
from droplet_graphene_analysis.util.graphene.sheet import generate_sheet

REF_C_C_DIST = 1.423
REF_O_H_DIST = 0.97
REF_H_O_H_ANG = 104.4

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--strain', type=float, default=0.0)
    args = parser.parse_args()

    root_folder = f's{args.strain:+.2f}'
    if not os.path.isdir(root_folder):
        os.mkdir(root_folder)

    cc_dist = REF_C_C_DIST * (1.0 + (args.strain / 100.0))
    
    system = generate_sheet(12.51 * cc_dist, 8.67 * cc_dist, interatomic_dist=cc_dist)
    cell_params = system.get_cell().cellpar()[0:3]
    carbons = np.mod(system.positions, cell_params)
    N_carbons = carbons.shape[0]
    x_O = 6.5 * cc_dist
    y_O = 1.5 * np.sqrt(3) * cc_dist
    OH_x = REF_O_H_DIST * np.sin((REF_H_O_H_ANG / 2.0) * np.pi / 180.0)
    OH_z = REF_O_H_DIST * np.cos((REF_H_O_H_ANG / 2.0) * np.pi / 180.0)

    for z_O in np.linspace(2.5, 4.5, 9):

        sub_folder = f'd{z_O:+.2f}'
        if not os.path.isdir(os.path.join(root_folder, sub_folder)):
            os.mkdir(os.path.join(root_folder, sub_folder))

        oxygens = np.array([[x_O, y_O, z_O],], dtype=float)
        hydrogens = np.array([[x_O + OH_x, y_O, z_O - OH_z], [x_O - OH_x, y_O, z_O - OH_z]], dtype=float)
        
        with open(os.path.join(root_folder, sub_folder, 'geometry.xyz'), 'w') as output_stream:
            output_stream.write(f'{N_carbons + 3}\n')
            output_stream.write(f'Lattice="{cell_params[0]} 0.0 0.0 0.0 {cell_params[1]} 0.0 0.0 0.0 {cell_params[2]}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
            for pos in carbons:
                output_stream.write(f'C       {pos[0]: 10.8f}      {pos[1]: 10.8f}      {pos[2]: 10.8f}\n')
            for pos in hydrogens:
                output_stream.write(f'H       {pos[0]: 10.8f}      {pos[1]: 10.8f}      {pos[2]: 10.8f}\n')
            for pos in oxygens:
                output_stream.write(f'O       {pos[0]: 10.8f}      {pos[1]: 10.8f}      {pos[2]: 10.8f}\n')

    