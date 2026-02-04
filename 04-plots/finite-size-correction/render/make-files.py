#! /usr/bin/env python

import os
import numpy as np
import ase.io
from ase import Atoms
from droplet_graphene_analysis.util import center_coordinates, read_droplet_trajectory
from droplet_graphene_analysis.util.droplet.coarse_grain import find_interface

SOURCE_FILE = '../../../01-unstrained/mace/4680-molecules/run_prod/nvt_prod_final.lammps-data'
SOURCE_TRAJ = '../../../01-unstrained/mace/4680-molecules/run_prod/nvt_prod_$$$$.lammpstrj'
N_PTS = 300

if (not os.path.isfile('render_carbons.xyz')) or (not os.path.isfile('render_water.xyz')):

    print('Generating xyz files for atomic coordinates...')

    src_atoms = ase.io.read(SOURCE_FILE)
    if np.array_equal(np.unique(src_atoms.numbers), [1, 2, 3]):
        src_atoms.numbers[src_atoms.numbers == 1] = 6
        src_atoms.numbers[src_atoms.numbers == 2] = 1
        src_atoms.numbers[src_atoms.numbers == 3] = 8
    cell_params = src_atoms.get_cell().cellpar()[0:3]
    src_oxy, src_car, src_hyd = center_coordinates(src_atoms, cell_params)

    car_copy0 = src_car + np.array([-cell_params[0], 0.0, 0.0])
    car_copy1 = src_car + np.array([cell_params[0], 0.0, 0.0])
    car_copy2 = src_car + np.array([0.0, cell_params[1], 0.0])
    car_copy3 = src_car + np.array([-cell_params[0], cell_params[1], 0.0])
    car_copy4 = src_car + np.array([cell_params[0], cell_params[1], 0.0])
    dest_car = np.concat((src_car, car_copy0, car_copy1, car_copy2, car_copy3, car_copy4), axis=0)
    dest_car = dest_car[dest_car[:,1] >= 0.0]

    carbons = Atoms(['C',] * dest_car.shape[0], dest_car)
    carbons.positions += cell_params / 2.0
    ase.io.write('render_carbons.xyz', carbons)

    dest_oxy = src_oxy[src_oxy[:,1] >= 0.0]
    pot_Hs = np.copy(src_hyd)
    dest_hyd = list()
    for pos in dest_oxy:
        dists = np.sum((pot_Hs - pos)**2, axis=-1)
        ind = np.argmin(dists)
        dest_hyd.append(pot_Hs[ind])
        pot_Hs = np.concat((pot_Hs[:ind], pot_Hs[ind+1:]), axis=0)
        dists = np.sum((pot_Hs - pos)**2, axis=-1)
        ind = np.argmin(dists)
        dest_hyd.append(pot_Hs[ind])
        pot_Hs = np.concat((pot_Hs[:ind], pot_Hs[ind+1:]), axis=0)
    dest_hyd = np.array(dest_hyd)

    hydrogens = Atoms(['H',] * dest_hyd.shape[0], dest_hyd)
    oxygens = Atoms(['O',] * dest_oxy.shape[0], dest_oxy)
    dest_atoms = hydrogens + oxygens
    dest_atoms.positions += cell_params / 2.0
    ase.io.write('render_water.xyz', dest_atoms)

    print('...generated xyz files for atomic coordinates.')

if not os.path.isfile('render_surface_ave.obj'):

    print('Generating obj file for time-averaged interface...')

    traj_files = list()
    for i in range(99):
        target = SOURCE_TRAJ.replace('$$$$', str(i))
        if os.path.isfile(target):
            traj_files.append(target)
    cell_params, waters, _, _ = read_droplet_trajectory(traj_files, index='::200')

    print('...read trajectory...')

    vertices = list()
    faces = list()

    polar_angles = np.linspace(0.0, np.pi, N_PTS + 2, endpoint=True)[1:-1]
    azimuthal_angles = np.linspace(0.0, np.pi, N_PTS, endpoint=True)

    CoM = np.mean(waters, axis=(0,1))
    vertices.append(find_interface(waters, CoM, axis=(0.0, 0.0, 1.0)))
    for theta in polar_angles:
        for phi in azimuthal_angles:
            axis = (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta))
            vertices.append(find_interface(waters, CoM, axis=axis))
    vertices.append(find_interface(waters, CoM, axis=(0.0, 0.0, -1.0)))

    for i in range(N_PTS - 1):
        faces.append([1, i + 2, i + 3])
    for i in range(N_PTS - 1):
        for j in range(N_PTS - 1):
            faces.append([(N_PTS * i) + 2 + j, (N_PTS * i) + 2 + N_PTS + j, (N_PTS * i) + 3 + j])
            faces.append([(N_PTS * i) + 3 + j, (N_PTS * i) + 2 + N_PTS + j, (N_PTS * i) + 3 + N_PTS + j])
    for i in range(N_PTS - 1):
        faces.append([(N_PTS**2) - N_PTS + 2 + i, (N_PTS**2) + 2, (N_PTS**2) - N_PTS + 3 + i])

    with open('render_surface_ave.obj', 'w') as output_stream:
        output_stream.write('# List of geometric vertices\n')
        for v in vertices:
            v += cell_params / 2.0
            output_stream.write(f'v {v[0]} {v[1]} {v[2]}\n')
        output_stream.write('\n# List of polygonal face elements\n')
        for f in faces:
            output_stream.write(f'f {f[0]} {f[1]} {f[2]}\n')

    print('...generated obj file for time-averaged interface.')

