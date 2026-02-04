#! /usr/bin/env python

import os
import numpy as np
import ase.io
import time
from scipy.optimize import curve_fit
from ase import Atoms
from droplet_graphene_analysis.util import read_droplet_trajectory, center_coordinates, elapsed_time
from droplet_graphene_analysis.util.interpolate import PeriodicGridInterpolator
from droplet_graphene_analysis.util.droplet.coarse_grain import find_interface
from droplet_graphene_analysis.util.droplet.contact_angle import find_spherical_cap
from droplet_graphene_analysis.util.graphene.sheet import raw_heightmap, CUTOFF_RADIUS
from droplet_graphene_analysis.util.graphene.angle import calc_inclination_angles

#==================================================================================================

SOURCE_TRAJ =     {'s+0.00': '../../01-unstrained/mace/1000-molecules/run_prod/nvt_prod_$$$$.lammpstrj',
                   's+0.20': '../../02-strained/s+0.20/run_prod/nvt_prod_$$$$.lammpstrj',
                   's+0.50': '../../02-strained/s+0.50/run_prod/nvt_prod_$$$$.lammpstrj',
                   's+2.00': '../../02-strained/s+2.00/run_prod/nvt_prod_$$$$.lammpstrj'}
SOURCE_SNAPSHOT = {'s+0.00': '../../01-unstrained/mace/1000-molecules/run_prod/nvt_prod_final.lammps-data',
                   's+0.20': '../../02-strained/s+0.20/run_prod/nvt_prod_final.lammps-data',
                   's+0.50': '../../02-strained/s+0.50/run_prod/nvt_prod_final.lammps-data',
                   's+2.00': '../../02-strained/s+2.00/run_prod/nvt_prod_final.lammps-data'}

N_PTS = 300
INDEX_OBJ = '::20'
SLICE_WIDTH = 2.0

#==================================================================================================

if not os.path.isdir('data'):
    os.mkdir('data')

def exp_curve(x, A, k, c):
    return A * np.exp(-k * x) + c

for key in SOURCE_TRAJ.keys():

    #----------------------------------------------------------------------------------------------

    output = os.path.join('data', f'sheet_{key}.npz')
    if not os.path.isfile(output):

        time_start_0 = time.time()
        print(f'Generating sheet npz files for {key} configuration...')

        traj_files = list()
        for i in range(99):
            target = SOURCE_TRAJ[key].replace('$$$$', str(i))
            if os.path.isfile(target):
                traj_files.append(target)
        cell_params, waters, carbons, _ = read_droplet_trajectory(traj_files, index=':')
        N_frames = waters.shape[0]

        print(f'    ...read trajectory in {elapsed_time(time_start_0)}...')
        time_start_1 = time.time()

        sheet_Nx = int(np.ceil(6.0 * cell_params[0] / CUTOFF_RADIUS))
        sheet_Ny = int(np.ceil(6.0 * cell_params[1] / CUTOFF_RADIUS))
        sheets = np.empty((N_frames, sheet_Nx, sheet_Ny), dtype=float)
        for f in range(N_frames):
            sheets[f] = raw_heightmap(carbons[f], cell_params[0:2], (sheet_Nx, sheet_Ny))
        mean_heightmap = PeriodicGridInterpolator(cell_params[0:2], np.mean(sheets, axis=0))
        std_heightmap = PeriodicGridInterpolator(cell_params[0:2], np.std(sheets, axis=0))
        incl_angles = calc_inclination_angles(carbons, cell_params[0:2], (sheet_Nx, sheet_Ny))
        incl_mean_map = PeriodicGridInterpolator(cell_params[0:2], np.mean(incl_angles, axis=0))
        incl_std_map = PeriodicGridInterpolator(cell_params[0:2], np.std(incl_angles, axis=0))

        autocorr = np.zeros((50, sheet_Nx, sheet_Ny), dtype=float)
        autocorr[0] = np.mean(np.square(incl_angles), axis=0)
        for tau in range(1, 50):
            autocorr[tau] = np.mean(incl_angles[:-tau] * incl_angles[tau:], axis=0)
        autocorr[:] /= autocorr[0]

        tau = np.array(range(50))
        incl_norm_autocorr = np.zeros((sheet_Nx, sheet_Ny), dtype=float)
        for i in range(sheet_Nx):
            for j in range(sheet_Ny):
                try:
                    popt, _ = curve_fit(exp_curve, tau, autocorr[:,i,j], p0=(0.2, 0.5, 0.8))
                    incl_norm_autocorr[i,j] = max(popt[-1], 0.0)
                except RuntimeError:
                    incl_norm_autocorr[i,j] = 0.0
        incl_norm_autocorr_map = PeriodicGridInterpolator(cell_params[0:2], incl_norm_autocorr)

        print(f'    ...calculated graphene sheet autocorrelations in {elapsed_time(time_start_1)}...')
        time_start_1 = time.time()

        spherical_cap = find_spherical_cap(waters, cell_params, mean_heightmap)
        sphere_a = spherical_cap['a']
        with open(os.path.join('data', f'log_{key}.txt'), 'w') as output_stream:
            azi = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False)
            test_pts = np.c_[sphere_a * np.cos(azi), sphere_a * np.sin(azi)]
            output_stream.write(f'{key}, a = {sphere_a} [A]\n')
            output_stream.write(f'{key}, <z>(a) = {np.mean(mean_heightmap(test_pts))} \u00b1 {np.mean(std_heightmap(test_pts))}\n')
            output_stream.write(f'{key}, <z>(2a) = {np.mean(mean_heightmap(2.0 * test_pts))} \u00b1 {np.mean(std_heightmap(2.0 * test_pts))}\n')
            output_stream.write(f'{key}, <\u0398>(a) = {np.mean(incl_mean_map(test_pts))} \u00b1 {np.mean(incl_std_map(test_pts))}\n')
            output_stream.write(f'{key}, <\u0398>(2a) = {np.mean(incl_mean_map(2.0 * test_pts))} \u00b1 {np.mean(incl_std_map(2.0 * test_pts))}\n')
            output_stream.write(f'{key}, C\u0398(a) = {np.mean(incl_norm_autocorr_map(test_pts))}\n')
            output_stream.write(f'{key}, C\u0398(2a) = {np.mean(incl_norm_autocorr_map(2.0 * test_pts))}\n')

        print(f'    ...calculated droplet radius in {elapsed_time(time_start_1)}...')

        np.savez_compressed(output, cell_params=cell_params, mean_sheet=np.mean(sheets, axis=0), incl_angles=np.mean(incl_angles, axis=0), incl_norm_autocorr=incl_norm_autocorr)

        print(f'...generated sheet npz files for {key} configuration in {elapsed_time(time_start_0)}.')
        time_start_1 = time.time()

        interface_points = list()
        CoM = np.mean(waters, axis=(0,1))
        CoM[1] = 0.0
        azi = np.linspace(0.0, np.pi, 180)
        axes = np.c_[-np.cos(azi), np.zeros((180,)), -np.sin(azi)]
        for axis in axes:
            interface_points.append(find_interface(waters, CoM, axis))
        np.save(os.path.join('data', f'interface_slice_{key}.npy'), np.array(interface_points))
        print(f'...generated sheet npz files for {key} configuration in {elapsed_time(time_start_1)}.\n')

    #----------------------------------------------------------------------------------------------
    
    output = os.path.join('data', f'interface_{key}.obj')
    if not os.path.isfile(output):

        time_start_0 = time.time()
        print(f'Generating obj file for time-averaged interface for {key} configuration...')

        traj_files = list()
        for i in range(99):
            target = SOURCE_TRAJ[key].replace('$$$$', str(i))
            if os.path.isfile(target):
                traj_files.append(target)
        cell_params, waters, _, _ = read_droplet_trajectory(traj_files, index=INDEX_OBJ)

        print(f'    ...read trajectory in {elapsed_time(time_start_0)}...')

        CoM = np.mean(waters, axis=(0,1))
        vertices = list()
        faces = list()

        polar_angles = np.linspace(0.0, np.pi, N_PTS + 2, endpoint=True)[1:-1]
        azimuthal_angles = np.linspace(0.0, 2.0 * np.pi, N_PTS, endpoint=False)

        vertices.append(find_interface(waters, CoM, axis=(0.0, 0.0, 1.0)))
        for theta in polar_angles:
            for phi in azimuthal_angles:
                axis = (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta))
                vertices.append(find_interface(waters, CoM, axis=axis))
        vertices.append(find_interface(waters, CoM, axis=(0.0, 0.0, -1.0)))

        for i in range(N_PTS - 1):
            faces.append([1, i + 2, i + 3])
        faces.append([1, N_PTS + 1, 2])
        for i in range(N_PTS - 1):
            for j in range(N_PTS - 1):
                faces.append([(N_PTS * i) + 2 + j, (N_PTS * i) + 2 + N_PTS + j, (N_PTS * i) + 3 + j])
                faces.append([(N_PTS * i) + 3 + j, (N_PTS * i) + 2 + N_PTS + j, (N_PTS * i) + 3 + N_PTS + j])
            faces.append([(N_PTS * i) + 1 + N_PTS, (N_PTS * (i + 2)) + 1, (N_PTS * i) + 2])
            faces.append([(N_PTS * i) + 2, (N_PTS * (i + 2)) + 1, (N_PTS * (i + 1)) + 2])
        for i in range(N_PTS - 1):
            faces.append([(N_PTS**2) - N_PTS + 2 + i, (N_PTS**2) + 2, (N_PTS**2) - N_PTS + 3 + i])
        faces.append([(N_PTS**2) + 1, (N_PTS**2) + 2, (N_PTS**2) - N_PTS + 2])

        with open(output, 'w') as output_stream:
            output_stream.write('# List of geometric vertices\n')
            for v in vertices:
                output_stream.write(f'v {v[0]} {v[1]} {v[2]}\n')
            output_stream.write('\n# List of polygonal face elements\n')
            for f in faces:
                output_stream.write(f'f {f[0]} {f[1]} {f[2]}\n')

        print(f'...generated obj file for time-averaged interface in {elapsed_time(time_start_0)}.\n')

    #----------------------------------------------------------------------------------------------

    output = os.path.join('data', f'water_slice_{key}.xyz')
    if not os.path.isfile(output):

        time_start_0 = time.time()
        print(f'Generating xyz files for water slice of {key} configuration...')

        system = ase.io.read(SOURCE_SNAPSHOT[key])
        cell_params = system.cell.cellpar()[0:3]
        src_O, _, src_H = center_coordinates(system, cell_params)

        dest_O = src_O[np.abs(src_O[:,1]) < SLICE_WIDTH]
        dest_H = list()
        for oxygen in dest_O:
            dists = np.sum((src_H - oxygen)**2, axis=1)
            idx = np.argmin(dists)
            dest_H.append(src_H[idx])
            src_H = np.concat((src_H[:idx], src_H[(idx+1):]))
            dists = np.sum((src_H - oxygen)**2, axis=1)
            idx = np.argmin(dists)
            dest_H.append(src_H[idx])
            src_H = np.concat((src_H[:idx], src_H[(idx+1):]))
        dest_H = np.array(dest_H)

        atoms_O = Atoms(['O',] * dest_O.shape[0], dest_O, cell=cell_params, pbc=True)
        atoms_H = Atoms(['H',] * dest_H.shape[0], dest_H, cell=cell_params, pbc=True)
        ase.io.write(output, atoms_O + atoms_H)
        print(f'...generated xyz files for water slice of {key} configuration in {elapsed_time(time_start_0)}.\n')


