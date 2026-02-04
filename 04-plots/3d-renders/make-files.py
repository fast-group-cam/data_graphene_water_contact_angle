#! /usr/bin/env python

import os
import time
import numpy as np
import ase.io
from droplet_graphene_analysis.util import elapsed_time, read_droplet_trajectory
from droplet_graphene_analysis.util.interpolate import PeriodicGridInterpolator
from droplet_graphene_analysis.util.droplet import find_interface
from droplet_graphene_analysis.util.graphene import raw_heightmap
from droplet_graphene_analysis.util.graphene.sheet import generate_sheet, C_C_DISTANCE
from droplet_graphene_analysis.util.graphene.sheet import CUTOFF_RADIUS as CARBON_RADIUS

SOURCES = ['../../01-unstrained/mace/1000-molecules/run_prod/nvt_prod_$$$$.lammpstrj',
           '../../02-strained/s-2.00/run_prod/nvt_prod_$$$$.lammpstrj',
           '../../02-strained/s+2.00/run_prod/nvt_prod_$$$$.lammpstrj']

DESTINATIONS = ['free', 'compressed', 'stretched']

N_FRAMES_TARGET = 100
N_PTS = 300

def transfer(src, dest):

    source_files = list()
    for i in range(99):
        if os.path.isfile(src.replace('$$$$', str(i))):
            source_files.append(src.replace('$$$$', str(i)))
        else:
            break
    
    if len(source_files) < 1:
        raise RuntimeError(f'No files found under pattern "{src}"!')
    
    slice_step = max(int(100 * len(source_files) / N_FRAMES_TARGET), 1)
    cell_params, waters, carbons, _ = read_droplet_trajectory(source_files, index=f'::{slice_step}')
    N_frames = waters.shape[0]

    sheet_Nx = int(np.ceil(6.0 * cell_params[0] / CARBON_RADIUS))
    sheet_Ny = int(np.ceil(6.0 * cell_params[1] / CARBON_RADIUS))
    sheets = np.empty((N_frames, sheet_Nx, sheet_Ny), dtype=float)
    for i in range(N_frames):
        sheets[i] = raw_heightmap(carbons[i], cell_params[0:2], (sheet_Nx, sheet_Ny))
    sheets = np.mean(sheets, axis=0)
    heightmap = PeriodicGridInterpolator(cell_params[0:2], sheets)

    atoms = generate_sheet(cell_params[0] + C_C_DISTANCE, cell_params[1] + C_C_DISTANCE)
    atoms.positions[:,2] = heightmap(atoms.positions[:,0:2])
    atoms.positions += cell_params / 2.0
    ase.io.write(dest + '-ave.xyz', atoms)

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

    with open(dest + '.obj', 'w') as output_stream:
        output_stream.write('# List of geometric vertices\n')
        for v in vertices:
            v += cell_params / 2.0
            output_stream.write(f'v {v[0]} {v[1]} {v[2]}\n')
        output_stream.write('\n# List of polygonal face elements\n')
        for f in faces:
            output_stream.write(f'f {f[0]} {f[1]} {f[2]}\n')

if __name__ == '__main__':
    for src, dest in zip(SOURCES, DESTINATIONS):
        print(f'Transferring from "{src}" to "{dest}"...', end='')
        time_start = time.time()
        transfer(src, dest)
        print(f'done in {elapsed_time(time_start)}.')
    print('All done.')
