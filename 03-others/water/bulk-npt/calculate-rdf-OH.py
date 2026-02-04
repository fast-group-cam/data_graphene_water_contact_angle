#! /usr/bin/env python

import sys
import os
import ase.io
import numpy as np
import matplotlib.pyplot as plt

FILENAME = 'run/npt_prod_snapshots.lammpstrj'
OUTPUT_FILE = 'rdf-mace/rdf_OH_mace.dat'
N_BINS = 300
MAX_DIST = 10

if __name__ == '__main__':

    if not os.path.isfile(FILENAME):
        raise RuntimeError(f'File "{FILENAME}" not found!')
    if os.path.isfile(OUTPUT_FILE):
        os.remove(OUTPUT_FILE)

    traj = ase.io.read(FILENAME, index=':')
    OH_distances = list()
    volumes = list()
    N_frames = 0
    N_O = None
    N_H = None
    for atoms in traj:
        
        N_frames += 1
        cell_params = atoms.cell.cellpar()[0:3]
        volumes.append(np.prod(cell_params))
        if np.array_equal(np.unique(atoms.numbers) , (1,2)):
            hydrogens = atoms.positions[atoms.numbers == 1]
            oxygens = atoms.positions[atoms.numbers == 2]
        else:
            hydrogens = atoms.positions[atoms.numbers == 1]
            oxygens = atoms.positions[atoms.numbers == 8]

        if N_O is None:
            N_H = hydrogens.shape[0]
            N_O = oxygens.shape[0]

        disp = oxygens[:,None] - hydrogens[None,:]
        disp -= cell_params * np.round(disp / cell_params)
        dist = np.sqrt(np.sum(disp**2, axis=-1))
        OH_distances.extend(dist[dist < MAX_DIST])

    OH_distances = np.array(OH_distances)
    mean_vol = np.mean(volumes)
    pair_density = N_O * N_H / mean_vol

    r_bin_edges = np.linspace(0, MAX_DIST, N_BINS + 1)
    r_range = (r_bin_edges[:-1] + r_bin_edges[1:]) / 2
    volumes = 4 * np.pi * ((r_bin_edges[1:]**3) - (r_bin_edges[:-1]**3)) / 3

    rdf = np.histogram(OH_distances, r_bin_edges)[0] / (volumes * N_frames * pair_density)

    with open(OUTPUT_FILE, 'w') as stream:
        stream.write('# O-H RDF for MACE\n')
        stream.write('# r g(r)')
        for i in range(N_BINS):
            stream.write(f'\n{r_range[i]:.10f} {rdf[i]:.10f}')

    sys.exit()
