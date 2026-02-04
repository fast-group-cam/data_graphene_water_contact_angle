#! /usr/bin/env python

import sys
import os
import ase.io
import numpy as np
import matplotlib.pyplot as plt

FILENAME = 'run/nvt_prod_snapshots.lammpstrj'
OUTPUT_FILE = 'rdf-mace/rdf_$$$_mace.dat'
N_BINS = 200
MAX_DIST = 6

if __name__ == '__main__':

    if not os.path.isfile(FILENAME):
        raise RuntimeError(f'File "{FILENAME}" not found!')
    if os.path.isfile(OUTPUT_FILE):
        os.remove(OUTPUT_FILE)

    traj = ase.io.read(FILENAME, index=':')
    OO_distances = list()
    OH_distances = list()
    HH_distances = list()
    cell_params = None
    need_to_reassign = None
    N_frames = 0
    N_O = None
    N_H = None
    for atoms in traj:
        
        N_frames += 1
        if cell_params is None:
            cell_params = atoms.cell.cellpar()[0:3]
        if need_to_reassign is None:
            need_to_reassign = np.array_equal(np.unique(atoms.numbers) , (1,2))
        
        hydrogens = atoms.positions[atoms.numbers == 1]
        oxygens = atoms.positions[atoms.numbers == (2 if need_to_reassign else 8)]

        if N_O is None:
            N_H = hydrogens.shape[0]
            N_O = oxygens.shape[0]

        for i, pos in enumerate(oxygens):
            disp = oxygens[i+1:] - pos
            disp -= cell_params * np.round(disp / cell_params)
            dist = np.sqrt(np.sum(disp**2, axis=-1))
            OO_distances.extend(dist[dist < MAX_DIST])
        
        disp = oxygens[:,None] - hydrogens[None,:]
        disp -= cell_params * np.round(disp / cell_params)
        dist = np.sqrt(np.sum(disp**2, axis=-1))
        OH_distances.extend(dist[dist < MAX_DIST])

        for i, pos in enumerate(hydrogens):
            disp = hydrogens[i+1:] - pos
            disp -= cell_params * np.round(disp / cell_params)
            dist = np.sqrt(np.sum(disp**2, axis=-1))
            HH_distances.extend(dist[dist < MAX_DIST])

    OO_distances = np.array(OO_distances)
    OH_distances = np.array(OH_distances)
    HH_distances = np.array(HH_distances)
    mean_vol = np.prod(cell_params)

    r_bin_edges = np.linspace(0, MAX_DIST, N_BINS + 1)
    r_range = (r_bin_edges[:-1] + r_bin_edges[1:]) / 2
    volumes = 4 * np.pi * ((r_bin_edges[1:]**3) - (r_bin_edges[:-1]**3)) / 3

    pair_density = N_O * (N_O - 1) / (2 * mean_vol)
    rdf = np.histogram(OO_distances, r_bin_edges)[0] / (volumes * N_frames * pair_density)
    with open(OUTPUT_FILE.replace('$$$', 'OO'), 'w') as stream:
        stream.write('# O-O RDF for MACE\n')
        stream.write('# r g(r)')
        for i in range(N_BINS):
            stream.write(f'\n{r_range[i]:.10f} {rdf[i]:.10f}')

    pair_density = N_O * N_H / mean_vol
    rdf = np.histogram(OH_distances, r_bin_edges)[0] / (volumes * N_frames * pair_density)
    with open(OUTPUT_FILE.replace('$$$', 'OH'), 'w') as stream:
        stream.write('# O-H RDF for MACE\n')
        stream.write('# r g(r)')
        for i in range(N_BINS):
            stream.write(f'\n{r_range[i]:.10f} {rdf[i]:.10f}')

    pair_density = N_H * (N_H - 1) / (2 * mean_vol)
    rdf = np.histogram(HH_distances, r_bin_edges)[0] / (volumes * N_frames * pair_density)
    with open(OUTPUT_FILE.replace('$$$', 'HH'), 'w') as stream:
        stream.write('# H-H RDF for MACE\n')
        stream.write('# r g(r)')
        for i in range(N_BINS):
            stream.write(f'\n{r_range[i]:.10f} {rdf[i]:.10f}')

    sys.exit()
