#! /usr/bin/env python

import sys
import warnings
import numpy as np
import ase.io

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from torch.cuda import is_available as cuda_is_available
from mace.calculators import MACECalculator

MACE_MODEL = '../../00-common/models/ch2o-dens-inv_swa.model'
STRAINS = [0.0, 1.0, 2.0]
DISTS = [2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5]

if __name__ == '__main__':

    device = 'cuda' if cuda_is_available() else 'cpu'
    mace_calc = MACECalculator(model_paths=MACE_MODEL, device=device)
    
    def read_cp2k_result(filename):
        result = None
        with open(filename, 'r') as input_stream:
            for line in input_stream:
                if line.startswith(' ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:'):
                    result = float(line.split()[-1]) * 27.2114
        if result is None:
            raise RuntimeError(f'Failed to read {filename}')
        return result
    
    def calculate_mace_result(filename):
        atoms = ase.io.read(filename)
        atoms.calc = mace_calc
        return atoms.get_potential_energy()

    cp2k_energies = np.empty((3, len(STRAINS), len(DISTS)), dtype=float)
    for i in range(3):
        for j, strain in enumerate(STRAINS):
            for k, dist in enumerate(DISTS):
                cp2k_energies[i,j,k] = read_cp2k_result(f'{i}-leg/s{strain:+.2f}/d{dist:+.2f}/cp2k_sub.out')

    mace_water_ref = calculate_mace_result('ref/water/geometry.xyz')
    mace_graphene_ref = np.empty(len(STRAINS), dtype=float)
    for i, strain in enumerate(STRAINS):
        mace_graphene_ref[i] = calculate_mace_result(f'ref/s{strain:+.2f}/geometry.xyz')
    mace_energies = np.empty((3, len(STRAINS), len(DISTS)), dtype=float)
    for i in range(3):
        for j, strain in enumerate(STRAINS):
            for k, dist in enumerate(DISTS):
                raw_energy = calculate_mace_result(f'{i}-leg/s{strain:+.2f}/d{dist:+.2f}/geometry.xyz')
                mace_energies[i,j,k] = raw_energy - mace_water_ref - mace_graphene_ref[j]

    np.savez_compressed('results.npz', strains=np.array(STRAINS), dists=np.array(DISTS), cp2k_energies=cp2k_energies, mace_energies=mace_energies)