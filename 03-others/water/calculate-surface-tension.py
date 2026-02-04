#! /usr/bin/env python

import sys
import os
import argparse
import yaml
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(prog='calculate-surface-tension.py')
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', default='surface-tension.txt')
    parser.add_argument('-L', '--L_z', dest='L_z', type=float, default=None)
    args = parser.parse_args()
    if args.L_z is None:
        raise RuntimeError('L_z must be specified!')

    for filename in args.filenames:
        if not os.path.isfile(filename):
            raise RuntimeError(f'File "{filename}" not found!')
    if os.path.isfile(args.output):
        os.remove(args.output)

    keywords = None
    data = None
    for filename in args.filenames:
        with open(filename) as stream:
            try:
                parsed = yaml.safe_load(stream)
                if keywords is None:
                    keywords = parsed['keywords']
                    data = parsed['data']
                else:
                    data.extend(parsed['data'])
            except yaml.YAMLError as exc:
                print(exc)
                sys.exit()
    data = np.array(data)
    N_data = data.shape[0]
    if N_data < 30:
        print('Too few data points!')
        sys.exit()

    fig, ax = plt.subplots(2, 3)
    indices = [1, 2, 5, 6, 7, 8]
    units = ['K', 'eV', 'bar', 'bar', 'bar', 'bar']
    for i in range(2):
        for j in range(3):
            idx = indices[3*i + j]
            ax[i][j].hist(data[:,idx], bins=50)
            ax[i][j].set_xlabel(keywords[idx] + f' [{units[3*i + j]}]')
            ax[i][j].set_ylabel('Count [a.u.]')
    fig.tight_layout()
    plt.show()

    with open(args.output, 'w') as output_stream:
        for i in range(6):
            name = keywords[indices[i]]
            mean_val = np.mean(data[:,indices[i]])
            err = np.std(data[:,indices[i]]) / np.sqrt(N_data - 1)
            output_stream.write(f'{name} = {mean_val} \u00b1 {err} [{units[i]}]\n')
        pxx = data[:,6] * 1e5
        pyy = data[:,7] * 1e5
        pzz = data[:,8] * 1e5
        gamma = args.L_z * 1e-7 * (pzz - ((pxx + pyy) / 2)) / 2
        deviations = gamma - np.mean(gamma)
        N_tau = min(100, int(N_data / 2))
        correlations = np.empty(N_tau, dtype=float)
        correlations[0] = np.mean(deviations**2)
        for tau in range(1, N_tau):
            correlations[tau] = np.mean(deviations[tau:] * deviations[:-tau]) / correlations[0]
        correlations[0] = 1.0
        N_eff = (N_data / np.sum(correlations)) - 1
        output_stream.write(f'Surface tension = {np.mean(gamma)} \u00b1 {np.std(gamma) / np.sqrt(N_eff)} [mN/m]\n')

    sys.exit()
