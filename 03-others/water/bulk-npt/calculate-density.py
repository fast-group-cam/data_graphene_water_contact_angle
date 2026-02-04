#! /usr/bin/env python

import sys
import os
import yaml
import numpy as np
import matplotlib.pyplot as plt

FILENAME = 'run/npt_prod.yaml'
OUTPUT_FILE = 'density.txt'
N_ATOMS = 520

if __name__ == '__main__':

    if not os.path.isfile(FILENAME):
        raise RuntimeError(f'File "{FILENAME}" not found!')
    if os.path.isfile(OUTPUT_FILE):
        os.remove(OUTPUT_FILE)

    keywords = None
    data = None
    with open(FILENAME) as stream:
        try:
            parsed = yaml.safe_load(stream)
            keywords = parsed['keywords']
            data = parsed['data']
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()
    data = np.array(data)
    N_data = data.shape[0]
    if N_data < 100:
        print('Too few data points!')
        sys.exit()

    fig, ax = plt.subplots(2, 3)
    indices = [1, 4, 7, 8, 9, 10]
    units = ['K', 'A', 'bar', 'bar', 'bar', 'bar']
    for i in range(2):
        for j in range(3):
            idx = indices[3*i + j]
            ax[i][j].plot(data[:100,0], data[:100,idx], 'b.')
            ax[i][j].set_xlabel(keywords[0])
            ax[i][j].set_ylabel(keywords[idx] + f' [{units[3*i + j]}]')
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots(2, 3)
    for i in range(2):
        for j in range(3):
            idx = indices[3*i + j]
            ax[i][j].hist(data[:,idx], bins=50)
            ax[i][j].set_xlabel(keywords[idx] + f' [{units[3*i + j]}]')
            ax[i][j].set_ylabel('Count [a.u.]')
    fig.tight_layout()
    plt.show()

    with open(OUTPUT_FILE, 'w') as output_stream:
        for i in range(6):
            name = keywords[indices[i]]
            mean_val = np.mean(data[:,indices[i]])
            err = np.std(data[:,indices[i]]) / np.sqrt(N_data - 1)
            output_stream.write(f'{name} = {mean_val} \u00b1 {err} [{units[i]}]\n')
        length = data[:,4]
        volume = length**3
        deviations = volume - np.mean(volume)
        N_tau = min(100, int(N_data / 2))
        correlations = np.empty(N_tau, dtype=float)
        correlations[0] = np.mean(deviations**2)
        for tau in range(1, N_tau):
            correlations[tau] = np.mean(deviations[tau:] * deviations[:-tau]) / correlations[0]
        correlations[0] = 1.0
        print(np.sum(correlations)) #DEBUG
        N_eff = (N_data / np.sum(correlations)) - 1
        output_stream.write(f'Volume = {np.mean(volume)} \u00b1 {np.std(volume) / np.sqrt(N_eff)} [\u212b\u00b3]\n')

        density_val = N_ATOMS / np.mean(volume)
        density_err = N_ATOMS * np.std(volume) / (np.sqrt(N_eff) * (np.mean(volume)**2))
        output_stream.write(f'Density = {density_val} \u00b1 {density_err} [\u212b\u207b\u00b3]\n')
        conversion_constant = 0.01801528 * (1e7 / (6.02214076))
        density_val *= conversion_constant
        density_err *= conversion_constant
        output_stream.write(f'Density = {density_val} \u00b1 {density_err} [kg/m\u00b3]\n')

    sys.exit()
