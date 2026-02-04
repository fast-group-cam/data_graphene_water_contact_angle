#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

NUMBERS = [520, 1170, 2080, 3250, 4680]
OUTPUT_FILE = 'surface-tension.png'

if __name__ == '__main__':

    L_x = list()
    L_y = list()
    gamma = list()
    gamma_err = list()
    for n in NUMBERS:
        filename = f'slab-{n}/start.lammps-data'
        with open(filename) as stream:
            for _ in range(5):
                stream.readline()
            L_x.append(float(stream.readline().split()[1]))
            L_y.append(float(stream.readline().split()[1]))
        filename = f'slab-{n}/surface-tension.txt'
        with open(filename) as stream:
            for line in stream:
                if line.startswith('Surface tension = '):
                    tokens = line.removeprefix('Surface tension = ').split()
                    gamma.append(float(tokens[0]))
                    gamma_err.append(float(tokens[2]))

    L_x = np.array(L_x)
    L_y = np.array(L_y)
    gamma = np.array(gamma)
    gamma_err = np.array(gamma_err)
    inv_L = (0.5 / L_x) + (0.5 / L_y)
    x_range = np.linspace(0.0, 1.1 * np.max(inv_L), 100)

    popt, cov = np.polyfit(inv_L, gamma, 1, full=False, cov=True)
    gamma_inf = popt[1]
    gamma_inf_err = 2.0 * np.sqrt(cov[1,1])

    fig, ax = plt.subplots()
    fig.set_size_inches(7.2, 4.8)
    ax.plot(x_range, ((popt[0] * x_range) + popt[1]), 'k--', label=r'Best-fit $1/L_{\parallel}$')
    ax.annotate(r'$' + f'{gamma_inf:.1f}' + r'\,\pm\,' + f'{gamma_inf_err:.1f}' + r'\, mN/m$',
                (0.02, 0.11), xycoords='axes fraction', ha='left', va='center')
    
    popt, cov = np.polyfit(inv_L**2, gamma, 1, full=False, cov=True)
    gamma_inf = popt[1]
    gamma_inf_err = 2.0 * np.sqrt(cov[1,1])
    ax.plot(x_range, ((popt[0] * x_range**2) + popt[1]), '-.', color='steelblue', label=r'Best-fit $\ln L_{z}/L_{\parallel}^{2}$')
    ax.annotate(r'$' + f'{gamma_inf:.1f}' + r'\,\pm\,' + f'{gamma_inf_err:.1f}' + r'\, mN/m$',
                (0.02, 0.33), xycoords='axes fraction', ha='left', va='center', color='steelblue')

    ax.errorbar(inv_L, gamma, yerr=(2*gamma_err), c='r', marker='s', ls='')
    ax.set_xlabel(r'$1/L_{\parallel}\; [\AA^{-1}]$')
    ax.set_ylabel(r'$\gamma\; [mN/m]$')
    ax.set_xlim(0.0, np.max(x_range))
    ax.set_ylim(70.5, 84.9)
    ax.legend(loc='lower right')
    fig.savefig(OUTPUT_FILE, dpi=(3*fig.dpi), bbox_inches='tight', pad_inches=0.05)

