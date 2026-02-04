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
        filename = f'../../03-others/water/slab-{n}/start.lammps-data'
        with open(filename) as stream:
            for _ in range(5):
                stream.readline()
            L_x.append(float(stream.readline().split()[1]))
            L_y.append(float(stream.readline().split()[1]))
        filename = f'../../03-others/water/slab-{n}/surface-tension.txt'
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
    L_par = 1.0 / inv_L
    L_range = np.linspace(13.0, 1.3 * np.max(L_par), 200)

    popt, cov = np.polyfit(inv_L**2, gamma, 1, full=False, cov=True)
    gamma_inf = popt[1]
    gamma_inf_err = 2.0 * np.sqrt(cov[1,1])

    fig, ax = plt.subplots()
    fig.set_size_inches(7.2, 4.8)
    
    ax.plot(L_range, ((popt[0] * np.power(L_range, -2)) + popt[1]), '--', color='steelblue', label='Best-fit correction')
    ax.annotate(r'$' + f'{gamma_inf:.1f}' + r'\,\pm\,' + f'{gamma_inf_err:.1f}' + r'\, mN/m$',
                (0.99, 0.28), xycoords='axes fraction', ha='right', va='top', color='steelblue')

    ax.errorbar(L_par, gamma, yerr=(1.5 * gamma_err), c='r', marker='s', ls='')

    for i, n in enumerate(NUMBERS):
        y = (gamma[i] + (1.5 * gamma_err[i]) + 0.3) if i < 2 else (gamma[i] - (1.5 * gamma_err[i]) - 0.3)
        va = 'bottom' if i < 2 else 'top'
        ax.annotate(f'{n}\nmolecules', (L_par[i], y), ha='center', va=va, color=(1.0, 0.0, 0.0, 0.5))

    ax.set_xlabel(r'$L_{\parallel}\; [\AA]$')
    ax.set_ylabel(r'$\gamma_{\text{lv}}\; [mN/m]$')
    ax.set_xlim(13.0, np.max(L_range))
    ax.set_ylim(70.5, 89.9)
    ax.legend(loc='lower left')
    fig.savefig(OUTPUT_FILE, dpi=(5*fig.dpi), bbox_inches='tight', pad_inches=0.05)

