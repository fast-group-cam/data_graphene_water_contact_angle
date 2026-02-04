#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

RESULTS_FILE = '../../03-others/interaction/results.npz'

REF_MORSE_PARAMS = np.array([[[90.1069, 95.5030,  84.4691], [3.10529, 3.12997, 3.08364], [1.30431, 1.21119, 1.39939]],
                             [[92.3613, 97.7956,  86.8458], [3.45671, 3.47248, 3.44256], [1.43449, 1.36350, 1.50699]],
                             [[98.9979, 105.0070, 92.8462], [3.37265, 3.38530, 3.36138], [1.34725, 1.29401, 1.40166]]])

def morse_potential(leg, variant, x):
    return REF_MORSE_PARAMS[leg, 0, variant] * ((1.0 - np.exp(-REF_MORSE_PARAMS[leg, 2, variant] * (x - REF_MORSE_PARAMS[leg, 1, variant])))**2 - 1.0)

def test_curve(x, D, x0, a, offset):
    return (D * ((1.0 - np.exp(-a * (x - x0)))**2 - 1.0)) + offset

results = np.load(RESULTS_FILE)
strains = results['strains']
dists = results['dists']
cp2k_energies = 1000 * results['cp2k_energies']
mace_energies = 1000 * results['mace_energies']

cp2k_D = np.zeros((3, 3), dtype=float)
cp2k_Derr = np.zeros((3, 3), dtype=float)
cp2k_popt = np.zeros((3, 3, 4), dtype=float)
cp2k_offset = np.zeros((3, 3), dtype=float)
mace_D = np.zeros((3, 3), dtype=float)
mace_Derr = np.zeros((3, 3), dtype=float)
mace_popt = np.zeros((3, 3, 4), dtype=float)
mace_offset = np.zeros((3, 3), dtype=float)
for i in range(3):
    for j in range(3):
        popt, pcov = curve_fit(test_curve, dists, cp2k_energies[i,j,:], p0=(*REF_MORSE_PARAMS[i,:,0], 58000000.0))
        cp2k_D[i,j] = popt[0]
        cp2k_Derr[i,j] = 2.0 * np.sqrt(np.diag(pcov))[0]
        cp2k_popt[i,j,:] = np.copy(popt)
        cp2k_offset[i,j] = popt[3]
        popt, pcov = curve_fit(test_curve, dists, mace_energies[i,j,:], p0=(*REF_MORSE_PARAMS[i,:,0], 0.0))
        mace_D[i,j] = popt[0]
        mace_Derr[i,j] = 2.0 * np.sqrt(np.diag(pcov))[0]
        mace_popt[i,j,:] = np.copy(popt)
        mace_offset[i,j] = popt[3]
cp2k_popt[:,:,3] = 0.0
mace_popt[:,:,3] = 0.0

fig, ax = plt.subplots(3, 3, sharex=True, sharey=True)
fig.set_size_inches(10, 6.3)
x = np.linspace(2.5, 4.5, 200)
for i in range(3):
    for j in range(3):
        if i == 0:
            ax[0,j].fill_between(x, morse_potential(j, 1, x), morse_potential(j, 2, x), fc=(0.8, 0.8, 0.8), ec='none', zorder=0.0)
            ax[0,j].plot(x, morse_potential(j, 0, x), 'k-', label='Ref. (DMC)', zorder=0.1)
            ax[0,j].set_title(f'{j}-leg')
        ax[i,j].plot(dists, cp2k_energies[j,i,:] - cp2k_offset[j,i], 's', mfc='none', mec='red', label='revPBE-D3', zorder=0.2)
        ax[i,j].plot(x, test_curve(x, *cp2k_popt[j,i]), '-', color='red', alpha=0.5, zorder=0.1)
        ax[i,j].plot(dists, mace_energies[j,i,:] - mace_offset[j,i], 'o', mfc='none', mec='steelblue', label='MACE', zorder=0.2)
        ax[i,j].plot(x, test_curve(x, *mace_popt[j,i]), '-', color='steelblue', alpha=0.5, zorder=0.1)
        if i == 2:
            ax[2,j].set_xlabel(r'Distance $d$ [$\AA$]')
        else:
            ax[i,j].tick_params(axis='x', bottom=False)
        if j == 0:
            ax[i,0].set_ylabel(r'$\epsilon = ' + f'{strains[i]:+.1f}' + r'\%$, $E_b$ [meV]')
            ax[i,0].set_ylim(-199.999, 399.999)
        else:
            ax[i,j].tick_params(axis='y', left=False)
ax[0,2].legend()
fig.tight_layout()
fig.savefig('interaction-energy-dft.png', dpi=(3.2*fig.dpi), bbox_inches='tight', pad_inches=0.05)

fig, ax = plt.subplots(1, 3, sharey=True)
fig.set_size_inches(10, 3.7)
for i in range(3):
    ax[i].errorbar((0.0,), (REF_MORSE_PARAMS[i, 0, 0],), yerr=((REF_MORSE_PARAMS[i, 0, 0] - REF_MORSE_PARAMS[i, 0, 2],), (REF_MORSE_PARAMS[i, 0, 1] - REF_MORSE_PARAMS[i, 0, 0],)), fmt='kD', label='Ref. (DMC)', zorder=0.1)
    ax[i].plot(strains, cp2k_D[i,:], '-', color='red', alpha=0.5, zorder=0.0)
    ax[i].errorbar(strains, cp2k_D[i,:], yerr=cp2k_Derr[i,:], fmt='s', c='red', mfc='none', mec='red', label='revPBE-D3', zorder=0.2)
    ax[i].plot(strains, mace_D[i,:], '-', color='steelblue', alpha=0.5, zorder=0.0)
    ax[i].errorbar(strains, mace_D[i,:], yerr=mace_Derr[i,:], fmt='o', c='steelblue', mfc='none', mec='steelblue', label='MACE', zorder=0.3)
    ax[i].set_xlim(-0.5, 2.5)
    ax[i].set_xticks((0.0, 1.0, 2.0), ('0.0', '1.0', '2.0'))
    ax[i].set_title(f'{i}-leg')
    ax[i].set_xlabel(r'Applied strain $\epsilon$ [%]')
ax[0].set_ylabel('Potential well depth [meV]')
ax[1].tick_params(axis='y', left=False)
ax[2].tick_params(axis='y', left=False)
ax[2].legend(loc='lower right')
fig.tight_layout()
fig.savefig('min-interaction-energy.png', dpi=(3.2*fig.dpi), bbox_inches='tight', pad_inches=0.05)

with open('min-interaction-energy-delta.txt', 'w') as output_stream:
    for i in range(3):
        output_stream.write('-------------------\n')
        output_stream.write(f'{i}-leg configuration\n')
        output_stream.write('-------------------\n\n')
        output_stream.write('DFT calculation (CP2K, revPBE-D3):\n')
        for j in range(3):
            output_stream.write(f'  - {strains[j]:.1f}% strain, D = {cp2k_D[i,j]:.1f} \u00b1 {cp2k_Derr[i,j]:.1f} meV\n')
        output_stream.write(f'  Relative change is {100.0 * ((cp2k_D[i,2] / cp2k_D[i,0]) - 1.0):+.1f}%\n\n')
        output_stream.write('MACE model:\n')
        for j in range(3):
            output_stream.write(f'  - {strains[j]:.1f}% strain, D = {mace_D[i,j]:.1f} \u00b1 {mace_Derr[i,j]:.1f} meV\n')
        output_stream.write(f'  Relative change is {100.0 * ((mace_D[i,2] / mace_D[i,0]) - 1.0):+.1f}%\n\n')
