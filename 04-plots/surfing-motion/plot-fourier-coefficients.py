#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

STRAINS = [-2.0, -1.5, -1.0, -0.7, -0.5, -0.4, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.2, 0.5, 1.0, 1.5, 2.0]
OUTPUT_FILE = 'fourier-coefficients.png'

if __name__ == '__main__':

    x_coeffs = list()
    y_coeffs = list()

    for strain in STRAINS:
        
        folder_name = 's+0.00' if strain == 0.0 else f's{strain:+.2f}'
        loaded = np.load(os.path.join('cache', f'fourier_{folder_name}.npz'))
        cx = loaded['cx']
        cy = loaded['cy']

        mags_x = np.abs(cx)
        mags_x_mean = np.mean(mags_x)
        mags_x_std = np.std(mags_x)
        mags_x = mags_x[mags_x < (mags_x_mean + (5 * mags_x_std))]
        mags_x = mags_x[mags_x > (mags_x_mean - (5 * mags_x_std))]

        mags_y = np.abs(cy)
        mags_y_mean = np.mean(mags_y)
        mags_y_std = np.std(mags_y)
        mags_y = mags_y[mags_y < (mags_y_mean + (5 * mags_y_std))]
        mags_y = mags_y[mags_y > (mags_y_mean - (5 * mags_y_std))]

        if mags_x_mean >= mags_y_mean:
            x_coeffs.append(list(mags_x))
            y_coeffs.append(list(mags_y))
        else:
            x_coeffs.append(list(mags_y))
            y_coeffs.append(list(mags_x))

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)
    
    parts = ax.violinplot(x_coeffs, showmeans=True, showextrema=False, side='low')
    for pc in parts['bodies']:
        pc.set_facecolor('red')
        pc.set_alpha(0.5)
    parts['cmeans'].set_color('red')
        
    parts = ax.violinplot(y_coeffs, showmeans=True, showextrema=False, side='high')
    for pc in parts['bodies']:
        pc.set_facecolor('steelblue')
        pc.set_alpha(0.5)
    parts['cmeans'].set_color('steelblue')

    ax.set_xlabel(r'Applied strain $\epsilon$  [%]')
    ax.set_ylabel(r'Distribution of $|c_{(m,n)}(t)|$ over time  [$\AA$]')
    ax.set_ylim(0.0, 2.6)
    ax.set_xticks([i + 1 for i in range(len(STRAINS))], [f'{s:+.2f}' for s in STRAINS])
    ax.tick_params(axis='x', labelrotation=20)

    labels = [(mpatches.Patch(color='red', alpha=0.6), r'1st x-direction, $|c_{(1,0)}|$'),
              (mpatches.Patch(color='steelblue', alpha=0.6), r'1st y-direction, $|c_{(0,1)}|$')]
    ax.legend(*zip(*labels))

    fig.savefig(OUTPUT_FILE, dpi=(2 * fig.dpi), bbox_inches='tight', pad_inches=0.05)

