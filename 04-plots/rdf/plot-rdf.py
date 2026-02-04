#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    def plot_from_file(axis, filename, options):
        with open(filename, 'r') as stream:
            x = list()
            y = list()
            stream.readline()
            stream.readline()
            for line in stream:
                tokens = line.split(' ')[:2]
                x.append(float(tokens[0]))
                y.append(float(tokens[1]))
            x = np.array(x)
            y = np.array(y)
            axis.plot(x, y, **options)
    
    expt_opt = {'linestyle': '--', 'color': 'firebrick', 'label': 'Expt.'}
    mace_npt_opt = {'linestyle': (0, (1.5, 1)), 'color': 'steelblue', 'label': 'MLP-MD, NpT'}
    aimd_opt = {'linestyle': '-.', 'color': 'peru', 'label': 'Ref. AIMD, NVT'}
    mace_nvt_opt = {'linestyle': (0, (1.5, 1)), 'color': 'steelblue', 'label': 'MLP-MD, NVT'}

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 3.5)
    plot_from_file(ax, '../../03-others/water/bulk-npt/rdf-ref-expt/rdf_OO_expt.dat', expt_opt)
    plot_from_file(ax, '../../03-others/water/bulk-npt/rdf-mace/rdf_OO_mace.dat', mace_npt_opt)
    ax.set_ylim(0, 2.9)
    ax.set_ylabel(r'$g_{OO}(r)$')
    ax.annotate('O-O', (0.01, 0.98), xycoords='axes fraction', ha='left', va='top')

    ax.set_xlabel(r'$r$ [$\AA$]')
    ax.set_xlim(0, 9.8)
    ax.legend()
    fig.savefig('rdf-npt.png', dpi=(3*fig.dpi), bbox_inches='tight', pad_inches=0.05)

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 3.5)
    plot_from_file(ax, '../../03-others/water/bulk-nvt/rdf-ref-aimd/rdf_OO_aimd.dat', aimd_opt)
    plot_from_file(ax, '../../03-others/water/bulk-nvt/rdf-mace/rdf_OO_mace.dat', mace_nvt_opt)
    ax.set_ylim(0, 2.9)
    ax.set_ylabel(r'$g_{OO}(r)$')
    ax.annotate('O-O', (0.01, 0.98), xycoords='axes fraction', ha='left', va='top')

    ax.set_xlabel(r'$r$ [$\AA$]')
    ax.set_xlim(2.1, 5.9)
    ax.legend()
    fig.savefig('rdf-nvt.png', dpi=(3*fig.dpi), bbox_inches='tight', pad_inches=0.05)

