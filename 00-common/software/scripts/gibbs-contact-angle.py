#! /usr/bin/env python

#==================================================================================================
# Python script for calculating the contact_angle using the Gibbs dividing surface definition (i.e.
# the standard method in the literature). Use as:
#
#     python gibbs-contact-angle.py <input_files(s)> [-o <output_dir>] [--options]
#
#==================================================================================================

import os
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from droplet_graphene_analysis.util import elapsed_time, read_droplet_trajectory
from droplet_graphene_analysis.util.droplet import find_interface
from droplet_graphene_analysis.util.droplet.coarse_grain import BULK_DENSITY

if __name__ == '__main__':

    N_Z_BINS = 8
    N_R_BINS = 50

    # Parse input arguments
    parser = argparse.ArgumentParser(prog='gibbs-contact-angle.py')
    parser.add_argument('input_file', nargs='+',
                        help='input file(s) to read data from')
    parser.add_argument('--index', default=':', dest='index',
                        help='index or slice of indices to take from each input file')
    parser.add_argument('-o', '--output', default='gibbs-contact-angle', dest='output_dir',
                        help='output folder to save log and graphical outputs to')
    parser.add_argument('--no-graphics', action='store_true', dest='no_graphics',
                        help='disables rendering of graphics (and speeds up the script)')
    args = parser.parse_args()

    for file in args.input_file:
        if not os.path.isfile(file):
            raise RuntimeError(f'File "{file}" not found.')
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    # Read input file and save coordinates
    file_msg = (f'"{args.input_file[0]}"' if len(args.input_file) == 1 else
                f'{len(args.input_file)} files')
    print(f'Reading {file_msg}...', end='')
    time_start_0 = time.time()
    cell_params, waters, _, _ = read_droplet_trajectory(args.input_file, index=args.index)
    N_frames, N_water, _ = waters.shape
    print(f'read {N_frames} frames from {file_msg} in {elapsed_time(time_start_0)}.')

    # Convert to cylindrical coordinates
    CoM = (0, 0, np.mean(waters[:,:,2]))
    z_floor = find_interface(waters, CoM, (0, 0, -1))[2]
    z_ceil = find_interface(waters, CoM, (0, 0, 1))[2]
    radial_coords = np.sqrt((waters[:,:,0]**2) + (waters[:,:,1]**2))

    # Fitting function
    def test_curve(x, A, x0, delta):
        return 0.5 * A * (np.tanh((x0 - x) / delta) + 1)

    # Split into bins and count
    print(f'Binning and counting...', end='')
    time_start_1 = time.time()
    z_bin_edges = np.linspace(z_floor, z_ceil, N_Z_BINS + 1)
    z_bin_centers = (z_bin_edges[:-1] + z_bin_edges[1:]) / 2
    z_bin_width = z_bin_edges[1] - z_bin_edges[0]
    radii = np.empty(N_Z_BINS, dtype=float)
    for i in range(N_Z_BINS):
        guess = find_interface(waters, (0, 0, z_bin_centers[i]), (1, 0, 0))[0]
        r_bin_edges = np.linspace(0, 1.5 * guess, N_R_BINS + 1)
        r_bin_centers = (r_bin_edges[:-1] + r_bin_edges[1:]) / 2
        r_bin_volumes = np.pi * z_bin_width * ((r_bin_edges[1:]**2) - (r_bin_edges[:-1]**2))
        sample = radial_coords[np.logical_and(waters[:,:,2] >= z_bin_edges[i], waters[:,:,2] < z_bin_edges[i+1])].flatten()
        density = np.histogram(sample, r_bin_edges)[0] / (N_frames * r_bin_volumes)
        popt, _ = curve_fit(test_curve, r_bin_centers, density, p0=(BULK_DENSITY, guess, 1))
        radii[i] = popt[1]
        if not args.no_graphics:
            fig, ax = plt.subplots()
            fig.set_size_inches(8, 6)
            plot_range = np.linspace(0, 1.5 * guess, N_R_BINS * 3)
            ax.plot(plot_range, test_curve(plot_range, *popt), 'k--')
            ax.plot(r_bin_centers, density, 'b.')
            ax.set_title(f'Density against r coord at z = {z_bin_centers[i]:.1f} ' + r'$\AA$')
            ax.set_xlabel(r'r [$\AA$]')
            ax.set_ylabel(r'$\rho$ [$\AA^{-3}$]')
            ax.annotate(r'$\rho_{0}\,=\,' + f'{popt[0]:.5f}' + r'\,\AA^{-3}$' + '\n' + r'$r(z)\,=\,' + f'{radii[i]:.1f}' + r'\,\AA$' +
                        '\n' + r'$\delta r\,=\,' + f'{popt[2]:.2f}' + r'\,\AA$', (0.99, 0.99), xycoords='axes fraction', ha='right', va='top')
            fig.savefig(os.path.join(args.output_dir, f'radial-fit_{i}.png'), dpi=(3*fig.dpi), bbox_inches='tight', pad_inches=0.05)
    print(f'done in {elapsed_time(time_start_1)}.')
    
    # Fitting function
    def sphere_func(z, R_sq, c):
        return R_sq - (z - c)**2

    # Find best-fit sphere
    popt, _ = curve_fit(sphere_func, z_bin_centers, radii**2, p0=(np.max(radii)**2, CoM[2]))
    sphere_R = np.sqrt(popt[0])
    sphere_c = popt[1]
    contact_angle = 90 + (np.arcsin(sphere_c / sphere_R) * 180 / np.pi)

    if not args.no_graphics:
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 6)
        ax.plot(radial_coords.flatten(), waters[:,:,2].flatten(), 'b,')
        for i in range(N_Z_BINS):
            ax.plot((radii[i], radii[i]), (z_bin_edges[i], z_bin_edges[i+1]), 'r-')
        theta = np.linspace(0, contact_angle * np.pi / 180, 200)
        ax.plot(sphere_R * np.sin(theta), sphere_c + sphere_R * np.cos(theta), '-', color=(1.0, 0.5, 0.0))
        ax.annotate(r'$R\,=\,' + f'{sphere_R:.1f}' + r'\,\AA$' + '\n' + r'$z_{c}\,=\,' + f'{sphere_c:.1f}' + r'\,\AA$' + '\n' +
                    r'$\theta\,=\,' + f'{contact_angle:.1f}' + r'\degree$', (0.99, 0.99), xycoords='axes fraction', ha='right', va='top')
        fig.savefig(os.path.join(args.output_dir, f'sphere-cap.png'), dpi=(3*fig.dpi), bbox_inches='tight', pad_inches=0.05)

    # Write to output
    results_file = open(os.path.join(args.output_dir, 'results.ini'), 'w', encoding='utf-8')
    results_file.write('[General]\n')
    results_file.write(f'No. of frames = {N_frames}\n')
    results_file.write(f'No. of water molecules = {N_water}\n')
    results_file.write(f'Droplet roof [A] = {z_ceil}\n')
    results_file.write(f'Droplet floor [A] = {z_floor}\n\n')
    results_file.write('[Gibbs Dividing Surface]\n')
    results_file.write(f'Contact angle [deg] = {contact_angle}\n')
    results_file.write(f'Three-phase line radius [A] = {np.sqrt(max(sphere_R**2 - sphere_c**2, 0.0))}\n')
    results_file.write(f'Best-fit sphere radius [A] = {sphere_R}\n')
    results_file.write(f'Best-fit sphere z-height [A] = {sphere_c}\n\n')

    final_elapsed_time = elapsed_time(time_start_0)
    results_file.write('[Misc]\n')
    results_file.write('Program type = gibbs_contact_angle\n')
    results_file.write(f'Program wall time = {final_elapsed_time}\n')
    results_file.close()
    print(f'Program completed in {final_elapsed_time}.')

