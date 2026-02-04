#! /usr/bin/env python

#==================================================================================================
# Python script for collating multiple outputs of the contact_angle script, and combining them as
# "manual" block averages. Use as:
#
#     python collate-blocks.py <prefix> -o <output_folder>
#
#==================================================================================================

import os
import argparse
import configparser
import numpy as np

if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(prog='collate-blocks.py')
    parser.add_argument('prefix', nargs='+')
    parser.add_argument('-o', '--output', default=None)
    args = parser.parse_args()
    if args.output is None:
        raise RuntimeError('Output file must be specified!')
    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # Detect valid targets
    targets = list()
    for folder_prefix in args.prefix:
        for i in range(100):
            folder_name = folder_prefix + f'_{i}'
            if os.path.isdir(folder_name) and os.path.isfile(os.path.join(folder_name, 'results.ini')):
                targets.append(os.path.join(folder_name, 'results.ini'))
    if len(targets) == 0:
        raise RuntimeError('No inputs found!')

    # Read targets
    N_data = len(targets)
    data = [None for _ in range(N_data)]
    for i in range(N_data):
        config = configparser.ConfigParser()
        config.read(targets[i])
        data[i] = config

    # Function that counts the number of header/key appearances
    def get_count(header, key):
        count = 0
        for i in range(N_data):
            if (header in data[i]) and (key in data[i][header]):
                count += 1
        return count
    
    # Function that performs block averaging
    def get_block_ave(header, key):
        values = list()
        for i in range(N_data):
            if (header in data[i]) and (key in data[i][header]):
                values.append(data[i].getfloat(header, key))
        return np.mean(values), (np.std(values) / np.sqrt(len(values) - 1))

    # Output to results file
    with open(os.path.join(args.output, 'results.ini'), 'w') as output_file:

        output_file.write('[General]\n')
        output_file.write(f'Collated from [files] = {N_data}\n')
        keys = ['Droplet roof [A]', 'Droplet floor [A]', 'Droplet CoM z-coordinate [A]',
                'Graphene sheet z at origin [A]', 'Droplet height [A]',
                'Nominal interfacial separation [A]']
        for key in keys:
            mean, unc = get_block_ave('General', key)
            output_file.write(f'{key}, mean of block means = {mean}\n')
            output_file.write(f'{key}, uncertainty = {unc}\n')
        output_file.write('\n')

        output_file.write('[Block-Averaged Interface]\n')
        keys = ['Contact angle [deg]', 'Three-phase line radius [A]', 'Best-fit sphere radius [A]',
                'Best-fit sphere z-height [A]']
        for key in keys:
            mean, unc = get_block_ave('Time-Averaged Interface', key)
            output_file.write(f'{key}, mean of block means = {mean}\n')
            output_file.write(f'{key}, uncertainty = {unc}\n')
        output_file.write('\n')

        if get_count('Gibbs Dividing Surface', 'Contact angle [deg]') > 0:

            output_file.write('[Block-Averaged Gibbs Dividing Surface]\n')
            keys = ['Contact angle [deg]', 'Three-phase line radius [A]',
                    'Best-fit sphere radius [A]', 'Best-fit sphere z-height [A]']
            for key in keys:
                mean, unc = get_block_ave('Gibbs Dividing Surface', key)
                output_file.write(f'{key}, mean of block means = {mean}\n')
                output_file.write(f'{key}, uncertainty = {unc}\n')
            output_file.write('\n')
    
