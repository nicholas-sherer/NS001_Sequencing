#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:20:11 2019

@author: nicholas
"""

import contextlib
import json
import numpy as np
import re
import os


# shamelessly stolen from StackOverflow
@contextlib.contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    yield
    os.chdir(previous_dir)


os.chdir('NS001_genomes_copy1')
directories = list(filter(os.path.isdir, os.listdir()))
directories.remove('44')

coverages = open('sample_coverages.csv', 'w')
coverage_arrays = {}


for directory in directories:
    if directory == 'Agg_over_ancestors_against_self':
        sample = 'Aggregate_NS001_Ancestors'
    elif directory == '21':
        sample = 'Ancestor_S1'
    elif directory == '22':
        sample = 'Ancestor_S2'
    elif directory == '43':
        sample = 'Ancestor_S3'
    else:
        with pushd(directory + '/output'):
            output_file = open('output.gd')
            line_with_samplename = output_file.readlines()[5]
            sample = re.findall('[a-zA-Z]+[1-4]t[1-2]',
                                line_with_samplename)[0]+'_S1'
    total_reads = 0
    total_reads_mapped = 0
    nbinom_mean = 0
    nbinom_dispersion = 0
    with pushd(directory+'/05_alignment_correction'):
        read_json = json.load(open('summary.json'))
        total_reads = read_json['total_reads']
        total_reads_mapped = read_json['total_reads_mapped_to_references']
    with pushd(directory+'/07_error_calibration'):
        read_json2 = json.load(open('summary.json'))
        nbinom_mean = read_json2['NC_000913']['nbinom_mean_parameter']
        nbinom_dispersion = read_json2['NC_000913']['nbinom_dispersion']
    coverages.write('{0}, {1}, {2}, {3}, {4}\n'
                    .format(sample, total_reads, total_reads_mapped,
                            nbinom_mean, nbinom_dispersion))
    with pushd(directory+'/08_mutation_identification'):
        coverage_arrays[sample] = np.loadtxt('NC_000913.coverage.tab',
                                             skiprows=1, usecols=(0, 1))
coverages.close()
np.savez_compressed('samples_coverage_arrays.npz', **coverage_arrays)
