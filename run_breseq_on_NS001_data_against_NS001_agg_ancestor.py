# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a one off script file for processing my sequencing data for evolved
replicates of NS001 Delta_cat with breseq
"""

import glob
import os

os.chdir('/home/nicholas/Documents/NS001_genomes_copy1')

subdirs_to_process = filter(os.path.isdir, os.listdir(os.curdir))

for subdir in subdirs_to_process:
    # move to subdirectory
    os.chdir(subdir)
    # check for output
    if 'output' not in os.listdir(os.curdir):
        # get fastq file names
        files = glob.glob('*.fastq')
        if len(files) > 0:
        # run breseq to find e. coli mutations
            os.system('breseq -r ~/Documents/NS001_genomes/Aggregate_over_Ancestors_21_22_and_43/output/NS001_Delta_cat.gff3 -j 6 -p ' + ' '.join(files))
    os.chdir('..')