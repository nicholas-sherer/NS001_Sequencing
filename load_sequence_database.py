#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 17:36:38 2019

@author: nicholas
"""

import contextlib
import os
import re
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import loadgdfiles as lgd

# shamelessly stolen from StackOverflow
@contextlib.contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    yield
    os.chdir(previous_dir)


engine = create_engine('sqlite:///NS001_evolved_mutations.db', echo=False)
session = sessionmaker(bind=engine)()

os.chdir('NS001_genomes_copy1')

# the directories with the output of breseq
directories = ['Agg_over_ancestors_against_self'] + \
    [str(i) for i in range(1, 21)] + [str(i) for i in range(23, 43)]

for directory in directories:
    if directory == 'Agg_over_ancestors_against_self':
        sample = 'Aggregate_NS001_Ancestors'
    else:
        with pushd(directory + '/output'):
            output_file = open('output.gd')
            line_with_samplename = output_file.readlines()[5]
            sample = re.findall('[a-zA-Z]+[1-4]t[1-2]',
                                line_with_samplename)[0]+'_S1'
    with pushd(directory+'/output/evidence'):
        anno = open('annotated.gd')
        gd_dicts = []
        for line in anno:
            try:
                gd_dicts.append(lgd.gd_line_to_dict(line))
            except ValueError:
                pass
        dictIDs_to_save = []
        for gd_dict in gd_dicts:
            if gd_dict['type'] == 'SNP':
                dictIDs_to_save.append(gd_dict['parent_id'])
        for gd_dict in gd_dicts:
            if gd_dict['id'] in dictIDs_to_save:
                gene, snp_mutation, snp_evidence = \
                    lgd.RA_dict_to_SNP_dicts(sample, gd_dict)
                lgd.commit_SNP_dicts(gene, snp_mutation, snp_evidence, session)
