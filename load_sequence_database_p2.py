#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 17:36:38 2019

@author: nicholas
"""

import contextlib
import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import loadgdfiles as lgd
import SequenceDataORM as SqD

# shamelessly stolen from StackOverflow
@contextlib.contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    yield
    os.chdir(previous_dir)


engine = create_engine('sqlite:///NS001_evolved_mutations.db', echo=False)
session = sessionmaker(bind=engine)()

session.add(SqD.DNA_Sample(name='Ancestor_S1', strain='NS001'))
session.add(SqD.DNA_Sample(name='Ancestor_S2', strain='NS001'))
session.add(SqD.DNA_Sample(name='Ancestor_S3', strain='NS001'))

os.chdir('NS001_genomes_copy1')

# the directories with the output of breseq
directories = [str(21), str(22), str(43)]

for directory in directories:
    if directory == '21':
        sample = 'Ancestor_S1'
    elif directory == '22':
        sample = 'Ancestor_S2'
    elif directory == '43':
        sample = 'Ancestor_S3'
    else:
        raise RuntimeError('somehow not in an anticipated directory')
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
