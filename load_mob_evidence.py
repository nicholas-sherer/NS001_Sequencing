#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 15:21:03 2019

@author: nicholas
"""

import contextlib
import os
import re
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import DBAPIError

import loadgdfiles as lgd

# shamelessly stolen from StackOverflow
@contextlib.contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    yield
    os.chdir(previous_dir)


engine = create_engine('sqlite:///NS001_evolved_mutations_copy4.db',
                       echo=False)
session = sessionmaker(bind=engine)()

os.chdir('NS001_genomes_copy1')

directories = os.listdir()

directories.remove('44')

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
    with pushd(directory+'/output/evidence'):
        anno = open('annotated.gd')
        for line in anno:
            try:
                gd_dict = lgd.gd_line_to_dict(line)
                if gd_dict['type'] == 'MOB':
                    mob_mut, mob_ev = \
                        lgd.gdMOB_dict_to_MOB_dicts(sample, gd_dict)
                    lgd.add_MOB_dicts(mob_mut, mob_ev, session)
            except ValueError:
                pass
    try:
        session.commit()
    except DBAPIError:
        print('failed in directory {0}'.format(directory))
        session.rollback()
