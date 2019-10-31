#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:41:07 2019

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


engine = create_engine('sqlite:///NS001_evolved_mutations_copy6.db',
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
        gd_dicts = []
        for line in anno:
            try:
                gd_dicts.append(lgd.gd_line_to_dict(line))
            except ValueError:
                pass
        for gd_dict in gd_dicts:
            if gd_dict['type'] == 'DEL':
                parent = gd_dict['parent_id'] - 1 # correct for empty first line
                if type(parent) is int:
                    gd_dict['parent_type'] = gd_dicts[parent]['type']
                else:
                    parents = [int(x) for x in parent.split(',')]
                    parent_types = set(gd_dicts[parent]['type']
                                       for parent in parents)
                    if len(parent_types) == 1:
                        gd_dict['parent_type'] = parent_types.pop()
                    else:
                        gd_dict['parent_type'] = parent_types
        possible_merged_dels = [[]]
        for gd_dict in gd_dicts:
            if gd_dict['type'] == 'DEL':
                if gd_dict['parent_type'] != 'RA':
                    del_mut, del_ev = lgd.gdDEL_dict_to_del_dicts(sample,
                                                                  gd_dict)
                    lgd.add_DEL_dicts(del_mut, del_ev, session)
                else:
                    parent = gd_dict['parent_id'] - 1 # correct for empty first line
                    for pmerges in possible_merged_dels:
                        if lgd.check_merge_RAdel_dicts(*(pmerges +
                                                         [gd_dicts[parent]])):
                            pmerges.append(gd_dicts[parent])
                            break
                    else:
                        possible_merged_dels.append([gd_dicts[parent]])
        for pmerges in possible_merged_dels:
            if len(pmerges) == 0:
                break
            del_mut, del_ev = lgd.merge_RAdel_dicts_to_del_dicts(sample,
                                                                 *pmerges)
            lgd.add_DEL_dicts(del_mut, del_ev, session)
    try:
        session.commit()
    except DBAPIError:
        print('failed in directory {0}'.format(directory))
        session.rollback()
