#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:16:45 2019

@author: nicholas
"""

import re
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import DBAPIError
import SequenceDataORM_updt as SqD


engine = create_engine('sqlite:///NS001_evolved_mutations_copy3.db',
                       echo=False)
session = sessionmaker(bind=engine)()


name_match = re.compile('Name=([^;]+);')
description_match = re.compile('Note=([^;]+)($|;)')

with open('NS001_Delta_cat.gff3') as f:
    for line in f:
        entries = line.split('\t')
        if len(entries) >= 2:
            entry_type = entries[2]
        else:
            entry_type = None
        gene_entry_types = ['CDS', 'rRNA', 'tRNA', 'fCDS', 'ncRNA']
        if entry_type in gene_entry_types:
            gene_dict = {}
            gene_dict['start'] = int(entries[3])
            gene_dict['end'] = int(entries[4])
            if entries[6] == '+':
                gene_dict['strand'] = True
            elif entries[6] == '-':
                gene_dict['strand'] = False
            else:
                gene_dict['strand'] = None
            gene_dict['name'] = name_match.search(entries[8]).group(1)
            gene_dict['description'] = description_match.search(entries[8]).group(1)
            if entry_type in gene_entry_types[0:2]:
                gene_dict['intergenic'] = False
            else:
                gene_dict['intergenic'] = True
            gene_matches = (session.query(SqD.Gene)
                                   .filter(SqD.Gene.name==gene_dict['name']))
            if gene_matches.count() == 0:
                session.add(SqD.Gene(**gene_dict))
            elif gene_matches.count() == 1:
                gene_indb = gene_matches.first()
                gene_indb.start = gene_dict['start']
                gene_indb.end = gene_dict['end']
                gene_indb.strand = gene_dict['strand']
            else:
                raise RuntimeError("uh-oh I didn't expect this was possible")
        else:
            pass
try:
    session.commit()
    session.close()
except DBAPIError:
    session.rollback()
    session.close()
    raise DBAPIError
