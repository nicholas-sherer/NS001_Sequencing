#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 13:09:18 2019

@author: nicholas
"""

from sqlalchemy import create_engine, event
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

import SequenceDataORM as SqD


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()


engine = create_engine('sqlite:///NS001_evolved_mutations.db', echo=True)
SqD.Base.metadata.create_all(engine)
session = sessionmaker(bind=engine)()

# Insert Wells into database
wells_seqd = ['A1', 'A2', 'A3', 'A4', 'B2', 'B3', 'B4', 'B5', 'C3', 'C4', 'C5',
              'C6', 'D4', 'D5', 'D6', 'D7', 'E5', 'E6', 'E7', 'E8']
session.add_all([SqD.Well(position=position) for position in wells_seqd])
# Insert Mutation Conditions into database
session.add_all([SqD.MutationCondition(name='Hi', mu_est=2.2*10**-7,
                                       mu_lowCI=1.6*10**-7,
                                       mu_hiCI=2.9*10**-7),
                 SqD.MutationCondition(name='HiMid', mu_est=4.1*10**-8,
                                       mu_lowCI=2.6*10**-8,
                                       mu_hiCI=5.9*10**-8),
                 SqD.MutationCondition(name='Mid', mu_est=1.4*10**-8,
                                       mu_lowCI=.64*10**-8,
                                       mu_hiCI=2.5*10**-8),
                 SqD.MutationCondition(name='LoMid', mu_est=3.8*10**-9,
                                       mu_lowCI=1.2*10**-9,
                                       mu_hiCI=7.4*10**-9),
                 SqD.MutationCondition(name='Lo', mu_est=1.7*10**-9,
                                       mu_lowCI=.41*10**-9,
                                       mu_hiCI=3.6*10**-9)])
# Put strains into the database
session.add(SqD.Strain(name='NS001', ancestor=None, day_in_platereader=None,
                            well=None, mutation_rate_condition=None))
session.add_all([SqD.Strain(name='Hi{0}t1'.format(n), ancestor='NS001',
                            day_in_platereader=24, well=wells_seqd[n-1],
                            mutation_rate_condition='Hi')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='HiMid{0}t1'.format(n), ancestor='NS001',
                            day_in_platereader=24, well=wells_seqd[n+3],
                            mutation_rate_condition='HiMid')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='Mid{0}t1'.format(n), ancestor='NS001',
                            day_in_platereader=24, well=wells_seqd[n+7],
                            mutation_rate_condition='Mid')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='LoMid{0}t1'.format(n), ancestor='NS001',
                            day_in_platereader=24, well=wells_seqd[n+11],
                            mutation_rate_condition='LoMid')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='Lo{0}t1'.format(n), ancestor='NS001',
                            day_in_platereader=24, well=wells_seqd[n+15],
                            mutation_rate_condition='Lo')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='Hi{0}t2'.format(n),
                            ancestor='Hi{0}t1'.format(n),
                            day_in_platereader=41, well=wells_seqd[n-1],
                            mutation_rate_condition='Hi')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='HiMid{0}t2'.format(n),
                            ancestor='HiMid{0}t1'.format(n),
                            day_in_platereader=41, well=wells_seqd[n+3],
                            mutation_rate_condition='HiMid')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='Mid{0}t2'.format(n),
                            ancestor='Mid{0}t1'.format(n),
                            day_in_platereader=41, well=wells_seqd[n+7],
                            mutation_rate_condition='Mid')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='LoMid{0}t2'.format(n),
                            ancestor='LoMid{0}t1'.format(n),
                            day_in_platereader=41, well=wells_seqd[n+11],
                            mutation_rate_condition='LoMid')
                 for n in range(1, 5)])
session.add_all([SqD.Strain(name='Lo{0}t2'.format(n),
                            ancestor='Lo{0}t1'.format(n),
                            day_in_platereader=41, well=wells_seqd[n+15],
                            mutation_rate_condition='Lo')
                 for n in range(1, 5)])
# Put samples into the database
Conditions = ['Hi', 'HiMid', 'Mid', 'LoMid', 'Lo']
session.add(SqD.DNA_Sample(name='Aggregate_NS001_Ancestors', strain='NS001'))
session.add_all([SqD.DNA_Sample(name='{0}{1}t{2}_S1'.format(condition, n, t),
                                strain='{0}{1}t{2}'.format(condition, n, t))
                 for condition in Conditions
                 for n in range(1, 5)
                 for t in [1, 2]
                 ])
session.commit()
