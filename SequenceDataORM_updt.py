#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 19:19:32 2019

@author: nicholas
"""

from sqlalchemy import Column, Boolean, Integer, Float, String,\
    ForeignKey, ForeignKeyConstraint
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship, scoped_session, sessionmaker, backref
from sqlalchemy.orm.exc import NoResultFound


Base = declarative_base()
maker = sessionmaker(autoflush=True, autocommit=False)
DBSession = scoped_session(maker)


class Well(Base):
    __tablename__ = 'wells'
    position = Column(String, primary_key=True)
    strains = relationship('Strain', back_populates='well_r')

    def __repr__(self):
        return "<Well(position='{0}')>".format(self.position)


class MutationCondition(Base):
    __tablename__ = 'mutation_conditions'
    name = Column(String, primary_key=True)
    mu_est = Column(Float, nullable=False)
    mu_lowCI = Column(Float, nullable=False)
    mu_hiCI = Column(Float, nullable=False)
    strains = relationship('Strain', back_populates='mutation_condition')

    def __repr__(self):
        return "<MutationCondition(name='{0}')>".format(self.name)


class Gene(Base):
    __tablename__ = 'genes'
    name = Column(String, primary_key=True)
    description = Column(String)
    intergenic = Column(Boolean, nullable=False)
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(Boolean)
    SNPs = relationship('SNP_Mutation', back_populates='gene_r')

    def __repr__(self):
        return ("<Gene(name='{0}', intergenic={1})>").format(self.name,
                                                             self.intergenic)


class Strain(Base):
    __tablename__ = 'strains'
    name = Column(String, primary_key=True)
    ancestor = Column(String, ForeignKey('strains.name'))
    day_in_platereader = Column(Integer)
    well = Column(String, ForeignKey('wells.position'))
    mutation_rate_condition = Column(String,
                                     ForeignKey('mutation_conditions.name'))
    children = relationship('Strain', backref=backref('ancestor_r',
                                                      remote_side=[name]))
    well_r = relationship('Well', back_populates='strains')
    mutation_condition = relationship('MutationCondition',
                                      back_populates='strains')
    dna_samples = relationship('DNA_Sample', back_populates='strain_r')

    def __repr__(self):
        return ("<Strain(name={0})>").format(self.name)


class SNP_Mutation(Base):
    __tablename__ = 'snp_mutations'
    chr_position = Column(Integer, primary_key=True)
    ref_base = Column(String(1), primary_key=True)
    new_base = Column(String(1), primary_key=True)
    gene = Column(String, ForeignKey('genes.name'), nullable=False)
    gene_pos = Column(String)
    ref_aa = Column(String)
    new_aa = Column(String)
    intergenic_left = Column(Integer)
    intergenic_right = Column(Integer)
    gene_r = relationship('Gene', back_populates='SNPs')
    SNPs_evidence = relationship('SNP_Evidence', back_populates='mutations')
    smpls = relationship('DNA_Sample',
                         secondary=lambda: SNP_Evidence.__table__,
                         viewonly=True)
    samples = association_proxy('smpls', 'name')

    @hybrid_property
    def mutation(self):
        return self.chr_position, self.ref_base, self.new_base

    @classmethod
    def by_mutation(cls, mutation):
        try:
            obj = DBSession.query(cls).filter_by(mutation=mutation).one()
        except NoResultFound:
            # obj = cls(iso_code=iso_code)
            # DBSession.add(obj)
            pass
        return obj

    @hybrid_property
    def synonymous(self):
        return self.ref_aa == self.new_aa

    def __repr__(self):
        return ("<SNP_Mutation(chr_position={0}, ref_base={1}, new_base={2}, "
                "gene={3})>").format(self.chr_position, self.ref_base,
                                     self.new_base, self.gene)


class DNA_Sample(Base):
    __tablename__ = 'dna_samples'
    name = Column(String, primary_key=True)
    strain = Column(String, ForeignKey('strains.name'), nullable=False)
    strain_r = relationship('Strain', back_populates='dna_samples')
    SNPs_evidence = relationship('SNP_Evidence', back_populates='samples')
    MCs_evidence = relationship('MC_Evidence', back_populates='samples')
    muts = relationship('SNP_Mutation',
                        secondary=lambda: SNP_Evidence.__table__,
                        viewonly=True)
    snp_mutations = association_proxy('muts', 'mutation',
                                      creator=SNP_Mutation.by_mutation)

    def snp_mutation_frequency(self, snp, session):
        if type(snp) is SNP_Mutation:
            chr_position = snp.chr_position
            ref_base = snp.ref_base
            new_base = snp.new_base
        elif type(snp) is tuple:
            if len(snp) == 3:
                chr_position = snp[0]
                ref_base = snp[1]
                new_base = snp[2]
            else:
                raise ValueError('snp must be an SNP_Mutation or a '
                                 'tuple of its 3 primary keys')
        else:
            raise ValueError('snp must be an SNP_Mutation or a '
                             'tuple of its 3 primary keys')
        if (chr_position, ref_base, new_base) in self.snp_mutations:
            evidence = (session.query(SNP_Evidence)
                               .filter(SNP_Evidence.sample == self.name)
                               .filter(SNP_Evidence.chr_position ==
                                       chr_position)
                               .filter(SNP_Evidence.ref_base == ref_base)
                               .filter(SNP_Evidence.new_base == new_base)
                               .first())
            return evidence.frequency
        else:
            return 0

    def __repr__(self):
        return "<DNA_Sample(name={0}, strain={1})>".format(self.name,
                                                           self.strain)


class SNP_Evidence(Base):
    __tablename__ = 'snp_evidence'
    sample = Column(String, ForeignKey('dna_samples.name'), primary_key=True)
    chr_position = Column(Integer, primary_key=True)
    ref_base = Column(String(1), primary_key=True)
    new_base = Column(String(1), primary_key=True)
    frequency = Column(Float)
    bias_e_value = Column(Float)
    bias_p_value = Column(Float)
    consensus_score = Column(Float)
    fisher_strand_p_value = Column(Float)
    ref_cov = Column(String)
    new_cov = Column(String)
    total_cov = Column(String)
    __table_args__ = (ForeignKeyConstraint([chr_position, ref_base, new_base],
                                           [SNP_Mutation.chr_position,
                                            SNP_Mutation.ref_base,
                                            SNP_Mutation.new_base]), {})

    samples = relationship('DNA_Sample', back_populates='SNPs_evidence')
    mutations = relationship('SNP_Mutation', back_populates='SNPs_evidence')

    @hybrid_property
    def bs_total_cov(self):
        return sum([int(x) for x in self.total_cov.split('/')])

    @hybrid_property
    def bs_ref_cov(self):
        return sum([int(x) for x in self.ref_cov.split('/')])

    @hybrid_property
    def bs_new_cov(self):
        return sum([int(x) for x in self.new_cov.split('/')])

    def __repr__(self):
        return ("<SNP_Evidence(sample={0}, chr_position={1}, ref_base={2}, "
                "new_base={3})>").format(self.sample, self.chr_position,
                                         self.ref_base, self.new_base)


class MC_Evidence(Base):
    __tablename__ = 'mc_evidence'
    sample = Column(String, ForeignKey('dna_samples.name'), primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)
    start_range = Column(Integer, nullable=False)
    end_range = Column(Integer, nullable=False)
    left_inside_cov = Column(Integer)
    left_outside_cov = Column(Integer)
    right_inside_cov = Column(Integer)
    right_outside_cov = Column(Integer)
    gene = Column(String)
    samples = relationship('DNA_Sample', back_populates='MCs_evidence')

    @hybrid_property
    def length(self):
        return self.end - self.start

    def __repr__(self):
        return ("<MC_Evidence(sample={0}, start={1}, "
                "end={2})>").format(self.sample, self.start, self.end)
