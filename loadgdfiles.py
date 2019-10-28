# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from collections import OrderedDict
import re
from sqlalchemy.exc import DBAPIError

import SequenceDataORM as SqD


def gd_line_to_dict(gd_line):
    if gd_line[0] == '#':
        raise ValueError("Comment lines aren't mutations or evidence for one")
    entries = gd_line.strip('\n').split('\t')
    gd_dict = OrderedDict()
    gd_dict['type'] = entries[0]
    gd_dict['id'] = int(entries[1])
    try:
        gd_dict['parent_id'] = int(entries[2])
    except ValueError:
        gd_dict['parent_id'] = entries[2]
    # types of genome diff mutations are SNP, SUB, DEL, INS, MOB, AMP, CON, INV
    # types of genome diff evidence are RA, MC, JC, UN
    # this doesn't handle validation types
    type_ = gd_dict['type']
    if type_ == 'SNP':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['new_seq'] = entries[5]
        next_idx = 6
    elif type_ == 'SUB':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        gd_dict['new_seq'] = entries[6]
        next_idx = 7
    elif type_ == 'DEL':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        next_idx = 6
    elif type_ == 'INS':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = entries[4]
        gd_dict['new_seq'] = entries[5]
        next_idx = 6
    elif type_ == 'MOB':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['repeat_name'] = entries[5]
        gd_dict['strand'] = entries[6]
        gd_dict['duplication_size'] = int(entries[7])
        next_idx = 8
    elif type_ == 'AMP':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        gd_dict['new_copy_number'] = int(entries[6])
        next_idx = 7
    elif type_ == 'CON':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        gd_dict['region'] = entries[6]
        next_idx = 7
    elif type_ == 'INV':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        next_idx = 6
    elif type_ == 'RA':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['insert_position'] = int(entries[5])
        gd_dict['ref_base'] = entries[6]
        gd_dict['new_base'] = entries[7]
        next_idx = 8
    elif type_ == 'MC':
        gd_dict['seq_id'] = entries[3]
        gd_dict['start'] = int(entries[4])
        gd_dict['end'] = int(entries[5])
        gd_dict['start_range'] = int(entries[6])
        gd_dict['end_range'] = int(entries[7])
        next_idx = 8
    elif type_ == 'JC':
        gd_dict['side_1_seq_id'] = entries[3]
        gd_dict['side_1_position'] = int(entries[4])
        gd_dict['side_1_strand'] = entries[5]
        gd_dict['side_2_seq_id'] = entries[6]
        gd_dict['side_2_position'] = int(entries[7])
        gd_dict['side_2_strand'] = entries[8]
        gd_dict['overlap'] = int(entries[9])
        next_idx = 10
    elif type_ == 'UN':
        gd_dict['seq_id'] = entries[3]
        gd_dict['start'] = int(entries[4])
        gd_dict['end'] = int(entries[5])
        next_idx = 6
    else:
        raise RuntimeError('unknown .gd entry type')
    for entry in entries[next_idx:]:
        try:
            key, value = entry.split('=')
            gd_dict[key] = value
        except ValueError:
            print('invalid key-value entry in .gd line')
    return gd_dict


def RA_dict_to_SNP_dicts(sample, gd_dict):
    # first unpack what we need to look up or make the entry in the genes table
    gene = {}
    gene['name'] = gd_dict.get('gene_name')
    if gene['name'] is None:
        return None
    gene['description'] = gd_dict.get('gene_product')
    gene['intergenic'] = (gd_dict['snp_type'] == 'intergenic')
    # next unpack what we need to look up or make the SNP_mutations entry
    snp_mutation = {}
    snp_mutation['chr_position'] = gd_dict.get('position')
    snp_mutation['ref_base'] = gd_dict.get('ref_base')
    snp_mutation['new_base'] = gd_dict.get('new_base')
    if (type(snp_mutation['chr_position']) is not int) or \
       (snp_mutation['ref_base'] is None) or \
       (snp_mutation['new_base'] is None):
        return None
    snp_mutation['gene'] = gene['name']
    snp_mutation['gene_pos'] = gd_dict.get('gene_position')
    snp_mutation['ref_aa'] = gd_dict.get('aa_ref_seq')
    snp_mutation['new_aa'] = gd_dict.get('aa_new_seq')
    if gene['intergenic'] and (snp_mutation['gene_pos'] is not None):
        nums = re.findall('[-+]?[0-9]+', snp_mutation['gene_pos'])
        if len(nums) == 2:
            snp_mutation['intergenic_left'] = int(nums[0])
            snp_mutation['intergenic_right'] = int(nums[1])
        else:
            snp_mutation['intergenic_left'] = None
            snp_mutation['intergenic_right'] = None
    else:
        snp_mutation['intergenic_left'] = None
        snp_mutation['intergenic_right'] = None
    # finally unpack what we need to look up or make the SNP_evidence entry
    snp_evidence = {}
    snp_evidence['sample'] = sample
    snp_evidence['chr_position'] = snp_mutation['chr_position']
    snp_evidence['ref_base'] = snp_mutation['ref_base']
    snp_evidence['new_base'] = snp_mutation['new_base']
    snp_evidence['frequency'] = gd_dict.get('frequency')
    snp_evidence['bias_e_value'] = gd_dict.get('bias_e_value')
    snp_evidence['bias_p_value'] = gd_dict.get('bias_p_value')
    snp_evidence['consensus_score'] = gd_dict.get('consensus_score')
    snp_evidence['fisher_strand_p_value'] = \
        gd_dict.get('fisher_strand_p_value')
    snp_evidence['ref_cov'] = gd_dict.get('ref_cov')
    snp_evidence['new_cov'] = gd_dict.get('new_cov')
    snp_evidence['total_cov'] = gd_dict.get('total_cov')
    return gene, snp_mutation, snp_evidence


def commit_SNP_dicts(gene, snp_mutation, snp_evidence, session):
    # check if the gene is already in the database. If not, add it.
    gene_matches = session.query(SqD.Gene).filter(SqD.Gene.name ==
                                                  gene['name'])
    if gene_matches.count() == 0:
        session.add(SqD.Gene(**gene))
    # check if the snp_mutation is already in the database. If not, add it.
    mu_query = session.query(SqD.SNP_Mutation)
    mu_matches = mu_query.filter(SqD.SNP_Mutation.chr_position ==
                                 snp_mutation['chr_position'])
    mu_matches = mu_matches.filter(SqD.SNP_Mutation.ref_base ==
                                   snp_mutation['ref_base'])
    mu_matches = mu_matches.filter(SqD.SNP_Mutation.new_base ==
                                   snp_mutation['new_base'])
    if mu_matches.count() == 0:
        session.add(SqD.SNP_Mutation(**snp_mutation))
    # check if the snp_evidence is already in the database. If not, add it.
    ev_query = session.query(SqD.SNP_Evidence)
    ev_matches = ev_query.filter(SqD.SNP_Evidence.sample ==
                                 snp_evidence['sample'])
    ev_matches = ev_matches.filter(SqD.SNP_Evidence.chr_position ==
                                   snp_evidence['chr_position'])
    ev_matches = ev_matches.filter(SqD.SNP_Evidence.ref_base ==
                                   snp_evidence['ref_base'])
    ev_matches = ev_matches.filter(SqD.SNP_Evidence.new_base ==
                                   snp_evidence['new_base'])
    if ev_matches.count() == 0:
        session.add(SqD.SNP_Evidence(**snp_evidence))
    # try committing the changes to the database. If it fails, rollback.
    try:
        session.commit()
    except DBAPIError:
        session.rollback()
